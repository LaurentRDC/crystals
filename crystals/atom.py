# -*- coding: utf-8 -*-

import numpy as np
from enum import Enum, auto, unique
from collections import namedtuple, OrderedDict
from copy import deepcopy

from .affine import change_of_basis, transform
from .atom_data import (
    ELEM_TO_MAGMOM,
    ELEM_TO_MASS,
    ELEM_TO_NAME,
    NAME_TO_ELEM,
    ELEM_TO_NUM,
    NUM_TO_ELEM,
    chemical_symbols,
    atomic_names,
)
from .lattice import Lattice


# TODO: make this class an enumeration?
class Element:
    """
    Class representing an abtract chemical element, but no particular atom.
    This class gives access to elemental properties, like atomic number,
    atomic mass, full element name, etc.

    Parameters
    ----------
    element : str, int, or :class:`Element`
        Elemental symbol (e.g. "He"), element name (e.g. "Helium"),
        atomic number, or another `Element` instance.

    Raises
    ------
    ValueError : if the element is not valid.
    """

    valid_symbols = frozenset(chemical_symbols)
    valid_names = frozenset(atomic_names)

    def __init__(self, element, *args, **kwargs):
        if isinstance(element, int):
            try:
                element = NUM_TO_ELEM[element]
            except KeyError:
                raise ValueError(f"Atomic number {element} not supported.")
        elif isinstance(element, Element):
            element = element.symbol

        # At this point, `element` is a string
        element = str(element).title()

        valid_string_inputs = self.valid_names.union(self.valid_symbols)
        if element not in valid_string_inputs:
            raise ValueError(f"Element {element} is not valid.")

        if element in self.valid_symbols:
            self.element = element
        elif element in self.valid_names:
            self.element = NAME_TO_ELEM[element]

    def __str__(self):
        return self.symbol

    def __repr__(self):
        return f"< {self.name} >"

    def __eq__(self, other):
        if isinstance(other, Element):
            return self.element == other.element
        return NotImplemented

    def __hash__(self):
        # Technically, if atomic_number is an int, hash(atomic_number) = atomic_number
        # However, just in case this won't be true in the future, we still use the hash() function
        return hash(self.atomic_number)

    @property
    def element_full(self):
        """Full element name, e.g. "Hydrogen" """
        return self.name

    @property
    def name(self):
        """Full element name, e.g. "Hydrogen" """
        return ELEM_TO_NAME[self.element]

    @property
    def symbol(self):
        """Elemental symbol, e.g. "He" """
        return self.element

    @property
    def atomic_number(self):
        """Atomic number"""
        return ELEM_TO_NUM[self.element]

    @property
    def mass(self):
        """Atomic mass [u]"""
        return ELEM_TO_MASS[self.element]

    @property
    def magnetic_moment_ground(self):
        """Ground state magnetic moment."""
        return ELEM_TO_MAGMOM.get(self.element, None)


class Atom(Element):
    """
    Container object for atomic data.

    Parameters
    ----------
    element : str or int
        Chemical element symbol or atomic number.
    coords : array-like, shape (3,)
        Coordinates of the atom in fractional form.
    lattice : Lattice or array-like, shape (3,3)
        Lattice on which the atom is positioned.
    displacement : array-like or None, optional
        Atomic maximum displacement [Angs].
    magmom : float, optional
        Magnetic moment. If None (default), the ground-state magnetic moment is used.
    occupancy : float, optional
        Fractional occupancy. If None (default), occupancy is set to 1.0.
    tag : int or None, optional
        Tag an atom with a unique identifier. Useful to keep track of atom order, for example
        in PWSCF output files. This is mostly for internal use.
    electronic_structure : ElectronicStructure or None, optional
        Electronic orbital structure for this atom. If `None` (default), the ground
        state for this element will be used.
    """

    # Because of the possibility of a large number of atoms (> 1e6), we use the __slots__
    # mechanism to make Atom objects as small as possible.
    __slots__ = (
        "element",
        "coords_fractional",
        "displacement",
        "magmom",
        "occupancy",
        "lattice",
        "electronic_structure",
    )

    def __init__(
        self,
        element,
        coords,
        lattice=None,
        displacement=None,
        magmom=None,
        occupancy=1.0,
        tag=None,
        electronic_structure=None,
        **kwargs,
    ):
        super().__init__(element=element)

        self.coords_fractional = np.asfarray(coords)
        self.lattice = lattice or Lattice(np.eye(3))
        self.displacement = np.asfarray(
            displacement if displacement is not None else (0, 0, 0)
        )
        self.magmom = magmom or self.magnetic_moment_ground
        self.occupancy = occupancy
        self.tag = tag

        # We distinguish between default (None) and other electronic structures
        # So that the atoms can be represented appropriately in __repr__
        self.electronic_structure = (
            ElectronicStructure.ground_state(self.element) or electronic_structure
        )

    def __repr__(self):
        x, y, z = tuple(self.coords_fractional)
        r = f"< Atom {self.element:<2} @ ({x:.2f}, {y:.2f}, {z:.2f})"
        # No point in polluting the representation if the electronic structure was not
        # set (i.e. default)
        if self.electronic_structure != ElectronicStructure.ground_state(self.element):
            r += f" | [{str(self.electronic_structure)}]"
        r += " >"
        return r

    def __eq__(self, other):
        if isinstance(other, Atom):
            return (
                super().__eq__(other)
                and (self.magmom == other.magmom)
                and distance_fractional(self, other) < 1e-3
                and (self.lattice == other.lattice)
                and np.allclose(self.displacement, other.displacement, atol=1e-3)
                and self.tag == other.tag
                and (self.electronic_structure == other.electronic_structure)
            )
        return False

    def __hash__(self):
        return hash(
            (
                super().__hash__(),
                self.magmom,
                tuple(np.round(self.coords_fractional, 4)),
                self.lattice,
                tuple(np.round(self.displacement, 3)),
                self.tag,
                hash(self.electronic_structure),
            )
        )

    def __lt__(self, other):
        # Here's a cool fact:
        # Sorting n-tuples works as follows:
        #   1. Sort by the first elements
        #   2. If the first elements are equal, sort by the second elements;
        #   3. and so on
        # This means that by sorting atoms by (a.element, c1, c2 ,c3),
        # we can have a stable order (i.e. sorted by fractional coordinates)
        # but with elements grouped together.
        # The part where sorting is stable is important because of doctests,
        # which are quite literal.
        return (self.element, *self.coords_fractional) < (
            other.element,
            *other.coords_fractional,
        )

    @classmethod
    def from_ase(cls, atom):
        """
        Returns an Atom instance from an ASE atom

        Parameters
        ----------
        atom : ase.Atom
        """
        lattice = np.eye(3)
        if atom.atoms is not None:
            lattice = np.array(atom.atoms.cell)

        return cls(
            element=atom.symbol,
            coords=frac_coords(atom.position, lattice),
            magmom=atom.magmom,
        )

    @property
    def coords_cartesian(self):
        """
        Real-space position of the atom on the lattice, in Angstroms.

        Returns
        -------
        pos : `~numpy.ndarray`, shape (3,)
            Atomic position

        Raises
        ------
        RuntimeError : if this atom is not place on a lattice
        """
        return real_coords(self.coords_fractional, self.lattice.lattice_vectors)

    def transform(self, *matrices):
        """
        Return an Atom with fractional coordinates transformed according to symmetry operators.

        Parameters
        ----------
        matrices : ndarrays, shape {(3,3), (4,4)}
            Transformation matrices.

        Returns
        -------
        atm : Atom
            Transformed atom. The original atom is left unchanged.
        """
        coords_fractional = np.array(self.coords_fractional, copy=True)

        for matrix in matrices:
            coords_fractional = transform(matrix, coords_fractional)

        new_atom = deepcopy(self)
        new_atom.coords_fractional[:] = coords_fractional

        return new_atom

    def __array__(self, *args, **kwargs):
        """Returns an array [Z, x, y, z]"""
        arr = np.empty(shape=(4,), *args, **kwargs)
        arr[0] = self.atomic_number
        arr[1::] = self.coords_fractional
        return arr


def real_coords(frac_coords, lattice_vectors):
    """
    Calculates the real-space coordinates of the atom from fractional coordinates and lattice vectors.

    Parameters
    ----------
    frac_coords : array-like, shape (3,)
        Fractional coordinates
    lattice_vectors : list of ndarrays, shape (3,)
        Lattice vectors of the crystal.

    Returns
    -------
    coords : ndarray, shape (3,)
    """
    COB = change_of_basis(np.array(lattice_vectors), np.eye(3))
    return transform(COB, frac_coords)


def frac_coords(real_coords, lattice_vectors):
    """
    Calculates the fractional coordinates of the atom from real-space coordinates and lattice vectors.

    Parameters
    ----------
    real_coords : array-like, shape (3,)
        Real-space coordinates
    lattice_vectors : list of ndarrays, shape (3,)
        Lattice vectors of the crystal.

    Returns
    -------
    coords : ndarray, shape (3,)
    """
    COB = change_of_basis(np.eye(3), np.array(lattice_vectors))
    return transform(COB, real_coords)


def distance_fractional(atm1, atm2):
    """
    Calculate the distance between two atoms in fractional coordinates.

    Parameters
    ----------
    atm1, atm2 : :class:`Atom`

    Returns
    -------
    dist : float
        Fractional distance between atoms.

    Raises
    ------
    RuntimeError : if atoms are not associated with the same lattice.
    """
    if atm1.lattice != atm2.lattice:
        raise RuntimeError(
            "Distance is undefined if atoms are sitting on different lattices."
        )

    return np.linalg.norm(atm1.coords_fractional - atm2.coords_fractional)


def distance_cartesian(atm1, atm2):
    """
    Calculate the distance between two atoms in cartesian coordinates.

    Parameters
    ----------
    atm1, atm2 : :class:`Atom`

    Returns
    -------
    dist : float
        Cartesian distance between atoms in Angstroms..

    Raises
    ------
    RuntimeError : if atoms are not associated with the same lattice.
    """
    if atm1.lattice != atm2.lattice:
        raise RuntimeError(
            "Distance is undefined if atoms are sitting on different lattices."
        )

    return np.linalg.norm(atm1.coords_cartesian - atm2.coords_cartesian)


def is_element(element):
    """
    Create a function that checks whether an atom is of a certain element.

    Parameters
    ----------
    element : str, int, or Element
        Elemental symbol (e.g. "He"), atomic number, or Element instance.

    Returns
    -------
    func : callable
        Returns a function that can be used to check whether a `crystals.Atom`
        instance is of a certain element.

    Examples
    --------
    >>> is_vanadium = is_element('V') # is_vanadium is a function
    >>> atm = Atom('V', [0,0,0])
    >>> is_vanadium(atm)
    True
    """
    if not isinstance(element, Element):
        element = Element(element)

    def _is_element(atm):
        return atm.atomic_number == element.atomic_number

    return _is_element


@unique
class Orbital(Enum):
    """
    Enumeration of electronic orbitals, used to described atomic
    orbital structure.

    We note that `Orbital` instances are ordered according to the Madelung rule.

    Examples
    --------
    >>> Orbital("1s")
    <Orbital.one_s: '1s'>
    >>> Orbital.four_p == Orbital("4p")
    True

    """

    # It is important that the Orbitals are listed in the order that they
    # are filled (Madelung rule) because this is how __lt__ is defined.
    one_s = "1s"
    two_s = "2s"
    two_p = "2p"
    three_s = "3s"
    three_p = "3p"
    four_s = "4s"
    three_d = "3d"
    four_p = "4p"
    five_s = "5s"
    four_d = "4d"
    five_p = "5p"
    six_s = "6s"
    four_f = "4f"
    five_d = "5d"
    six_p = "6p"
    seven_s = "7s"
    five_f = "5f"
    six_d = "6d"
    seven_p = "7p"

    def __lt__(self, other):
        madelung_rule = [v for k, v in Orbital.__members__.items()]
        return madelung_rule.index(self) < madelung_rule.index(other)

    @classmethod
    def maximum_electrons(cls, shell):
        """
        Maximum number of electrons that can be placed in an orbital.

        Parameters
        ----------
        shell : Orbital or str

        Returns
        -------
        max : int
        """
        shell = Orbital(shell)
        maxima = {
            "s": 2,
            "p": 6,
            "d": 10,
            "f": 14,
        }
        return maxima[shell.value[-1]]


# To print superscript in electronic structures.
# Note that this requires UTF-8 output.
superscript_map = {
    "0": "⁰",
    "1": "¹",
    "2": "²",
    "3": "³",
    "4": "⁴",
    "5": "⁵",
    "6": "⁶",
    "7": "⁷",
    "8": "⁸",
    "9": "⁹",
}
superscript_trans = str.maketrans(
    "".join(superscript_map.keys()), "".join(superscript_map.values())
)


class ElectronicStructure:
    """
    Description of the atomic orbital structure.

    Parameters
    ----------
    shells : dict[Orbital,int]
        Dictionary containing the number of electrons in each Orbital, e.g. `{"1s": 2}`.

    Raises
    ------
    ValueError : if the electronic structure is not representable

    Examples
    --------
    Electronic structures can be specified by hand:

    >>> ElectronicStructure({"1s":2, "2s":2, "2p":2})
    < ElectronicStructure: 1s²2s²2p² >

    A shortcut exists for atomic ground states:

    >>> ElectronicStructure.ground_state("Ti")
    < ElectronicStructure: 1s²2s²2p⁶3s²3p⁶4s²3d² >

    Notes
    -----
    Shells are allowed to be filled our of order deliberately, given that unusual
    electronic structures can arise from ultrafast photoexcitation.
    """

    def __init__(self, shells):
        shells = {Orbital(k): v for k, v in shells.items()}

        # Subclassing OrderedDict causes problems with pickling
        # Instead, we dress this class on top of an OrderedDict property.
        self._structure = OrderedDict([])
        for k, v in shells.items():
            self.__setitem__(k, v)

    def __setitem__(self, key, value):
        # We check that the number of electrons in each Orbital does not
        # go above maximum possible.
        shell = Orbital(key)
        maximum_allowed_electrons = Orbital.maximum_electrons(shell)
        if value > maximum_allowed_electrons:
            raise ValueError(
                f"There cannot be {value} electrons in orbital {shell.value}"
            )
        self._structure.__setitem__(shell, value)

        # We ensure that orbital order is maintained
        # by re-creating the OrderedDict, filled
        self._structure = OrderedDict(sorted(self._structure.items()))

    def __getitem__(self, key):
        # In case the key doesn't exist, we return 0
        # (i.e. 0 electrons in this orbital) because this allows
        # to add electrons in-place, e.g.:
        # >>> struct = ElectronicStructure({"1s":2})
        # >>> struct["2p"] += 1
        # even though there were no electrons there.
        key = Orbital(key)
        try:
            return self._structure.__getitem__(key)
        except KeyError:
            return 0

    def __str__(self):
        result = ""
        for shell, occ in self._structure.items():
            result += shell.value + f"{occ}".translate(superscript_trans)
        return result

    def __repr__(self):
        return f"< ElectronicStructure: {str(self)} >"

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

    @classmethod
    def ground_state(cls, element):
        """
        Ground state electronic structure for a particular element.

        Parameters
        ----------
        element : Element, str, or int
            Element, Symbol, or atomic number.

        Returns
        -------
        structure : ElectronicStructure
            Return the ground state electronic structure for a particular element.

        Examples
        --------
        >>> ElectronicStructure.ground_state("Ne")
        < ElectronicStructure: 1s²2s²2p⁶ >
        """
        element = Element(element)
        num_elec = element.atomic_number

        structure = dict()
        for shell in Orbital:
            shell_elec = min([Orbital.maximum_electrons(shell), num_elec])
            structure[shell] = shell_elec
            num_elec -= shell_elec

            if num_elec == 0:
                break

        return cls(structure)

    @property
    def outer_shell(self):
        """The outermost shell, or valence orbital."""
        shells = list(self._structure.keys())
        return shells[-1]
