# -*- coding: utf-8 -*-

from functools import lru_cache

import numpy as np

from .affine import change_of_basis, transform
from .atom_data import (
    ELEM_TO_MAGMOM,
    ELEM_TO_MASS,
    ELEM_TO_NAME,
    ELEM_TO_NUM,
    NUM_TO_ELEM,
)
from .lattice import Lattice


# TODO: store atomic data as class attributes?
class Atom(object):
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
    """

    __slots__ = (
        "element",
        "coords_fractional",
        "displacement",
        "magmom",
        "occupancy",
        "lattice",
    )

    def __init__(
        self,
        element,
        coords,
        lattice=None,
        displacement=None,
        magmom=None,
        occupancy=1.0,
        **kwargs
    ):
        if isinstance(element, int):
            element = NUM_TO_ELEM[element]
        elif element not in ELEM_TO_NUM:
            raise ValueError("Invalid chemical element {}".format(element))

        if magmom is None:
            magmom = ELEM_TO_MAGMOM[element]

        self.element = element.title()
        self.coords_fractional = np.asfarray(coords)
        self.lattice = lattice or Lattice(np.eye(3))
        self.displacement = np.asfarray(
            displacement if displacement is not None else (0, 0, 0)
        )
        self.magmom = magmom
        self.occupancy = occupancy

    def __repr__(self):
        return "< Atom {:<2} @ ({:.2f}, {:.2f}, {:.2f}) >".format(
            self.element, *tuple(self.coords_fractional)
        )

    def __eq__(self, other):
        if type(other) is type(self):
            return (
                (self.element == other.element)
                and (self.magmom == other.magmom)
                and np.allclose(
                    self.coords_fractional, other.coords_fractional, atol=1e-3
                )
                # Lattice information is encoded in coords_cartesian
                # We don't compare lattice directly because that might cause a recursive check if
                # self.lattice is something else, like Supercell or Crystal
                and np.allclose(
                    self.coords_cartesian, other.coords_cartesian, atol=1e-3
                )
                and np.allclose(self.displacement, other.displacement, atol=1e-3)
            )
        return NotImplemented

    def __hash__(self):
        return hash(
            (
                self.element,
                self.magmom,
                tuple(np.round(self.coords_fractional, 3)),
                tuple(np.round(self.coords_cartesian, 3)),
                tuple(np.round(self.displacement, 3)),
            )
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
    def atomic_number(self):
        return ELEM_TO_NUM[self.element]

    @property
    def mass(self):
        return ELEM_TO_MASS[self.element]

    @property
    def coords_cartesian(self):
        """ 
        Real-space position of the atom on the lattice.
                    
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

        return Atom(
            element=self.element,
            coords=coords_fractional,
            lattice=self.lattice,
            displacement=self.displacement,
            magmom=self.magmom,
            occupancy=self.occupancy,
        )

    def __array__(self, *args, **kwargs):
        """ Returns an array [Z, x, y, z] """
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
    Calculates and sets the real-space coordinates of the atom from fractional coordinates and lattice vectors.
    Only valid for inorganic compounds.
    
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
    return np.mod(transform(COB, real_coords), 1)


def distance_fractional(atm1, atm2):
    """ 
    Calculate the distance between two atoms in fractional coordinates.
    
    Parameters
    ----------
    atm1, atm2 : ``crystals.Atom``

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
    atm1, atm2 : ``crystals.Atom``

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
