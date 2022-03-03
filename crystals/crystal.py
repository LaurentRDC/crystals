# -*- coding: utf-8 -*-
from copy import deepcopy
from enum import Enum, unique
from operator import attrgetter
from functools import lru_cache
from itertools import islice, product, chain, combinations
from pathlib import Path

import numpy as np
from spglib import (
    get_error_message,
    get_spacegroup_type,
    get_symmetry,
    get_symmetry_dataset,
    standardize_cell,
)

from .affine import affine_map, change_of_basis
from .atom import Atom
from .atom_data import chemical_symbols
from .base import AtomicStructure
from .lattice import Lattice
from .parsers import (
    CIFParser,
    CODParser,
    MPJParser,
    PDBParser,
    POSCARParser,
    PWSCFParser,
)
from .spg_data import Hall2HM
from .writers import write_cif, write_poscar, write_xyz, ase_atoms

CIF_ENTRIES = frozenset((Path(__file__).parent / "cifs").glob("*.cif"))
LONGEST_CHEMICAL_SYMBOL = max(len(s) for s in chemical_symbols)

is_atom = lambda a: isinstance(a, Atom)
is_structure = lambda s: isinstance(s, AtomicStructure)


class Crystal(AtomicStructure, Lattice):
    """
    The :class:`Crystal` class is a set-like container that represent
    crystalline structures. In addition to constructing the :class:`Crystal`
    object yourself, other constructors are also available
    (and preferred):

    * :meth:`Crystal.from_cif`: create an instance from a CIF file;

    * :meth:`Crystal.from_pdb`: create an instance from a Protein Data Bank entry;

    * :meth:`Crystal.from_database`: create an instance from the internal database of CIF files;

    * :meth:`Crystal.from_cod`: create an instance from a Crystallography Open Database entry;

    * :meth:`Crystal.from_mp`: create an instance from the Materials Project database;

    * :meth:`Crystal.from_pwscf`: create an instance from the output of the PWSCF program;

    * :meth:`Crystal.from_ase`: create an instance from an ``ase.Atoms`` instance;

    * :meth:`Crystal.from_poscar`: create an instance from VASP POSCAR files.

    Parameters
    ----------
    unitcell : iterable of :class:`Atom` or :class:`AtomicStructure`
        Unit cell atoms or substructures. It is assumed that the atoms are
        in fractional coordinates.
    lattice_vectors : iterable of array_like
        Lattice vectors. If ``lattice_vectors`` is provided as a 3x3 array, it
        is assumed that each lattice vector is a row.
    source : str or None, optional
        Provenance, e.g. filename. Only used for bookkeeping.
    """

    builtins = frozenset(map(attrgetter("stem"), CIF_ENTRIES))

    def __init__(self, unitcell, lattice_vectors, source=None, **kwargs):
        unitcell = list(unitcell)
        # Atoms need to be modified BEFORE they are fed to the constructor
        # of AtomicStructure
        for atom in filter(is_atom, unitcell):
            atom.lattice = Lattice(lattice_vectors)

        super().__init__(
            atoms=filter(is_atom, unitcell),
            substructures=filter(is_structure, unitcell),
            lattice_vectors=lattice_vectors,
            **kwargs,
        )

        self.source = source

    def __eq__(self, other):
        if isinstance(other, Crystal):
            # The explicit comparison of lattices sidesteps the problem
            # of evaluating NotImplemented in a bool context, which has been deprecated
            # as of Python 3.9
            return (
                Lattice(self.lattice_vectors) == Lattice(other.lattice_vectors)
            ) and super(AtomicStructure, self).__eq__(other)
        return NotImplemented

    def __hash__(self):
        return super(AtomicStructure, self).__hash__() | super(Lattice, self).__hash__()

    @property
    def unitcell(self):
        """Generator of atoms forming the crystal unit cell."""
        return super().__iter__()

    @lru_cache(
        maxsize=1
    )  # This operation is very heavy; better to cache it than recalculate.
    def asymmetric_cell(self):
        """
        Calculates the asymmetric cell that generates the crystal unit cell.

        .. versionadded:: 1.2.0
        """
        return symmetry_reduction(self.unitcell, self.symmetry_operations())

    @classmethod
    @lru_cache(maxsize=len(builtins))
    def from_cif(cls, path, **kwargs):
        """
        Returns a Crystal object created from a CIF 1.0, 1.1 or 2.0 file.
        Keyword arguments are passed to the Crystal constructor.

        Parameters
        ----------
        path : path-like
            File path
        """
        with CIFParser(filename=path) as parser:
            return cls(
                unitcell=symmetry_expansion(
                    parser.atoms(), parser.symmetry_operators()
                ),
                lattice_vectors=parser.lattice_vectors(),
                source=str(path),
                **kwargs,
            )

    @classmethod
    def from_database(cls, name, **kwargs):
        """
        Returns a Crystal object create from the internal CIF database.
        Keyword arguments are passed to the class constructor.

        Parameters
        ----------
        name : path-like
            Name of the database entry. Available items can be retrieved from `Crystal.builtins`
        """
        if name not in cls.builtins:
            raise ValueError(
                f"Entry {name} is not available in the database. See `Crystal.builtins` for valid entries."
            )

        path = Path(__file__).parent / "cifs" / (name + ".cif")
        return cls.from_cif(path, **kwargs)

    @classmethod
    def from_cod(cls, num, revision=None, download_dir=None, overwrite=False, **kwargs):
        """
        Returns a Crystal object built from the Crystallography Open Database.
        Keyword arguments are passed to the class constructor.

        Parameters
        ----------
        num : int
            COD identification number.
        revision : int or None, optional
            Revision number. If None (default), the latest revision is used.
        download_dir : path-like object, optional
            Directory where to save the CIF file. Default is a local folder in the current directory
        overwrite : bool, optional
            Whether or not to overwrite files in cache if they exist. If no revision
            number is provided, files will always be overwritten.
        """
        with CODParser(num, revision, download_dir, overwrite) as parser:
            return cls(
                unitcell=symmetry_expansion(
                    parser.atoms(), parser.symmetry_operators()
                ),
                lattice_vectors=parser.lattice_vectors(),
                source=f"COD num:{num} rev:{revision}",
                **kwargs,
            )

    @classmethod
    def from_mp(cls, query, api_key=None, download_dir=None, overwrite=False, **kwargs):
        """
        Returns a Crystal object built from the Materials Project.
        Keyword arguments are passed to the class constructor.

        Parameters
        ----------
        query : str
            The query can be a Materials Project material id (e.g., `"mp-1234"`), a
            formula, e.g. (`"Fe2O3"`), or a chemical system ("-" separated list of elemments,
            e.g., `"Li-Fe-O"`).
        api_key : str or None, optional
            An API key for accessing the Materials Project REST interface.
            Please obtain your API key at https://www.materialsproject.org/dashboard.
            If `None` (default), ``crystals`` will look for your API key in the
            `MATERIALS_PROJECT_API_KEY` environment variable.
        download_dir : path-like object or None, optional
            Directory where to save the CIF file. This is used for caching.
        overwrite : bool, optional
            Whether or not to overwrite files in cache if they exist. If True,
            a new file will be downloaded, possibly overwriting previously-downloaded file.
        """
        with MPJParser(
            query=query,
            api_key=api_key,
            download_dir=download_dir,
            overwrite=overwrite,
            **kwargs,
        ) as parser:
            return cls(
                unitcell=symmetry_expansion(
                    parser.atoms(), parser.symmetry_operators()
                ),
                lattice_vectors=parser.lattice_vectors(),
                source=f"Materials Project ID: {parser.material_id}",
                **kwargs,
            )

    @classmethod
    def from_pdb(cls, ID, download_dir=None, overwrite=False, **kwargs):
        """
        Returns a Crystal object created from a Protein DataBank entry.
        Keyword arguments are passed to the class constructor.

        Parameters
        ----------
        ID : str
            Protein DataBank identification. The correct .pdb file will be downloaded,
            cached and parsed.
        download_dir : path-like object, optional
            Directory where to save the PDB file.
        overwrite : bool, optional
            Whether or not to overwrite files in cache if they exist. If no revision
            number is provided, files will always be overwritten.
        """
        with PDBParser(ID=ID, download_dir=download_dir) as parser:
            return cls(
                unitcell=symmetry_expansion(
                    parser.residues(), parser.symmetry_operators()
                ),
                lattice_vectors=parser.lattice_vectors(),
                source=parser.filename,
                **kwargs,
            )

    @classmethod
    def from_pwscf(cls, path, **kwargs):
        """
        Returns a Crystal object created from an output file of PWSCF.
        Keyword arguments are passed to the class constructor.

        Parameters
        ----------
        path : path-like
            File path
        """
        with PWSCFParser(path) as parser:
            return cls(
                unitcell=parser.atoms(),
                lattice_vectors=parser.lattice_vectors(),
                source=parser.filename,
                **kwargs,
            )

    @classmethod
    def from_ase(cls, atoms, **kwargs):
        """
        Returns a Crystal object created from an ASE Atoms object.
        Keyword arguments are passed to the class constructor.

        Parameters
        ----------
        atoms : ase.Atoms
            Atoms group.
        """
        lattice_vectors = atoms.get_cell()

        return cls(
            unitcell=[Atom.from_ase(atm) for atm in atoms],
            lattice_vectors=lattice_vectors,
            **kwargs,
        )

    @classmethod
    def from_poscar(cls, path, **kwargs):
        """
        Returns a Crystal object created from a VASP's POSCAR file.
        Keyword arguments are passed to the class constructor.

        .. versionadded:: 1.4.0

        Parameters
        ----------
        path : path-like
            File path
        """
        with POSCARParser(path) as parser:
            return cls(
                unitcell=parser.atoms(),
                lattice_vectors=parser.lattice_vectors(),
                source=parser.filename,
                **kwargs,
            )

    def _spglib_cell(self):
        """Returns an array in spglib's cell format."""
        # To get symmetry information, we only give spglib the unit cell atoms
        # This way, all spglib-related methods (like symmetry()) will act on the unit cell only.
        # This distinction is important for Crystal subclasses, like Supercell.
        unitcell = np.stack([np.asarray(atm) for atm in self.unitcell])
        return np.array(self.lattice_vectors), unitcell[:, 1:], unitcell[:, 0]

    @classmethod
    def _from_spglib_cell(cls, lattice_vectors, scaled_positions, numbers, **kwargs):
        """
        Build a Crystal object from the return value of many SPGLIB routines.

        Parameters
        ----------
        lattice_vectors : ndarray, shape (3,3)
            Lattice vectors
        scaled_positions : ndarray, shape (N, 3)
            Fractional atomic positions, row-wise.
        numbers : ndarray, shape (N,)
            Atomic numbers, associated with each row of `scaled_positions`
        """
        atoms = [
            Atom(int(Z), coords=coords) for Z, coords in zip(numbers, scaled_positions)
        ]
        # Preserve whatever subclass this object already is
        # This is important because some properties can be extracted from
        # source files (e.g. PWSCF output files)
        return cls(unitcell=atoms, lattice_vectors=lattice_vectors, **kwargs)

    def primitive(self, symprec=1e-2):
        """
        Returns a Crystal object in the primitive unit cell.

        Parameters
        ----------
        symprec : float, optional
            Symmetry-search distance tolerance in Cartesian coordinates [Angstroms].

        Returns
        -------
        primitive : Crystal
            Crystal with primitive cell. Even if the crystal already has a primitive
            cell, a new crystal is returned.

        Raises
        ------
        RuntimeError : If primitive cell could not be found.

        Notes
        -----
        Optional atomic properties (e.g magnetic moment) might be lost in the reduction.
        """
        search = standardize_cell(
            self._spglib_cell(), to_primitive=True, no_idealize=True, symprec=symprec
        )
        if search is None:
            raise RuntimeError("Primitive cell could not be found.")

        return self._from_spglib_cell(*search, source=self.source)

    def ideal(self, symprec=1e-2):
        """
        Returns a Crystal object with an idealized unit cell.

        Parameters
        ----------
        symprec : float, optional
            Symmetry-search distance tolerance in Cartesian coordinates [Angstroms].

        Returns
        -------
        ideal : Crystal
            Crystal with idealized cell.

        Raises
        ------
        RuntimeError : If an ideal cell could not be found.

        Notes
        -----
        Optional atomic properties (e.g magnetic moment) might be lost in the symmetrization.
        """
        search = standardize_cell(
            self._spglib_cell(), to_primitive=True, no_idealize=False, symprec=symprec
        )
        if search is None:
            raise RuntimeError("Ideal cell could not be found.")

        return self._from_spglib_cell(*search, source=self.source)

    def supercell(self, n1, n2, n3):
        """
        Create a supercell from this crystal, i.e. an atomic structure where the crystal unit cell
        is duplicated along lattice vectors.

        Parameters
        ----------
        n1, n2, n3 : int
            Repeat along the `a1`, `a2`, and `a3` lattice vectors. For example,
            ``(1, 1, 1)`` represents the trivial supercell.

        Returns
        -------
        cell : AtomicStructure
            Iterable of `crystals.Atom` objects following the supercell dimensions.
        """
        multicell = list()
        for atm in self:
            for factors in product(range(n1), range(n2), range(n3)):
                fractional_offset = np.asarray(factors)
                newatom = deepcopy(atm)
                newatom.coords_fractional += fractional_offset
                multicell.append(newatom)

        return Supercell(unitcell=multicell, lattice_vectors=self.lattice_vectors)

    def symmetry(self, symprec=1e-2, angle_tolerance=-1.0):
        """
        Returns a dictionary containing space-group information. This information
        is computed from the crystal unit cell.

        Parameters
        ----------
        symprec : float, optional
            Symmetry-search distance tolerance in Cartesian coordinates [Angstroms].
        angle_tolerance: float, optional
            Symmetry-search tolerance in degrees. If the value is negative (default),
            an internally optimized routine is used to judge symmetry.

        Returns
        -------
        info : dict
            Dictionary of space-group information. The following keys are available:

            * ``'international_symbol'``: International Tables of Crystallography
              space-group symbol (short);

            * ``'international_full'``: International Tables of
              Crystallography space-group full symbol;

            * ``'hall_symbol'`` : Hall symbol;

            * ``'hm_symbol'`` : Hermann-Mauguin symbol;

            *``'centering'``: Centering-type ("P", "F", etc.);

            * ``'pointgroup'`` : International Tables of
              Crystallography point-group;

            * ``'international_number'`` : International Tables of
              Crystallography space-group number (between 1 and 230);

            * ``'hall_number'`` : Hall number (between 1 and 531).

        Raises
        ------
        RuntimeError : if symmetry-determination has not succeeded.

        Notes
        -----
        Note that crystals generated from the Protein Data Bank are often incomplete;
        in such cases the space-group information will be incorrect.
        """
        dataset = get_symmetry_dataset(
            cell=self._spglib_cell(), symprec=symprec, angle_tolerance=angle_tolerance
        )

        if dataset is None:
            raise RuntimeError("[SPGLIB] Symmetry-determination has not found a match.")

        spg_type = get_spacegroup_type(dataset["hall_number"])
        hm_symbol = Hall2HM[dataset["hall"]]

        # We do not distinguish between base-centered "A", "B", and "C"
        # "A" and "B" are translated to "C"
        centering = CenteringType(
            hm_symbol[0] if hm_symbol[0] not in {"A", "B"} else "C"
        )
        info = {
            "international_symbol": dataset["international"],
            "hall_symbol": dataset["hall"],
            "hm_symbol": hm_symbol,
            "centering": centering,
            "international_number": dataset["number"],
            "hall_number": dataset["hall_number"],
            "international_full": spg_type["international_full"],
            "pointgroup": spg_type["pointgroup_international"],
        }

        err_msg = get_error_message()
        if err_msg != "no error":
            raise RuntimeError(
                "[SPGLIB] Symmetry-determination has returned the following error: {err_msg}"
            )

        return info

    def symmetry_operations(self, symprec=1e-2):
        """
        Get the symmetry operations that the crystal unit cell respects. These symmetry operations
        are expressed in fractional coordinates.

        Parameters
        ----------
        symprec : float, optional
            Symmetry-search distance tolerance in Cartesian coordinates [Angstroms].

        Returns
        -------
        sym_ops : iterable of array_like, shapes (4,4)
            Iterable of affine transforms, where ``m[:3,:3]`` is the rotation part,
            while ``m[:3,-1]`` is the translation.

        Raises
        ------
        RuntimeError : if symmetry-determination has not succeeded.

        See also
        --------
        Crystal.reciprocal_symmetry_operations : symmetry operations in reciprocal basis
        """
        dataset = get_symmetry(cell=self._spglib_cell(), symprec=symprec)

        def _to_affine(r, t):
            """Convert rotation and translation into single 4x4 affine transformation"""
            m = np.eye(4)
            m[:3, :3] = r
            m[:3, -1] = t
            return m

        return [
            _to_affine(r, t)
            for r, t in zip(dataset["rotations"], dataset["translations"])
        ]

    def reciprocal_symmetry_operations(self, symprec=1e-2):
        """
        Get the symmetry operations that the reciprocal unit cell respects. These symmetry operations
        are expressed in reciprocal fractional coordinates.

        Parameters
        ----------
        symprec : float, optional
            Symmetry-search distance tolerance in Cartesian coordinates [Angstroms].

        Returns
        -------
        sym_ops : iterable of array_like, shapes (4,4)
            Iterable of affine transforms, where ``m[:3,:3]`` is the rotation part,
            while ``m[:3,-1]`` is the translation.

        Raises
        ------
        RuntimeError : if symmetry-determination has not succeeded.

        See also
        --------
        Crystal.symmetry_operations : symmetry operations in lattice basis
        """
        transformations = self.symmetry_operations(symprec=symprec)

        # Change of basis matrices allow to express
        # transformations in other bases
        to_reciprocal = change_of_basis(
            np.array(self.lattice_vectors), np.array(self.reciprocal_vectors)
        )
        from_reciprocal = np.linalg.inv(to_reciprocal)

        cast = lambda m: to_reciprocal @ m @ from_reciprocal

        # Pack and unpack rotation and translation from/to affine transform
        # TODO: this is ugly
        def pack(r, t):
            m = np.eye(4)
            m[:3, :3] = r
            m[:3, -1] = t
            return m

        ops = list()
        for m in transformations:
            r, t = m[:3, :3], m[:3, -1]
            ops.append(pack(cast(r), cast(t)))
        return ops

    @property
    def international_symbol(self):
        """International Tables of Crystallography space-group short symbol."""
        return self.symmetry()["international_symbol"]

    @property
    def international_full(self):
        """International Tables of Crystallography space-group full symbol."""
        return self.symmetry()["international_full"]

    @property
    def hall_symbol(self):
        """Hall symbol."""
        return self.symmetry()["hall_symbol"]

    @property
    def hm_symbol(self):
        """Hermann-Mauguin symbol."""
        return self.symmetry()["hm_symbol"]

    @property
    def pointgroup(self):
        """International Tables of Crystallography point-group."""
        return self.symmetry()["pointgroup"]

    @property
    def international_number(self):
        """International Tables of Crystallography space-group number (between 1 and 230)."""
        return self.symmetry()["international_number"]

    @property
    def hall_number(self):
        """Hall number (between 1 and 531)."""
        return self.symmetry()["hall_number"]

    @property
    def centering(self):
        """Centering type of this crystals."""
        return self.symmetry()["centering"]

    def indexed_by(self, lattice):
        """
        Return a crystal structure, indexed by another lattice/crystal structure.

        Parameters
        ----------
        lattice : :class:`Crystal`, :class:`Lattice`, or array_like, shape (3,3)
            Lattice or Crystal by which to index the structure. ``lattice`` can also be a 3x3 array,
            where every row is a lattice vector.

        Returns
        -------
        crystal : :class:`Crystal`
            New crystal, indexed by ``lattice``.
        """
        # TODO: example
        if not isinstance(lattice, (Crystal, Lattice)):
            lattice = Lattice(lattice)

        old_basis = self.lattice_vectors
        new_basis = lattice.lattice_vectors

        cob = change_of_basis(old_basis, new_basis)
        return self.__class__(
            unitcell=(atm.transform(cob) for atm in self),
            lattice_vectors=new_basis,
            source=self.source,
        )

    def __str__(self):
        """String representation of this instance. Atoms may be omitted."""
        return self._to_string(natoms=10)

    def __repr__(self):
        """Verbose string representation of this instance."""
        return self._to_string(natoms=len(self))

    def _to_string(self, natoms):
        """Generate a string representation of this Crystal. Only include
        a maximum of `natoms` if provided."""

        # Note : Crystal subclasses need not override this method
        # since the class name is dynamically determined
        rep = f"< {self.__class__.__name__} object with following unit cell:"
        atoms = islice(sorted(self), natoms)

        # Note that repr(Atom(...)) includes these '< ... >'
        # We remove those for cleaner string representation
        rep += "".join(
            "\n    " + repr(atm).replace("<", "").replace(">", "").strip()
            for atm in atoms
        )

        num_omitted_atms = len(self) - natoms
        if num_omitted_atms > 0:
            rep += f"\n      ... omitting {num_omitted_atms:d} atoms ..."
            rep += "\n      ... use repr() to show the full cell ... "

        # Lattice parameters are split between lengths and angles
        a, b, c, alpha, beta, gamma = self.lattice_parameters
        rep += "\nLattice parameters:"
        rep += f"\n    a={a:.3f}Å, b={b:.3f}Å, c={c:.3f}Å"
        rep += f"\n    α={alpha:.3f}°, β={beta:.3f}°, γ={gamma:.3f}°"

        # Show stochiometric information
        # It is important that the chemical symbols are presented in order
        # so that doctests are reproducible
        rep += "\nChemical composition:"
        for chem_symbol, composition in sorted(
            self.chemical_composition.items(), key=lambda t: t[0]
        ):
            rep += f"\n    {chem_symbol.rjust(LONGEST_CHEMICAL_SYMBOL)}: {100 * composition:.3f}%"
        rep += " >"
        return rep

    def to_cif(self, filename):
        """
        Convert this :class:`Crystal` instance to a CIF file.

        Note that some information may be lost in the translation. However, we guarantee that
        reading a structure from a file, and then writing back to the same format is idempotent.

        Parameters
        ----------
        filename : path-like
            Path to a file. If the file already exists, it will be overwritten.

        See also
        --------
        Crystal.to_xyz : write a structure to a `.xyz` file.
        Crystal.to_ase : convert a structure into an ``ase.Atoms`` object.
        """
        write_cif(self, filename)

    def to_xyz(self, filename):
        """
        Convert this :class:`Crystal` instance to a XYZ file.

        Note that some information may be lost in the translation. However, we guarantee that
        reading a structure from a file, and then writing back to the same format is idempotent.

        Parameters
        ----------
        filename : path-like
            Path to a file. If the file already exists, it will be overwritten.

        See also
        --------
        Crystal.to_cif : write a structure to a `.cif` file.
        Crystal.to_ase : convert a structure into an ``ase.Atoms`` object.
        """
        write_xyz(self, filename)

    def to_ase(self, **kwargs):
        """
        Convert a into an :class:`ase.Atoms` object. Keyword arguments are passed
        to :class:`ase.Atoms` constructor.

        Note that some information may be lost in the translation. However, we guarantee that
        reading a structure from a file, and then writing back to the same format is idempotent.

        Returns
        -------
        atoms : ase.Atoms
            Group of atoms ready for ASE's routines.

        Raises
        ------
        ImportError : If ASE is not installed

        See also
        --------
        Crystal.to_cif : write a structure to a `.cif` file.
        Crystal.to_xyz : write a structure to a `.xyz` file.
        """
        return ase_atoms(self)

    def to_poscar(self, filename, **kwargs):
        """
        Convert this :class:`Crystal` instance to a POSCAR file.
        Keyword arguments are passed to :meth:`writers.write_poscar`.

        Note that some information may be lost in the translation. However, we guarantee that
        reading a structure from a file, and then writing back to the same format is idempotent.

        Parameters
        ----------
        filename : path-like
            Path to a file. If the file already exists, it will be overwritten.
        kwargs:
            Keyword arguments are passed to :meth:`writers.write_poscar`.
        """
        write_poscar(self, filename, **kwargs)


class Supercell(Crystal):
    """
    The :class:`Supercell` class is a set-like container that represents a
    supercell of crystalline structures.

    It is strongly recommended that you do not instantiate a :class:`Supercell` by hand, but rather
    create a :class:`Crystal` object and use the :meth:`Crystal.supercell` method.

    To iterate over all atoms in the supercell, use this object as an iterable.
    """

    pass


@unique
class CenteringType(Enum):
    """
    Enumeration of possible centering types. Together with the lattice system,
    these centering types defined all 14 Bravais lattices in 3D.

    The possible centering types are:

    * ``'P'`` : Primitive

    * ``'I'`` : Body-centered

    * ``'F'`` : Face-centered

    * ``'C'`` : Base-centered

    * ``'R'`` : Rhombohedral in hexagonal setting.
    """

    primitive = "P"
    base_centered = "C"
    body_centered = "I"
    face_centered = "F"
    rhombohedral = "R"


def symmetry_expansion(atoms, symmetry_operators):
    """
    Generate a set of unique atoms and structures from an asymmetric cell
    and symmetry operators.

    Parameters
    ----------
    atoms : iterable of :class:`Atom` and/or :class:`AtomicStructure`
        Asymmetric unit cell atoms. It is assumed that the atomic
        coordinates are in fractional form. Transformations work
        the same way for :class:`Atom` objects and :class:`AtomicStructure`
        objects: a copy is made and moved to the symmetric location.
    symmetry_operators : iterable of array_like
        Symmetry operators that generate the full unit cell.

    Yields
    ------
    it : :class:`Atom` and/or :class:`AtomicStructure`
        Appropriately-transformed object. Original objects are left untouched.

    See Also
    --------
    symmetry_reduction : Determine the asymmetric cell that can generate a unit cell.
    """
    # TODO: provide ability to reduce to primitive, niggli_reduce, etc.
    #       using spglib?

    is_identity = lambda op: np.allclose(op, np.eye(4))
    symmetry_operators = [
        m for m in map(affine_map, symmetry_operators) if not is_identity(m)
    ]

    unique_atoms = set([])
    for atm in filter(is_atom, atoms):
        # At least one atom is kept
        unique_atoms.add(atm)
        for sym_op in symmetry_operators:
            new = atm.transform(sym_op)
            new.coords_fractional[:] = np.mod(new.coords_fractional, 1)
            unique_atoms.add(new)

    def _normalize_fractional_coordinates(struct):
        """Ensure that all atoms in this structure have fractional coordinates that are valid, i.e.
        in [0, 1). Atoms are changed in-place."""
        for atm in struct.atoms:
            atm.coords_fractional[:] = np.mod(atm.coords_fractional, 1)

        for substruct in struct.substructures:
            _normalize_fractional_coordinates(substruct)

    unique_structures = set([])
    for structure in filter(is_structure, atoms):
        unique_structures.add(structure)
        for sym_op in symmetry_operators:
            new = structure.transform(sym_op)
            unique_structures.add(_normalize_fractional_coordinates(new))

    yield from unique_atoms
    yield from unique_structures


def symmetry_reduction(unitcell, symmetry_operators):
    """
    Determine the asymmetric cell that generates `unitcell` when combined with
    symmetry operations. Effectively, this function is the reciprocal operation
    to `symmetry_expansion`.

    Note that the resulting asymmetric cell is not unique.

    .. versionadded:: 1.2.0

    Parameters
    ----------
    unitcell : iterable of :class:`Atom` and/or :class:`AtomicStructure`
        Unit cell. It is assumed that the atomic
        coordinates are in fractional form.
    symmetry_operators : iterable of array_like
        Symmetry operators that generate the full unit cell.

    Returns
    -------
    asym_cell: set of :class:`Atom` and/or :class:`AtomicStructure`
        Asymmetric cell.

    See Also
    --------
    symmetry_expansion : Expand the asymmetric cell into a unit cell.
    """
    # NOTE: This function does the dumb thing of trying to recreate the input
    #       by searching through all the possibilities.
    #       By nature, it is very slow; scaling like n^2 (where n is the unit cell length)
    #       However, other approaches that I tried were not reliable. For example,
    #       using inverse symmetry operators would result in numerical instabilities.
    #       If you are reading this and have a better idea, please submit a pull request :).

    unitcell = set(unitcell)

    # The powerset of the unit cell is the set of all possible sets
    # made from elements in `unitcell`
    powerset = chain.from_iterable(
        combinations(unitcell, r) for r in range(1, len(unitcell) + 1)
    )

    # We sort by length so that the resulting asymmetric cell is the smallest
    # possible one.
    for asym_cell in sorted(powerset, key=len):
        reconstituted = set(symmetry_expansion(asym_cell, symmetry_operators))
        if reconstituted == unitcell:
            return set(asym_cell)

    # Base case
    return unitcell
