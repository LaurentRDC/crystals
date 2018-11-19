# -*- coding: utf-8 -*-

from collections import Counter, OrderedDict
from copy import deepcopy as copy
from functools import reduce
from itertools import chain
from math import gcd

import numpy as np

from .affine import affine_map


class Base:
    """ 
    Base class that overrides the builtin behavior:
    >>> object() == object()
    False

    This allows for transparent multiple inheritance of subclasses.
    """

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return True
        return NotImplemented

    def __hash__(self):
        return 0


class AtomicStructure(Base):
    """
    Base class for atomic structures. These structures can be made
    out of :class:`Atom` objects, or other AtomicStructure subclasses.

    The AtomicStructure class provides an abstraction over structure with and without
    substructure. Subclasses can be iterated over as an iterable of atoms. Order of iteration
    is not guaranteed.
    
    Hierarchical containership of AtomicStructures is implemented, as well as 
    containership of :class:`Atom` instances.

    Parameters
    ----------
    atoms : iterable, optional
        Atoms not attached to a substructure.
    substructures : iterable, optional
        Atomic substructures. For example, different atom chains could form
        a secondary structure in a protein.
    """

    def __init__(self, atoms=None, substructures=None, **kwargs):
        self.atoms = frozenset(atoms or {})
        self.substructures = frozenset(substructures or {})
        super().__init__(**kwargs)

    def __iter__(self):
        """ Yields :class:`Atom` instances from the structure and substructures 
        recursively. Order is not guaranteed. """
        yield from iter(self.atoms)
        yield from chain.from_iterable(self.substructures)

    def itersorted(self, *, key=None, reverse=False):
        """ 
        Yields :class:`Atom` in sorted order. By default, atoms are sorted by element. 

        Parameters
        ----------
        key : callable or None, optional 
            Function of one argument that is used to extract a comparison key for an :class:`Atom` instance.
        reverse : bool, optional
            If True, elements are yielded in reverse order.
        
        Yields
        ------
        atm : `Atom` 
        """
        if key is None:
            key = lambda atm: atm.element
        yield from sorted(iter(self), key=key, reverse=reverse)

    @property
    def chemical_composition(self):
        """ 
        Chemical composition of this structure as an ordered dictionary. Keys are elemental symbols. 
        Elements are in descending order of prevalence.
        """
        # We can't use a Counter directly since Counter values are integer by default
        number_atoms = len(self)
        counter = Counter(atm.element for atm in self)
        # For display reasons, it makes more sense to sort elements in descending percentage
        # i.e. more prevalent elemts are inserted in the dictionary first
        sorted_by_percentage = sorted(
            counter.items(), key=lambda tup: tup[-1], reverse=True
        )
        return OrderedDict((k, v / number_atoms) for k, v in sorted_by_percentage)

    @property
    def chemical_formula(self):
        """ 
        Empirical chemical formula for this structure based on the chemical symbols. The string is returned 
        in Hill notation: symbols are alphabetically ordered except for carbon (C) and hydrogen (H), which are put first. 
        """
        symbols_count = Counter(atm.element for atm in self)

        # Special case: only one distinct element
        if len(symbols_count) == 1:
            return next(iter(symbols_count.keys()))

        # Calculate greatest common divisor for all elements
        # Note that we are guaranteed to have at least two distinct elements
        counts = list(symbols_count.values())
        common_divisor = reduce(gcd, counts[1:], counts[0])
        symbols_count_empirical = {
            k: int(v / common_divisor) for k, v in symbols_count.items()
        }

        # Build a list of elements and their counts
        # We put carbon and hydrogen first because convention
        elements = []
        for elem in ("C", "H"):
            if elem in symbols_count_empirical:
                elements.append((elem, symbols_count_empirical.pop(elem)))

        # Then, elements are in alphabetical number
        elements.extend(
            [(e, symbols_count_empirical[e]) for e in sorted(symbols_count_empirical)]
        )
        return "".join(
            f"{symbol}{count}" if count > 1 else symbol for symbol, count in elements
        )

    def transform(self, *operators):
        """
        Return a transformed AtomicStructure based on symmetry operators.

        operators : iterable of array_like
            Symmetry operators, either 3x3 or 4x4.
        
        Returns
        -------
        transformed : AtomicStructure
            New structure.
        """
        operators = tuple(map(affine_map, operators))

        transformed_atoms = set()
        for atom in self.atoms:
            for operator in operators:
                transformed_atoms.add(atom.transform(operator))

        # Substructures are recursively transformed
        transformed_substructures = set()
        for substructure in self.substructures:
            transformed_substructures.add(substructure.transform(*operators))

        # We defer construction to the current class. Therefore, subclasses of AtomicStructure
        # will transform into their own class
        return self.__class__(
            atoms=transformed_atoms, substructures=transformed_substructures
        )

    def __contains__(self, item):
        """ Check containership of :class:`Atom` instances or :class:`AtomicStructure` substructures recursively."""
        if isinstance(item, AtomicStructure):
            return item in self.substructures

        # Either the item is an orphan atom or
        # it is in one of the substructures
        # Checking containership of sets is faster than iterating
        if item in self.atoms:
            return True
        else:
            return any((item in struct) for struct in self.substructures)

    def __len__(self):
        """ Number of :class:`Atom` instances present in the structure and substructures """
        return len(self.atoms) + sum(len(struct) for struct in self.substructures)

    def __hash__(self):
        return hash((self.atoms, self.substructures, super().__hash__()))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (
                set(self.atoms) == set(other.atoms)
                and set(self.substructures) == set(other.substructures)
                # Subclasses (like Crystal via Lattice) might have extra requirements for equality
                and super().__eq__(other)
            )
        return NotImplemented

    def __repr__(self):
        """ Verbose string representation of this instance. """
        # AtomicStructure subclasses need not override this method
        # since the class name is dynamically determined
        rep = "< {clsname} object with following orphan atoms:".format(
            clsname=self.__class__.__name__
        )

        # Note that repr(Atom(...)) includes these '< ... >'
        # We remove those for cleaner string representation
        for atm in self.itersorted():
            rep += "\n    " + repr(atm).replace("<", "").replace(">", "").strip()

        if self.substructures:
            rep += "and the following substructures:"
            for struct in self.substructures:
                rep += "\n" + repr(struct)

        return rep + " >"

    def __array__(self, *args, **kwargs):
        """ Returns an array in which each row represents an :class:`Atom` instance. Atoms are ordered by atomic number """
        arr = np.empty(shape=(len(self), 4), *args, **kwargs)
        atoms = self.itersorted(key=lambda atm: atm.atomic_number)
        for row, atm in enumerate(atoms):
            arr[row, :] = np.array(atm, *args, **kwargs)
        return arr


class Residue(AtomicStructure):
    """
    Container representing a protein residue.

    Parameters
    ----------
    atoms : iterable, optional
        Atoms not attached to a substructure.
    name : str
        Residue name, e.g. "GLY"
    sequence_number : int
        Sequence number within a protein.
    """

    def __init__(self, atoms, name, sequence_number, **kwargs):
        super().__init__(atoms=atoms, **kwargs)
        self.name = name
        self.sequence_number = int(sequence_number)

    def transform(self, *operators):
        """
        Return a transformed Residue based on symmetry operators.

        operators : iterable of array_like
            Symmetry operators, either 3x3 or 4x4.
        
        Returns
        -------
        transformed : Residue
            New structure.
        """
        operators = tuple(map(affine_map, operators))

        transformed_atoms = set()
        for atom in self.atoms:
            for operator in operators:
                transformed_atoms.add(atom.transform(operator))

        # We defer construction to the current class. Therefore, subclasses of AtomicStructure
        # will transform into their own class
        return self.__class__(
            atoms=transformed_atoms,
            name=self.name,
            sequence_number=self.sequence_number,
        )


# Helper class to not repeat code betwee Helix and Sheet
class SecondaryStructure(AtomicStructure):
    """
    Abstract container representing protein secondary structures (i.e. helices, sheets, etc.)

    Parameters
    ----------
    residues : iterable, optional
        Residues making this structure.
    sequence_number : int
        Sequence number within a protein.
    """
    
    def __init__(self, residues, sequence_number, **kwargs):
        super().__init__(substructures=residues, **kwargs)
        self.sequence_number = int(sequence_number)

    def transform(self, *operators):
        """
        Return a transformed SecondaryStructure based on symmetry operators.

        operators : iterable of array_like
            Symmetry operators, either 3x3 or 4x4.
        
        Returns
        -------
        transformed : SecondaryStructure
            New structure.
        """
        operators = tuple(map(affine_map, operators))

        # Substructures are recursively transformed
        transformed_substructures = set()
        for substructure in self.substructures:
            transformed_substructures.add(substructure.transform(*operators))

        # We defer construction to the current class. Therefore, subclasses of AtomicStructure
        # will transform into their own class
        return self.__class__(
            residues=transformed_substructures, sequence_number=self.sequence_number
        )


class Helix(SecondaryStructure):
    """
    Container representing a protein helix.

    Parameters
    ----------
    residues : iterable, optional
        Residues making this helix.
    sequence_number : int
        Sequence number within a protein.
    """

    def transform(self, *operators):
        """
        Return a transformed helix based on symmetry operators.

        operators : iterable of array_like
            Symmetry operators, either 3x3 or 4x4.
        
        Returns
        -------
        transformed : Helix
            New structure.
        """
        return super().transform(*operators)


class Sheet(SecondaryStructure):
    """
    Container representing a protein sheet.

    Parameters
    ----------
    residues : iterable, optional
        Residues making this sheet.
    sequence_number : int
        Sequence number within a protein.
    """

    def transform(self, *operators):
        """
        Return a transformed sheet based on symmetry operators.

        operators : iterable of array_like
            Symmetry operators, either 3x3 or 4x4.
        
        Returns
        -------
        transformed : Sheet
            New structure.
        """
        return super().transform(*operators)
