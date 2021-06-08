# -*- coding: utf-8 -*-
from .affine import affine_map
from .base import AtomicStructure


class Residue(AtomicStructure):
    """
    Fundamental building block of a protein: amino acid residue.

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

    @property
    def residues(self):
        """Residues making this structure."""
        return self.substructures

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
