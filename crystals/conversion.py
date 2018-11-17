# -*- coding: utf-8 -*-
"""
Conversion between ``crystals`` data structures and other modules.
"""
import numpy as np

from .atom import Atom
from .crystal import Crystal
from .lattice import Lattice

try:
    import ase
except ImportError:
    WITH_ASE = False
else:
    WITH_ASE = True

# TODO: test against known XYZ file
def write_xyz(crystal, fname, comment=None):
    """
    Generate an atomic coordinates .xyz file from a crystal structure.

    Parameters
    ----------
    crystal : crystals.Crystal
        Crystal to be converted.
    fname : path-like
        The XYZ file will be written to this file. If the file already exists,
        it will be overwritten.
    comment : str or None, optional
        Comment to include at the second line of ``fname``.
    """
    # Format is specified here:
    #   http://openbabel.org/wiki/XYZ_%28format%29
    comment = comment or ""
    atom_format_str = "  {:<2}       {:10.5f}       {:10.5f}       {:10.5f}"

    with open(fname, "wt", encoding="ascii") as file:
        # First two lines are:
        #   1. Number of atoms described in the file
        #   2. Optional comment
        file.write(str(len(crystal)) + "\n")
        file.write(comment + "\n")

        # Write atomic data row-by-row
        # For easier human readability, atoms are sorted
        # by element
        for atom in crystal.itersorted():
            row = atom_format_str.format(atom.element, *atom.coords_cartesian)
            file.write(row + "\n")


def ase_atoms(crystal, **kwargs):
    """ 
    Convert a ``crystals.Crystal`` object into an ``ase.Atoms`` object. 
    Keyword arguments are passed to ``ase.Atoms`` constructor.
    
    Parameters
    ----------
    crystal : crystals.Crystal
        Crystal to be converted.
    
    Returns
    -------
    atoms : ase.Atoms
        Group of atoms ready for ASE's routines.
    
    Raises
    ------
    ImportError : If ASE is not installed
    """
    if not WITH_ASE:
        raise ImportError("ASE is not installed/importable.")

    return ase.Atoms(
        symbols=[
            ase.Atom(
                symbol=atom.element,
                position=atom.coords_cartesian,
                magmom=atom.magmom,
                mass=atom.mass,
            )
            for atom in crystal
        ],
        cell=np.array(crystal.lattice_vectors),
        **kwargs
    )
