# -*- coding: utf-8 -*-
"""
Conversion between ``crystals`` data structures and other modules.
"""
import numpy as np

from .atom import Atom
from .lattice import Lattice
from .crystal import Crystal

try:
    import ase
except ImportError:
    WITH_ASE = False
else:
    WITH_ASE = True


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
