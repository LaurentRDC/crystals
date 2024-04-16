# -*- coding: utf-8 -*-
"""
This package allows for manipulation and modelling of atomic structures in crystalline form.
"""
__author__ = "Laurent P. René de Cotret"
__email__ = "laurent.renedecotret@mail.mcgill.ca"
__license__ = "GPLv3"
__version__ = "1.7.0a1"

from .atom import (
    Atom,
    ElectronicStructure,
    Element,
    Orbital,
    distance_cartesian,
    distance_fractional,
    frac_coords,
    is_element,
    real_coords,
)
from .atom_data import (
    ELEM_TO_MAGMOM,
    ELEM_TO_MASS,
    ELEM_TO_NAME,
    ELEM_TO_NUM,
    NUM_TO_ELEM,
)
from .base import AtomicStructure
from .biological import Helix, Residue, Sheet
from .crystal import (
    CenteringType,
    Crystal,
    Supercell,
    symmetry_expansion,
    symmetry_reduction,
)
from .indexing import IndexingError, index_dirax
from .lattice import Lattice, LatticeSystem, lattice_system
from .parsers import CIFParser, CODParser, MPJParser, ParseError, PDBParser, PWSCFParser
from .writers import ase_atoms
