# -*- coding: utf-8 -*-
"""
This package allows for manipulation and modelling of atomic structures in crystalline form.
"""
__author__ = 'Laurent P. René de Cotret'
__email__ = 'laurent.renedecotret@mail.mcgill.ca'
__license__ = 'MIT'
__version__ = '0.0.1'
    
from .atom import Atom, real_coords, frac_coords
from .atom_data import ELEM_TO_MAGMOM, ELEM_TO_MASS, ELEM_TO_NAME, ELEM_TO_NUM, NUM_TO_ELEM
from .base import AtomicStructure, Base
from .lattice import Lattice, LatticeSystem, lattice_system
from .parsers import CIFParser, PDBParser, CODParser, ParseError
from .crystal import Crystal, symmetry_expansion
