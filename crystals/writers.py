# -*- coding: utf-8 -*-
"""
Conversion between ``crystals`` data structures and other modules.

These functions are not expected to be used on their own; see the associated
 `Crystal` methods instead, like `Crystal.to_cif`.
"""
from abc import abstractmethod
from contextlib import AbstractContextManager, redirect_stdout
from io import StringIO

import numpy as np
from CifFile import CifBlock, CifFile

from . import __version__

try:
    import ase
except ImportError:
    WITH_ASE = False
else:
    WITH_ASE = True

CIF_HEADER = """
#\#CIF_2.0
#
# File generated by `crystals` {__version__}, documented at https://crystals.rtfd.io
# Reference:  L. P. René de Cotret et al, An open-source software ecosystem for the interactive exploration
#             of ultrafast electron scattering data, Advanced Structural and Chemical Imaging 4:11 (2018) DOI: 10.1186/s40679-018-0060-y.
#
# For more information on this type of file, please visit https://www.iucr.org/resources/cif
"""

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
            x, y, z = atom.coords_cartesian
            row = f"  {atom.element:<2}       {x:10.5f}       {y:10.5f}       {z:10.5f}"
            file.write(row + "\n")


def write_cif(crystal, fname):
    """
    Generate an atomic coordinates .cif file from a crystal structure.

    Parameters
    ----------
    crystal : crystals.Crystal
        Crystal to be converted.
    fname : path-like
        The CIF file will be written to this file. If the file already exists,
        it will be overwritten.
    comment : str or None, optional
        Comment to include at the second line of ``fname``.
    """
    cf = CifFile(strict=False)
    a, b, c, alpha, beta, gamma = crystal.lattice_parameters
    lattice_items = {
        "_cell_length_a": a,
        "_cell_length_b": b,
        "_cell_length_c": c,
        "_cell_angle_alpha": alpha,
        "_cell_angle_beta": beta,
        "_cell_angle_gamma": gamma,
    }

    sym = crystal.symmetry()
    symmetry_items = {
        "_symmetry_Int_Tables_number": sym["international_number"],
        "_symmetry_space_group_name_Hall": sym["hall_symbol"],
    }

    block = CifBlock()
    for key, val in lattice_items.items():
        block[key] = val

    for key, val in symmetry_items.items():
        block[key] = val

    # Note that we are using all atoms in the unit-cell,
    # and not the asymmetric unit cell + symmetry operators
    # This is valid CIF! And it is much simpler to implement
    # TODO: how to determine asymmetric cell + symmetry operations?
    atoms = list(crystal.primitive().unitcell)
    symbols = [atm.symbol for atm in atoms]
    xf = [atm.coords_fractional[0] for atm in atoms]
    yf = [atm.coords_fractional[1] for atm in atoms]
    zf = [atm.coords_fractional[2] for atm in atoms]

    block.CreateLoop(
        datanames=[
            "_atom_site_type_symbol",
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
        ],
        length_check=False,
    )
    block["_atom_site_type_symbol"] = symbols
    block["_atom_site_fract_x"] = xf
    block["_atom_site_fract_y"] = yf
    block["_atom_site_fract_z"] = zf

    # Name of the block cannot be empty!
    block_name = crystal.chemical_formula.replace(" ", "_")
    cf[block_name] = block

    # Converting to string writes to stdout for some reason
    with redirect_stdout(StringIO()):
        lines = str(cf).splitlines()

    with open(fname, "w", encoding="utf-8") as f:
        f.write(CIF_HEADER)
        f.write("\n".join(lines[13::]))  # Skip the fixed header


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
        **kwargs,
    )
