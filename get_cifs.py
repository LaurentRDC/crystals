"""
Update the Crystal Information Files (CIF) built into `crystals`. These CIFs
are provided by the Crystallography Open Database (COD).

This script exists because COD heavily curates CIF entries. Therefore, updates
to the CIFs are to be expected.
"""

import csv
import os
from pathlib import Path
from multiprocessing.dummy import (
    Pool,
)  # multiprocessing API using threads rather than processes

from crystals import CODParser

# Ultimate download dir
DOWNLOAD_DIR = Path(".") / "crystals" / "cifs"

BUILTINS = {
    "Er": 9008496,
    "vo2-rutile": 1537412,
    "Nb": 9008546,
    "Mg": 9008506,
    "Ca": 9008464,
    "Cs": 9008532,
    "Pu-epsilon": 9008548,
    "C": 9011577,
    "diamond": 9008564,
    "Ba": 9008528,
    "As": 9008574,
    "La": 9008525,
    "Gd": 9008498,
    "alpha-Mn": 9008589,
    "Pd": 9008478,
    "SiC": 1010995,
    "Pu-delta": 9008481,
    "I": 9008595,
    "Np": 9008585,
    "In": 2100456,
    "Am": 9008524,
    "Cu": 9008468,
    "Tb": 9010992,
    "K": 9008539,
    "gamma-Pu": 9008588,
    "Ac": 9008458,
    "Ho": 9010994,
    "Rh": 9008482,
    "Zr": 9008523,
    "Sm": 9010999,
    "V": 9008557,
    "Sb": 9008575,
    "Lu": 9010998,
    "Zn": 9008522,
    "BaTiO3_cubic": 2100862,
    "Ti": 9008517,
    "Si": 9008565,
    "Be": 9008488,
    "Rb": 9008549,
    "Na": 9008544,
    "Au": 9008463,
    "Nd": 9008526,
    "Ce-gamma": 9008465,
    "Ru": 9008513,
    "Al": 9008460,
    "Dy": 9008494,
    "Eu": 9008534,
    "Pr": 9008527,
    "Li": 9008541,
    "Cd": 9008490,
    "Y": 9008521,
    "vo2-m1": 9009089,
    "Tl": 9008518,
    "Pu-gamma": 9008588,
    "Yb": 9010997,
    "Bi": 5000215,
    "Re": 9008512,
    "Mn": 9011108,
    "Pt": 9008480,
    "Ni": 9008476,
    "P": 9008572,
    "Os": 9008510,
    "Sn-alpha": 9008568,
    "Pu-alpha": 9008587,
    "Ag": 9008459,
    "Mo": 9008543,
    "Br": 9008594,
    "Fe": 9008536,
    "Ir": 9008470,
    "Sn-beta": 9008570,
    "B": 9008561,
    "Ge": 9008567,
    "Ta": 9008552,
    "Pb": 9008477,
    "FeAs": 9007668,
    "Sr": 9008484,
    "Ga": 9008562,
    "Cr": 9008531,
    "Sc": 9008514,
    "Hf": 9008501,
    "Tm": 9010996,
    "S": 9008577,
    "Th": 9008485,
    "U": 9008584,
    "W": 9008558,
    "Co": 9008492,
    "GaAs": 9008845,
}


def download(name, cod_id):
    """
    Parameters
    ----------
    name : str
        Name of the structure to download
    cod_id : str or int
        Crystallography Open Database ID of the structure.
    """
    path = CODParser.download_cif(
        DOWNLOAD_DIR, num=int(cod_id), revision=None, overwrite=True
    )
    os.replace(path, path.parent / f"{name}.cif")
    print(f"Downloaded ID: {cod_id}".ljust(20), f"| {name} ")


if __name__ == "__main__":
    Pool(None).starmap(download, BUILTINS.items())
