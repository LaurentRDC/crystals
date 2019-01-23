"""
Update the Crystal Information Files (CIF) built into `crystals`. These CIFs
are provided by the Crystallography Open Database (COD).

This script exists because COD heavily curates CIF entries. Therefore, updates
to the CIFs are to be expected.
"""

import csv
import os
from pathlib import Path
from multiprocessing.dummy import Pool # multiprocessing API using threads rather than processes

from crystals import CODParser

# Ultimate download dir
DOWNLOAD_DIR = Path(".") / "crystals" / "cifs"

# By defining a function to download a specific structure,
# we can use it in conjunction with a thread pool.
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
    print(f"Downloaded ID: {cod_id}".ljust(20),  f"| {name} ")

if __name__ == "__main__":

    with open(Path(".") / "builtins_cifs.csv") as built_in_cifs:
        reader = list(csv.reader(built_in_cifs, delimiter=","))
    
    Pool(None).starmap(download, reader)
