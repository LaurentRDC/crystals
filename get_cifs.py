"""
Update the Crystal Information Files (CIF) built into `crystals`. These CIFs
are provided by the Crystallography Open Database (COD).

This script exists because COD heavily curates CIF entries. Therefore, updates
to the CIFs are to be expected.
"""

import csv
import os
from pathlib import Path

from crystals import CODParser

if __name__ == "__main__":

    # Ultimate download dir
    download_dir = Path(".") / "crystals" / "cifs"

    with open(Path(".") / "builtins_cifs.csv") as built_in_cifs:
        reader = csv.reader(built_in_cifs, delimiter=",")
        for (name, cod_id) in reader:
            # The file will be downloaded to the appropriate folder but it's name
            # will not be very evocative. We therefore rename to something more appropriate
            path = CODParser.download_cif(
                download_dir, num=int(cod_id), revision=None, overwrite=True
            )
            os.replace(path, path.parent / f"{name}.cif")

            print(f"Downloaded ID: {cod_id} [ {name} ] ")
