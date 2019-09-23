# -*- coding: utf-8 -*-
"""
crystals command-line utilities
"""

import argparse
import sys
from pathlib import Path
from . import Crystal, __version__

constructors = {
    "database": Crystal.from_database,
    "cif": Crystal.from_cif,
    "pwscf": Crystal.from_pwscf,
    "cod": Crystal.from_cod,
    "pdb": Crystal.from_pdb,
}

parser = argparse.ArgumentParser(
    prog="crystals", description=f"crystals {__version__} command-line utilities"
)

subparsers = parser.add_subparsers(title="command", dest="command")

info_parser = subparsers.add_parser(
    "info", help="Display the crystallographic information related to a crystal file or database input"
)
info_parser.add_argument("input", nargs=1, help="Path to a file, or database input.")

info_parser.add_argument(
    "--type",
    help="Type of input. By default, `crystals` will try to guess.",
    type=str,
    choices=[f"{itype}" for itype in constructors],
    required=False,
)


def guess_constructor(i):
    """ 
    Try to guess the appropriate Crystal constructor
    
    Returns
    -------
    cryst : Crystal
    
    Raises
    ------
    RuntimeError: if the input type could not be determined.
    """
    # Some easy heuristics first
    if str(i).endswith(".cif"):
        return Crystal.from_cif(i)
    elif str(i).endswith(".pwscf"):
        return Crystal.from_pwscf(i)

    # We iterate over constructors that don't require network access first
    cryst = None
    for constructor in (
        Crystal.from_cif,
        Crystal.from_pwscf,
        Crystal.from_database,
        Crystal.from_cod,
        Crystal.from_pdb,
    ):
        try:
            cryst = constructor(i)
        except:
            continue
        else:
            break
    if cryst is None:
        raise RuntimeError

    return cryst


def main():
    args = parser.parse_args()

    if args.command == "info":
        if args.type is not None:
            cryst = constructors[args.type](args.input[0])
        else:
            try:
                cryst = guess_constructor(args.input[0])
            except RuntimeError:
                print(
                    "Input type could not be determined. Try to use the --type option."
                )
                sys.exit(1)

        unitcell = (
            repr(cryst).replace("< ", "").replace(" >", "")
        )  # This is from python's object representation
        sym = f"""Symmetry information:
    International symbol 
                (short) ..... {cryst.international_symbol}
                (full) ...... {cryst.international_full}
    International number .... {cryst.international_number}
    Hermann-Mauguin symbol .. {cryst.hm_symbol}
    Pointgroup .............. {cryst.pointgroup}
    Hall Number ............. {cryst.hall_number}
    Centering ............... {cryst.centering}
    """

        print(unitcell)
        print(sym)
        sys.exit(0)


if __name__ == "__main__":
    main()
