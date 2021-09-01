# -*- coding: utf-8 -*-
import os
import socket
import tempfile
import pytest
from collections import Counter, namedtuple
from contextlib import suppress
from functools import lru_cache
from itertools import chain
from pathlib import Path
from tempfile import gettempdir
from warnings import catch_warnings, filterwarnings

import numpy as np
from crystals import CIFParser, Crystal, MPJParser, PDBParser, frac_coords, is_element
from crystals.affine import transform
from crystals.parsers import STRUCTURE_CACHE, PWSCFParser
from crystals.spg_data import Hall2Number
from spglib import get_symmetry_dataset

try:
    import Bio.PDB as biopdb
except ImportError:
    WITH_BIOPYTHON = False
else:
    WITH_BIOPYTHON = True

# API key to test Materials Project-related things
MPJ_API_KEY = os.environ.get("MATERIALS_PROJECT_API_KEY", None)

CIF_FILES = [Crystal.from_database(name).source for name in Crystal.builtins]
PWSCF_FILES = [
    Path(__file__).parent / "data" / stem
    for stem in ("pwscf_tise2.out", "pwscf_snse.out", "pwscf_graphite.out")
]

# Used to compare crystals.Atom instances and Bio.PDB.Atom instances
GenericAtom = namedtuple("GenericAtom", ["element", "coords"])

filterwarnings("ignore", category=UserWarning)


@lru_cache(maxsize=1)
def web_avail():
    """Returns whether or not an internet connection is available"""
    with suppress(OSError), socket.create_connection(("www.google.com", 80)):
        return True
    return False


@pytest.mark.skipif(not web_avail(), reason="Internet connection is required.")
def test_pdb_parser_fractional_atoms():
    """Test the PDBParser returns fractional atomic coordinates."""
    with tempfile.TemporaryDirectory() as temp_dir:
        with PDBParser("1fbb", download_dir=temp_dir) as parser:
            for atm in parser.atoms():
                pass  # TODO: find a test for this.


@pytest.mark.skipif(not web_avail(), reason="Internet connection is required.")
def test_pdb_parser_symmetry_operators():
    """Test that the non-translation part of the symmetry_operators is an invertible
    matrix of determinant 1 | -1"""
    with tempfile.TemporaryDirectory() as temp_dir:
        with PDBParser("1fbb", download_dir=temp_dir) as parser:
            for sym_op in parser.symmetry_operators():
                t = sym_op[:3, :3]
                assert round(abs(abs(np.linalg.det(t)) - 1), 5) == 0


@pytest.mark.skipif(not web_avail(), reason="Internet connection is required.")
def test_pdb_parser_residues():
    """Test the parsing of residues for 1fbb"""
    with tempfile.TemporaryDirectory() as temp_dir:
        with PDBParser("1fbb", download_dir=temp_dir) as parser:
            residues = list(parser.residues())
            atoms = list(parser.atoms())

    # Flatten residues into a list of atoms
    residue_atoms = list(chain.from_iterable(residues))
    assert len(residue_atoms) == len(atoms)


@pytest.mark.skipif(not web_avail(), reason="Internet connection is required.")
def test_pdb_parser_default_download_dir():
    """Test that the file is saved in the correct temporary directory by default"""
    filename = PDBParser.download_pdb_file("1fbb")

    assert filename.exists()
    assert filename.parent == STRUCTURE_CACHE


@pytest.mark.skipif(not WITH_BIOPYTHON, reason="Biopython is not installed/importable.")
@pytest.mark.skipif(not web_avail(), reason="Internet connection is required.")
@pytest.mark.parametrize("protein_id", ("1fbb", "1fat", "1gzx"))
def test_chemical_composition(protein_id):
    """Test crystals.PDBParser returns the same chemical composition as BIO.PDB.PDBParser implementation,
    i.e. the same elements in the right proportions."""
    pdb_list = biopdb.PDBList(verbose=False, obsolete_pdb=gettempdir())
    biopdb_parser = biopdb.PDBParser()

    with catch_warnings():
        filterwarnings("ignore", category=biopdb.PDBExceptions.PDBConstructionWarning)
        with tempfile.TemporaryDirectory() as temp_dir:
            with PDBParser(protein_id, download_dir=temp_dir) as parser:
                fname = pdb_list.retrieve_pdb_file(
                    protein_id, pdir=temp_dir, file_format="pdb"
                )

                # Note: Bio.PDB atoms store element as uppercase strings. Thus, they must be changed to titlecase
                crystals_chemical_composition = Counter(
                    [atm.element for atm in parser.atoms()]
                )
                biopdb_chemical_composition = Counter(
                    [
                        atm.element.title()
                        for atm in biopdb_parser.get_structure(
                            protein_id, fname
                        ).get_atoms()
                    ]
                )

                assert biopdb_chemical_composition == crystals_chemical_composition


@pytest.mark.parametrize("file", CIF_FILES)
def test_cif_compatibility(file):
    """Test the CIFParser on all CIF files stored herein to check build errors"""
    Crystal.from_cif(file)


def test_cif_uncertainties_issue7():
    """Test that the CIF parser handles atomic positions with uncertainties (issue #7)."""
    path = Path(__file__).parent / "data" / "issue7.cif"
    Crystal.from_cif(path)


@pytest.mark.parametrize("file", CIF_FILES)
def test_cif_fractional_atoms(file):
    """Test the CIFParser returns fractional atomic coordinates."""
    with CIFParser(file) as p:
        for atm in p.atoms():
            assert atm.coords_fractional.max() <= 1
            assert atm.coords_fractional.min() >= 0


@pytest.mark.parametrize("file", CIF_FILES)
def test_cif_symmetry_operators(file):
    """Test that the non-translation part of the symmetry_operators is an invertible
    matrix of determinant 1 | -1"""
    with CIFParser(file) as p:
        for sym_op in p.symmetry_operators():
            t = sym_op[:3, :3]
            assert round(abs(abs(np.linalg.det(t)) - 1), 7) == 0


@pytest.mark.parametrize("file", CIF_FILES)
def test_cif_international_number(file):
    """Test that the international space group number  found by
    CIFParser is the same as spglib's"""
    with CIFParser(file) as p:
        from_parser = Hall2Number[p.hall_symbol()]

        crystal = Crystal.from_cif(file)
        from_spglib = crystal.international_number
        assert from_parser == from_spglib


def test_cif_silicon():
    """Test CIFParser on Si.cif (diamond structure)"""
    Si_path = Path(__file__).parent.parent / "cifs" / "Si.cif"
    si = Crystal.from_cif(Si_path)

    assert len(si) == 8


def test_cif_issue_5():
    """Test that the parsing of MgSiO3 is working. See Issue 5 for details on why this was a problem."""
    c = Crystal.from_cif(Path(__file__).parent / "data" / "issue5_MgSiO3.cif")
    assert len(c) == 20


def test_cif_vo2():
    """Test CIFParser on vo2.cif (monoclinic M1)"""
    VO2_path = Path(__file__).parent.parent / "cifs" / "vo2-m1.cif"
    vo2 = Crystal.from_cif(VO2_path)

    assert len(vo2) == 12
    assert np.allclose(vo2.lattice_parameters, (5.743, 4.517, 5.375, 90.0, 122.6, 90.0))
    assert round(abs(vo2.volume - 117.466_153_0), 7) == 0  # from cif2cell


def test_cif_site_occupancy():
    """Test that atom site occupancy is correctly parsed from CIF files."""
    path = Path(__file__).parent / "data" / "SiC_partial_site_occ.cif"
    with CIFParser(path) as parser:
        atoms = list(parser.atoms())

    for atm in filter(is_element("Si"), atoms):
        assert atm.occupancy == 0.75
    for atm in filter(is_element("C"), atoms):
        assert atm.occupancy == 0.85


@pytest.mark.xfail(reason="Bad SSL certificates")
@pytest.mark.skipif(
    MPJ_API_KEY is None,
    reason="Materials Project API key not defined in the environment.",
)
@pytest.mark.skipif(not web_avail(), reason="Internet connection is required.")
def test_mpj_example():
    """Test that the API example given on the Materials Project website is working as expected."""
    with tempfile.TemporaryDirectory() as temp_dir:
        cryst = Crystal.from_mp(
            api_key=MPJ_API_KEY,
            query="Fe2O3",
            download_dir=temp_dir,
            overwrite=True,
        )
        assert isinstance(cryst, Crystal)


@pytest.mark.xfail(reason="Bad SSL certificates")
@pytest.mark.skipif(
    MPJ_API_KEY is None,
    reason="Materials Project API key not defined in the environment.",
)
@pytest.mark.skipif(not web_avail(), reason="Internet connection is required.")
def test_material_id():
    """Test that that material ID for Fe2O3 is as expected."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # We use the Zr5Te6 formula because it has a single Materials ID
        with MPJParser(
            api_key=MPJ_API_KEY,
            query="Zr5Te6",
            download_dir=temp_dir,
            overwrite=True,
        ) as parser:
            assert parser.material_id in "mp-29957"


@pytest.mark.parametrize("file,alat", zip(PWSCF_FILES, [6.6764, 21.4862, 4.6563]))
def test_pwscf_alat(file, alat):
    """Test the parsing of the lattice parameter (alat)"""
    parser = PWSCFParser(file)

    assert parser.alat == alat

    assert round(abs(parser.alat - parser.celldm[1]), 4) == 0


@pytest.mark.parametrize("file,natoms", zip(PWSCF_FILES, [3, 8, 4]))
def test_pwscf_natoms(file, natoms):
    """Test the parsing of the number of unit cell atoms"""
    parser = PWSCFParser(file)
    assert parser.natoms == natoms


def test_pwscf_lattice_vectors_alat():
    """Test the parsing of lattice vectors in alat units"""
    parser_tise2 = PWSCFParser(Path(__file__).parent / "data" / "pwscf_tise2.out")
    a1, a2, a3 = parser_tise2.lattice_vectors_alat()
    assert np.allclose(a1, np.array([0.989_891, -0.001_560, -0.001_990]))
    assert np.allclose(a2, np.array([-0.496_296, 0.856_491, 0.001_990]))
    assert np.allclose(a3, np.array([-0.003_315, 0.001_914, 1.669_621]))


def test_pwscf_reciprocal_vectors_alat():
    """Test the parsing of reciprocal vectors in alat units"""
    parser_tise2 = PWSCFParser(Path(__file__).parent / "data" / "pwscf_tise2.out")
    b1, b2, b3 = parser_tise2.reciprocal_vectors_alat()
    assert np.allclose(b1, np.array([1.011_138, 0.585_904, 0.001_336]))
    assert np.allclose(b2, np.array([0.001_839, 1.168_623, -0.001_336]))
    assert np.allclose(b3, np.array([0.001_203, -0.000_695, 0.598_942]))


@pytest.mark.parametrize("file", PWSCF_FILES)
def test_pwscf_atoms(file):
    parser = PWSCFParser(file)
    atoms = parser.atoms()
    assert len(atoms) == parser.natoms


@pytest.mark.parametrize("file,expected_spg", zip(PWSCF_FILES, [164, 62, 194]))
def test_pwscf_crystal_instance(file, expected_spg):
    """Test the construction of Crystal instances, and check against expected symmetry properties"""

    parser = PWSCFParser(file)
    crystal = Crystal.from_pwscf(file)

    # COmparison of lattice vectors
    assert np.allclose(
        np.array(crystal.lattice_vectors),
        np.array(parser.lattice_vectors()),
    )

    # Comparison requires slighly relaxed precision
    assert crystal.symmetry(symprec=1e-1)["international_number"] == expected_spg


@pytest.mark.parametrize(
    "poscar", ["issue5_MgSiO3.cartesian.POSCAR", "issue5_MgSiO3.direct.POSCAR"]
)
def test_pwscf_lattice_vectors_alat(poscar):
    """Test POSCAR parser against CIF parser"""
    c0 = Crystal.from_cif(Path(__file__).parent / "data" / "issue5_MgSiO3.cif")
    c1 = Crystal.from_poscar(Path(__file__).parent / "data" / poscar)
    assert c0 == c1
