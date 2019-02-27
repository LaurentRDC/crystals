# -*- coding: utf-8 -*-
import os
import socket
import tempfile
import unittest
from collections import Counter, namedtuple
from contextlib import suppress
from itertools import chain
from pathlib import Path
from tempfile import gettempdir
from warnings import catch_warnings, filterwarnings

import numpy as np
from spglib import get_symmetry_dataset

from crystals import CIFParser, Crystal, PDBParser, frac_coords
from crystals.affine import transform
from crystals.parsers import STRUCTURE_CACHE, PWSCFParser
from crystals.spg_data import Hall2Number

try:
    import Bio.PDB as biopdb
except ImportError:
    WITH_BIOPYTHON = False
else:
    WITH_BIOPYTHON = True

# Used to compare crystals.Atom instances and Bio.PDB.Atom instances
GenericAtom = namedtuple("GenericAtom", ["element", "coords"])

filterwarnings("ignore", category=UserWarning)


def connection_available():
    """ Returns whether or not an internet connection is available """
    with suppress(OSError):
        try:
            socket.create_connection(("www.google.com", 80))
        except:
            return False
        else:
            return True
    return False


@unittest.skipUnless(connection_available(), "Internet connection is required.")
class TestPDBParser(unittest.TestCase):
    def test_fractional_atoms(self):
        """ Test the PDBParser returns fractional atomic coordinates. """
        with tempfile.TemporaryDirectory() as temp_dir:
            with PDBParser("1fbb", download_dir=temp_dir) as parser:
                for atm in parser.atoms():
                    self.assertLessEqual(atm.coords_fractional.max(), 1)
                    self.assertGreaterEqual(atm.coords_fractional.min(), 0)

    def test_symmetry_operators(self):
        """ Test that the non-translation part of the symmetry_operators is an invertible
        matrix of determinant 1 | -1 """
        with tempfile.TemporaryDirectory() as temp_dir:
            with PDBParser("1fbb", download_dir=temp_dir) as parser:
                for sym_op in parser.symmetry_operators():
                    t = sym_op[:3, :3]
                    self.assertAlmostEqual(abs(np.linalg.det(t)), 1, places=5)

    def test_residues(self):
        """ Test the parsing of residues for 1fbb """
        with tempfile.TemporaryDirectory() as temp_dir:
            with PDBParser("1fbb", download_dir=temp_dir) as parser:
                residues = list(parser.residues())
                atoms = list(parser.atoms())

        # Flatten residues into a list of atoms
        residue_atoms = list(chain.from_iterable(residues))
        self.assertEqual(len(residue_atoms), len(atoms))

    def test_default_download_dir(self):
        """ Test that the file is saved in the correct temporary directory by default """
        filename = PDBParser.download_pdb_file("1fbb")

        self.assertTrue(filename.exists())
        self.assertEqual(filename.parent, STRUCTURE_CACHE)


@unittest.skipUnless(WITH_BIOPYTHON, "Biopython is not installed/importable.")
@unittest.skipUnless(connection_available(), "Internet connection is required.")
class TestPDBParserAgainstBioPython(unittest.TestCase):

    # Each test will be performed on the following structures
    test_ids = ("1fbb", "1fat", "1gzx")

    def setUp(self):
        self.pdb_list = biopdb.PDBList(verbose=False, obsolete_pdb=gettempdir())
        self.biopdb_parser = biopdb.PDBParser()

    def test_chemical_composition(self):
        """ Test crystals.PDBParser returns the same chemical composition as BIO.PDB.PDBParser implementation,
        i.e. the same elements in the right proportions. """
        with catch_warnings():
            filterwarnings(
                "ignore", category=biopdb.PDBExceptions.PDBConstructionWarning
            )
            with tempfile.TemporaryDirectory() as temp_dir:
                for protein_id in self.test_ids:
                    with self.subTest(f"Protein ID: {protein_id}"):
                        with PDBParser(protein_id, download_dir=temp_dir) as parser:
                            fname = self.pdb_list.retrieve_pdb_file(
                                protein_id, pdir=temp_dir, file_format="pdb"
                            )

                            # Note: Bio.PDB atoms store element as uppercase strings. Thus, they must be changed to titlecase
                            crystals_chemical_composition = Counter(
                                [atm.element for atm in parser.atoms()]
                            )
                            biopdb_chemical_composition = Counter(
                                [
                                    atm.element.title()
                                    for atm in self.biopdb_parser.get_structure(
                                        protein_id, fname
                                    ).get_atoms()
                                ]
                            )

                            self.assertDictEqual(
                                biopdb_chemical_composition,
                                crystals_chemical_composition,
                            )

    @unittest.skip("")
    def test_atomic_positions(self):
        """ Test crystals.PDBParser returns atoms in the same position as the BIO.PDB.PDBParser implementation """
        with catch_warnings():
            filterwarnings(
                "ignore", category=biopdb.PDBExceptions.PDBConstructionWarning
            )

            with tempfile.TemporaryDirectory() as temp_dir:
                for protein_id in ("1fbb", "1fat", "1gzx"):
                    with self.subTest(f"Protein ID: {protein_id}"):
                        with PDBParser(protein_id, download_dir=temp_dir) as parser:
                            fname = self.pdb_list.retrieve_pdb_file(
                                protein_id, pdir=temp_dir, file_format="pdb"
                            )

                            biopdb_atoms = self.biopdb_parser.get_structure(
                                protein_id, fname
                            ).get_atoms()
                            crystals_atoms = parser.atoms()

                            # To compare atom positions, we build "generic" atoms (just tuples (elem, coords))
                            # Note: Bio.PDB atoms store element as uppercase strings. Thus, they must be changed to titlecase
                            # Since numpy arrays are unhashable, they are converted to tuples
                            # crystals.PDBParser returns atoms in fractional coordinates, so we must also do the same with Bio.PDB atoms
                            bio_pdb_generic_atoms = set()
                            for atm in biopdb_atoms:
                                coords = np.round(
                                    frac_coords(
                                        atm.coord_fractional, parser.lattice_vectors()
                                    ),
                                    3,
                                )
                                elem = atm.element.title()
                                bio_pdb_generic_atoms.add(
                                    GenericAtom(elem, tuple(coords))
                                )

                            crystals_generic_atoms = set()
                            for atm in crystals_atoms:
                                coords = np.round(atm.coords_fractional, 3)
                                crystals_generic_atoms.add(
                                    GenericAtom(atm.element, tuple(coords))
                                )
                            self.assertEqual(
                                bio_pdb_generic_atoms, crystals_generic_atoms
                            )


class TestCIFParser(unittest.TestCase):
    """ Test the CIFParser on all CIF files stored herein """

    def _cif_files(self):
        """ Yields cif files included in crystals """
        for root, _, files in os.walk(os.path.join("crystals", "cifs")):
            for name in filter(lambda path: path.endswith(".cif"), files):
                yield os.path.join(root, name)

    def test_compatibility(self):
        """ Test the CIFParser on all CIF files stored herein to check build errors"""
        for name in self._cif_files():
            with self.subTest(name.split("\\")[-1]):
                Crystal.from_cif(name)

    def test_fractional_atoms(self):
        """ Test the CIFParser returns fractional atomic coordinates. """
        for name in self._cif_files():
            with self.subTest(name.split("\\")[-1]):
                with CIFParser(name) as p:
                    for atm in p.atoms():
                        self.assertLessEqual(atm.coords_fractional.max(), 1)
                        self.assertGreaterEqual(atm.coords_fractional.min(), 0)

    def test_symmetry_operators(self):
        """ Test that the non-translation part of the symmetry_operators is an invertible
        matrix of determinant 1 | -1 """
        for name in self._cif_files():
            with self.subTest(name.split("\\")[-1]):
                with CIFParser(name) as p:
                    for sym_op in p.symmetry_operators():
                        t = sym_op[:3, :3]
                        self.assertAlmostEqual(abs(np.linalg.det(t)), 1)

    def test_international_number(self):
        """ Test that the international space group number  found by 
        CIFParser is the same as spglib's """
        for name in self._cif_files():
            with self.subTest(name.split("\\")[-1]):
                with CIFParser(name) as p:
                    from_parser = Hall2Number[p.hall_symbol()]

                    crystal = Crystal.from_cif(name)
                    from_spglib = crystal.international_number
                    self.assertEqual(from_parser, from_spglib)

    def test_silicon(self):
        """ Test CIFParser on Si.cif (diamond structure) """
        Si_path = os.path.join("crystals", "cifs", "Si.cif")
        si = Crystal.from_cif(Si_path)

        self.assertEqual(len(si), 8)

    def test_vo2(self):
        """ Test CIFParser on vo2.cif (monoclinic M1) """
        VO2_path = os.path.join("crystals", "cifs", "vo2-m1.cif")
        vo2 = Crystal.from_cif(VO2_path)

        self.assertEqual(len(vo2), 12)
        self.assertSequenceEqual(
            vo2.lattice_parameters,
            (
                5.743_000_000_000_000_3,
                4.517_000_000_000_000_3,
                5.375,
                90.0,
                122.600_000_000_000_01,
                90.0,
            ),
        )
        self.assertAlmostEqual(vo2.volume, 117.466_153_0)  # from cif2cell


class TestPWSCFParser(unittest.TestCase):
    def setUp(self):
        self.parser_tise2 = PWSCFParser(
            Path(".") / "tests" / "data" / "pwscf_tise2.out"
        )
        self.parser_snse = PWSCFParser(Path(".") / "tests" / "data" / "pwscf_snse.out")

    def test_alat(self):
        """ Test the parsing of the lattice parameter (alat) """
        self.assertEqual(self.parser_tise2.alat, 6.6764)
        self.assertEqual(self.parser_snse.alat, 21.4862)

    def test_natoms(self):
        """ Test the parsing of the number of unit cell atoms """
        self.assertEqual(self.parser_tise2.natoms, 3)
        self.assertEqual(self.parser_snse.natoms, 8)

    def test_lattice_vectors_alat(self):
        """ Test the parsing of lattice vectors in alat units """
        a1, a2, a3 = self.parser_tise2.lattice_vectors_alat()
        self.assertTrue(np.allclose(a1, np.array([0.989_891, -0.001_560, -0.001_990])))
        self.assertTrue(np.allclose(a2, np.array([-0.496_296, 0.856_491, 0.001_990])))
        self.assertTrue(np.allclose(a3, np.array([-0.003_315, 0.001_914, 1.669_621])))

    def test_reciprocal_vectors_alat(self):
        """ Test the parsing of reciprocal vectors in alat units """
        b1, b2, b3 = self.parser_tise2.reciprocal_vectors_alat()
        self.assertTrue(np.allclose(b1, np.array([1.011_138, 0.585_904, 0.001_336])))
        self.assertTrue(np.allclose(b2, np.array([0.001_839, 1.168_623, -0.001_336])))
        self.assertTrue(np.allclose(b3, np.array([0.001_203, -0.000_695, 0.598_942])))

    def test_atoms(self):
        for parser in (self.parser_tise2, self.parser_snse):
            atoms = parser.atoms()
            self.assertEqual(len(atoms), parser.natoms)

    def test_crystal_instance(self):
        """ Test the construction of Crystal instances, and check against expected symmetry properties """

        for parser, expected_spg in zip(
            (self.parser_tise2, self.parser_snse), (164, 62)
        ):
            with self.subTest(parser.filename):
                crystal = Crystal.from_pwscf(parser.filename)

                # COmparison of lattice vectors
                self.assertTrue(
                    np.allclose(
                        np.array(crystal.lattice_vectors),
                        np.array(parser.lattice_vectors()),
                    )
                )

                # Comparison requires slighly relaxed precision
                self.assertEqual(
                    crystal.symmetry(symprec=1e-1)["international_number"], expected_spg
                )


if __name__ == "__main__":
    unittest.main()
