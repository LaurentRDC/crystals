# -*- coding: utf-8 -*-
import pickle
import socket
import tempfile
import unittest
from contextlib import suppress
from copy import copy, deepcopy
from math import radians
from itertools import islice
from pathlib import Path

import numpy as np

from crystals import Atom, AtomicStructure, Crystal, Lattice, CenteringType
from crystals.crystal import SymmetryOperation
from crystals.affine import rotation_matrix, transform


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


class TestSpglibMethods(unittest.TestCase):
    def test_symmetry_graphite(self):
        """ Test that Crystal.symmetry() works correctly for graphite """
        c = Crystal.from_database("C")
        info = c.symmetry()

        supposed = {
            "international_number": 194,
            "hall_number": 488,
            "international_symbol": "P6_3/mmc",
            "international_full": "P 6_3/m 2/m 2/c",
            "hall_symbol": "-P 6c 2c",
            "hm_symbol": "P63/mmc",
            "centering": CenteringType("P"),
            "pointgroup": "D6h",
        }

        self.assertDictEqual(info, supposed)

    def test_primitive_for_builtins(self):
        """ Test that all built-in crystal have a primitive cell """
        for name in Crystal.builtins:
            with self.subTest(name):
                c = Crystal.from_database(name)
                prim = c.primitive(symprec=0.1)
                self.assertLessEqual(len(prim), len(c))

    def test_primitive_preserves_subclass(self):
        """ Check that Crystal returned by Crystal.primitive() preserve subclass """

        class TestCrystal(Crystal):
            pass

        c = TestCrystal.from_database("C")
        prim = c.primitive()
        self.assertEqual(TestCrystal, type(prim))

    def test_ideal_for_builtins(self):
        """ Test that all built-in crystal have an ideal cell """
        for name in Crystal.builtins:
            with self.subTest(name):
                # This will raise an error if no idealized cell is found
                c = Crystal.from_database(name).ideal()

    def test_ideal_preserves_subclass(self):
        """ Check that Crystal returned by Crystal.ideal() preserve subclass """

        class TestCrystal(Crystal):
            pass

        c = TestCrystal.from_database("C")
        ideal = c.ideal()
        self.assertEqual(TestCrystal, type(ideal))

    def test_symmetry_operations(self):
        """ Test that the symmetry operations output makes sense """
        identity = SymmetryOperation(np.eye(3, dtype=np.int32), np.zeros((3,)))

        for name in Crystal.builtins:
            with self.subTest(name):
                c = Crystal.from_database(name)
                symops = c.symmetry_operations()

                self.assertTrue(np.allclose(identity.rotation, symops[0].rotation))

                self.assertTrue(
                    np.allclose(identity.translation, symops[0].translation)
                )

    def test_reciprocal_symmetry_operations(self):
        """ Test that the reciprocal symmetry operations output makes sense """
        identity = SymmetryOperation(np.eye(3, dtype=np.int32), np.zeros((3,)))

        for name in Crystal.builtins:
            with self.subTest(name):
                c = Crystal.from_database(name)
                symops = c.reciprocal_symmetry_operations()

                self.assertTrue(np.allclose(identity.rotation, symops[0].rotation))

                self.assertTrue(
                    np.allclose(identity.translation, symops[0].translation)
                )


class TestCrystalSpecialMethods(unittest.TestCase):
    def test_str_vs_repr(self):
        """ Test that str and repr are workign as expected """
        for name in Crystal.builtins:
            with self.subTest(name):
                c = Crystal.from_database(name)

                # If small crystal, repr and str should be the same
                if len(c) <= 10:
                    self.assertEqual(repr(c), str(c))
                else:
                    self.assertNotEqual(repr(c), str(c))

    @unittest.skip("Not ready for primetime")
    def test_equality(self):
        """ Test that __eq__ works as expected """
        c1 = Crystal.from_database("Pu-alpha")
        c2 = deepcopy(c1)
        self.assertEqual(c1, c2)


class TestCrystalConstructors(unittest.TestCase):
    def test_builtins(self):
        """ Test that all names in Crystal.builtins build without errors,
        and that Crystal.source is correctly recorded. """
        for name in Crystal.builtins:
            with self.subTest(name):
                c = Crystal.from_database(name)

                self.assertIn(name, c.source)

    def test_builtins_wrong_name(self):
        """ Test that a name not in Crystal.builtins will raise a ValueError """
        with self.assertRaises(ValueError):
            Crystal.from_database("___")

    def test_substructure_preservation(self):
        """ Test that initializing a crystal with substructures preserves the substructures """
        atoms = [Atom("Ag", [0, 0, 0]), Atom("Ag", [1, 1, 1])]
        substructures = [AtomicStructure(atoms=[Atom("U", [0, 0, 0])])]
        c = Crystal(unitcell=atoms + substructures, lattice_vectors=np.eye(3))

        self.assertEqual(len(c), 3)
        self.assertIn(substructures[0], c.substructures)

    @unittest.skipUnless(connection_available(), "Internet connection is required.")
    def test_from_pdb(self):
        """ Test Crystal.from_pdb constructor """
        c = Crystal.from_pdb("1fbb")
        self.assertIn("1fbb", c.source)

    @unittest.skipUnless(connection_available(), "Internet connection is required.")
    def test_from_cod(self):
        """ Test building a Crystal object from the COD """
        # revision = None and latest revision should give the same Crystal
        c = Crystal.from_cod(1521124)
        c2 = Crystal.from_cod(1521124, revision=176429)

        self.assertEqual(c, c2)

    @unittest.skipUnless(connection_available(), "Internet connection is required.")
    def test_from_cod_new_dir(self):
        """ Test that a cache dir is created by Crystal.from_cod """
        with tempfile.TemporaryDirectory() as temp_dir:
            download_dir = Path(temp_dir) / "test_cod"
            self.assertFalse(download_dir.exists())
            c = Crystal.from_cod(1521124, download_dir=download_dir)
            self.assertTrue(download_dir.exists())


class TestSupercell(unittest.TestCase):
    def test_constructors(self):
        """ Test Supercell constructors for varyous 'builtin' structures """
        for name in islice(Crystal.builtins, 20):
            with self.subTest(name):
                s = Crystal.from_database(name).supercell(2, 2, 2)

                self.assertEqual(len(s), 8 * len(s.crystal))


if __name__ == "__main__":
    unittest.main()
