# -*- coding: utf-8 -*-

import unittest
from itertools import islice

import numpy as np

from crystals import Atom, Crystal, Supercell


class TestSupercell(unittest.TestCase):
    def test_constructors(self):
        """ Test Supercell constructors for varyous 'builtin' structures """
        for name in islice(Crystal.builtins, 20):
            with self.subTest(name):
                s = Supercell.from_database(name, dimensions=(2, 2, 2))

                self.assertEqual(len(s), 8 * len(list(s.unitcell)))

    def test_conversion_to_crystal(self):
        """ Test conversion of crystal to supercell, back to crystal """
        c = Crystal(unitcell=[Atom("Ag", [0, 0, 0])], lattice_vectors=2 * np.eye(3))
        s = c.supercell(2, 2, 2)

        self.assertEqual(c, s.crystal)

    def test_symmetry(self):
        """ Test that the symmetry information is the same between a Crystal and associated Supercell """
        for name in islice(Crystal.builtins, 10):
            with self.subTest(name):
                c = Crystal.from_database(name)
                s = c.supercell(2, 2, 2)

                self.assertEqual(c.symmetry(), s.symmetry())

    def test_string_representation(self):
        """ Test that the string representation of a Supercell is similar to Crystals """
        for name in islice(Crystal.builtins, 20):
            with self.subTest(name):
                cryst_lines = str(Crystal.from_database(name)).split("\n")
                super_lines = str(
                    Supercell.from_database(name, dimensions=(1, 1, 1))
                ).split("\n")

                self.assertEqual(len(cryst_lines) + 2, len(super_lines))


if __name__ == "__main__":
    unittest.main()
