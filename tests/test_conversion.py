# -*- coding: utf-8 -*-

import unittest
from random import choice
from random import randint

import numpy as np

from crystals import Crystal
from crystals import ase_atoms

try:
    import ase
except ImportError:
    WITH_ASE = False
else:
    WITH_ASE = True


@unittest.skipIf(not WITH_ASE, "ASE not installed or importable")
class TestAseAtoms(unittest.TestCase):
    def setUp(self):
        name = choice(list(Crystal.builtins))
        self.crystal = Crystal.from_database(name)

    def test_construction(self):
        """ Test that ase_atoms returns without error """
        to_ase = ase_atoms(self.crystal)
        self.assertEqual(len(self.crystal), len(to_ase))

    def test_back_and_forth(self):
        """ Test conversion to and from ase Atoms """
        to_ase = ase_atoms(self.crystal)
        crystal2 = Crystal.from_ase(to_ase)

        # ase has different handling of coordinates which can lead to
        # rounding beyond 1e-3. Therefore, we cannot compare directly sets
        # self.assertSetEqual(set(self.crystal), set(crystal2))
        self.assertEqual(len(self.crystal), len(crystal2))


if __name__ == "__main__":
    unittest.main()
