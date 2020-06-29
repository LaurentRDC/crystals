# -*- coding: utf-8 -*-

import tempfile
import unittest
from random import choice, randint

import numpy as np
from pathlib import Path
from crystals import Crystal, ase_atoms

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


class TestCifWriter(unittest.TestCase):

    def test_idempotence(self):
        """ Test that conversion to CIF of a structure loaded from CIF is idempotent. """
        # Testing on all built-in structure assures us that corner cases
        # are taken care of.
        for name in Crystal.builtins:
            with self.subTest(f"CIF idempotence {name}"):
                
                cryst = Crystal.from_database(name)
                with tempfile.TemporaryDirectory() as temp_dir:
                    f = Path(temp_dir) / "temp.cif"
                    cryst.to_cif(f)
                    cryst2 = Crystal.from_cif(f)
                    self.assertEqual(cryst, cryst2)



if __name__ == "__main__":
    unittest.main()
