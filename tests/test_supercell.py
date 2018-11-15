# -*- coding: utf-8 -*-

import unittest
import numpy as np
from crystals import Crystal, Supercell, Atom


class TestSupercell(unittest.TestCase):
    def test_length(self):
        """ Test that supercell length matches dimensions """
        c = Crystal.from_database("C")
        s = c.supercell(2, 2, 2)

        self.assertEqual(len(s), 8 * len(c))

    def test_conversion_to_crystal(self):
        """ Test conversion of crystal to supercell, back to crystal """
        c = Crystal(unitcell=[Atom('Ag', [0,0,0])], lattice_vectors=2*np.eye(3))
        s = c.supercell(2, 2, 2)

        self.assertEqual(
            c, s.crystal
        )


if __name__ == "__main__":
    unittest.main()
