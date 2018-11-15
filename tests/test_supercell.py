# -*- coding: utf-8 -*-

import unittest

from crystals import Crystal, Supercell


class TestSupercell(unittest.TestCase):
    def test_length(self):
        """ Test that supercell length matches dimensions """
        c = Crystal.from_database("C")
        s = c.supercell(2, 2, 2)

        self.assertEqual(len(s), 8 * len(c))

    @unittest.skip("Not ready for primetime")
    def test_conversion_to_crystal(self):
        """ Test conversion of crystal to supercell, back to crystal """
        c = Crystal.from_database("C")
        s = c.supercell(2, 2, 2)
        c2 = s.crystal

        self.assertEqual(s.crystal, c)


if __name__ == "__main__":
    unittest.main()
