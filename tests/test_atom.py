# -*- coding: utf-8 -*-
import unittest
from copy import deepcopy
from random import randint
from random import random
from random import seed

import numpy as np

from crystals import Atom
from crystals import Lattice
from crystals.affine import rotation_matrix

seed(23)
np.random.seed(23)


def random_transform():
    return rotation_matrix(random(), axis=np.random.random((3,)))


class TestAtom(unittest.TestCase):
    def setUp(self):
        self.atom = Atom(randint(1, 103), coords=np.random.random((3,)))

    def test_init(self):
        """ Test that Atom can be instantiated with an element str or atomic number """
        by_element = Atom("C", coords=(0, 0, 0))
        by_number = Atom(6, coords=(0, 0, 0))
        self.assertEqual(by_element, by_number)

    def test_equality(self):
        """ Test __eq__ for atoms """
        other = deepcopy(self.atom)
        self.assertEqual(self.atom, self.atom)
        self.assertEqual(self.atom, other)

        other.coords_fractional = self.atom.coords_fractional + 1
        self.assertNotEqual(self.atom, other)

    def test_trivial_transform(self):
        """ Test Atom.transform() with the identity """
        before = np.array(self.atom.coords_fractional, copy=True)
        self.atom.transform(np.eye(3))
        after = np.array(self.atom.coords_fractional, copy=True)

        self.assertSequenceEqual(tuple(before), tuple(after))

    def test_transform_back_and_forth(self):
        """ Test Atom.transform() with a random transformation back and forth """
        before = np.array(self.atom.coords_fractional, copy=True)

        transf = random_transform()
        self.atom.transform(transf)
        self.atom.transform(np.linalg.inv(transf))
        after = np.array(self.atom.coords_fractional, copy=True)

        # No assert sequence almost equal
        for x1, x2 in zip(tuple(before), tuple(after)):
            self.assertAlmostEqual(x1, x2)

    def test_atom_array(self):
        """ Test that numpy.array(Atom(...)) works as expected """
        arr = np.array(self.atom)
        self.assertTupleEqual(arr.shape, (4,))
        self.assertEqual(arr[0], self.atom.atomic_number)
        self.assertTrue(np.allclose(arr[1::], self.atom.coords_fractional))

    def test_distance(self):
        """ Test the distance between atoms """
        atm1 = Atom("He", [0, 0, 0])
        atm2 = Atom("He", [1, 0, 0])
        self.assertEqual(atm1 - atm2, 1)


if __name__ == "__main__":
    unittest.main()
