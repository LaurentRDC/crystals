# -*- coding: utf-8 -*-
import collections.abc as abc
import pickle
import unittest
from copy import deepcopy
from itertools import islice

import numpy as np

from crystals import Atom, AtomicStructure, Crystal

np.random.seed(23)


class TestAtomicStructure(unittest.TestCase):
    def setUp(self):
        self.substructure = AtomicStructure(atoms=[Atom("U", [0, 0, 0])])
        self.structure = AtomicStructure(
            atoms=[Atom("Ag", [0, 0, 0]), Atom("Ag", [1, 1, 1])],
            substructures=[self.substructure],
        )

    def test_iteration(self):
        """ Test iteration of AtomicStructure yields from orphan atoms and substructure atoms alike """
        elements = [atm.element for atm in self.structure]
        self.assertTrue(len(elements), 3)

    def test_addition_trivial(self):
        """ Test that the addition of two AtomicStructures, one being empty, works as expected """
        addition = self.structure + AtomicStructure()
        self.assertEqual(addition, self.structure)
        self.assertIsNot(addition, self.structure)

    def test_addition_uniqueness(self):
        """ Test that the addition of two AtomicStructures, works as expected regarding unique atoms """
        self.assertEqual(self.structure + self.structure, self.structure)

    def test_addition(self):
        """ Test the addition of two different AtomicStructures works as expected. """
        new_struct = AtomicStructure(
            atoms=[Atom("U", [0, 1, 0])],
            substructures=[
                AtomicStructure(
                    atoms=[Atom("Ag", [0.5, 0, 0]), Atom("Ag", [1, 0.3, 1])]
                )
            ],
        )

        addition = self.structure + new_struct

        self.assertEqual(len(new_struct) + len(self.structure), len(addition))
        self.assertEqual(
            len(new_struct.atoms) + len(self.structure.atoms), len(addition.atoms)
        )
        self.assertEqual(
            len(new_struct.substructures) + len(self.structure.substructures),
            len(addition.substructures),
        )

    def test_addition_subclasses(self):
        """ Test that the addition of two subclass of AtomicStructures is preserved under addition. """

        class NewAtomicStructure(AtomicStructure):
            pass

        addition = NewAtomicStructure() + NewAtomicStructure()
        self.assertIs(type(addition), NewAtomicStructure)

    def test_truthiness(self):
        """ Test that empty AtomicStructures are falsey, and truthy otherwise. """
        empty_structure = AtomicStructure()
        self.assertFalse(empty_structure)

        self.assertTrue(self.structure)

    def test_trivial_transformation(self):
        """ Test that the identity transformation of an AtomicStructure works as expected. """
        transformed = self.structure.transform(np.eye(3))

        # transformed structure should be different, but equal, to original structure
        self.assertIsNot(transformed, self.structure)
        self.assertEqual(transformed, self.structure)

    def test_transformations_inversions(self):
        """ Test that symmetry operations work as expected when inverted. """
        operator = np.random.random(size=(3, 3))
        inv_op = np.linalg.inv(operator)

        transformed1 = self.structure.transform(operator)
        transformed2 = transformed1.transform(inv_op)

        # transformed2 structure should be different, but equal, to original structure
        self.assertIsNot(transformed2, self.structure)
        self.assertEqual(transformed2, self.structure)

    def test_transform_subclass(self):
        """ Test that the object returned by the transform() method is the 
        same class as the method caller. """

        class NewAtomicStructure(AtomicStructure):
            pass

        structure = NewAtomicStructure(
            atoms=[Atom("Ag", [0, 0, 0]), Atom("Ag", [1, 1, 1])]
        )
        transformed = structure.transform(np.eye(3))

        self.assertIs(type(transformed), type(structure))

    def test_transformations_correctness(self):
        """ Test that AtomicStructure.transform() works as expected. """
        operator = 2 * np.eye(3)
        transformed = self.structure.transform(operator)

        expected_atoms = [atm.transform(operator) for atm in self.structure]

        for atm in expected_atoms:
            self.assertIn(atm, transformed)

    def test_itersorted(self):
        """ Test that AtomicStructure.itersorted() works as expected """
        sorted_from_structure = list(self.structure.itersorted())
        sorted_from_list = list(sorted(self.structure, key=lambda a: a.element))

        self.assertListEqual(sorted_from_structure, sorted_from_list)

    def test_chemical_composition_trivial(self):
        """ Test that AtomicStructure.chemical_composition works as expected """
        expected = {"U": 1 / 3, "Ag": 2 / 3}
        self.assertDictEqual(self.structure.chemical_composition, expected)

    def test_chemical_composition_add_to_unity(self):
        """ Test that AtomicStructure.chemical_composition always adds up to 1 """
        # Faster to create a large atomic structure from a Crystal object
        # Testing for 10 crystal structures only
        for name in islice(Crystal.builtins, 10):
            with self.subTest("Chemical composition: " + name):
                structure = AtomicStructure(atoms=Crystal.from_database(name))
                self.assertAlmostEqual(sum(structure.chemical_composition.values()), 1)

    def test_chemical_formula(self):
        """ Test that AtomicStructure.chemical_formula is working as expected. """
        self.assertEqual(self.structure.chemical_formula, "Ag2 U")

    def test_chemical_formula_hill_notation(self):
        """ Test that the Hill notation, where elements are alphabetically ordered except C and H, which are first. """
        structure = AtomicStructure(
            atoms=[
                Atom("Ag", [0, 1, 0]),
                Atom("C", [0, 0, 0]),
                Atom("H", [0, 1, 0]),
                Atom("U", [1, 1, 1]),
            ]
        )
        self.assertEqual(structure.chemical_formula, "C H Ag U")

    def test_length(self):
        """ Test the __len__ methods """
        self.assertTrue(len(self.structure), 3)

    def test_containership_substructures(self):
        """ Test that containership works on substructure and atoms separately """
        self.assertIn(self.substructure, self.structure)
        self.assertNotIn(self.structure, self.substructure)

    def test_containership_atoms(self):
        """ Test that atom containership testing is working, even in substructures """
        atm = next(iter(self.substructure))
        self.assertIn(atm, self.structure)

    def test_equality(self):
        """ Test that AtomicStructure is equal to itself but not others """
        self.assertEqual(self.structure, self.structure)
        self.assertEqual(self.structure, deepcopy(self.structure))
        self.assertNotEqual(self.structure, self.substructure)

        # Special case: make structures from Crystals
        c1 = Crystal.from_database("vo2-m1")
        c2 = deepcopy(c1)
        s1 = AtomicStructure(atoms=c1)
        s2 = AtomicStructure(atoms=c2.atoms)
        self.assertEqual(s1, s2)

    def test_array(self):
        """ Test AtomicStructure.__array__ """
        arr = np.array(self.structure)
        self.assertSequenceEqual(arr.shape, (len(self.structure), 4))

    def test_picklable(self):
        """ Test that Crystal instances can be pickled, and that the unpickled instance
        is identical to the source """
        pickled = pickle.dumps(self.structure)
        unpickled = pickle.loads(pickled)
        self.assertEqual(self.structure, unpickled)

    def test_abstract_base_classes(self):
        """ Test that AtomicStructure fits with collections.abc module """
        for abstract_base_class in (abc.Hashable, abc.Iterable, abc.Sized):
            self.assertIsInstance(self.structure, abstract_base_class)

    def test_satisfying(self):
        """ Test the AtomicStructure.satisfying method """
        uranium = self.structure.satisfying(lambda a: a.element == "U")
        silver = self.structure.satisfying(lambda a: a.element == "Ag")

        self.assertEqual(len(uranium), 1)
        self.assertEqual(len(silver), 2)


if __name__ == "__main__":
    unittest.main()
