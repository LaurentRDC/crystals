# -*- coding: utf-8 -*-
import collections.abc as abc
import pickle
import pytest
from copy import deepcopy
from itertools import islice

import numpy as np

from crystals import Atom, AtomicStructure, Crystal

np.random.seed(23)


@pytest.fixture
def structure():
    substructure = AtomicStructure(atoms=[Atom("U", [0, 0, 0])])
    return AtomicStructure(
        atoms=[Atom("Ag", [0, 0, 0]), Atom("Ag", [1, 1, 1])],
        substructures=[substructure],
    )


def test_iteration(structure):
    """Test iteration of AtomicStructure yields from orphan atoms and substructure atoms alike"""
    elements = [atm.element for atm in structure]
    assert len(elements), 3


def test_addition_trivial(structure):
    """Test that the addition of two AtomicStructures, one being empty, works as expected"""
    addition = structure + AtomicStructure()
    assert addition == structure
    assert addition is not structure


def test_addition_uniqueness(structure):
    """Test that the addition of two AtomicStructures, works as expected regarding unique atoms"""
    assert structure + structure == structure


def test_addition(structure):
    """Test the addition of two different AtomicStructures works as expected."""
    new_struct = AtomicStructure(
        atoms=[Atom("U", [0, 1, 0])],
        substructures=[
            AtomicStructure(atoms=[Atom("Ag", [0.5, 0, 0]), Atom("Ag", [1, 0.3, 1])])
        ],
    )

    addition = structure + new_struct

    assert len(new_struct) + len(structure) == len(addition)
    assert len(new_struct.atoms) + len(structure.atoms) == len(addition.atoms)
    assert len(new_struct.substructures) + len(structure.substructures) == len(
        addition.substructures
    )


def test_addition_subclasses(structure):
    """Test that the addition of two subclass of AtomicStructures is preserved under addition."""

    class NewAtomicStructure(AtomicStructure):
        pass

    addition = NewAtomicStructure() + NewAtomicStructure()
    assert type(addition) is NewAtomicStructure


def test_truthiness(structure):
    """Test that empty AtomicStructures are falsey, and truthy otherwise."""
    empty_structure = AtomicStructure()
    assert not empty_structure

    assert structure


def test_trivial_transformation(structure):
    """Test that the identity transformation of an AtomicStructure works as expected."""
    transformed = structure.transform(np.eye(3))

    # transformed structure should be different, but equal, to original structure
    assert transformed is not structure
    assert transformed == structure


def test_transformations_inversions(structure):
    """Test that symmetry operations work as expected when inverted."""
    operator = np.random.random(size=(3, 3))
    inv_op = np.linalg.inv(operator)

    transformed1 = structure.transform(operator)
    transformed2 = transformed1.transform(inv_op)

    # transformed2 structure should be different, but equal, to original structure
    assert transformed2 is not structure
    assert transformed2 == structure


def test_transform_subclass(structure):
    """Test that the object returned by the transform() method is the
    same class as the method caller."""

    class NewAtomicStructure(AtomicStructure):
        pass

    structure = NewAtomicStructure(atoms=[Atom("Ag", [0, 0, 0]), Atom("Ag", [1, 1, 1])])
    transformed = structure.transform(np.eye(3))

    assert type(transformed) is type(structure)


def test_transformations_correctness(structure):
    """Test that AtomicStructure.transform() works as expected."""
    operator = 2 * np.eye(3)
    transformed = structure.transform(operator)

    expected_atoms = [atm.transform(operator) for atm in structure]

    for atm in expected_atoms:
        assert atm in transformed


def test_itersorted(structure):
    """Test that AtomicStructure.itersorted() works as expected"""
    sorted_from_structure = list(structure.itersorted())
    sorted_from_list = list(sorted(structure))

    assert sorted_from_structure == sorted_from_list


def test_chemical_composition_trivial(structure):
    """Test that AtomicStructure.chemical_composition works as expected"""
    expected = {"U": 1 / 3, "Ag": 2 / 3}
    assert structure.chemical_composition == expected


@pytest.mark.parametrize("name", islice(Crystal.builtins, 10))
def test_chemical_composition_add_to_unity(name):
    """Test that AtomicStructure.chemical_composition always adds up to 1"""
    # Faster to create a large atomic structure from a Crystal object
    # Testing for 10 crystal structures only
    structure = AtomicStructure(atoms=Crystal.from_database(name))
    assert round(abs(sum(structure.chemical_composition.values()) - 1), 7) == 0


def test_chemical_formula(structure):
    """Test that AtomicStructure.chemical_formula is working as expected."""
    assert structure.chemical_formula == "Ag2 U"


def test_chemical_formula_hill_notation(structure):
    """Test that the Hill notation, where elements are alphabetically ordered except C and H, which are first."""
    structure = AtomicStructure(
        atoms=[
            Atom("Ag", [0, 1, 0]),
            Atom("C", [0, 0, 0]),
            Atom("H", [0, 1, 0]),
            Atom("U", [1, 1, 1]),
        ]
    )
    assert structure.chemical_formula == "C H Ag U"


def test_length(structure):
    """Test the __len__ methods"""
    assert len(structure), 3


def test_containership_substructures(structure):
    """Test that containership works on substructure and atoms separately"""
    substructure = next(iter(structure.substructures))
    assert substructure in structure
    assert structure not in substructure


def test_containership_atoms(structure):
    """Test that atom containership testing is working, even in substructures"""
    atm = next(iter(structure))
    assert atm in structure


def test_equality(structure):
    """Test that AtomicStructure is equal to itself but not others"""
    substructure = next(iter(structure.substructures))
    assert structure == structure
    assert structure == deepcopy(structure)
    assert structure != substructure

    # Special case: make structures from Crystals
    c1 = Crystal.from_database("vo2-m1")
    c2 = deepcopy(c1)
    s1 = AtomicStructure(atoms=c1)
    s2 = AtomicStructure(atoms=c2.atoms)
    assert s1 == s2


def test_array(structure):
    """Test AtomicStructure.__array__"""
    arr = np.array(structure)
    assert arr.shape == (len(structure), 4)


def test_picklable(structure):
    """Test that Crystal instances can be pickled, and that the unpickled instance
    is identical to the source"""
    pickled = pickle.dumps(structure)
    unpickled = pickle.loads(pickled)
    assert structure == unpickled


def test_abstract_base_classes(structure):
    """Test that AtomicStructure fits with collections.abc module"""
    for abstract_base_class in (abc.Hashable, abc.Iterable, abc.Sized):
        assert isinstance(structure, abstract_base_class)


def test_satisfying(structure):
    """Test the AtomicStructure.satisfying method"""
    uranium = structure.satisfying(lambda a: a.element == "U")
    silver = structure.satisfying(lambda a: a.element == "Ag")

    assert len(uranium) == 1
    assert len(silver) == 2
