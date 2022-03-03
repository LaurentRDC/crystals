# -*- coding: utf-8 -*-
from copy import deepcopy
from random import random, seed
import pickle

import numpy as np
from itertools import islice

from crystals import (
    Atom,
    Orbital,
    ElectronicStructure,
    Element,
    Lattice,
    distance_cartesian,
    distance_fractional,
    is_element,
)
from crystals.atom_data import chemical_symbols
from crystals.affine import rotation_matrix
import pytest

seed(23)
np.random.seed(23)

# Madelung rule
MADELUNG = [
    "1s",
    "2s",
    "2p",
    "3s",
    "3p",
    "4s",
    "3d",
    "4p",
    "5s",
    "4d",
    "5p",
    "6s",
    "4f",
    "5d",
    "6p",
    "7s",
    "5f",
    "6d",
    "7p",
]


def random_transform():
    """Create a rotation matrix, around a random axis, of a random amount"""
    return rotation_matrix(random(), axis=np.random.random((3,)))


@pytest.mark.parametrize("symbol", chemical_symbols)
def test_element_build(symbol):
    """Test that all valid chemical symbols can be used to create an Element instance"""
    Element(symbol)


def test_element_invalid_element():
    """Test that an invalid chemical symbol will raise an error"""
    with pytest.raises(ValueError):
        Element("montreal")


@pytest.mark.parametrize("symbol", chemical_symbols)
def test_element_init_with_atomic_number(symbol):
    """Test that constructing an Element from a symbol or atomic number results in the same element"""
    from_symbol = Element(symbol)
    from_number = Element(from_symbol.atomic_number)
    assert from_number == from_symbol


def test_element_atomic_number_out_of_range():
    """Test that constructing an Element from an atomic number out-of-range raises a ValueError"""
    with pytest.raises(ValueError):
        Element(256)


def test_element_input_equivalence():
    """Test that the constructor for `Element` supports atomic numbers,
    strings, and other `Element`s"""
    expected = Element("H")

    assert Element("H") == expected
    assert Element("Hydrogen") == expected
    assert Element("hydrogen") == expected  # not case sensitive
    assert Element(1) == expected
    assert Element(expected) == expected


def test_atom_init():
    """Test that Atom can be instantiated with an element str or atomic number"""
    by_element = Atom("C", coords=(0, 0, 0))
    by_number = Atom(6, coords=(0, 0, 0))
    assert by_element == by_number


@pytest.mark.parametrize("tag", [None, 1, 2])
@pytest.mark.parametrize("symbol", islice(chemical_symbols, 50))
def test_atom_equality(symbol, tag):
    """Test __eq__ for atoms"""
    atom = Atom(symbol, coords=[0, 0, 0], tag=tag)
    other = deepcopy(atom)
    assert atom == atom
    assert atom == other
    assert hash(atom) == hash(other)

    other.coords_fractional = atom.coords_fractional + 1
    assert atom != other
    assert hash(atom) != hash(other)

    other = deepcopy(atom)
    other.tag = 5
    assert atom != other
    assert hash(atom) != hash(other)


@pytest.mark.parametrize(
    "atom", map(lambda s: Atom(s, [0, 0, 0]), islice(chemical_symbols, 50))
)
def test_atom_trivial_transform(atom):
    """Test Atom.transform() with the identity"""
    transformed = atom.transform(np.eye(3))

    assert transformed is not atom
    assert transformed == atom


def test_atom_subclasses_transform_preserved():
    """Test transform() preserves subclasses"""

    class NewAtom(Atom):
        pass

    atm = NewAtom("He", [0, 0, 0])
    transformed = atm.transform(np.eye(3))

    assert type(atm) is type(transformed)


@pytest.mark.parametrize(
    "atom", map(lambda s: Atom(s, [0, 0, 0]), islice(chemical_symbols, 50))
)
def test_atom_transform_back_and_forth(atom):
    """Test Atom.transform() with a random transformation back and forth"""

    transf = random_transform()
    transformed1 = atom.transform(transf)
    transformed2 = transformed1.transform(np.linalg.inv(transf))

    assert transformed2 is not atom
    assert transformed2 == atom


@pytest.mark.parametrize("atom", map(lambda s: Atom(s, [0, 0, 0]), chemical_symbols))
def test_atom_atom_array(atom):
    """Test that numpy.array(Atom(...)) works as expected"""
    arr = np.array(atom)
    assert arr.shape == (4,)
    assert arr[0] == atom.atomic_number
    assert np.allclose(arr[1::], atom.coords_fractional)


def test_atomic_distance_fractional():
    """Test the fractional distance between atoms"""
    atm1 = Atom("He", [0, 0, 0])
    atm2 = Atom("He", [1, 0, 0])
    assert distance_fractional(atm1, atm2) == 1
    assert distance_fractional(atm1, atm2) == distance_fractional(atm2, atm1)


def test_atomic_distance_cartesian():
    """Test the cartesian distance between atom"""
    lattice = Lattice(4 * np.eye(3))  # Cubic lattice side length 4 angs

    atm1 = Atom("He", [0, 0, 0], lattice=lattice)
    atm2 = Atom("He", [1, 0, 0], lattice=lattice)

    assert distance_cartesian(atm1, atm2) == 4.0


def test_atomic_distance_different_lattice():
    """Test that fractional and cartesian distances
    between atoms in different lattices raises an error."""
    lattice1 = Lattice(np.eye(3))
    lattice2 = Lattice(2 * np.eye(3))

    atm1 = Atom("He", [0, 0, 0], lattice=lattice1)
    atm2 = Atom("He", [1, 0, 0], lattice=lattice2)

    with pytest.raises(RuntimeError):
        distance_fractional(atm1, atm2)

    with pytest.raises(RuntimeError):
        distance_cartesian(atm1, atm2)


def test_is_element_input_types():
    """Test that is_element() works as expected for all
    supported input types"""
    atm = Atom("V", [0, 0, 0])

    assert is_element(atm.element)(atm)
    assert is_element(atm.atomic_number)(atm)
    assert is_element(atm)(atm)


def test_orbital_madelung_rule():
    """Test that the orbitals are listed in the Madelung rule order,
    which is the filling order."""

    enumeration = [shell.value for shell in Orbital]
    assert MADELUNG == enumeration


def test_orbital_maximum_electrons():
    """That that the maximum number of electrons per Orbital is as expected."""
    maxima = {"s": 2, "p": 6, "d": 10, "f": 14}
    for shell in Orbital:
        assert maxima[shell.value[-1]] == Orbital.maximum_electrons(shell)


def test_electronic_structure_maximum_electrons():
    """Test that an error is raised for impossible electronic structures."""
    with pytest.raises(ValueError):
        ElectronicStructure({"1s": 3})


def test_electronic_structure_pickable():
    """Test that ElectronicStructure instances are pickable."""
    structure = ElectronicStructure.ground_state("W")
    assert structure == pickle.loads(pickle.dumps(structure))


def test_electronic_structure_missing_orbital():
    """Test that 'missing' orbitals will return 0 electrons"""
    struct = ElectronicStructure.ground_state("He")
    assert struct["2p"] == 0


def test_electronic_structure_orbital_keys():
    """Test that orbital occupancies can be accessed either with strings or Orbital"""
    struct = ElectronicStructure.ground_state("He")
    assert struct["1s"] == struct[Orbital.one_s]


def test_electronic_structure_modification_in_place():
    """Test that ElectronicStructure instances can be modified in-place"""
    struct = ElectronicStructure({"1s": 2, "2s": 2, "2p": 4})
    expected = ElectronicStructure({"1s": 2, "2s": 2, "2p": 3, "4s": 1})
    struct["2p"] -= 1
    struct["4s"] += 1
    assert struct == expected


def test_electronic_structure_valence_shell():
    """Test that the outermost shell is as expected"""
    struct = ElectronicStructure({"1s": 2})
    assert struct.outer_shell == Orbital("1s")

    # Intentionally omitting 2s
    struct = ElectronicStructure({"1s": 2, "2p": 4})
    assert struct.outer_shell == Orbital("2p")
