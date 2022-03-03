# -*- coding: utf-8 -*-
import tempfile
from copy import deepcopy
from functools import lru_cache
from itertools import islice
from pathlib import Path
import socket
from contextlib import suppress

import numpy as np
from crystals import Atom, AtomicStructure, CenteringType, Crystal
from crystals.affine import translation_matrix
from crystals.crystal import symmetry_expansion, symmetry_reduction
import pytest

np.random.seed(23)


@lru_cache(maxsize=1)
def connection_available():
    """Returns whether or not an internet connection is available"""
    with suppress(OSError), socket.create_connection(("www.google.com", 80)):
        return True
    return False


def test_symmetry_graphite():
    """Test that Crystal.symmetry() works correctly for graphite"""
    c = Crystal.from_database("C")
    info = c.symmetry(1e-1)

    supposed = {
        "international_number": 194,
        "hall_number": 488,
        "international_symbol": "P6_3/mmc",
        "international_full": "P 6_3/m 2/m 2/c",
        "hall_symbol": "-P 6c 2c",
        "hm_symbol": "P63/mmc",
        "centering": CenteringType("P"),
        "pointgroup": "6/mmm",
    }

    assert info == supposed


@pytest.mark.parametrize("name", Crystal.builtins)
def test_primitive_for_builtins(name):
    """Test that all built-in crystal have a primitive cell"""
    c = Crystal.from_database(name)
    prim = c.primitive(symprec=0.1)
    assert len(prim) <= len(c)


def test_primitive_preserves_subclass():
    """Check that Crystal returned by Crystal.primitive() preserve subclass"""

    class TestCrystal(Crystal):
        pass

    c = TestCrystal.from_database("C")
    prim = c.primitive()
    assert TestCrystal == type(prim)


@pytest.mark.parametrize("name", Crystal.builtins)
def test_ideal_for_builtins(name):
    """Test that all built-in crystal have an ideal cell"""
    # This will raise an error if no idealized cell is found
    c = Crystal.from_database(name).ideal()


def test_ideal_preserves_subclass():
    """Check that Crystal returned by Crystal.ideal() preserve subclass"""

    class TestCrystal(Crystal):
        pass

    c = TestCrystal.from_database("C")
    ideal = c.ideal()
    assert TestCrystal == type(ideal)


@pytest.mark.parametrize("name", Crystal.builtins)
def test_symmetry_operations(name):
    """Test that the symmetry operations output makes sense"""
    identity = np.eye(4)

    c = Crystal.from_database(name)
    symops = c.symmetry_operations()
    assert np.allclose(identity, symops[0])


@pytest.mark.parametrize("name", Crystal.builtins)
def test_reciprocal_symmetry_operations(name):
    """Test that the reciprocal symmetry operations output makes sense"""
    identity = np.eye(4)

    c = Crystal.from_database(name)
    symops = c.reciprocal_symmetry_operations()
    assert np.allclose(identity, symops[0])


@pytest.mark.parametrize("name", Crystal.builtins)
def test_str_vs_repr(name):
    """Test that str and repr are workign as expected"""
    c = Crystal.from_database(name)

    # If small crystal, repr and str should be the same
    if len(c) <= 10:
        assert repr(c) == str(c)
    else:
        assert repr(c) != str(c)


def test_equality():
    """Test that __eq__ works as expected"""
    c1 = Crystal.from_database("Pu-alpha")
    c2 = deepcopy(c1)
    assert c1 == c2

    c3 = Crystal.from_database("Pu-epsilon")
    assert c1 != c3


def test_containership():
    """Test that __contains__ works as expected"""
    graphite = Crystal.from_database("C")
    carbon = next(iter(graphite))  # Pick any atom from the crystal

    assert carbon in graphite


@pytest.mark.parametrize("name", Crystal.builtins)
def test_builtins(name):
    """Test that all names in Crystal.builtins build without errors,
    and that Crystal.source is correctly recorded."""
    c = Crystal.from_database(name)

    assert name in c.source


def test_builtins_wrong_name():
    """Test that a name not in Crystal.builtins will raise a ValueError"""
    with pytest.raises(ValueError):
        Crystal.from_database("___")


def test_substructure_preservation():
    """Test that initializing a crystal with substructures preserves the substructures"""
    atoms = [Atom("Ag", [0, 0, 0]), Atom("Ag", [1, 1, 1])]
    substructures = [AtomicStructure(atoms=[Atom("U", [0, 0, 0])])]
    c = Crystal(unitcell=atoms + substructures, lattice_vectors=np.eye(3))

    assert len(c) == 3
    assert substructures[0] in c.substructures


@pytest.mark.skipif(
    not connection_available(), reason="Internet connection is required."
)
def test_from_pdb():
    """Test Crystal.from_pdb constructor"""
    c = Crystal.from_pdb("1fbb")
    assert "1fbb" in c.source


@pytest.mark.skipif(
    not connection_available(), reason="Internet connection is required."
)
def test_from_cod():
    """Test building a Crystal object from the COD"""
    # revision = None and latest revision should give the same Crystal
    c = Crystal.from_cod(1521124)
    c2 = Crystal.from_cod(1521124, revision=176429)

    assert c == c2


@pytest.mark.skipif(
    not connection_available(), reason="Internet connection is required."
)
def test_from_cod_new_dir():
    """Test that a cache dir is created by Crystal.from_cod"""
    with tempfile.TemporaryDirectory() as temp_dir:
        download_dir = Path(temp_dir) / "test_cod"
        assert not download_dir.exists()
        c = Crystal.from_cod(1521124, download_dir=download_dir)
        assert download_dir.exists()


@pytest.mark.parametrize("name", islice(Crystal.builtins, 20))
def test_supercell_constructors(name):
    """Test Supercell constructors for varyous 'builtin' structures"""
    c = Crystal.from_database(name)
    s = c.supercell(2, 2, 2)

    assert len(s) == 8 * len(c)


@pytest.mark.parametrize("tag", [None, 1, 2])
@pytest.mark.parametrize("occupancy", [1.0, 1.5, 2])
@pytest.mark.parametrize("symbol", ["H", "He", "C"])
def test_supercell_preserve_attributes(tag, symbol, occupancy):
    """Test that the `Crystal.supercell` method preserves atom attributes (see issue #9)"""
    atom = Atom(element=symbol, coords=[0, 0, 0], tag=tag, occupancy=occupancy)
    c = Crystal(unitcell=[atom], lattice_vectors=3 * np.eye(3))
    s = c.supercell(2, 2, 2)

    for supercell_atom in s:
        assert supercell_atom.tag == atom.tag
        assert supercell_atom.symbol == atom.symbol


def test_indexed_by_trivial_reindexing():
    """Test re-indexing a crystal by itself"""
    c1 = Crystal.from_database("Pu-gamma")
    c2 = c1.indexed_by(c1)

    assert c1 == c2


def test_symmetry_reduction_trivial():
    """Test that the symmetry_reduction function returns the unit cell when
    there is only one possibility."""
    ucell = set([Atom("H", coords=[0, 0, 0])])
    symops = [np.eye(4)]
    asym_cell = symmetry_reduction(ucell, symops)
    assert ucell == asym_cell


def test_symmetry_reduction_simple_translation():
    """Test that symmetry_reduction works on a unitcell where two atoms are
    linked by a translation"""
    symops = [np.eye(4), translation_matrix([0, 0, 1 / 3])]
    asym_cell = set([Atom("H", coords=[0, 0, 0])])
    unitcell = set(symmetry_expansion(asym_cell, symmetry_operators=symops))

    asym_cell2 = symmetry_reduction(unitcell, symops)
    assert asym_cell == asym_cell2


@pytest.mark.parametrize("name", ["vo2-m1", "Os", "Na"])
def test_symmetry_reduction_reciprocity_with_symmetry_expansion(name):
    """Test that symmetry_reduction is reciprocal to symmetry_expansion"""
    cryst = Crystal.from_database(name)
    asym_cell = symmetry_reduction(cryst.unitcell, cryst.symmetry_operations())
    ucell = set(symmetry_expansion(asym_cell, cryst.symmetry_operations()))
    assert set(cryst.unitcell) == ucell
