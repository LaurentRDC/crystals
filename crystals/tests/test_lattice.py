# -*- coding: utf-8 -*-
from copy import deepcopy
from math import degrees, radians
from itertools import islice

import numpy as np
from numpy.linalg import norm

from crystals import Crystal, Lattice, LatticeSystem
from crystals.lattice import lattice_vectors_from_parameters
from crystals.affine import rotation_matrix
import pytest

np.random.seed(23)


def angle_between(v1, v2):
    """Returns the angle in degrees between vectors `v1` and `v2`"""
    v1 /= norm(v1)
    v2 /= norm(v2)
    rad = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
    return degrees(rad)


def test_euclidian_lattice():
    lattice = Lattice(np.eye(3))
    assert np.array_equal(lattice.a1, [1, 0, 0])
    assert np.array_equal(lattice.a2, [0, 1, 0])
    assert np.array_equal(lattice.a3, [0, 0, 1])

    assert lattice.volume == 1


def test_euclidian_lattice_equality():
    """Test equality between identical Lattice instances and copies"""
    lattice = Lattice(np.eye(3))
    assert lattice == lattice
    assert lattice == deepcopy(lattice)

    assert lattice != Lattice(2 * np.eye(3))


def test_lattice_array_shape():
    """Test that array(Lattice(...)) is always 3x3"""
    arr = np.random.random(size=(3, 3))
    lattice = Lattice(arr)
    assert np.allclose(arr, np.array(lattice))


def test_lattice_array_dtype():
    """Test that the data-type of array(Lattice(...)) is respected"""
    arr = np.random.random(size=(3, 3))
    lattice = Lattice(arr)
    assert np.array(lattice, dtype=np.int32).dtype, np.int32


def test_lattice_frac_mesh():
    """Test that Lattice.frac_mesh is working as expected compared to numpy.meshgrid"""
    lattice = Lattice(np.eye(3))
    x = np.linspace(0, 1, num=8)

    for out, n in zip(lattice.frac_mesh(x), np.meshgrid(x, x, x)):
        assert np.allclose(out, n)


def test_lattice_frac_mesh_two_arr():
    """Test that Lattice.frac_mesh is raising an exception for two input vectors"""
    lattice = Lattice(np.eye(3))
    x = np.linspace(0, 1, num=2)

    with pytest.raises(ValueError):
        lattice.frac_mesh(x, x)


def test_lattice_real_mesh_trivial():
    """Test that Lattice.mesh works identically to Lattice.frac_mesh for trivial lattice"""
    lattice = Lattice(np.eye(3))
    x = np.linspace(0, 1, num=8)

    for frac, real in zip(lattice.frac_mesh(x), lattice.mesh(x)):
        assert np.allclose(frac, real)


def test_lattice_real_mesh():
    """Test that Lattice.mesh works as expected"""
    lattice = Lattice(2 * np.eye(3))
    x = np.linspace(0, 1, num=8)

    # since lattice is a stretched euclidian lattice, we expect
    # a maximum length of 2
    assert np.max(lattice.mesh(x)[0]) == 2


def test_lattice_parameters_orthorombic():
    """alpha = beta = gamma = 90"""
    a1, a2, a3 = Lattice.from_parameters(2, 1, 5, 90, 90, 90).lattice_vectors
    assert np.allclose(a1, [2, 0, 0])
    assert np.allclose(a2, [0, 1, 0])
    assert np.allclose(a3, [0, 0, 5])


def test_lattice_parameters_monoclinic():
    """beta =/= 90"""
    a1, a2, a3 = Lattice.from_parameters(1, 2, 3, 90, 120, 90).lattice_vectors

    assert round(abs(norm(a1) - 1), 7) == 0
    assert round(abs(norm(a2) - 2), 7) == 0
    assert round(abs(norm(a3) - 3), 7) == 0

    assert round(abs(angle_between(a2, a3) - 90), 7) == 0
    assert round(abs(angle_between(a3, a1) - 120), 7) == 0
    assert round(abs(angle_between(a1, a2) - 90), 7) == 0


def test_lattice_parameters_triclinic():
    """alpha, beta, gama =/= 90"""
    a1, a2, a3 = Lattice.from_parameters(1, 2, 3, 75, 40, 81).lattice_vectors

    assert round(abs(norm(a1) - 1), 7) == 0
    assert round(abs(norm(a2) - 2), 7) == 0
    assert round(abs(norm(a3) - 3), 7) == 0

    assert round(abs(angle_between(a2, a3) - 75), 7) == 0
    assert round(abs(angle_between(a3, a1) - 40), 7) == 0
    assert round(abs(angle_between(a1, a2) - 81), 7) == 0


def test_lattice_parameters_reciprocal_and_back():
    """Create lattice from parameters, take reciprocal twice,
    and see if the parameters have changed."""
    triclinic = (3, 4, 20, 45, 90, 126)
    triclinic2 = Lattice.from_parameters(
        *triclinic
    ).reciprocal.reciprocal.lattice_parameters
    assert np.allclose(triclinic, triclinic2)


@pytest.mark.parametrize("name", islice(Crystal.builtins, 20))
def test_lattice_parameters_back_and_forth(name):
    """Test that the conversion between lattice vectors and lattice parameters
    is working"""

    c = Crystal.from_database(name)
    lv1 = c.lattice_vectors
    params = c.lattice_parameters
    lv2 = lattice_vectors_from_parameters(*params)

    assert np.allclose(lv1, lv2)


def test_cubic_lattice_system():
    """Test that the Lattice.lattice_system attribute is working properly for cubic lattice"""
    assert Lattice(2 * np.eye(3)).lattice_system == LatticeSystem.cubic


def test_tetragonal_lattice_system():
    """Test that the Lattice.lattice_system attribute is working properly for tetragonal lattice"""
    parameters = (2, 2, 3, 90, 90, 90)
    l = Lattice.from_parameters(*parameters)
    assert l.lattice_system == LatticeSystem.tetragonal


def test_rhombohedral_lattice_system():
    """Test that the Lattice.lattice_system attribute is working properly for rhombohedral lattice"""
    parameters = (1, 1, 1, 87, 87, 87)
    l = Lattice.from_parameters(*parameters)
    assert l.lattice_system == LatticeSystem.rhombohedral


def test_monoclinic_lattice_system():
    """Test that the Lattice.lattice_system attribute is working properly for monoclinic lattice
    including all possible permutations."""
    parameters = (1, 2, 3, 90, 115, 90)
    l = Lattice.from_parameters(*parameters)
    assert l.lattice_system == LatticeSystem.monoclinic

    parameters = (2, 3, 1, 115, 90, 90)
    l = Lattice.from_parameters(*parameters)
    assert l.lattice_system == LatticeSystem.monoclinic

    parameters = (3, 1, 2, 90, 90, 115)
    l = Lattice.from_parameters(*parameters)
    assert l.lattice_system == LatticeSystem.monoclinic


def test_hexagonal_lattice_system():
    """Test that the Lattice.lattice_system attribute is working properly for hexagonal lattice,
    including all possible permutations of lattice parameters."""
    parameters = (2, 2, 3, 90, 90, 120)
    l = Lattice.from_parameters(*parameters)
    assert l.lattice_system == LatticeSystem.hexagonal

    parameters = (3, 2, 2, 120, 90, 90)
    l = Lattice.from_parameters(*parameters)
    assert l.lattice_system == LatticeSystem.hexagonal

    parameters = (2, 3, 2, 90, 120, 90)
    l = Lattice.from_parameters(*parameters)
    assert l.lattice_system == LatticeSystem.hexagonal

    parameters = (2, 2, 2, 90, 120, 90)
    l = Lattice.from_parameters(*parameters)
    assert l.lattice_system == LatticeSystem.hexagonal


def test_triclinic_lattice_system():
    """Test that the Lattice.lattice_system attribute is working properly for triclinic lattice"""
    l = Lattice.from_parameters(1, 2, 3, 75, 40, 81)
    assert l.lattice_system == LatticeSystem.triclinic


def test_lattice_system_graphite():
    """Test that the builtin Crystal for graphite has a hexagonal lattice system"""
    graphite = Crystal.from_database("C")
    assert graphite.lattice_system == LatticeSystem.hexagonal


def test_lattice_system_lead():
    """Test that the builtin Crystal for lead has a cubic lattice system"""
    pb = Crystal.from_database("Pb")
    assert pb.lattice_system == LatticeSystem.cubic


def test_lattice_system_vo2():
    """Test that the builtin Crystal for monoclinic M1 VO2 has a monoclinic lattice system"""
    vo2 = Crystal.from_database("vo2-m1")
    assert vo2.lattice_system == LatticeSystem.monoclinic


def test_scattering_vector_trivial():
    """Test that Lattice.scattering_vectors is working"""
    lattice = Lattice(np.random.random((3, 3)))
    assert np.allclose(lattice.scattering_vector((0, 0, 0)), ((0, 0, 0)))


def test_scattering_vector_table():
    """Test that Lattice.scattering_vector is working on tables of reflections"""
    lattice = Lattice(2 * np.pi * np.eye(3))
    vectors = np.random.random(size=(10, 3))
    assert np.allclose(lattice.scattering_vector(vectors), vectors)


def test_miller_indices_trivial():
    """Test that Lattice.miller_indices is working"""
    lattice = Lattice(np.random.random((3, 3)))
    assert np.allclose(lattice.miller_indices((0, 0, 0)), ((0, 0, 0)))


def test_miller_indices_table():
    """Test that Lattice.miller_indices is working on tables of vectors"""
    lattice = Lattice(2 * np.pi * np.eye(3))
    vectors = np.random.random(size=(10, 3))
    assert np.allclose(lattice.miller_indices(vectors), vectors)


def test_back_and_forth():
    """Test that Lattice.miller_indices and Lattice.scattering_vector are
    reciprocal to each other"""
    lattice = Lattice(np.random.random((3, 3)))
    h, k, l = np.random.randint(-10, 10, size=(3,))
    vector = lattice.scattering_vector((h, k, l))
    hp, kp, lp = lattice.miller_indices(vector)

    assert round(abs(h - float(hp)), 7) == 0
    assert round(abs(k - float(kp)), 7) == 0
    assert round(abs(l - float(lp)), 7) == 0


def test_bounded_reflections_negative():
    """Test that negative reflection bounds raise an Exception.
    Otherwise, an infinite number of reflections will be generated"""
    crystal = Crystal.from_database("C")
    with pytest.raises(ValueError):
        next(crystal.bounded_reflections(-1))


def test_bounded_reflections_zero():
    """Check that bounded_reflections returns (000) for a zero bound"""
    crystal = Crystal.from_database("C")
    reflections = set(crystal.bounded_reflections(0))
    assert (0, 0, 0) in reflections
    assert len(reflections), 1


@pytest.mark.parametrize(
    "crystal", map(Crystal.from_database, islice(Crystal.builtins, 3))
)
def test_bounded_reflections_all_within_bounds(crystal):
    """Check that every reflection is within the bound"""
    bound = 10
    vectors = (
        crystal.scattering_vector(refl) for refl in crystal.bounded_reflections(bound)
    )
    for vector in vectors:
        assert np.linalg.norm(vector) <= bound
