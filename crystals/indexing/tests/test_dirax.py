# -*- coding: utf-8 -*-
import numpy as np

from crystals import Crystal, Lattice, IndexingError, index_dirax
import pytest

np.random.seed(2021)


# The structures to test have been chosen so that it doesn't take too long.
@pytest.mark.parametrize(
    "name,bound", zip(["Pu-epsilon", "C", "vo2-m1", "BaTiO3_cubic"], [2, 3, 2, 2])
)
def test_dirax_indexing_ideal(name, bound):
    """
    Test that indexing always succeeds with an ideal
    set of reflections: no noise, no alien reflections.
    """
    cryst = Crystal.from_database(name)
    refls = [cryst.scattering_vector(r) for r in cryst.bounded_reflections(bound=bound)]
    lat, hkls = index_dirax(refls)

    assert np.allclose(hkls - np.rint(hkls), 0, atol=0.001)


@pytest.mark.parametrize("name", ["Pu-epsilon", "C", "vo2-m1", "BaTiO3_cubic"])
def test_dirax_indexing_initial_guess(name):
    """
    Test that indexing succeeds with an initial guess and very few reflections.
    """
    cryst = Crystal.from_database(name)
    # We restrict the number of reflections to a single (0,0,0); with
    # the initial guess, indexing should still succeed!
    lat, _ = index_dirax([(0, 0, 0)], initial=cryst)

    assert np.allclose(lat.lattice_parameters, cryst.lattice_parameters, atol=1)


@pytest.mark.parametrize(
    "name,bound", zip(["Pu-epsilon", "C", "vo2-m1", "BaTiO3_cubic"], [2, 3, 2, 2])
)
def test_dirax_indexing_alien_reflections(name, bound):
    """
    Test that indexing always succeeds even with 20% of alien reflection.
    """
    cryst = Crystal.from_database(name)
    refls = [cryst.scattering_vector(r) for r in cryst.bounded_reflections(bound=bound)]
    num_aliens = len(refls) // 5
    aliens = [
        cryst.lattice_parameters[0] * np.random.random(size=(3,))
        for _ in range(num_aliens)
    ]
    lat, hkls = index_dirax(refls + aliens)
    # The alien reflections will not be indexed correctly, of course
    hkls = hkls[:num_aliens]
    assert np.allclose(hkls - np.rint(hkls), 0, atol=0.01)


@pytest.mark.parametrize(
    "name,bound", zip(["Pu-epsilon", "C", "vo2-m1", "BaTiO3_cubic"], [2, 3, 2, 2])
)
def test_dirax_indexing_noise(name, bound):
    """
    Test that indexing always succeeds despite noise in reflection positions.
    """
    np.random.seed(0)

    cryst = Crystal.from_database(name)
    refls = [
        cryst.scattering_vector(r) + np.random.normal(0, scale=0.01, size=(3,))
        for r in cryst.bounded_reflections(bound=bound)
    ]
    lat, hkls = index_dirax(refls)
    assert np.allclose(hkls - np.rint(hkls), 0, atol=0.1)


def test_dirax_indexing_length_bound():
    """
    Test that indexing fails as expected if the lattice length bounds are too restrictive.
    """
    cryst = Crystal.from_database("C")
    refls = [cryst.scattering_vector(r) for r in cryst.bounded_reflections(bound=2)]
    # We restrict the number of reflections to a single (0,0,0); with
    # the initial guess, indexing should still succeed!
    with pytest.raises(IndexingError):
        index_dirax(
            refls,
            length_bounds=(0.01, 2),
        )
