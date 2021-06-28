# -*- coding: utf-8 -*-
"""
Crystal structure indexing with the DirAx algorithm.
"""
from itertools import product
import numpy as np
from ..lattice import Lattice
from .common import IndexingError, row_echelon_form


def index_dirax(reflections, initial=None, length_bounds=(2, 20)):
    """
    Find the lattice associated with a list of
    reflections using the DirAx algorithm.

    .. versionadded:: 1.3.0

    Parameters
    ----------
    reflections : iterable of 3-tuple or ndarray, shape (N, 3)
        Iterable of reflections with their three-dimensional reciprocal space
        coordinates, or ndarray where each row is a reflections. Coordinates are
        in inverse Angstroms.
    initial : :class:`Lattice` or :class:`Crystal`, optional
        Initial guess for a lattice. The DirAx algorithm does not need an initial guess, but
        it can certainly help when many reflections are missing.
    length_bounds : 2-tuple of floats, optional
        Minimum and maximum lattice vector lengths to consider, in Angstrom.

    Returns
    -------
    indexed : :class:`Lattice`
        Lattice that best indexes `reflections`.
    indices : ndarray, shape (N, 3)
        Miller indices associated with each vector in `reflections`, in order.

    Raises
    ------
    IndexingError
        The indexing has failed.

    Notes
    -----
    This indexing routine is based on the reference below. The algorithm is well-suited to
    situations where reflections might be missing, or situations where the list of reflections
    contains "alien" reflections not associated with the lattice.

    Examples
    --------
    We generate reflections from a crystal structure and re-index it for demonstration
    purposes.

    >>> from crystals import index_dirax, Crystal
    >>> import numpy as np
    >>> graphite = Crystal.from_database('C')
    >>> # The list of reflections `qs` might be experimental measurements from either
    >>> # x-ray or electron diffraction.
    >>> qs = [graphite.scattering_vector(r) for r in graphite.bounded_reflections(bound=3.5)]
    >>> lattice, hkls = index_dirax(qs)
    >>> lattice # doctest: +SKIP
    < Lattice object with parameters 2.464Å, 2.464Å, 6.711Å, 90.00°, 90.00°, 120.00° >

    References
    ----------
    A. J. M. Duisenberg, Indexing in Single-Crystal Diffractometry with an Obstinate List of Reflections (1992),
    J. Appl. Cryst vol. 25 pp. 92 - 96
    """
    reflections = np.asfarray(reflections)

    length_min, length_max = sorted(length_bounds)
    d_max = 2 * np.pi / length_min
    d_min = 2 * np.pi / length_max

    ns = np.arange(start=1, stop=13, step=1).reshape((-1, 1))

    potential_direct_vectors = set()

    # If a "guess" lattice is provided, we include its
    # lattice vectors as being high-priority.
    if initial is not None:
        for v in initial.lattice_vectors:
            potential_direct_vectors.add(LatVec(nf=len(reflections), vector=v))

    points = [np.squeeze(a) for a in np.vsplit(reflections, reflections.shape[0])]
    for a1, a2, a3 in product(points, repeat=3):
        normal = np.cross(a2 - a1, a3 - a1)
        if np.allclose(normal, 0, atol=1e-2):
            continue
        normal /= np.linalg.norm(normal)

        proj = np.sort(np.sum(reflections * normal[None, :], axis=1))
        d_star_start = np.diff(proj).max()
        if d_star_start <= 0:
            continue

        frac = proj / d_star_start / ns
        residues = np.sum(np.abs(frac - np.rint(frac)), axis=1)
        best_n = np.argmin(residues) + 1

        d_star = d_star_start / best_n
        if (d_star < d_min) or (d_star > d_max):
            continue

        frac_dist = proj / d_star
        nf = np.sum(np.isclose(frac_dist, np.rint(frac_dist), atol=1 / 24))

        t = 2 * np.pi * normal / d_star
        potential_direct_vectors.add(LatVec(nf, t))

    if len(potential_direct_vectors) == 0:
        raise IndexingError("No candidate lattice vectors could be determined.")

    max_nf = max((v.nf for v in potential_direct_vectors))
    acceptance_level = 0.8 * max_nf
    potential_direct_vectors = list(
        filter(lambda v: v.nf >= acceptance_level, potential_direct_vectors)
    )

    lattice = _find_basis(
        (v.vector for v in potential_direct_vectors), reflections=reflections
    )
    return lattice, lattice.miller_indices(reflections)


class LatVec:
    """
    Potential lattice vector, including the number indexed reflections `nf`
    """

    def __init__(self, nf, vector):
        self.nf = nf
        self.vector = np.around(vector, decimals=3)
        if np.sum(vector) < 0:
            self.vector *= -1

    def __hash__(self):
        return hash((self.nf, *tuple(self.vector)))

    def __repr__(self):
        return f"< LatVec: nf={self.nf}, vector={self.vector} >"


def _find_basis(vectors, reflections):
    """Find the shorted three linearly-independent vectors from a list."""
    vectors = sorted(vectors, key=np.linalg.norm)
    a1 = vectors.pop(0)
    try:
        index, a2 = next(
            (i, v)
            for i, v in enumerate(vectors)
            if np.linalg.norm(np.cross(a1, v)) >= 1
        )
    except StopIteration:
        raise IndexingError(
            "No set of three linearly-independent lattice vectors were found in the candidates."
        )
    vectors.pop(index)

    # The best way to find the last vector is to iterate through the remaining
    # vectors, from shortest to longest, until a matrix with the candidate
    # basis as rows has rank 3
    m = np.empty(shape=(3, 3), dtype=float)
    m[0, :] = a1
    m[1, :] = a2
    for v in vectors:
        m[2, :] = v
        if np.linalg.matrix_rank(m, tol=0.1) == 3:
            return Lattice(row_echelon_form([a1, a2, v]))

    raise IndexingError(
        "No set of three linearly-independent lattice vectors were found in the candidates."
    )
