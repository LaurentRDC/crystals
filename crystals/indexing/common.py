# -*- coding: utf-8 -*-
"""
Basic definitions common to all indexers.
"""
import numpy as np


class IndexingError(RuntimeError):
    """
    Indexing has failed.

    .. versionadded :: 1.3.0
    """

    pass


def ratio_indexed(lattice, reflections):
    """
    Determine the proportion of `reflections` that are indexed by `lattice`, i.e.
    that have integer Miller indices.

    Parameters
    ----------
    lattice : :class:`Lattice` or :class:`Crystal`

    reflections : iterable of 3-tuple or ndarray, shape (N, 3)
        Iterable of reflections with their three-dimensional reciprocal space
        coordinates, or ndarray where each row is a reflections. Coordinates are
        in inverse Angstroms.

    Returns
    -------
    ratio : float
        Ratio of reflections that are indexed correctly, in the [0, 1] range.
    """
    reflections = np.asarray(reflections, dtype=float)

    hkl = lattice.miller_indices(reflections)
    return (
        np.sum(np.all(np.isclose(hkl - np.rint(hkl), 0, atol=0.1), axis=1))
        / reflections.shape[0]
    )


def row_echelon_form(A):
    """
    Transform the matrix `A` in row echelon form.

    .. versionadded :: 1.3.0
    """
    # Algorithm adapted from:
    #   https://johnfoster.pge.utexas.edu/numerical-methods-book/LinearAlgebra_DirectSolvers.html
    A = np.array(A, copy=True)

    for i in range(A.shape[0]):
        pivot = i + np.abs(A[i:, i]).argmax()

        # Swapping rows to make the maximal entry the
        # pivot (if needed).
        if pivot != i:
            A[[pivot, i]] = A[[i, pivot]]

        # Eliminate all entries below the pivot
        factor = A[i + 1 :, i] / A[i, i]
        A[i + 1 :] -= factor[:, None] * A[i]

    return A
