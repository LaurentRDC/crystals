# -*- coding: utf-8 -*-
"""
Crystal structure indexing with the pinkindexer algorithm.
"""
from typing import Tuple
from warnings import warn

import numpy as np

from . import _pinkindexer
from .common import IndexingError


def index_pink(
    peaks,
    intensities,
    detector_distance,
    beam_energy,
    divergence_angle,
    non_monochromaticity,
    detector_radius,
    tolerance: float,
    reflection_radius: float,
    initial_reciprocal_guess,
    num_threads: int = 6,
) -> Tuple[np.ndarray, int]:
    """
    Index reflections using pinkindexer, and indexing routine that can be
    used in a variety of contexts including measurements made with a
    monochromatic radiation source, a polychromatic source or with radiation
    of very short wavelength.

    .. versionadded:: 1.3.1

    Parameters
    ----------
    peaks : ndarray, shape (N, 2)
        Each row represents a peak location [x, y] on the detector [m]
    intensities : ndarray, shape (N,)
        Scattering intensity for each row in `peaks` [a.u.]
    detector_distance : float
        Distance between the sample and detector [m]
    beam_energy : float
        Radiation energy [eV]
    divergence_angle : float
        Divergence angle [deg]
    non_monochromaticity: float
        I don't know what that is yet.
    detector_radius : float
        Detector radius [m]
    tolerance : float
        I don't know what that is yet.
    reflection_radius : float
        I don't know what that is yet [1/A]
    initial_reciprocal_guess : ndarray, shape (3,3)
        Initial guess of the reciprocal lattice vectors of the crystal
        structure. These vectors must representa] a *primitive* lattice.
    num_threads : int, optional
        Number of threads to use for indexing.

    Returns
    -------
    indexed_lattice : :class:`Lattice`
        Lattice that indexes `peaks` best.
    num_indexed : int
        Number of reflections successfully indexed.

    Raises
    ------
    ValueError
        If the number of peaks does not match the intensities provided.
        If the dimensions of the reciprocal lattice are not adequate.

    References
    ----------
    Y. Gevorkov, A. Barty, W. Brehm, T. A. White, A. Tolstikova,
    M. O. Wiedorn, A. Meents R.-R. Grigat, H. N. Chapman and O. Yefanova,
    pinkIndexer - a universal indexer for pink-beam X-ray and electron
    diffraction snapshots (2020). Acta Cryst. A, vol 76, pages 121-132.
    """
    intensities = np.asfarray(intensities)
    peaks = np.asfarray(peaks)

    if intensities.shape[0] != peaks.shape[0]:
        raise ValueError(
            f"Number of peaks ({peaks.shape[0]}) and intensities ({intensities.shape[0]}) does not match."
        )

    if peaks.shape[1] != 2:
        raise ValueError(
            f"Expected peaks to be an (N,2) iterable, but got {peaks.shape}"
        )

    try:
        recip, num_indexed = _pinkindexer.index_pink(
            intensities=intensities,
            peaks=peaks,
            detector_distance=float(detector_distance),
            beam_energy=float(beam_energy),
            divergence_angle=float(divergence_angle),
            non_monochromaticity=float(non_monochromaticity),
            detector_radius=float(detector_radius),
            tolerance=tolerance,
            reflectionRadius_1_per_A=reflection_radius,
            reciprocal_lattice=initial_reciprocal_guess,
            num_threads=num_threads,
        )
    except _pinkindexer.PinkIndexerError:
        warn("Indexing has failed; returning the initial guess")
        raise IndexingError()

    if num_indexed == 0:
        raise IndexingError("Indexing has failed; no peaks were successfully indexed.")

    return np.asarray(recip), num_indexed
