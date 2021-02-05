# -*- coding: utf-8 -*-
"""
Crystal structure indexing with the pinkindexer algorithm.
"""
import numpy as np
from . import _pinkindexer
from .common import IndexingError
from ..lattice import Lattice
from ..crystal import Crystal


def index_pink(
    peaks,
    intensities,
    detector_distance,
    beam_energy,
    divergence_angle,
    non_monochromaticity,
    detector_radius,
    initial,
):
    """
    Index reflections using pinkindexer, and indexing routine that can be
    used in a variety of contexts including measurements made with a
    monochromatic radiation source, a polychromatic source or with radiation
    of very short wavelength.

    Parameters
    ----------
    peaks : ndarray, shape (N, 2)
        Peak locations on detector [m]
    intensities : ndarray, shape (N,)
        Scatterign intensity for each peak in `peaks` [a.u.]
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
    initial : Lattice or Crystal
        Initial guess for the crystal structure. If a :class:`Lattice`
        instance is provided, it is assumed that it is the **primitive**
        lattice. If a :class:`Crystal` instances is provided,
        the primitive lattice will be calculated.

    Returns
    -------
    indexed_lattice : :class:`Lattice`
        Lattice that indexes `peaks` best.

    Raises
    ------
    ValueError
        If the number of peaks does not match the intensities provided.
        If the dimensions of the reciprocal lattice are not adequate.

    References
    ----------
    Y. Gevorkov, A. Barty, W. Brehm, T. A. White, A. Tolstikova,
    M. O. Wiedorn, A. Meents R.-R. Grigat, H. N. Chapman and O. Yefanova,
    pinkIndexer – a universal indexer for pink-beam X-ray and electron
    diffraction snapshots (2020). Acta Cryst. A, vol 76, pages 121-132.
    """
    intensities = np.asarray(intensities)
    peaks = np.asarray(peaks)

    if intensities.shape[0] != peaks.shape[0]:
        raise ValueError(
            f"Number of peaks ({peaks.shape[0]}) and intensities ({intensities.shape[0]}) does not match."
        )

    if peaks.shape[1] != 2:
        raise ValueError(
            f"Expected peaks to be an (N,2) iterable, but got {peaks.shape}"
        )

    # Yaroslav Gevorkov has mentionned that the indexing routine
    # only works provided an initial lattice that is PRIMITIVE
    if isinstance(initial, Crystal):
        initial = initial.primitive()

    try:
        lat = _pinkindexer.index_pink(
            intensities=intensities,
            peaks=peaks,
            detector_distance=detector_distance,
            beam_energy=beam_energy,
            divergence_angle=divergence_angle,
            non_monochromaticity=non_monochromaticity,
            detector_radius=detector_radius,
            reciprocal_lattice=np.array(initial.reciprocal_vectors),
        )
    except _pinkindexer.PinkIndexerError:
        return initial
    return Lattice(np.asarray(lat))
