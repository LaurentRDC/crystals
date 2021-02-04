# -*- coding: utf-8 -*-
"""
Crystal structure indexing with the pinkindexer algorithm.
"""
import numpy as np
from . import _pinkindexer
from ..lattice import Lattice


def index_pink(
    peaks,
    intensities,
    detector_distance,
    beam_energy,
    divergence_angle,
    non_monochromaticity,
    detector_radius,
    reciprocal_lattice,
):
    """
    Index reflections using pinkindexer.

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
    reciprocal_lattice : ndarray, shape (3,3)
        Initial guess for the eciprocal lattice [1/A]

    Raises
    ------
    ValueError
        If the number of peaks does not match the intensities provided.
        If the dimensions of the reciprocal lattice are not adequate.
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

    lat = _pinkindexer.index_pink(
        intensities=intensities,
        peaks=peaks,
        detector_distance=detector_distance,
        beam_energy=beam_energy,
        divergence_angle=divergence_angle,
        non_monochromaticity=non_monochromaticity,
        detector_radius=detector_radius,
        reciprocal_lattice=reciprocal_lattice,
    )
    return Lattice(np.asarray(lat))
