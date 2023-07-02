# -*- coding: utf-8 -*-
"""
Crystal structure indexing with the pinkindexer algorithm.
"""
import faulthandler
from dataclasses import dataclass
from enum import IntEnum, unique
from typing import Tuple
from warnings import warn

import numpy as np

from ..lattice import Lattice
from . import _pinkindexer
from .common import IndexingError


@unique
class ConsideredPeaksCount(IntEnum):
    """This enumeration specifies the number of found Bragg spots that are used to
    compute the initial indexing in `index_pink` solution from the maximum of the rotogram.
    Regardless of this enumeration, all Bragg spots are considered during refinement.
    """

    very_few = 0
    few = 1
    standard = 2
    many = 3
    manyMany = 4


@unique
class AngleResolution(IntEnum):
    """This enumeration controls the resolution of the rotogram in terms of number of
    voxels spanning :math:`-\\arctan \\pi/4` to :math:`\\arctan \\pi/4`

    Choosing larger voxels (lower resolution) leads to a faster calculation but
    lower precision in the initial step of determining the orientation from the rotogram.
    """

    extremely_loose = 0
    loose = 1
    standard = 2
    dense = 3
    extremely_dense = 4


@unique
class RefinementType(IntEnum):
    """Refinement can be performed by a gradient descent method, fitting
    all parameters of the lattice or keeping the cell parameters constant
    and just refining the orientation.
    """

    none = 0
    fixedLattice_parameters = 1
    variable_lattice_parameters = 2
    first_fixed_then_variable_lattice_parameters = 3
    first_fixed_then_variable_lattice_parameters_multi_seed_lengths = 4
    first_fixed_then_variable_lattice_parameters_multi_seed = 5
    first_fixed_then_variable_lattice_parameters_center_adjustment_multi_seed = 6


@dataclass(frozen=True)
class Geometry:
    """Parameters from a geometry file (*.geom)."""

    camera_length: float
    camera_offset: float
    resolution: float
    detector_half_width: int

    @property
    def pixel_length(self) -> float:
        """Length of a pixel in meters."""
        return 1 / self.resolution

    @property
    def detector_distance(self) -> float:
        """Returns the detector distance in meters."""
        return self.camera_length + self.camera_offset

    @property
    def detector_radius(self) -> float:
        """Returns the detector radius (or half-width) in meters."""
        return self.detector_half_width * self.pixel_length


def index_pink(
    peaks,
    intensities,
    beam_energy,
    divergence_angle,
    non_monochromaticity,
    detector_distance: float,
    detector_radius: float,
    tolerance: float,
    reflection_radius: float,
    initial_guess: Lattice,
    considered_peaks_count: ConsideredPeaksCount = ConsideredPeaksCount.standard,
    angle_resolution: AngleResolution = AngleResolution.standard,
    refinement_type: RefinementType = RefinementType.first_fixed_then_variable_lattice_parameters,
    num_threads: int = 1,
) -> Tuple[Lattice, int]:
    """
    Index reflections using `pinkindexer`, an indexing routine that can be
    used in a variety of contexts including measurements made with a
    monochromatic radiation source, a polychromatic source or with radiation
    of very short wavelength.

    Some of the inputs below may not be known. In case, you may use the `Geometry` class
    to derive some of these inputs. See the user guide for an example.

    .. versionadded:: 1.7.0

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
        Bandwidth of the beam energy, as a fraction of the beam energy [a.u.].
    detector_radius : float
        Detector radius, or half-width. [m]
    tolerance : float
        Fractional tolerance for reflections to be considered indexed during refinement [a.u.].
    reflection_radius : float
        Radius of a typical reflection [1/A]
    initial_guess : Lattice
        Initial guess of the real-space lattice vectors of the crystal
        structure. This must represent a **primitive** lattice.
    considered_peaks_count: :class:`crystals.indexing.pinkindexer.ConsideredPeaksCount`, optional
        Controls the number of Bragg spots which are used to compute the indexing solution.
    angle_resolution: :class:`crystals.indexing.pinkindexer.AngleResolution`, optional
        Set the resolution of the rotogram in terms of number of voxels
        spanning :math:`-\\arctan \\pi/4` to :math:`\\arctan \\pi/4`
    refinement_type: :class:`crystals.indexing.pinkindexer.RefinementType`, optional
        Determines which type of refinement to perform after initial indexing.
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

    # The unit convention for reciprocal lattice vectors in `crystals`
    # is that reciprocal lattice vectors have dimensions (rad/angstrom)
    # while pinkindexer expects (1/angstrom)
    guess = np.asarray(initial_guess.reciprocal) / (2 * np.pi)

    # We enable `faulthandler` in order to
    # catch segmentation faults in pinkindexer as exceptions.
    faulthandler.enable()
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
            reciprocal_lattice=guess,
            considered_peaks_count=int(considered_peaks_count),
            angle_resolution=int(angle_resolution),
            refinement_type=int(refinement_type),
            num_threads=num_threads,
        )
    except _pinkindexer.PinkIndexerError:
        warn("Indexing has failed; returning the initial guess")
        raise IndexingError()
    finally:
        faulthandler.disable()

    if num_indexed == 0:
        raise IndexingError("Indexing has failed; no peaks were successfully indexed.")

    return Lattice(2 * np.pi * np.asarray(recip)).reciprocal, num_indexed
