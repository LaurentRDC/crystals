from pathlib import Path

import numpy as np

from crystals import Lattice
from crystals.indexing import index_pink

DATADIR = Path(__file__).parent / "data"


def test_index_pink_desy_data():
    """Test index_pink with data provided by DESY's Yaroslav Gevorkov"""
    # The following corresponds to the pinkindexer tests (tests.cpp),
    # specifically the `testPatternPrediction` test
    # using `getExperimentSettingLysPink`
    intensities = np.loadtxt(DATADIR / "intensities_new")
    peaks = np.loadtxt(DATADIR / "peaksOnDetector_m_new").T
    peaks = np.array(peaks[:, ::-1])

    assert peaks.shape[0] == intensities.shape[0]
    assert peaks.shape[1] == 2

    # The guess reciprocal lattice has been given to me by Y. Gevorkov,
    # but for a better user experience we need to provide the real-space
    # lattice as a guess
    guess = Lattice(
        2
        * np.pi
        * np.diag([0.0126422250316056, 0.0126422250316056, 0.0263157894736842])
    ).reciprocal

    indexed, num_indexed = index_pink(
        peaks=peaks,
        intensities=intensities,
        detector_distance=0.2500,
        beam_energy=8.0010e03,
        divergence_angle=0.001745329251994,
        non_monochromaticity=0.25,
        detector_radius=88.6e-6 * 1300,
        tolerance=0.02,
        reflection_radius=2.528445006321113e-04,
        initial_guess=guess,
    )
    assert num_indexed > 100

    assert np.allclose(
        np.linalg.norm(indexed, axis=0), np.linalg.norm(guess, axis=0), atol=1
    )
