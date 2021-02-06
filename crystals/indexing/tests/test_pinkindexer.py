from crystals.indexing.pinkindexer import index_pink
import crystals.indexing._pinkindexer as _pinkindexer
from crystals import Lattice
import numpy as np
import pytest
from pathlib import Path

DATADIR = Path(__file__).parent / "data"


def test_index_pink_desy_data():
    """ Test index_pink with data provided by DESY's Yaroslav Gevorkov """
    intensities = np.loadtxt(DATADIR / "intensities_new")
    peaks = np.loadtxt(DATADIR / "peaksOnDetector_m_new").T
    peaks = np.array(peaks[:, ::-1])

    assert peaks.shape[0] == intensities.shape[0]
    assert peaks.shape[1] == 2

    basis = np.diag([0.0126422250316056, 0.0126422250316056, 0.0263157894736842])

    indexed, num_indexed = _pinkindexer.index_pink(
        peaks=peaks,
        intensities=intensities,
        detector_distance=float(0.2500),
        beam_energy=float(8.0010e03),
        divergence_angle=float(0.1 * np.pi / 180),
        non_monochromaticity=float(0.25),
        detector_radius=float(88.6e-6 * 1300),
        reciprocal_lattice=basis,
    )

    # TODO: this is a very low bar to set...
    assert num_indexed > 0
