from crystals.indexing.pinkindexer import index_pink
from crystals import Lattice
import numpy as np
import pytest


def test_index_pink_trivial():
    """ Test a bogus indexing to see if things don't crash. """
    intensities = np.array([0])
    peaks = np.array([[0, 0]])
    r = index_pink(
        peaks, intensities, 1, 1, 1, 1, detector_radius=1, initial=Lattice(np.eye(3))
    )
    assert type(r) == Lattice
