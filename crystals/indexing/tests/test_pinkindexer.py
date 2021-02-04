
from crystals.indexing._pinkindexer import index_pink
import numpy as np
import pytest

def test_index_pink_trivial():
    """ Test a bogus indexing to see if things don't crash. """
    r = index_pink(1,1,1,1,detector_radius=1, reciprocal_lattice=np.eye(3))
    assert r is None