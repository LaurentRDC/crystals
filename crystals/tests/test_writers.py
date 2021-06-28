# -*- coding: utf-8 -*-

import tempfile
from itertools import islice
from pathlib import Path
from random import choice, randint

import numpy as np
import pytest
from crystals import Crystal

try:
    import ase
except ImportError:
    WITH_ASE = False
else:
    WITH_ASE = True


@pytest.mark.skipif(not WITH_ASE, reason="ASE not installed or importable")
@pytest.mark.parametrize("name", islice(Crystal.builtins, 10))
def test_ase_atoms_construction(name):
    """Test that ase_atoms returns without error"""
    crystal = Crystal.from_database(name)
    to_ase = crystal.to_ase()
    assert len(crystal) == len(to_ase)


@pytest.mark.skipif(not WITH_ASE, reason="ASE not installed or importable")
@pytest.mark.parametrize("name", islice(Crystal.builtins, 10))
def test_ase_atoms_back_and_forth(name):
    """Test conversion to and from ase Atoms"""
    crystal = Crystal.from_database(name)
    to_ase = crystal.to_ase()
    crystal2 = Crystal.from_ase(to_ase)

    # ase has different handling of coordinates which can lead to
    # rounding beyond 1e-3. Therefore, we cannot compare directly sets
    # assertSetEqual(set(crystal), set(crystal2))
    assert len(crystal) == len(crystal2)


@pytest.mark.parametrize("name", Crystal.builtins)
def test_cif_writer_idempotence(name):
    """Test that conversion to CIF of a structure loaded from CIF is idempotent."""
    # Testing on all built-in structure assures us that corner cases
    # are taken care of.
    cryst = Crystal.from_database(name)
    with tempfile.TemporaryDirectory() as temp_dir:
        f = Path(temp_dir) / "temp.cif"
        cryst.to_cif(f)
        cryst2 = Crystal.from_cif(f)
        assert cryst == cryst2


@pytest.mark.parametrize("name", Crystal.builtins)
def test_vasp_writer_idempotence(name):
    """Test that conversion to VASP of a structure loaded from VASP is idempotent."""
    # Testing on all built-in structure assures us that corner cases
    # are taken care of.
    cryst = Crystal.from_database(name)
    with tempfile.TemporaryDirectory() as temp_dir:
        f = Path(temp_dir) / "temp.POSCAR"
        cryst.to_poscar(f)
        cryst2 = Crystal.from_poscar(f)
        assert cryst == cryst2
