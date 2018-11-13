.. include:: references.txt

.. _api:

*************
Reference/API
*************

.. currentmodule:: crystals

You will find the detailed documentation below. Every exposed data structure and function is listed herein.
You can also take a look at the short :ref:`user guide <user_guide>` for some usage examples.

Data Structures
---------------

Structure manipulation is done through the following classes:

.. autosummary::
    :toctree: classes/
    :nosignatures:

    Crystal
    Atom

Bases classes
-------------

The :class:`Lattice` class allows for manipulating lattice information separately from
atomic information.

.. autosummary::
    :toctree: classes/
    :nosignatures:

    Lattice
    LatticeSystem
    AtomicStructure

Utilities
---------

To help with fleshing out unit cell atoms from symmetry operators:

.. autosummary::
    :toctree: functions/
    :nosignatures:

    symmetry_expansion
    lattice_system

.. autosummary::
    :toctree: functions/
    :nosignatures:

    distance_fractional
    distance_cartesian
    
Conversion between ``Crystal`` and other packages

.. autosummary::
    :toctree: functions/
    :nosignatures:

    ase_atoms
    Crystal.from_ase