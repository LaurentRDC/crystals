.. include:: references.txt

.. _api:

*************
Reference/API
*************

.. currentmodule:: crystals

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

Parsers
-------

Structure parsers are used to build :class:`Crystal` instances, mostly through :class:`Crystal` class methods.

.. autosummary::
    :toctree: classes/
    :nosignatures:

    CIFParser
    CODParser
    PDBParser

Affine Transforms
-----------------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    affine_map
    transform
    change_of_basis
    change_basis_mesh
    is_basis
    is_rotation_matrix
    minimum_image_distance
    rotation_matrix
    translation_matrix
    translation_rotation_matrix