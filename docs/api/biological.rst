.. _biological:

**************************
Biological data structures
**************************

.. currentmodule:: crystals

For crystals created from the Protein Data Bank, the atoms will be stored in substructures. These substructures are described below.

===========================
Primary Structure: Residues
===========================

.. autoclass:: Residue
    :show-inheritance:

===================
Secondary Structure
===================

Secondary structures are really all the same, with different names:

.. currentmodule:: crystals.base

.. autoclass:: SecondaryStructure
    :show-inheritance:

.. currentmodule:: crystals

-----
Helix
-----

.. autoclass:: Helix
    :show-inheritance:

-----
Sheet
-----

.. autoclass:: Sheet
    :show-inheritance: