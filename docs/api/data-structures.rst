.. _data_structures:

***********************
Primary Data Structures
***********************

.. currentmodule:: crystals

Here are the primary data structures available in :mod:`crystals`. You will work directly with these structures, 
maybe even either creating them directly.

-------
Crystal
-------

.. autoclass:: Crystal
    :show-inheritance:
    :exclude-members: from_parameters

---------
Supercell
---------

.. autoclass:: Supercell
    :show-inheritance:
    :exclude-members: from_parameters

----
Atom
----

To deal with atoms with coordinates, take a look at the :class:`Atom` class:

.. autoclass:: Atom
    :show-inheritance:

-------
Element
-------

To access atomic data only, take a look at the :class:`Element` class:

.. autoclass:: Element
    :show-inheritance: