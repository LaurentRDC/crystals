.. _indexing:

********************
Indexing reflections
********************

.. currentmodule:: crystals

-----------
Error types
-----------

.. autoclass:: IndexingError
    :show-inheritance:

------------------------------------------------
Single-crystal indexing with the DirAx algorithm
------------------------------------------------

.. autofunction:: index_dirax

---------------------------------------
Polychromatic indexing with pinkIndexer
---------------------------------------

.. autofunction:: index_pink

Extra options
`````````````

The following enumerations may be used to control the indexing performed by `index_pink`. See the `pinkindexer` reference
on what the various options do.

.. autoclass:: crystals.indexing.pinkindexer.ConsideredPeaksCount

.. autoclass:: crystals.indexing.pinkindexer.AngleResolution

.. autoclass:: crystals.indexing.pinkindexer.RefinementType
