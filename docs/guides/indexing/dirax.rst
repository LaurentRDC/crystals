.. currentmodule:: crystals

General purpose direct indexing with DirAx
==========================================

The simplest indexing routine, :func:`index_dirax`, is based on the DirAx algorithm presented in:

`A. J. M. Duisenberg, Indexing in Single-Crystal Diffractometry with an Obstinate List of Reflections (1992),
J. Appl. Cryst vol. 25 pp. 92 - 96. <https://doi.org/10.1107/S0021889891010634>`_

For this example, we will need a list of reflections oriented in reciprocal space. Let's start with a very simple example using the simplest crystal structure, cubic polonium:

    >>> from crystals import Crystal
    >>> import numpy as np
    >>>
    >>> indices = [ (0,0,0), (1,0,0), (0,1,0), (0,0,1) ]
    >>> 
    >>> cryst = Crystal.from_database('Pu-epsilon') # cubic polonium
    >>> peaks = [cryst.scattering_vector(hkl) for hkl in indices]
    
The list of ``peaks`` are peak positions in three-dimensional reciprocal space. To index:

.. testsetup::

    from crystals import Crystal
    import numpy as np
    indices = [ (0,0,0), (1,0,0), (0,1,0), (0,0,1) ]
    cryst = Crystal.from_database('Pu-epsilon') # cubic polonium
    peaks = [cryst.scattering_vector(hkl) for hkl in indices]

.. doctest:: 

    >>> from crystals import index_dirax
    >>> lattice, hkls = index_dirax(peaks)
    >>> lattice
    < Lattice object with parameters 3.638Å, 3.638Å, 3.638Å, 90.00°, 90.00°, 90.00° >
    >>> hkls.astype(int)
    array([[0, 0, 0],
           [1, 0, 0],
           [0, 1, 0],
           [0, 0, 1]])

Here, ``hkls`` are the indices of the reflections in ``peaks`` according to the indexing lattice, ``lattice``. 
As expected, the Miller indices in ``hkls`` match the ones we created in ``indices``.

Indexing a cubic lattice with perfect peak placement is pretty easy: we only needed four reflections. 
The DirAx algorithm is robust against missing reflections, but also *alien* reflections that do not belong there.
In the next example, we will introduce 10% of non-fitting reflections:


.. testsetup::

    import numpy as np
    np.random.seed(0) # ensure reproducible example

.. doctest::

    >>> from crystals import Crystal, index_dirax
    >>> import numpy as np
    >>>
    >>> cryst = Crystal.from_database('Pu-epsilon') # cubic polonium
    >>> indices = list(cryst.bounded_reflections(3))
    >>> num_aliens = len(indices) // 10 # 10% alien reflections
    >>> aliens  = [ np.random.random(size=(3,)) for _ in range(num_aliens) ]
    >>> 
    >>> peaks = [cryst.scattering_vector(hkl) for hkl in indices + aliens] 
    >>> lattice, hkls = index_dirax(peaks)
    >>> lattice
    < Lattice object with parameters 3.638Å, 3.638Å, 3.638Å, 90.00°, 90.00°, 90.00° >
    >>> hkls.round(decimals=3)
    array([[-1.   , -1.   , -1.   ],
           [-1.   , -1.   ,  0.   ],
           [-1.   , -1.   ,  1.   ],
           [-1.   ,  0.   , -1.   ],
           [-1.   ,  0.   ,  0.   ],
           [-1.   , -0.   ,  1.   ],
           [-1.   ,  1.   , -1.   ],
           [-1.   ,  1.   ,  0.   ],
           [-1.   ,  1.   ,  1.   ],
           [ 0.   , -1.   , -1.   ],
           [ 0.   , -1.   ,  0.   ],
           [ 0.   , -1.   ,  1.   ],
           [ 0.   ,  0.   , -1.   ],
           [ 0.   ,  0.   ,  0.   ],
           [-0.   , -0.   ,  1.   ],
           [-0.   ,  1.   , -1.   ],
           [-0.   ,  1.   ,  0.   ],
           [-0.   ,  1.   ,  1.   ],
           [ 1.   , -1.   , -1.   ],
           [ 1.   , -1.   ,  0.   ],
           [ 1.   , -1.   ,  1.   ],
           [ 1.   ,  0.   , -1.   ],
           [ 1.   ,  0.   ,  0.   ],
           [ 1.   , -0.   ,  1.   ],
           [ 1.   ,  1.   , -1.   ],
           [ 1.   ,  1.   ,  0.   ],
           [ 1.   ,  1.   ,  1.   ],
           [ 0.549,  0.715,  0.603],
           [ 0.545,  0.424,  0.646]])

As you can see, the last indexed reflections in ``hkls`` don't have integer 
Miller indices; those are the non-fitting (alien) reflections!

Of course, the DirAx algorithm is also resistant to noise. Here is an example where
we add noise to our perfect peaks:

.. testsetup::

    import numpy as np
    np.random.seed(0)

.. doctest::

    >>> from crystals import Crystal, Lattice, index_dirax
    >>> import numpy as np
    >>> 
    >>> cryst = Crystal.from_database("Pu-epsilon")
    >>> indices = list(cryst.bounded_reflections(2))
    >>> 
    >>> peaks = [
    ...     cryst.scattering_vector(hkl) + np.random.normal(0, scale=0.01, size=(3,))
    ...     for hkl in indices
    ... ]
    >>> lattice, hkls = index_dirax(peaks)
    >>> hkls.round(decimals=2)
    array([[-0.98, -0.  , -0.  ],
           [ 0.  , -1.  ,  0.01],
           [ 0.  , -0.  , -1.  ],
           [ 0.  ,  0.01, -0.  ],
           [ 0.  , -0.  ,  0.99],
           [ 0.  ,  1.01, -0.  ],
           [ 1.  , -0.01,  0.  ]])

Dealing with a restricted set of reflections
--------------------------------------------

:func:`index_dirax` can also deal with a restricted set of reflections. For example, in 
electron diffraction measurements, it might occur that the only reflections measured are 
of the form :math:`\{ (hk0) \}`. This happens, for example, when measuring diffraction from layered compounds
like graphite.

Let's try to index **without** prior knowledge first:

    >>> from crystals import Crystal, index_dirax
    >>> import numpy as np
    >>> from itertools import product
    >>>
    >>> cryst = Crystal.from_database('C') # graphite
    >>> indices = product(
    ...     [-2,-1,0,1,2], #-3 < h < 3
    ...     [-2,-1,0,1,2], #-3 < k < 3
    ...     [0]            #     l = 0
    ... )
    >>> peaks = [cryst.scattering_vector(hkl) for hkl in indices] 
    >>> lattice, hkls = index_dirax(peaks)
    Traceback (most recent call last):
    crystals.indexing.common.IndexingError: No candidate lattice vectors could be determined.

As you can see, there not enough information to index. Let's now use an initial guess:

.. testsetup::

    from crystals import Crystal, index_dirax
    from itertools import product
    cryst = Crystal.from_database('C') # graphite
    indices = product([-2,-1,0,1,2],[-2,-1,0,1,2],[0])
    peaks = [cryst.scattering_vector(hkl) for hkl in indices] 

.. doctest:: 

    >>> lattice, hkls = index_dirax(peaks, initial=cryst) # Lattice or Crystal works here
    >>> lattice
    < Lattice object with parameters 2.464Å, 2.464Å, 6.711Å, 90.00°, 90.00°, 120.00° >

Great!

Other indexing methods
----------------------

More specialized and performant indexing methods are planned to be included in :mod:`crystals`. 
If there's a particular method you would like to see included, please do not hesitate to 
`raise an issue <https://github.com/LaurentRDC/crystals/issues/new>`_!