.. currentmodule:: crystals

Polychromatic indexing with pinkindexer
=======================================

.. versionadded:: 1.7.0

`pinkindexer` is an indexing and refinement algorithm presented in:

`Y. Gevorkov, A. Barty, W. Brehm, T. A. White, A. Tolstikova, M. O. Wiedorn, A. Meents R.-R. Grigat, 
H. N. Chapman and O. Yefanova, pinkIndexer - a universal indexer for pink-beam X-ray and electron 
diffraction snapshots (2020). Acta Cryst. A, vol 76, pages 121-132. <https://doi.org/10.1107/S2053273319015559>`_

As described in the attached reference

    (pinkindexer) can be used in a variety of contexts including measurements made with a monochromatic radiation source, 
    a polychromatic source or with radiation of very short wavelength. As such, the algorithm is particularly suited to 
    automated data processing for two emerging measurement techniques for macromolecular structure determination: serial pink-beam 
    X-ray crystallography and serial electron crystallography

Let's work through a real-world example, with data provided by Vukica Srajer from the University of Chicago.

We consider X-ray diffraction data, where diffraction peaks have already been extracted in some way. The data looks like this:

.. code-block::

    fs/px   ss/px (1/d)/nm^-1    Intensity
    2111.91  620.29       4.04       82.48
    2249.23  636.26       4.05       65.91
    2241.82  754.15       3.74       99.43
    2226.25  787.50       3.65       65.06
    1600.35  858.21       3.57      171.23
    2098.81  877.50       3.36      138.99
    (...)    (...)        (...)     (...)

This data was acquired at a beamline which provides a geometry file. We can specify the geometry using the :class:`crystals.indexing.Geometry` class:

.. doctest::

    >>> from crystals.indexing import Geometry
    >>> geom = Geometry(
    ...     camera_length=0.2511,          # Camera length in meters
    ...     camera_offset=0,               # Camera offset in meters
    ...     detector_half_width=1972,      # Half of the camera width in pixels
    ...     resolution=11287.47795414462,  # Detector resolution in px/m
    ... )

Based on this geometry, we extract the peak locations and intensities::

.. doctest::

    >>> import numpy as np
    >>> data = np.loadtxt("docs/guides/indexing/data/pinkindexer-example.txt", skiprows=1, usecols=[0, 1, 3])
    >>> peaklocs = data[:, 0:2]
    >>> # Centering peaks and conversion from pixels to meters
    >>> peaklocs[:, 0] -= 1987.061423
    >>> peaklocs[:, 1] -= 1972.143086
    >>> peaklocs /= geom.resolution
    >>> intensities = data[:, 2]
    >>> print(peaklocs[0:5, :])
    [[ 0.0110608  -0.11976573]
     [ 0.0232265  -0.11835089]
     [ 0.02257002 -0.10790657]
     [ 0.02119061 -0.10495197]
     [-0.03426022 -0.09868751]]

We need to provide an initial guess for the real-space lattice being indexed. In this case, we expect a hexagonal lattice::

    >>> from crystals import lattice
    >>> guess = Lattice.from_parameters(a=76.03, b=76.03, c=76.14, alpha=90.00, beta=90.00, gamma=120.00)

Finally, we're ready to index using :func:`index_pink`::

    >>> from crystals.indexing import index_pink
    >>> indexed, num_indexed = index_pink(
    ...     peaks=peaklocs,
    ...     intensities=intensities,
    ...     detector_distance=geom.detector_distance,
    ...     beam_energy=11300.0,                        # Beam energy in eV
    ...     divergence_angle=0.0,                       # Beam divergence angle in degrees
    ...     non_monochromaticity=4e-2,                  # Non-monochromaticity, or energy bandwidth, in fraction of the beam energy
    ...     detector_radius=geom.detector_radius,       
    ...     tolerance=0.03,                             # Fractional tolerance to consider a peak being successfully indexed as part of the refinement procedure
    ...     reflection_radius=2.66e-10,                 # Average radius of the peaks in meters
    ...     initial_guess=guess,
    ... )

The first output `indexed` is a real-space lattice, while the second output `num_indexed` is the number of reflections which are being indexed (up to the provided tolerance).
In case no reflections are indexed, an :class:`IndexingError` exception is raised; therefore, if no exception is thrown, we always expect that `num_indexed` > 0.

:func:`index_pink` accepts additional parameters which are used to control both the indexing and refinement steps. These options include:

* :class:`crystals.indexing.pinkindexer.ConsideredPeaksCount`
* :class:`crystals.indexing.pinkindexer.AngleResolution`
* :class:`crystals.indexing.pinkindexer.RefinementType`
