
*********************************************
crystals : data structure for crystallography
*********************************************

.. image:: https://readthedocs.org/projects/crystals/badge/?version=master
    :target: https://readthedocs.org/projects/crystals
    :alt: Documentation Status

.. image:: https://ci.appveyor.com/api/projects/status/github/LaurentRDC/crystals?branch=master&svg=true
    :alt: AppVeyor Build Status
    :target: https://ci.appveyor.com/project/LaurentRDC/crystals

.. image:: https://img.shields.io/pypi/v/crystals.svg
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/crystals

.. image:: https://img.shields.io/pypi/pyversions/crystals.svg
    :alt: Supported Python versions
    :target: https://pypi.python.org/pypi/crystals

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :alt: Code formatting style
    :target: https://github.com/ambv/black

``crystals`` is a library of data structure and algorithms to manipulate abstract crystals. ``crystals`` helps with reading crystallographic 
files (like .cif and .pdb), provides access to atomic positions, and allows for space-group determination. Although ``crystals`` can be used on its own,
it was made to be integrated into larger projects (like `scikit-ued <https://github.com/LaurentRDC/scikit-ued>`_).

Table of content
================

.. toctree::
    :maxdepth: 2
    
    installation
    guide
    reference
    contributing

Usage example
=============

``crystals`` is all about constructing crystals and getting information about the resulting object. Crystals can be built from a variety of sources:

* From files on disk, such as Crystallography Information Files (CIF);
* From the internal database of over 90 structure files (mostly elemental crystals);
* From online databases, such as the `RCSB Protein DataBank <http://www.rcsb.org/>`_ or the 
  `Crystallography Open Database <http://www.crystallography.net/cod/>`_.

Here's a quick example of building a crystal from the internal database::

    >>> from crystals import Crystal
    >>>
    >>> vo2 = Crystal.from_database('vo2-m1')
    >>> print(vo2)	   # Short string representation
    < Crystal object with following unit cell:
        Atom O  @ (0.90, 0.79, 0.80)
        Atom O  @ (0.90, 0.71, 0.30)
        Atom O  @ (0.61, 0.31, 0.71)
        Atom O  @ (0.39, 0.69, 0.29)
        Atom O  @ (0.61, 0.19, 0.21)
        Atom O  @ (0.10, 0.29, 0.70)
        Atom O  @ (0.10, 0.21, 0.20)
        Atom O  @ (0.39, 0.81, 0.79)
        Atom V  @ (0.76, 0.03, 0.97)
        Atom V  @ (0.76, 0.48, 0.47)
        ... omitting 2 atoms ...
    Lattice parameters:
        a=5.743Å, b=4.517Å, c=5.375Å
        α=90.000°, β=122.600°, γ=90.000°
    Chemical composition:
        O: 66.667%
        V: 33.333%
    Source:
        (...omitted...)\crystals\cifs\vo2-m1.cif >

Symmetry information is also readily available::

    >>> print(vo2.symmetry())
    {'international_symbol': 'P2_1/c', 
     'hall_symbol': '-P 2ybc', 
     'international_number': 14, 
     'hall_number': 81, 
     'international_full': 'P 1 2_1/c 1', 
     'pointgroup': 'C2h'}
    
Aknowledgements
===============

This package depends on the work of some amazing people. Of note are the `spglib contributors <https://github.com/atztogo/spglib>`_

Citations
=========

As this package is a spinoff from ``scikit-ued``, please consider citing the following publication if you find ``crystals`` useful:

.. [#] L. P. René de Cotret, M. R. Otto, M. J. Stern. and B. J. Siwick, *An open-source software ecosystem for the interactive 
       exploration of ultrafast electron scattering data*, Advanced Structural and Chemical Imaging **4**:11 (2018) DOI: 10.1186/s40679-018-0060-y

Underlying algorithms provided by ``spglib`` are described in the following publication:

.. [#] A. Togo and I. Tanaka, *spglib: a software library for crystal symmetry search*. https://arxiv.org/abs/1808.01590 (written at version 1.10.4).

Structure parsing from CIF files has been tested for correctness against CIF2CELL, detailed here:

.. [#] Torbjorn Bjorkman, *CIF2Cell: Generating geometries for electronic structure programs*, 
       Computer Physics Communications **182**, 1183-1186 (2011) doi: 10.1016/j.cpc.2011.01.013

Structure parsing from PDB files has been tested for correctness against ``Bio.PDB``, detailed here:

.. [#] Hamelryck, T., Manderick, B. *PDB parser and structure class implemented in Python*. Bioinformatics **19**: 2308–2310 (2003)

Support / Report Issues
=======================

All support requests and issue reports should be `filed on Github as an issue <https://github.com/LaurentRDC/crystals/issues>`_.

License
=======

``crystals`` is made available under the BSD 3-clause license. For more details, see `LICENSE.txt <https://github.com/LaurentRDC/crystals/blob/master/LICENSE.txt>`_.

Related projects
================

- Streaming operations on NumPy arrays are available in the `npstreams package <https://pypi.org/pypi/npstreams>`_.

- Interactive exploration of ultrafast electron diffraction data with the `iris-ued package <https://pypi.org/project/iris-ued/>`_.

- Data structures and algorithms to handle ultrafast electron scattering data in the `scikit-ued package <https://pypi.org/project/scikit-ued>`_.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`