=====
Index
=====

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

``crystals`` is a library of data structure and algorithms to manipulate abstract crystals. ``crystals`` helps with reading crystallographic 
files (like .cif and .pdb), provides access to atomic positions, and allows for space-group determination. Although ``crystals`` can be used on its own,
it was made to be integrated into larger projects (like `scikit-ued <https://github.com/LaurentRDC/scikit-ued>`_).

Take a look at the `documentation <https://crystals.readthedocs.io/>`_ for more information.

.. toctree::
    :maxdepth: 2
    
    index
    installation
    guide
    reference
    contributing

Citations
=========

If you find this software useful, please consider citing the following publication:

.. [#] L. P. René de Cotret, M. R. Otto, M. J. Stern. and B. J. Siwick, *An open-source software ecosystem for the interactive 
       exploration of ultrafast electron scattering data*, Advanced Structural and Chemical Imaging **4**:11 (2018) DOI: 10.1186/s40679-018-0060-y

Underlying algorithms provided by ``spglib`` are described in the following publication:

.. [#] A. Togo and I. Tanaka, *spglib: a software library for crystal symmetry search*. https://arxiv.org/abs/1808.01590 (written at version 1.10.4).

Aknowledgements
===============

This package depends on the work of some amazing people. Of note are the `spglib contributors <https://github.com/atztogo/spglib>`_

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