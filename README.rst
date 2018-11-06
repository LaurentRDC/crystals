========
crystals
========

.. image:: https://readthedocs.org/projects/crystals/badge/?style=flat
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
files (like .cif and .pdb), provides access to atomic positions, and allows for space-group determination. Take a look at the `documentation <https://crystals.readthedocs.io/>`_
for more information.

Installation
============

``crystals`` is available on the Python Package Index::

    pip install crystals

``crystals`` can also be installed from source::

    git clone https://github.com/LaurentRDC/crystals.git
    cd crystals
    python setup.py install

You can install the latest development version using ``pip`` as well::

    python -m pip install git+git://github.com/LaurentRDC/crystals.git

To build documentation, you will need a few more packages, listed in ``dev-requirements.txt``. For example, to build documentation from source::

    git clone https://github.com/LaurentRDC/crystals.git
    cd crystals
    pip install -r dev-requirements.txt
    python setup.py build_sphinx

Documentation
=============

The documentation, including a user guide as well as detailed reference, is available here: https://crystals.readthedocs.io/

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

Development
===========

Tests can be run with the standard library's `unittest` module:: 

    python -m unittest discover

Some optional tests might be skipped if dependencies are not installed, e.g. `ASE`.

Related projects
================

- Streaming operations on NumPy arrays are available in the `npstreams package <https://pypi.org/pypi/npstreams>`_.

- Interactive exploration of ultrafast electron diffraction data with the `iris-ued package <https://pypi.org/project/iris-ued/>`_.

- Data structures and algorithms to handle ultrafast electron scattering data in the `scikit-ued package <https://pypi.org/project/scikit-ued>`_.