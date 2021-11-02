crystals
========
[![Documentation Status](https://readthedocs.org/projects/crystals/badge/?version=master)](https://crystals.readthedocs.io/) [![PyPI Package latest release](https://img.shields.io/pypi/v/crystals.svg)](https://pypi.python.org/pypi/crystals) [![Conda-forge Version](https://img.shields.io/conda/vn/conda-forge/crystals.svg)](https://anaconda.org/conda-forge/crystals) [![DOI badge](https://img.shields.io/badge/DOI-10.1186%2Fs40679--018--0060--y-blue)](https://doi.org/10.1186/s40679-018-0060-y)

`crystals` is a library of data structure and algorithms to manipulate
abstract crystals in a Pythonic way. `crystals` helps with reading
crystallographic files (like .cif and .pdb), provides access to atomic
positions, scattering utilities, allows for symmetry determination, and
indexing of diffraction peaks. Although `crystals` can be used on its own, 
it was made to be integrated into larger projects (like
[scikit-ued](https://github.com/LaurentRDC/scikit-ued)).

Take a look at the [documentation](https://crystals.readthedocs.io/) for
more information and examples.

Installation
------------

`crystals` is available on the Python Package Index:

    pip install crystals

For users of the conda package manager, `crystals` is also available
from the conda-forge channel:

    conda install -c conda-forge crystals

### From source

`crystals` can also be installed from source:

    git clone https://github.com/LaurentRDC/crystals.git
    cd crystals
    python setup.py install

You can install the latest development version using `pip` as well:

    python -m pip install git+git://github.com/LaurentRDC/crystals.git

To build documentation, you will need a few more packages, listed in
`dev-requirements.txt`. For example, to build documentation from source:

    git clone https://github.com/LaurentRDC/crystals.git
    cd crystals
    pip install -r dev-requirements.txt
    python setup.py build_sphinx

Documentation
-------------

The documentation, including user guides as well as detailed reference,
is available here: <https://crystals.readthedocs.io/>

Development
-----------

Tests can be run with the `pytest` package:

    python -m pytest --pyargs crystals

Some optional tests might be skipped if dependencies are not installed,
e.g. [ASE](https://wiki.fysik.dtu.dk/ase/).

Citations
---------

As this package is a spinoff from `scikit-ued`, please consider citing
the following publication if you find `crystals` useful:

> L. P. René de Cotret, M. R. Otto, M. J. Stern. and B. J. Siwick, *An open-source software ecosystem for the interactive exploration of ultrafast electron scattering data*, Advanced Structural and Chemical Imaging 4:11 (2018) [DOI: 10.1186/s40679-018-0060-y.](https://ascimaging.springeropen.com/articles/10.1186/s40679-018-0060-y)

Underlying algorithms provided by `spglib` are described in the
following publication:

> A. Togo and I. Tanaka, *spglib: a software library for crystal symmetry search*. [https://arxiv.org/abs/1808.01590](https://arxiv.org/abs/1808.01590) (written at version 1.10.4).

Structure parsing from CIF files has been tested for correctness against
CIF2CELL, detailed here:

> Torbjorn Bjorkman, *CIF2Cell: Generating geometries for electronic structure programs*, Computer Physics Communications 182, 1183-1186 (2011) [DOI: 10.1016/j.cpc.2011.01.013](https://doi.org/10.1016/j.cpc.2011.01.013)

Structure parsing from PDB files has been tested for correctness against
`Bio.PDB`, detailed here:

> Hamelryck, T., Manderick, B. *PDB parser and structure class implemented in Python*. Bioinformatics 19: 2308–2310 (2003)

Atomic weights are reported in the following publication:

> Meija, J., Coplen, T., Berglund, M., et al. (2016). *Atomic weights of the elements 2013* (IUPAC Technical Report). Pure and Applied Chemistry, 88(3), pp. 265-291. Retrieved 30 Nov. 2016, [DOI:10.1515/pac-2015-0305](https://doi.org/10.1515/pac-2015-0305)

Covalent radii are reported in the following article:

> Cordero, B. et al. (2008). *Covalent radii revisited*. Dalton Transactions, issue 21, pp. 2832-2838. The Royal Society of Chemistry. [DOI: 10.1039/B801115j](https://dx.doi.org/10.1039/B801115J)

Support / Report Issues
-----------------------

All support requests and issue reports should be [filed on Github as an
issue](https://github.com/LaurentRDC/crystals/issues).

License
-------

`crystals` is made available under the GPLv3 license. For more
details, see [LICENSE](https://github.com/LaurentRDC/crystals/blob/master/LICENSE).

Related projects
----------------

-   Streaming operations on NumPy arrays are available in the [npstreams
    package](https://pypi.org/pypi/npstreams).
-   Interactive exploration of ultrafast electron diffraction data with
    the [iris-ued package](https://pypi.org/project/iris-ued/).
-   Data structures and algorithms to handle ultrafast electron
    scattering data in the [scikit-ued
    package](https://pypi.org/project/scikit-ued).
