============
Installation
============

``crystals`` is available on the Python Package Index::

    pip install crystals

From source
-----------

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