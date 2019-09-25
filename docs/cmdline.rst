.. _cmdline:

**********************
Command-line utilities
**********************

Starting with version 0.6.6, ``crystals`` includes command-line utilities for repetitive tasks. To see the available commands
on your system:

.. code:: bash

    crystals --help

Crystal information
-------------------

The most important command-line utility is the crystallographic information script, ``crystals info``. 
To show information about a crystal file (``vo2.cif``, for example):

.. code:: bash

    crystals info vo2.cif

``crystals`` is smart about inputs. You can let it guess. For example, the Crystallographic Open Database
entry for Bismuth is 5000215. You can look at the crystallographic information as follows:

.. code:: bash

    crystals info 5000215

or, more precisely:

.. code:: bash

    crystals info 5000215 --type cod

See ``crystals info --help`` for more details.
