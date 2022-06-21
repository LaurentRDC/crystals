
What's new
==========

.. currentmodule:: crystals

Release 1.6.0
-------------

* Added the :meth:`Crystal.groupby` method to group unit cell atoms by site-symmetry (#12).

Release 1.5.0
-------------

* Added some typing information.
* Added the :attr:`Supercell.scaled_lattice_vectors` property and associated documentation (#11).
* Protein Data Bank downloads are now done through HTTPS rather than FTP, which is recommended by the RCSB data bank.
* Fixed some documentation formatting.

Release 1.4.1
-------------

* Fixed an issue with the `tag` attribute of `Atom` not being propagated properly (#9).

Release 1.4.0
-------------

* Added the ability to read and write POSCAR files from the Vienna Ab initio Simulation Package (VASP). Contributed by Chenxing Luo (#8).

Release 1.3.2
-------------

* Releases are now automatically performed using Github Actions
* Fixed an issue where uncertainties in atom site occupancy in CIF files would not be parsed correctly (#7).

Release 1.3.1
-------------

* The distinction between :class:`Supercell` and :class:`Crystal` no longer exists; :class:`Supercell` objects can be used everywhere a :class:`Crystal` is expected.

1.3.0
-----

* Starting with this version, ``crystals`` is licensed under GPLv3.
* General purpose single-crystal structure indexing with the DirAx algorithm has been added: :func:`index_dirax`.
* :meth:`Lattice.scattering_vector` and :meth:`Lattice.miller_indices` now accept tables of reflections/scattering vectors. This calculation is vectorized using NumPy.
* Migration of testing infrastructure to pytest.
* `Support for Python 3.6 and NumPy<1.17 has been dropped <https://numpy.org/neps/nep-0029-deprecation_policy.html>`_

1.2.2
-----

* The sorting of `AtomicStructure` objects is now stable.
* Fixed an issue where checking the containership of atoms did not work properly.
* Fixed an issue where downloading structures from the Materials Project failed on macOS/Linux.
* Code snippets in documentation are now tested for correctness.
* Tests are now included in the package itself.

1.2.1
-----

* Fixed deprecations that came with Python 3.9, involving comparison of crystal structures.

1.2.0
-----

* Added the ability to write crystal structures to CIF/XYZ files with the `Crystal.to_cif` and `Crystal.to_xyz` methods. Also, structures can be converted to ASE's ``Atoms`` class with `Crystal.to_ase`. This can be used to convert crystal structures from one format to another!
* Added the `symmetry_reduction` function, and associated method `Crystal.asymmetric_cell()`, to determine the smallest asymmetric cell that generates a unit cell.
* The method `Lattice.bounded_reflections` now takes in an additional parameter, `min_bound`, to find reflections between a lower and upper bound.
* Fixed an issue where in certain cases, atoms related by symmetry were not pruned appropriately (#5)
* Official support for Python 3.9.

1.1.2
-----

* Atom site occupancies are now parsed from CIF files `(#3) <https://github.com/LaurentRDC/crystals/issues/3>`_.

1.1.1
-----

* ``Orbital`` instances are now orderable, i.e. ``Orbital("1s") < Orbital("2p")``.
* The ``Element`` class can now be initialized from another ``Element``, and also from full atomic names (e.g. ``"hydrogen"``). This simplifies the normalization for input types.
* Added the ``ElectronicStructure.outer_shell`` method to quickly get the outermost shell.

1.1.0
-----

* Added the ability to describe the electronic structure of atoms using the ``ElectronicStructure`` class.

1.0.0
-----

* Added the ability to set a Materials Project API key as environment variable ``MATERIALS_PROJECT_API_KEY``.
* Added the Materials Project constructor to the ``crystals info`` script.

0.6.7
-----

* Added the ``Crystal.from_mp`` constructor to create crystals through the Materials Project API.
* Added the method ``Crystal.indexed_by``, which allows to index reflections of a crystal by using another lattice.
* Re-organized the ``Crystal`` class hierarchy to be more specific. This change should not affect anyone.

0.6.6
-----

* Added command-line utilities, including the `crystals info` script.

0.6.5
-----

* New algorithm to determine the lattice vectors from lattice parameters. Results should be the same except in edge cases.
* Space-group and symmetry information has been inserted into immutable dictionaries
* Fixed an issue where the parsing of CIF files without Hall symbol would fail because the International Table number was parsed as a string, not an integer.
* Updated minimum dependency bound on PyCifRW to 4.4.1 because of potential issues with Python 3.7.
* Refactoring of internal string handling to use f-strings.

0.6.4
-----

* ``Crystal.ideal`` and ``Crystal.primitive`` will now preserve subclasses.
* ``Crystal.primitive`` will now always return a new ``Crystal``, even if the ``Crystal`` already has a primitive cell.
* ``Supercell`` is no longer a subclass of ``Crystal``. The only recommended way of building supercells is now ``Crystal.supercell``.
* ``Lattice.scattering_vector`` and ``Lattice.miller_indices`` now work on single vectors of shape (3,).
* Added the ``Lattice.bounded_reflections`` generator. 


0.6.3
-----

* Added the constructor ``Crystal.from_pwscf`` to create crystal instances from output files generated by the `Plane-Wave Self-Consistent Field (PWSCF) program <https://www.quantum-espresso.org/Doc/pw_user_guide/>`_.  
* Added the ``Crystal.ideal`` method to create idealized crystal structures.
* Added the ``Crystal.reciprocal_symmetry_operations`` method to determine crystal symmetry operations in reciprocal basis.
* Symmetry-determination via ``Crystal.symmetry()`` and related methods now raises a ``RuntimeError`` if symmetry-determination fails, rather than returning ``None``
* Added the ``Atom.tag`` property, which can be used to tag ``Atom`` instances with extra information. This is for internal use. So far, it is used to label the order of Atoms from a PWSCF output file.

0.6.2
-----

* Added the ``Crystal.symmetry_operations`` method to determine crystal symmetry operations.

0.6.1
-----

* Fixed an issue where cache directories on MacOS could not be used.
* `CODParser` will now try to download files from three Crystallography Open Database mirrors.
* `AtomicStructure` and subclasses now support "truthiness", i.e. they are considered `False` if empty, and `True` otherwise.
* Added the `AtomicStructure.satisfying` method, to extract atoms satisfying a predicate from structures
* Added the `is_element` function. It can be used to make `AtomicStructure.satisfying` more readable.