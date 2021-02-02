
.. currentmodule:: crystals

User guide: atoms and elements
==============================
The basis of structure manipulations is to manipulate atoms. :class:`Atom` objects are in the
category of `Transformable` objects, meaning that their coordinates can be transformed
according to any affine transform.

To create an atom, simply provide its element and coordinates:
    
    >>> from crystals import Atom
    >>> 
    >>> copper = Atom(element = 'Cu', coords = [0,0,0])

Optional information can be given, such as magnetic moment and mean-squared displacement. For users of :mod:`ase`, 
another possibility is to instantiate an :class:`Atom` from an :class:`ase.Atom` using the :meth:`Atom.from_ase` 
constructor.

:class:`Atom` instances are hashable; they can be used as ``dict`` keys or stored in a ``set``.

Since we are most concerned with atoms in crystals, the coordinates here are assumed to be fractional.
If the atom was created as part of a structure, the real-space position with respect to its parent (:class:`Crystal` 
or :class:`Lattice`) can be accessed using the :meth:`Atom.coords_cartesian` method:

    >>> from crystals import Crystal
    >>> graphite = Crystal.from_database('C')
    >>> 
    >>> carbon = sorted(graphite)[-1] 
    >>> carbon.coords_fractional
    array([0.66667, 0.33334, 0.75   ])
    >>> carbon.coords_cartesian
    array([1.42259818e+00, 1.23200000e-05, 5.03325000e+00])

Atomic distances
----------------

The fractional/cartesian distance between two atoms sitting *on the same lattice* is possible:

    >>> from crystals import distance_fractional, distance_cartesian
    >>> graphite = Crystal.from_database('C')
    >>> 
    >>> carbon1, carbon2, *_ = tuple(sorted(graphite))
    >>> carbon1
    < Atom C  @ (0.00, 0.00, 0.25) >
    >>> carbon2
    < Atom C  @ (0.00, 0.00, 0.75) >
    >>> distance_fractional(carbon1, carbon2)
    0.5
    >>> distance_cartesian(carbon1, carbon2)    # in Angstroms
    3.3555000000000006

If atoms are not sitting on the same lattice, calculating the distance should not be defined. In this case, an exception is raised:

    >>> gold = Crystal.from_database('Au')
    >>> silver = Crystal.from_database('Ag')
    >>>
    >>> gold1, *_ = tuple(sorted(gold))
    >>> silver1, *_ = tuple(sorted(silver))
    >>>
    >>> distance_cartesian(gold1, silver1) # doctest: +SKIP
    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "(...omitted...)\crystals\atom.py",  in distance_cartesian
        "Cartesian distance is undefined if atoms are sitting on different lattices."
    RuntimeError: Distance is undefined if atoms are sitting on different lattices.

The :class:`Element` class
--------------------------

If all you want is access to elemental information, like atomic weights, you can instantiate an :class:`Element` instead of an :class:`Atom`:

    >>> from crystals import Element
    >>> Element("H")
    < Hydrogen >
    >>> Element("Hydrogen")
    < Hydrogen >

You can specify elements by their atomic number as well:

    >>> Element(29)
    < Copper >
    >>> Element(29).symbol
    'Cu'
    >>> Element(29).name
    'Copper'

:class:`Element` instances give access to atomic properties:

    >>> Element("Cu").mass # Atomic mass in [u]
    63.546
    >>> Element("Cu").atomic_number
    29

Since :class:`Atom` is a subclass of :class:`Element`, all of the above examples also work for :class:`Atom`:

    >>> from crystals import Atom
    >>> Atom("Cu", coords = [0,0,0]).mass
    63.546

Handling the atomic orbital structure of atoms
----------------------------------------------

For certain applications, access to the electronic structure of orbitals is important. To that end, :class:`Atom` instances carry this information in instances of :class:`ElectronicStructure`.

You can create an electronic structure by hand:

    >>> from crystals import ElectronicStructure
    >>> ElectronicStructure({"1s": 2, "2s": 2, "2p": 2})
    < ElectronicStructure: 1s²2s²2p² >

It is much more ergonomic to start from the ground state of an element. For example:

    >>> ElectronicStructure.ground_state("C")
    < ElectronicStructure: 1s²2s²2p² >

Once you have a starting point, the electronic structure can be modified for your application:

    >>> structure = ElectronicStructure.ground_state("C")
    >>> structure["2p"] -= 1
    >>> structure["3d"] += 1
    >>> structure
    < ElectronicStructure: 1s²2s²2p¹3d¹ >

You can always check which orbital is the outermost orbital:

    >>> structure = ElectronicStructure.ground_state("Ar")
    >>> structure.outer_shell
    <Orbital.three_p: '3p'>

Note that you cannot create impossible electronic structures, however:

    >>> structure = ElectronicStructure.ground_state("C")
    >>> structure["2s"] = 3 # doctest: +SKIP
    ValueError: There cannot be 3 electrons in orbital 2s

Finally, you can modify atomic electronic structures on a particular atom. By default, 
the electronic structure of atoms is set to the ground state. Let's move an electron up from 
the "2p" to "3d" orbital in one atom of graphite:

    >>> graphite = Crystal.from_database('C')
    >>> graphite
    < Crystal object with following unit cell:
        Atom C  @ (0.00, 0.00, 0.25)
        Atom C  @ (0.00, 0.00, 0.75)
        Atom C  @ (0.33, 0.67, 0.25)
        Atom C  @ (0.67, 0.33, 0.75)
    Lattice parameters:
        a=2.464Å, b=2.464Å, c=6.711Å
        α=90.000°, β=90.000°, γ=120.000°
    Chemical composition:
        C: 100.000% >
    >>> atom, *_ = sorted(graphite)
    >>> atom.electronic_structure["2p"] -= 1
    >>> atom.electronic_structure["3d"] += 1
    >>> graphite
    < Crystal object with following unit cell:
        Atom C  @ (0.00, 0.00, 0.25) | [1s²2s²2p¹3d¹]
        Atom C  @ (0.00, 0.00, 0.75)
        Atom C  @ (0.33, 0.67, 0.25)
        Atom C  @ (0.67, 0.33, 0.75)
    Lattice parameters:
        a=2.464Å, b=2.464Å, c=6.711Å
        α=90.000°, β=90.000°, γ=120.000°
    Chemical composition:
        C: 100.000% >
    
Note that atoms with ground-state electronic structure don't show it explicitly. 
You could entirely replace the electronic structure of an atom:

    >>> graphite = Crystal.from_database('C')
    >>> atom, *_ = sorted(graphite)
    >>> atom.electronic_structure = ElectronicStructure.ground_state("Ti") # This is just an example.
    >>> graphite
    < Crystal object with following unit cell:
        Atom C  @ (0.00, 0.00, 0.25) | [1s²2s²2p⁶3s²3p⁶4s²3d²]
        Atom C  @ (0.00, 0.00, 0.75)
        Atom C  @ (0.33, 0.67, 0.25)
        Atom C  @ (0.67, 0.33, 0.75)
    Lattice parameters:
        a=2.464Å, b=2.464Å, c=6.711Å
        α=90.000°, β=90.000°, γ=120.000°
    Chemical composition:
        C: 100.000% >