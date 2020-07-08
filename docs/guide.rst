.. _user_guide:

**********
User guide
**********

This page highlights some of the features of ``crystals``. In case anything is unclear, 
`raising an issue <https://github.com/LaurentRDC/crystals/issues/>`_ would be appreciated.

The :class:`Crystal` class
==========================

Handling crystal models is the main feature of the :mod:`crystals` library. This is done through the :class:`Crystal` class, a representation of a crystal including 
unit cell atoms, lattice parameters, and other goodies.

Since working with :class:`Crystal` instances is so important, there are many ways to construct them.

Constructing a :class:`Crystal` object
--------------------------------------

Creating a :class:`Crystal` object can be done most easily from a Crystal Information File (CIF, .cif)::
    
    >>> from crystals import Crystal
    >>> TiSe2 = Crystal.from_cif('tise2.cif')

:mod:`crystals` also has an internal database of CIF files. Valid names are stored in :attr:`Crystal.builtins` and can be
constructed like so::

    >>> Crystal.builtins
    frozenset({'Ac',
               'Ag',
               'Al',
               (...omitted...)
               'alpha-Mn',
               'b-Pu',
               'diamond',
               'gamma-Pu'
               })
    >>> 'Au' in Crystal.builtins
    True
    >>> Crystal.from_database('Au')
    < Crystal object with following unit cell:
        Atom Au @ (0.00, 0.00, 0.00)
        Atom Au @ (0.00, 0.50, 0.50)
        Atom Au @ (0.50, 0.00, 0.50)
        Atom Au @ (0.50, 0.50, 0.00)
    Lattice parameters:
        a=4.078Å, b=4.078Å, c=4.078Å
        α=90.000°, β=90.000°, γ=90.000°
    Chemical composition:
        Au: 100.000%
    Source:
        (...omitted...)crystals\cifs\Au.cif >

`RCSB Protein DataBank <http://www.rcsb.org/>`_ files are even easier to handle; simply provide the 4-letter identification code
and the structure file will be taken care of by :mod:`crystals`::
    
    >>> hemoglobin = Crystal.from_pdb('1gzx')
    >>> print(hemoglobin)
    < Crystal object with following unit cell:
        Atom C  @ (0.85, 0.79, 0.96)
        Atom C  @ (0.86, 0.85, 0.94)
        Atom C  @ (0.92, 0.99, 0.63)
        Atom C  @ (0.38, 0.68, 0.09)
        Atom C  @ (0.14, 0.62, 0.88)
        Atom C  @ (0.53, 0.84, 0.07)
        Atom C  @ (0.11, 0.77, 0.40)
        Atom C  @ (0.70, 0.23, 0.66)
        Atom C  @ (0.60, 0.10, 0.83)
        Atom C  @ (0.92, 0.10, 0.49)
        ... omitting 19066 atoms ...
    Lattice parameters:
        a=97.050Å, b=99.500Å, c=66.110Å
        α=90.000°, β=90.000°, γ=90.000°
    Chemical composition:
        C: 61.942%
        O: 21.367%
        N: 16.356%
        S: 0.252%
        Fe: 0.084%
    Source:
        (...omitted...)\crystals_cache\pdb1gzx.ent >

Another convenient way to construct a :class:`Crystal` is through the `Crystallography Open Database <http://www.crystallography.net/cod/>`_::

    >>> # Default is the latest revision
    >>> vo2 = Crystal.from_cod(1521124)
    >>> # Revisions are accessible as well
    >>> old_vo2 = Crystal.from_cod(1521124, revision = 140771)

the `Materials Project <https://materialsproject.org/>`_ provides another avenue where to get crystal structures.
You will need an API key from your `account dashboard <https://materialsproject.org/dashboard>`_:: 

    >>> fe2o3 = Crystal.from_mp(api_key="xxxxxxxxxxxxxxxx", query = "Fe2O3")
    >>> print(fe2o3)
    < Crystal object with following unit cell:
        Atom Fe @ (0.89, 0.37, 0.67)
        Atom Fe @ (0.99, 0.75, 0.96)
        Atom Fe @ (0.49, 0.75, 0.79)
        Atom Fe @ (0.25, 0.51, 0.71)
        Atom Fe @ (0.13, 0.39, 0.92)
        Atom Fe @ (0.88, 0.12, 0.59)
        Atom Fe @ (0.38, 0.38, 0.84)
        Atom Fe @ (0.62, 0.62, 0.66)
        Atom Fe @ (0.61, 0.87, 0.58)
        Atom Fe @ (0.38, 0.38, 0.16)
        ... omitting 150 atoms ...
        ... use repr() to show the full cell ...
    Lattice parameters:
        a=8.525Å, b=8.525Å, c=25.593Å
        α=90.000°, β=90.000°, γ=90.000°
    Chemical composition:
        O: 60.000%
        Fe: 40.000%
    Source:
        Materials Project: mp-1456 >

Other constructors are supported. See the reference for the :class:`Crystal` class for more details.

Constructing a :class:`Crystal` object by hand
----------------------------------------------
If you don't have a file on hand, or want to create an idealized crystal, consider building a :class:`Crystal`
object by hand.

To do this, you need:

1. iterable of :class:`Atom` objects, with coordinates. These atoms must be the full unit cell;
2. three lattice vectors;

As an example, let's create the simplest crystal structure known: 
`alpha-Polonium (simple cubic) <https://en.wikipedia.org/wiki/Polonium#Solid_state_form>`_::
    
    >>> from crystals import Crystal, Atom
    >>> import numpy as np
    >>>
    >>> lattice_vectors = 3.35 * np.eye(3)
    >>> unitcell = [Atom('Po', coords = [0,0,0])]
    >>>
    >>> polonium = Crystal(unitcell, lattice_vectors)

In the case where atoms are given as an asymmetric unit cell and a set of symmetry operators, you can use the
:func:`symmetry_expansion` function to generate a set of *unique* atoms (even if some symmetry operators might be redundant).
The generated set of atoms can be passed to the constructor of :class:`Crystal`.

Converting a Crystal to other formats
-------------------------------------

You can use the `crystals` package to convert crystal structures from one format to another. Currently, you can write a structure 
either to an `.xyz` (:meth:``Crystal.to_xyz`) file, a Crystallography Information Framework `.cif` (:meth:``Crystal.to_cif`), 
or to an ``ase.Atoms`` structure (:meth:``Crystal.to_ase``). Here is an example::

    >>> from crystals import Crystal
    >>> import numpy as np
    >>> # Create a crystal structure by hand
    >>> lattice_vectors = 3.35 * np.eye(3)
    >>> unitcell = [Atom('Po', coords = [0,0,0])]
    >>> polonium = Crystal(unitcell, lattice_vectors)
    >>> # Convert to CIF
    >>> polonium.to_cif('polonium.cif')

Would you like to convert to another format that is not supported yet? Please `raise an issue <https://github.com/LaurentRDC/crystals/issues/new>`_!

Crystal attributes
------------------
The :class:`Crystal` object provides some interfaces for easy structure manipulation. First, a :class:`Crystal` is an iterable::

    >>> from crystals import Crystal
    >>> graphite = Crystal.from_database('C')
    >>> 
    >>> for atm in graphite:
    ...	    print(atm)
    ...		
    < Atom C  @ (0.00, 0.00, 0.25) >
    < Atom C  @ (0.33, 0.67, 0.25) >
    < Atom C  @ (0.00, 0.00, 0.75) >
    < Atom C  @ (0.67, 0.33, 0.75) >
    
The :func:`len` of a :class:`Crystal` is the unit cell size (in number of atoms)::

    >>> hemoglobin = Crystal.from_pdb('1gzx')
    >>> len(hemoglobin)
    17536

The :class:`Crystal` class is a set-like container; checking containership (with the builtin ``in`` statement) is very fast::

    >>> graphite = Crystal.from_database('C')
    >>> carbon = next(iter(graphite))			# Pick any atom from the crystal
    >>> 
    >>> carbon in graphite 
    True

:class:`Crystal` instances can be equated to each other::

    >>> gold = Crystal.from_database('Au')
    >>> silver = Crystal.from_database('Ag')
    >>>
    >>> gold == silver
    False

Just like lists and other container types, a :class:`Crystal` instance is `False` if empty, and `True` otherwise::

    >>> mycrystal = Crystal.from_database('Cu')
    >>> if mycrystal: # equivalent to: if len(mycrystal) > 0:
    ...     pass

Structures can be extracted from a :class:`Crystal` instance by making use of its superclass, :class:`AtomicStructure`. For example, 
all atoms satisfying a certain condition can be found using `Crystal.satisfying`::

    >>> vo2 = Crystal.from_database('vo2-m1')
    >>> 
    >>> vo2.satisfying( lambda atom: atom.element == 'V' )
    < AtomicStructure object with following orphan atoms:
        Atom V  @ (0.24, 0.97, 0.03)
        Atom V  @ (0.76, 0.48, 0.47)
        Atom V  @ (0.76, 0.03, 0.97)
        Atom V  @ (0.24, 0.53, 0.53) >

To make it easier, take a look at the :func:`is_element` function:

    >>> vo2.satisfying( is_element('O') )
    < AtomicStructure object with following orphan atoms:
        Atom O  @ (0.61, 0.31, 0.71)
        Atom O  @ (0.10, 0.21, 0.20)
        Atom O  @ (0.10, 0.29, 0.70)
        Atom O  @ (0.39, 0.81, 0.79)
        Atom O  @ (0.90, 0.79, 0.80)
        Atom O  @ (0.61, 0.19, 0.21)
        Atom O  @ (0.39, 0.69, 0.29)
        Atom O  @ (0.90, 0.71, 0.30) >

If a :class:`Crystal` was generated from a file, the path to its file can be retrieved
from the :attr:`source` attribute::

    >>> hemoglobin = Crystal.from_pdb('1gzx')
    >>> print(hemoglobin.source)
    '(...omitted...)\pdb1gzx.ent'

:class:`Crystal` instances have a nice string representation (with :func:`str`), and a complete string representation (with :func:`repr`):

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
    

:class:`Crystal` instances can be converted to NumPy arrays as well::

    >>> import numpy
    >>> numpy.array(Crystal.from_database('Si'))
    array([[14.  ,  0.  ,  0.5 ,  0.5 ],
           [14.  ,  0.5 ,  0.5 ,  0.  ],
           [14.  ,  0.  ,  0.  ,  0.  ],
           [14.  ,  0.75,  0.75,  0.25],
           [14.  ,  0.5 ,  0.  ,  0.5 ],
           [14.  ,  0.75,  0.25,  0.75],
           [14.  ,  0.25,  0.75,  0.75],
           [14.  ,  0.25,  0.25,  0.25]])

:data:`arr` will contain one row per unit cell atom:

You can calculate what the asymmetric cell of a :class:`Class` is with the :meth:`Crystal.asymmetric_cell()`
method::

    >>> Crystal.from_database('C').asymmetric_cell()
    {< Atom C  @ (0.00, 0.00, 0.25) | [1s²2s²2p²] >,
    < Atom C  @ (0.67, 0.33, 0.75) | [1s²2s²2p²] >}


.. table::
    :widths: grid

    +--------------+-----------------------+
    |Atomic Number |Fractional coordinates |
    +==============+=======+========+======+
    |      Z1      |  x1   |   y1   |  z1  |
    +--------------+-------+--------+------+
    |      Z2      |  x2   |   y2   |  z2  |
    +--------------+-------+--------+------+
    |      Z3      |  x3   |   y3   |  z3  |
    +--------------+-------+--------+------+
    |            ...                       |
    +--------------------------------------+

Lattice vectors and reciprocal space
====================================
Once a :class:`Crystal` object is ready, you can manipulate the lattice parameters via the :class:`Lattice`
super-class. Let's use the built-in example of graphite::

    >>> from crystals import Crystal
    >>> import numpy as np
    >>> 
    >>> graphite = Crystal.from_database('C')
    >>> 	
    >>> a1, a2, a3 = graphite.lattice_vectors
    >>> a1
    array([2.464, 0.   , 0.   ])
    >>> b1, b2, b3 = graphite.reciprocal_vectors
    >>> b1
    array([ 2.54999404,  1.47223974, 0.    ])
    >>>
    >>> np.dot(a1, b1)	# 2 pi
    6.283185307179586

The standard `three lengths and angles` description of a lattice is also accessible::

    >>> a, b, c, alpha, beta, gamma = graphite.lattice_parameters

The unit cell volume (and by extensions, density) is also accessible::

    >>> vol = graphite.volume	# Angstroms cubed
    >>> density = vol/len(graphite)

As a lattice can be fully described by its basis vectors, a NumPy array can be created from a :class:`Lattice` instance.

Space-group Information
=======================
The `lattice system <https://en.wikipedia.org/wiki/Bravais_lattice#Bravais_lattices_in_3_dimensions>`_ of a Lattice or Crystal instance is also available via the :attr:`lattice_system` attribute::

    >>> vo2 = Crystal.from_database('vo2-m1') # Monoclinic M1 VO2
    >>> vo2.lattice_system
    <LatticeSystem.monoclinic: 2>

Better control on length tolerances is available via the :func:`lattice_system` function.

Thanks to `spglib <http://atztogo.github.io/spglib/>`_, we can get space-group information from a :class:`Crystal` instance::

    >>> from crystals import Crystal
    >>>
    >>> gold = Crystal.from_database('Au')
    >>> gold.symmetry()
    {'international_symbol': 'Fm-3m', 
     'hall_symbol': '-F 4 2 3',
     'hm_symbol': 'Fm-3m', 
     'international_number': 225, 
     'hall_number': 523, 
     'international_full': 'F 4/m -3 2/m', 
     'pointgroup': 'Oh'}

In the above example, :data:`spg_info` is a dictionary with the following keys:

* ``'international_symbol'``: International Tables of Crystallography space-group symbol (short);

* ``'international_full'``: International Tables of Crystallography space-group full symbol;

* ``'hall_symbol'`` : Hall symbol;

* ``'hm_symbol'`` : Hermann-Mauguin symbol;

* ``'pointgroup'`` : International Tables of Crystallography point-group;

* ``'international_number'`` : International Tables of Crystallography space-group number (between 1 and 230);

* ``'hall_number'`` : Hall number (between 1 and 531).

Each of those items are also available directly from the :class:`Crystal` instance. The Hall number of a crystal structure 
is located in the ``crystal.hall_number`` attribute, the short international symbol is located in the ``crystal.international_symbol`` 
attribute, and so on.

Symmetry operations
-------------------
You can get the matrix symmetry operations directly from the :class:`Crystal` class:

    >>> cryst = Crystal.from_database('C')
    >>> first_symop = cryst.symmetry_operations()[0]
    >>> 
    >>> print(first_symop)
    array([[1, 0, 0, 0],
           [0, 1, 0, 0],
           [0, 0, 1, 0]
           [0, 0, 0, 1]], dtype=int32), 
    
Symmetry operations are described using 4x4 affine matrices, where the rotation is the top 3x3 block,
and the translation is the right-most column. Example with iteration::

    >>> for m in cryst.symmetry_operations():
    ...     rotation = m[:3, :3]
    ...     translation = m[:3, -1]
    ...
    >>>

Symmetry operations in reciprocal space are also made available via ``Crystal.reciprocal_symmetry_operations``

Cell refinements and manipulations
----------------------------------

Again, through `spglib <http://atztogo.github.io/spglib/>`_, we can create different versions of ``Crystal`` instances. For example,
a primitive ``Crystal`` can be created using the ``Crystal.primitive`` method::

    >>> gold = Crystal.from_database('Au')
    >>> print(gold)
    < Crystal object with following unit cell:
        Atom Au @ (0.00, 0.00, 0.00)
        Atom Au @ (0.00, 0.50, 0.50)
        Atom Au @ (0.50, 0.00, 0.50)
        Atom Au @ (0.50, 0.50, 0.00)
    Lattice parameters:
        a=4.078Å, b=4.078Å, c=4.078Å
        α=90.000°, β=90.000°, γ=90.000°
    Chemical composition:
        Au: 100.000%
    Source:
        (...omitted...)\crystals\cifs\Au.cif >
    >>>
    >>> primitive_gold = gold.primitive() # this is a whole new Crystal instance
    >>> print(primitive_gold)
    < Crystal object with following unit cell:
        Atom Au @ (0.00, 0.00, 0.00)
    Lattice parameters:
        a=2.884Å, b=2.884Å, c=2.884Å
        α=60.000°, β=60.000°, γ=60.000°
    Chemical composition:
        Au: 100.000%
    Source:
        (...omitted...)\crystals\cifs\Au.cif >

Notice how the primitive structure is much simpler.

Idealized versions of ``Crystal``s are also made available. Let's take the example of iron arsenide::

    >>> iron_arsenide = Crystal.from_database('FeAs')
    >>> print(iron_arsenide)
    < Crystal object with following unit cell:
        Atom As @ (0.30, 0.08, 0.75)
        Atom As @ (0.70, 0.92, 0.25)
        Atom As @ (0.20, 0.58, 0.25)
        Atom As @ (0.80, 0.42, 0.75)
        Atom Fe @ (0.00, 0.20, 0.25)
        Atom Fe @ (1.00, 0.80, 0.75)
        Atom Fe @ (0.50, 0.30, 0.25)
        Atom Fe @ (0.50, 0.70, 0.75)
    Lattice parameters:
        a=5.440Å, b=6.026Å, c=3.371Å
        α=90.000°, β=90.000°, γ=90.000°
    Chemical composition:
        Fe: 50.000%
        As: 50.000%
    Source:
        (...omitted...)\crystals\cifs\FeAs.cif >
    >>>
    >>> idealized = iron_arsenide.ideal() # this is a whole new Crystal instance
    >>> print(idealized)
    < Crystal object with following unit cell:
        Atom As @ (0.70, 0.25, 0.08)
        Atom As @ (0.20, 0.25, 0.42)
        Atom As @ (0.80, 0.75, 0.58)
        Atom As @ (0.30, 0.75, 0.92)
        Atom Fe @ (0.00, 0.25, 0.80)
        Atom Fe @ (0.50, 0.25, 0.70)
        Atom Fe @ (0.50, 0.75, 0.30)
        Atom Fe @ (1.00, 0.75, 0.20)
    Lattice parameters:
        a=5.440Å, b=3.371Å, c=6.026Å
        α=90.000°, β=90.000°, γ=90.000°
    Chemical composition:
        Fe: 50.000%
        As: 50.000%
    Source:
        (...omitted...)\crystals\cifs\FeAs.cif >

The atomic coordinates have been swapped to follow conventions of the appropriate space-group.

Scattering utilities
====================
:class:`Lattice` objects have a few methods that make life easier when dealing with scattering data and modeling.

The conversion between Miller indices and scattering vectors is available:: 

    >>> from crystals import Crystal
    >>> graphite = Crystal.from_database('C')
    >>>
    >>> # Behavior inherited from Lattice superclass
    >>> e1 = (1, 0, 0)
    >>> G = graphite.scattering_vector(e1)
    >>> graphite.miller_indices(G)
    array([ 1.00000000e+00, -1.11022302e-16,  4.93038066e-32])

Supercells
==========
For various reasons, creating a supercell from a crystal might be desirable. The process is very easy. 
Let's imagine we want to create a 2x2x2 supercell of graphite::

    >>> from crystals import Crystal
    >>> graphite = Crystal.from_database('C')
    >>>
    >>> for atm in graphite: # Iterate over the unitcell
    ...     print(atm)
    ...
    < Atom C  @ (0.00, 0.00, 0.25) >
    < Atom C  @ (0.67, 0.33, 0.75) >
    < Atom C  @ (0.00, 0.00, 0.75) >
    < Atom C  @ (0.33, 0.67, 0.25) >
    >>>
    >>> for atm in graphite.supercell(2,2,2):   # Iterate over the supercell
    ...    print(atm)
    ...
    < Atom C  @ (0.67, 0.33, 0.75) >
    < Atom C  @ (0.67, 0.33, 1.75) >
    < Atom C  @ (0.67, 1.33, 0.75) >
    < Atom C  @ (0.67, 1.33, 1.75) >
    (... omitting 24 atoms ...)
    < Atom C  @ (1.00, 0.00, 0.75) >
    < Atom C  @ (1.00, 0.00, 1.75) >
    < Atom C  @ (1.00, 1.00, 0.75) >
    < Atom C  @ (1.00, 1.00, 1.75) >

Supercells are generated by copying unit cells along crystallographic axes. For example, a 2x3x5 supercell will be extended by 2x along 
the ``a1`` lattice vector, extended by 3x along the ``a2`` lattice vector, and extended by 5x along the ``a3`` lattice vector.

Users can recover the ``Crystal`` object underlying a supercell via the ``crystal`` attribute::

    >>> spcl = Crystal.from_database('C').supercell(2,3,4)
    >>> graphite = spcl.crystal

Compatibility with ASE
======================
The `Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase/index.html>`_ is a powerful tool. You can harness its power and convert
between :class:`ase.Atoms` and :class:`crystals.Crystal` at will.

To create an :class:`ase.Atoms` object from a :class:`Crystal`, use the :meth:`Crystal.ase_atoms` method::

    >>> from ase.calculators.abinit import Abinit
    >>> from crystals import Crystal
    >>>
    >>> gold = Crystal.from_database('Au')
    >>> ase_gold = gold.to_ase(calculator = Abinit(...))

All keywords of the :class:`ase.Atoms` constructor are supported. To get back to a :class:`Crystal` instance::

    >>> gold2 = Crystal.from_ase(ase_gold)

The :class:`Atom` Class
=======================
The basis of structure manipulations is to manipulate atoms. :class:`Atom` objects are in the
category of `Transformable` objects, meaning that their coordinates can be transformed
according to any affine transform.

To create an atom, simply provide its element and coordinates::
    
    >>> from crystals import Atom
    >>> 
    >>> copper = Atom(element = 'Cu', coords = [0,0,0])

Optional information can be given, such as magnetic moment and mean-squared displacement. For users of :mod:`ase`, 
another possibility is to instantiate an :class:`Atom` from an :class:`ase.Atom` using the :meth:`Atom.from_ase` 
constructor.

:class:`Atom` instances are hashable; they can be used as ``dict`` keys or stored in a ``set``.

Since we are most concerned with atoms in crystals, the coordinates here are assumed to be fractional.
If the atom was created as part of a structure, the real-space position with respect to its parent (:class:`Crystal` 
or :class:`Lattice`) can be accessed using the :meth:`Atom.coords_cartesian` method::

    >>> from crystals import Crystal
    >>> graphite = Crystal.from_database('C')
    >>> 
    >>> carbon = list(graphite)[-1]
    >>> carbon.coords_fractional
    array([0.  , 0.  , 0.75])
    >>> carbon.coords_cartesian
    array([0.  , 0.  , 5.033])

Atomic distances
----------------

The fractional/cartesian distance between two atoms sitting *on the same lattice* is possible::

    >>> from crystals import Crystal, distance_fractional, distance_cartesian
    >>> graphite = Crystal.from_database('C')
    >>> 
    >>> carbon1, carbon2, *_ = tuple(graphite)
    >>> carbon1
    < Atom C  @ (0.00, 0.00, 0.25) >
    >>> carbon2
    < Atom C  @ (0.67, 0.33, 0.75) >
    >>> distance_fractional(carbon1, carbon2)
    0.8975324197487241
    >>> distance_cartesian(carbon1, carbon2)    # in Angstroms
    3.6446077732986644

If atoms are not sitting on the same lattice, calculating the distance should not be defined. In this case, an exception is raised::

    >>> from crystals import Crystal
    >>> gold = Crystal.from_database('Au')
    >>> silver = Crystal.from_database('Ag')
    >>>
    >>> gold1, *_ = tuple(gold)
    >>> silver1, *_ = tuple(silver)
    >>>
    >>> distance_cartesian(gold1, silver1)
    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "(...omitted...)\crystals\atom.py", line 181, in distance_cartesian
        "Cartesian distance is undefined if atoms are sitting on different lattices."
    RuntimeError: Distance is undefined if atoms are sitting on different lattices.

The :class:`Element` class
==========================

If all you want is access to elemental information, like atomic weights, you can instantiate an :class:`Element` instead of an :class:`Atom`::

    >>> from crystals import Element
    >>> Element("H")
    < Hydrogen >
    >>> Element("Hydrogen")
    < Hydrogen >

You can specify elements by their atomic number as well::

    >>> Element(29)
    < Copper >
    >>> Element(29).symbol
    'Cu'
    >>> Element(29).name
    'Copper'

:class:`Element` instances give access to atomic properties::

    >>> Element("Cu").mass # Atomic mass in [u]
    63.546
    >>> Element("Cu").atomic_number
    29

Since :class:`Atom` is a subclass of :class:`Element`, all of the above examples also work for :class:`Atom`:

    >>> from crystals import Atom
    >>> Atom("Cu", coords = [0,0,0]).mass
    63.546

Handling the atomic orbital structure of atoms
==============================================

For certain applications, access to the electronic structure of orbitals is important. To that end, :class:`Atom` instances carry this information in instances of :class:`ElectronicStructure`.

You can create an electronic structure by hand::

    >>> from crystals import ElectronicStructure
    >>> ElectronicStructure({"1s": 2, "2s": 2, "2p": 2})
    < ElectronicStructure: 1s²2s²2p² >

It is much more ergonomic to start from the ground state of an element. For example::

    >>> ElectronicStructure.ground_state("C")
    < ElectronicStructure: 1s²2s²2p² >

Once you have a starting point, the electronic structure can be modified for your application:

    >>> structure = ElectronicStructure.ground_state("C")
    >>> structure["2p"] -= 1
    >>> structure["3d]  += 1
    >>> structure
    < ElectronicStructure: 1s²2s²2p¹3d¹ >

You can always check which orbital is the outermost orbital:: 

    >>> structure = ElectronicStructure.ground_state("Ar")
    >>> structure.outer_shell
    <Orbital.three_p: '3p'>

Note that you cannot create impossible electronic structures, however!

    >>> structure = ElectronicStructure.ground_state("C")
    >>> structure["2s"] = 3
    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "(...omitted...)\crystals\crystals\atom.py", line 489, in __setitem__
        f"There cannot be {value} electrons in orbital {shell.value}"
    ValueError: There cannot be 3 electrons in orbital 2s

Finally, you can modify atomic electronic structures on a particular atom. By default, the electronic structure of atoms is set to the ground state. Let's move an electron up from the "2p" to "3d" orbital in one atom of graphite::

    >>> graphite = Crystal.from_database('C')
    >>> graphite
    < Crystal object with following unit cell:
        Atom C  @ (0.00, 0.00, 0.25) | [1s²2s²2p²]
        Atom C  @ (0.67, 0.33, 0.75) | [1s²2s²2p²]
        Atom C  @ (0.33, 0.67, 0.25) | [1s²2s²2p²]
        Atom C  @ (0.00, 0.00, 0.75) | [1s²2s²2p²]
    Lattice parameters:
        a=2.464Å, b=2.464Å, c=6.711Å
        α=90.000°, β=90.000°, γ=120.000°
    Chemical composition:
        C: 100.000%
    Source:
        (...omitted...)\crystals\crystals\cifs\C.cif >
    >>> atom, *_ = list(graphite)
    >>> atom.electronic_structure["2p"] -= 1
    >>> atom.electronic_structure["3d"] += 1
    >>> graphite
    < Crystal object with following unit cell:
        Atom C  @ (0.00, 0.00, 0.25) | [1s²2s²2p¹3d¹]
        Atom C  @ (0.67, 0.33, 0.75) | [1s²2s²2p²]
        Atom C  @ (0.33, 0.67, 0.25) | [1s²2s²2p²]
        Atom C  @ (0.00, 0.00, 0.75) | [1s²2s²2p²]
    Lattice parameters:
        a=2.464Å, b=2.464Å, c=6.711Å
        α=90.000°, β=90.000°, γ=120.000°
    Chemical composition:
        C: 100.000%
    Source:
        (...omitted...)\crystals\crystals\cifs\C.cif >
    
You could entirely replace the electronic structure of an atom::

    >>> graphite = Crystal.from_database('C')
    >>> graphite
    < Crystal object with following unit cell:
        Atom C  @ (0.00, 0.00, 0.25) | [1s²2s²2p²]
        Atom C  @ (0.67, 0.33, 0.75) | [1s²2s²2p²]
        Atom C  @ (0.33, 0.67, 0.25) | [1s²2s²2p²]
        Atom C  @ (0.00, 0.00, 0.75) | [1s²2s²2p²]
    Lattice parameters:
        a=2.464Å, b=2.464Å, c=6.711Å
        α=90.000°, β=90.000°, γ=120.000°
    Chemical composition:
        C: 100.000%
    Source:
        (...omitted...)\crystals\crystals\cifs\C.cif >
    >>> atom, *_ = list(graphite)
    >>> atom.electronic_structure = ElectronicStructure.ground_state("Ti") # This is just an example.
    >>> graphite
    < Crystal object with following unit cell:
        Atom C  @ (0.00, 0.00, 0.25) | [1s²2s²2p⁶3s²3p⁶4s²3d²]
        Atom C  @ (0.00, 0.00, 0.75) | [1s²2s²2p²]
        Atom C  @ (0.33, 0.67, 0.25) | [1s²2s²2p²]
        Atom C  @ (0.67, 0.33, 0.75) | [1s²2s²2p²]
    Lattice parameters:
        a=2.464Å, b=2.464Å, c=6.711Å
        α=90.000°, β=90.000°, γ=120.000°
    Chemical composition:
        C: 100.000%
    Source:
        (...omitted...)\crystals\crystals\cifs\C.cif >


:ref:`Return to Top <user_guide>`