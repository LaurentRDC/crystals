# -*- coding: utf-8 -*-
from enum import Enum, unique
from functools import partial, wraps
from itertools import count, product, repeat, takewhile
from math import cos, isclose, radians, sin, sqrt
from typing import Any, Generator, Iterator, Iterable, Tuple, Union

from numpy.typing import ArrayLike
import numpy as np
from numpy.linalg import norm

from .affine import change_basis_mesh, change_of_basis


def primed(gen: Generator) -> Generator:
    """
    Decorator that primes a generator function, i.e. runs the function
    until the first ``yield`` statement. Useful in cases where there
    are preliminary checks when creating the generator.
    """

    @wraps(gen)
    def primed_gen(*args, **kwargs):
        generator = gen(*args, **kwargs)
        next(generator)
        return generator

    return primed_gen


def matmulrow(matrix: np.ndarray, arr: np.ndarray) -> np.ndarray:
    """Row-wise matrix multiplication."""
    return np.transpose(matrix @ arr.T)


@unique
class LatticeSystem(Enum):
    """
    Lattice system enumeration.

    Equivalent to Crystal family
    except that the hexagonal crystal family is split between
    the rhombohedral system and hexagonal system.
    """

    triclinic = 1
    monoclinic = 2
    orthorhombic = 3
    tetragonal = 4
    rhombohedral = 5
    hexagonal = 6
    cubic = 7


class Lattice:
    """
    Container class for lattice information and manipulations.

    Instances can also be create from the standard 'three lengths and angles'
    parameters via ``Lattice.from_parameters``:

    Parameters
    ----------
    lattice_vectors: iterable of `~numpy.ndarray`, shape (3,)
        Lattice vectors.
    """

    def __init__(self, lattice_vectors: ArrayLike, **kwargs):
        a1, a2, a3 = lattice_vectors
        self.a1 = np.asarray(a1, dtype=float)
        self.a2 = np.asarray(a2, dtype=float)
        self.a3 = np.asarray(a3, dtype=float)
        super().__init__(**kwargs)

    def __repr__(self) -> str:
        a, b, c, alpha, beta, gamma = self.lattice_parameters
        return f"< Lattice object with parameters {a:.3f}Å, {b:.3f}Å, {c:.3f}Å, {alpha:.2f}°, {beta:.2f}°, {gamma:.2f}° >"

    def __hash__(self) -> int:
        return hash(self.lattice_parameters)

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, Lattice):
            return np.allclose(self.lattice_vectors, other.lattice_vectors, atol=1e-3)
        return NotImplemented

    def __array__(self, *args, **kwargs) -> np.ndarray:
        """Returns a 3x3 float array in which each row is a lattice vector"""
        return np.array(self.lattice_vectors, *args, **kwargs)

    @classmethod
    def from_parameters(
        cls, a: float, b: float, c: float, alpha: float, beta: float, gamma: float
    ):
        """
        Create a Lattice instance from three lengths and angles.

        Parameters
        ----------
        a, b, c : floats
            Lattice vectors lengths [Å]
        alpha, beta, gamma : floats
            Angles between lattice vectors [deg]

        Raises
        ------
        ValueError : if lattice parameters are invalid.
        """
        return cls(lattice_vectors_from_parameters(a, b, c, alpha, beta, gamma))

    @property
    def lattice_parameters(self) -> Tuple[float, float, float, float, float, float]:
        """Lattice parameters as three lengths [Å] and three angles [degrees]."""
        a, b, c = norm(self.a1), norm(self.a2), norm(self.a3)
        alpha = np.arccos(np.vdot(self.a2, self.a3) / (b * c))
        beta = np.arccos(np.vdot(self.a1, self.a3) / (a * c))
        gamma = np.arccos(np.vdot(self.a1, self.a2) / (a * b))
        return a, b, c, np.rad2deg(alpha), np.rad2deg(beta), np.rad2deg(gamma)

    @property
    def lattice_system(self) -> "LatticeSystem":
        """One of the seven lattice system, returned in the form of the :class:`LatticeSystem` enumeration."""
        return lattice_system(*self.lattice_parameters, atol=5e-2)

    @property
    def volume(self) -> float:
        """Lattice cell volume Angtroms cubed"""
        return np.dot(self.a1, np.cross(self.a2, self.a3))

    @property
    def lattice_vectors(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Iterable of lattice vectors"""
        return self.a1, self.a2, self.a3

    @lattice_vectors.setter
    def lattice_vectors(self, vectors: ArrayLike):
        self.a1, self.a2, self.a3 = vectors

    @property
    def reciprocal(self) -> "Lattice":
        """Reciprocal lattice"""
        return Lattice(lattice_vectors=self.reciprocal_vectors)

    @property
    def reciprocal_vectors(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Reciprocal lattice vectors, defined as:

        .. math::

            b_i = 2 \\pi \\frac{a_j \\times a_k}{v}

        For :math:`v` the unit cell volume.
        """
        cell_volume = self.volume
        b1 = 2 * np.pi * np.cross(self.a2, self.a3) / cell_volume
        b2 = 2 * np.pi * np.cross(self.a3, self.a1) / cell_volume
        b3 = 2 * np.pi * np.cross(self.a1, self.a2) / cell_volume
        return b1, b2, b3

    @property
    def periodicity(self) -> Tuple[float, float, float]:
        """
        Crystal periodicity in x, y and z direction from the lattice constants.
        This is effectively a bounding cube for the unit cell, which is itself a unit cell.

        Returns
        -------
        x, y, z : float
            Periodicity in Angstroms along the x, y, and z direction, respectively.
        """
        # Add the absolute value of the component of every lattice vector
        # along the three euclidian vectors, which is effectively the sum of
        # absolutes of columns
        lv = np.abs(np.array(self.lattice_vectors))
        return tuple(lv.sum(axis=0))

    def scattering_vector(self, reflection: ArrayLike) -> np.ndarray:
        """
        Scattering vector from Miller indices.

        .. versionchanged:: 1.3
           Can now operate on tables of reflections, where every
           reflection is a row.

        Parameters
        ----------
        reflection : array_like, shape (3,) or (N, 3)
            Miller indices.

        Returns
        -------
        G : ndarray, shape (3,) or (N, 3)
            Scattering vector in :math:`A^{-1}`.
        """
        COB = change_of_basis(basis1=self.reciprocal_vectors, basis2=np.eye(3))
        return matmulrow(COB, np.asarray(reflection))

    def miller_indices(self, scattering_vector: ArrayLike) -> np.ndarray:
        """
        Miller indices from scattering vector components.

        .. versionchanged:: 1.3
           Can now operate on tables of vectors, where every
           vector is a row.

        Parameters
        ----------
        scattering_vector : array_like, shape (3,) or (N, 3)
            Scattering vector in :math:`A^{-1}`.

        Returns
        -------
        reflection : `~numpy.ndarray`, shape (3,) or (N, 3)
            Miller indices.
        """
        COB = change_of_basis(basis1=np.eye(3), basis2=self.reciprocal_vectors)
        return matmulrow(COB, np.asarray(scattering_vector))

    @staticmethod
    def frac_mesh(
        *xi,
        indexing: str = "xy",  # Not using Literal type because it landed in Python 3.8
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Coordinate arrays for fractional coordinates.

        Parameters
        ----------
        x1, x2, x3 : `~numpy.ndarray`, shape (N,)
            1d coordinate vectors. If only ``x1`` is provided, it is assumed
            that ``x1 = x2 = x3``. Otherwise, three coordinate vectors are expected.
        indexing : str, {'ij', 'xy'}
            Cartesian (‘xy’, default) or matrix (‘ij’) indexing of output.

        Returns
        -------
        out1, out2, out3 : `~numpy.ndarray`
            Fractional coordinate arrays.

        Raises
        ------
        ValueError : if number of input vectors is neither 1 nor 3.

        See Also
        --------
        numpy.meshgrid : Coordinate arrays from coordinate vectors
        Lattice.mesh : Real-space coordinate arrays from fractional coordinate vectors
        """
        if len(xi) == 1:
            xi = tuple(repeat(xi[0], times=3))
        elif len(xi) != 3:
            raise ValueError(
                f"1 or 3 coordinate arrays are required, but received {len(xi)}"
            )

        return np.meshgrid(*xi, indexing=indexing)

    def mesh(
        self,
        *xi: np.ndarray,
        indexing: str = "xy",  # Not using Literal type because it landed in Python 3.8
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Cartesian coordinate arrays from fractional coordinate vectors.

        Parameters
        ----------
        x1, x2, x3 : `~numpy.ndarray`, shape (N,)
            1d coordinate vectors in fractional coordinates.
            If only ``x1`` is provided, it is assumed that ``x1 = x2 = x3``.
            Otherwise, three coordinate vectors are expected.
        indexing : str, {'ij', 'xy'}
            Cartesian (‘xy’, default) or matrix (‘ij’) indexing of output.

        Returns
        -------
        out1, out2, out3 : `~numpy.ndarray`
            Real-space (cartesian) coordinate arrays.

        Raises
        ------
        ValueError : if number of input vectors is neither 1 nor 3.

        See Also
        --------
        numpy.meshgrid : Coordinate arrays from coordinate vectors
        Lattice.frac_mesh : Coordinate arrays for fractional coordinates
        """
        return change_basis_mesh(
            *self.frac_mesh(*xi, indexing=indexing),
            basis1=np.array(self.lattice_vectors),
            basis2=np.eye(3),
        )

    # Primed generators allows for checks on creation
    # All lines of code before the first 'yield' will run at first
    @primed
    def bounded_reflections(
        self, bound: float, min_bound: float = 0
    ) -> Iterator[Tuple[int, int, int]]:
        """
        Generates reflections (hkl) with norm(G) <= bound

        Parameters
        ----------
        bound : float
            Maximal scattering vector norm [:math:`A^{-1}`].
        min_bound : float, optional
            Minimal scattering vector norm [:math:`A^{-1}`].

            .. versionadded:: 1.2.0

        Yields
        ------
        reflection : 3-tuple of ints
            Miller indices of a bounded reflection.

        Raises
        ------
        ValueError : if `bound` is smaller than `min_bound`.

        Examples
        --------
        >>> cryst = Crystal.from_database('C')
        >>> refls = cryst.bounded_reflections(1.5) # 1.5 inverse Angstroms
        >>> list(refls)
        [(0, 0, -1), (0, 0, 0), (0, 0, 1)]
        """
        if bound < min_bound:
            raise ValueError(f"Bound {bound} is smaller than the minimum {min_bound}.")

        # Determine the maximum index such that (i00) family is still within data limits
        # This provides a (large) upper bound so that we are sure that the overall filtering will terminate
        bounded = lambda i: any(
            [np.linalg.norm(i * b) <= bound for b in self.reciprocal_vectors]
        )

        max_index = max(takewhile(bounded, count(0)))
        extent = range(-max_index, max_index + 1)
        refls = product(extent, repeat=3)

        yield  # Priming the generator ends here

        # Generalized hypotenuse
        def _hypot(*args):
            return sqrt(sum(map(lambda i: i**2, args)))

        # The above bound was only a first pass. We can refine further
        in_bounds = (
            lambda refl: min_bound <= _hypot(*self.scattering_vector(refl)) <= bound
        )
        yield from filter(in_bounds, refls)


def lattice_vectors_from_parameters(
    a: float, b: float, c: float, alpha: float, beta: float, gamma: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """

    Parameters
    ----------
    a, b, c : floats
        Lattice vectors lengths [Å]
    alpha, beta, gamma : floats
        Angles between lattice vectors [deg]
    """
    # The algorithm below is a refactoring of the
    # one made available in Diffpy.Structure
    alpha, beta, gamma = map(radians, (alpha, beta, gamma))

    # The unit cell volume when a = b = c = 1.
    unit_volume = sqrt(
        1.0
        + 2.0 * cos(alpha) * cos(beta) * cos(gamma)
        - cos(alpha) ** 2
        - cos(beta) ** 2
        - cos(gamma) ** 2
    )

    # reciprocal lattice
    a_recip = sin(alpha) / (a * unit_volume)
    cos_gamma_recip = (cos(alpha) * cos(beta) - cos(gamma)) / (sin(alpha) * sin(beta))
    sin_gamma_recip = sqrt(1 - cos_gamma_recip**2)

    a1 = np.asfarray(
        [1 / a_recip, -cos_gamma_recip / sin_gamma_recip / a_recip, cos(beta) * a]
    )
    a2 = np.asfarray([0, b * sin(alpha), b * cos(alpha)])
    a3 = np.asfarray([0, 0, c])
    return a1, a2, a3


# TODO: also determine body-centered, primitive, face-centered, etc.
#       https://en.wikipedia.org/wiki/Bravais_lattice#Bravais_lattices_in_3_dimensions
def lattice_system(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
    atol: float = 1e-2,
) -> LatticeSystem:
    """
    Determine the lattice system. All cyclic permutations are checked,
    so that no convention on ordering of lattice parameters is assumed.

    Parameters
    ----------
    a, b, c : floats
        Lattice vectors lengths [Å]
    alpha, beta, gamma : floats
        Angles between lattice vectors [deg]
    atol : float, optional
        Absolute tolerance (in Angstroms)

    Returns
    -------
    system : LatticeSystem
        One of the seven lattice system.
    """
    lengths, angles = (a, b, c), (alpha, beta, gamma)

    angleclose = partial(isclose, abs_tol=1)
    lengthclose = partial(isclose, abs_tol=atol)

    lengths_equal = all(lengthclose(length, a) for length in lengths)
    angles_equal = all(angleclose(angle, alpha) for angle in angles)

    # Checking for monoclinic lattice system is generalized
    # to the case where (a, b, c) can be cycled
    # i.e. a != c and beta != 90
    #   or b != c and alpha != 90
    #   or a != b and gamma != 90
    for clengths, cangles in zip(cyclic(lengths), cyclic(angles)):
        (l1, _, l3), (a1, a2, a3) = clengths, cangles
        if (
            (not lengthclose(l1, l3))
            and angleclose(a1, 90)
            and angleclose(a3, 90)
            and (not angleclose(a2, 90))
        ):
            return LatticeSystem["monoclinic"]

    if lengths_equal and angles_equal:
        if angleclose(alpha, 90):
            return LatticeSystem["cubic"]
        else:
            return LatticeSystem["rhombohedral"]

    # Special note : technically, a hexagonal lattice system
    # could have all three lengths equal
    elif lengths_equal and (not angles_equal):
        if any(isclose(angle, 120) for angle in angles) and (
            sum(isclose(i, 90) for i in angles) == 2
        ):
            return LatticeSystem["hexagonal"]

    # At this point, two lengths are equal at most
    elif _two_equal(lengths, atol=atol):
        if angles_equal and angleclose(alpha, 90):
            return LatticeSystem["tetragonal"]

        elif any(isclose(angle, 120) for angle in angles) and (
            sum(isclose(i, 90) for i in angles) == 2
        ):
            return LatticeSystem["hexagonal"]

    # At this point, all lengths are unequal
    elif angles_equal and angleclose(alpha, 90):
        return LatticeSystem["orthorhombic"]

    else:
        return LatticeSystem["triclinic"]


def _two_equal(iterable, atol):
    """Returns True if and only if two items from an iterable are equal"""
    iterable = tuple(iterable)
    for i in iterable:
        if sum(isclose(i, l, abs_tol=atol) for l in iterable) == 2:
            return True
    return False


def cyclic(iterable: Iterable[Any]) -> Iterable[Iterable[Any]]:
    """
    Yields cyclic permutations of an iterable.

    Examples
    --------
    >>> list(cyclic((1,2,3)))
    [(1,2,3), (2,3,1), (3,1,2)]
    """
    iterable = tuple(iterable)
    n = len(iterable)
    yield from (tuple(iterable[i - j] for i in range(n)) for j in range(n))
