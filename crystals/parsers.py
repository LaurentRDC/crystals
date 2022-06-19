# -*- coding: utf-8 -*-
"""
Atomic structure parsers.
"""
import re
from typing import Any, Iterable, Optional, Tuple, Union
import warnings
from abc import abstractmethod
from contextlib import AbstractContextManager, suppress
from functools import lru_cache
from itertools import repeat
from os import PathLike, environ
from pathlib import Path
from platform import system
from string import digits, punctuation
from tempfile import gettempdir
from urllib.error import URLError
from urllib.request import urlretrieve
from warnings import warn

import numpy as np
import requests
from CifFile import ReadCif
import CifFile as Cif
from numpy.linalg import inv

from . import __version__
from .affine import affine_map, transform
from .atom import Atom, frac_coords
from .biological import Helix, Residue, SecondaryStructure, Sheet
from .lattice import Lattice
from .spg_data import HM2Hall, Number2Hall, SymOpsHall

# Temporary directory in which to cache crystal structure files
# downloaded from the internet.
# Note : temporary directory on MacOS systems (Darwin) do not work with gettempdir
#        In this case, we use "/tmp"
STRUCTURE_CACHE = (
    Path("/tmp" if system() == "Darwin" else gettempdir()) / "crystals_cache"
)

MISSING_MP_API_KEY_MESSAGE = """
You have not provided a Materials Project API key. Either provide it, or set it 
as the `MATERIALS_PROJECT_API_KEY` environment variable.
"""


def get_number_with_esd(x: Any) -> Tuple[float, float]:
    """pycifrw's version cannot handle floats, only strings."""
    try:
        return Cif.get_number_with_esd(x)
    except TypeError:
        return Cif.get_number_with_esd(str(x))


class ParseError(IOError):
    """General parsing error type."""

    pass


class AbstractStructureParser(AbstractContextManager):
    """
    Abstract base class for structure parsers. The preferred method
    of using this object is as a context manager.

    Parameters
    ----------
    filename : str or path-like
        Location of the CIF file.
    """

    @abstractmethod
    def __init__(self, *args, **kwargs):
        pass

    @abstractmethod
    def lattice_vectors(self) -> Iterable[np.ndarray]:
        """
        Returns the lattice vectors associated to a structure.

        Returns
        -------
        lv : iterable of ndarrays, shape (3,)
        """
        pass

    def symmetry_operators(self) -> Iterable[np.ndarray]:
        """
        Returns the symmetry operators that map the fractional atomic positions in a
        structure to the crystal *conventional* unit cell.

        Returns
        -------
        sym_ops : iterable of ndarray, shape (4,4)
            Transformation matrices. Since translations and rotation are combined,
            the transformation matrices are 4x4.
        """
        # Default implementation is the trivial transformation
        return [np.eye(4)]

    @abstractmethod
    def atoms(self):
        """
        Asymmetric unit cell. Combine with AbstractStructureParser.symmetry_operators()
        for a full unit cell.

        Returns
        -------
        atoms : iterable of Atom instance
        """
        pass


class PDBParser(AbstractStructureParser):
    """
    Collection of methods that parses Protein DataBank (PDB) files. This object should be used as a context manager.

    Parameters
    ----------
    ID : str
        Protein DataBank identification. The correct .pdb file will be downloaded,
        cached and parsed.
    download_dir : path-like object or None, optional
        Directory where to save the PDB file.
    overwrite : bool, optional
        Whether or not to overwrite files in cache if they exist. If no revision
        number is provided, files will always be overwritten.
    """

    def __init__(
        self, ID: str, download_dir: Optional[PathLike] = None, overwrite: bool = False
    ):
        if download_dir is None:
            download_dir = STRUCTURE_CACHE

        filename = self.download_pdb_file(
            pdb_code=ID, download_dir=download_dir, overwrite=overwrite
        )
        self._handle = open(filename, "r")

    def __exit__(self, *args, **kwargs):
        self._handle.close()

    @staticmethod
    def download_pdb_file(
        pdb_code: str,
        download_dir: Optional[PathLike] = None,
        server: str = "https://files.rcsb.org",
        overwrite: bool = False,
    ) -> Path:
        """
        Retrieves a PDB structure file from a PDB server and
        stores it in a local file tree.

        Parameters
        ----------
        pdf_code : str, len 4
            PDB ID code
        download_dir : path-like object
            Directory where to save the PDB file. Default is a local folder in the current directory
        server : str, optional
            Root address of the server from which to download the PDB file. Default is the main server.
        overwrite : bool, optional
            If True, existing PDB file with the same structure will be overwritten. Default is False.

        Returns
        -------
        file : pathlib.Path
            Pointer to the downloaded file
        """
        if download_dir is None:
            path = STRUCTURE_CACHE
        else:
            path = Path(download_dir)

        path.mkdir(exist_ok=True)

        final_file = path / f"pdb{pdb_code.lower()}.ent"  # (decompressed)

        # Skip download if the file already exists
        if (not overwrite) and (final_file.exists()):
            return final_file

        resp = requests.get(server + f"/download/{pdb_code.upper()}.pdb")
        resp.raise_for_status()

        with open(final_file, "wb") as out:
            out.write(resp.content)

        return Path(final_file)

    @property
    def filename(self) -> str:
        return self._handle.name

    @lru_cache(maxsize=1)
    def lattice_vectors(self) -> Iterable[np.ndarray]:
        """
        Returns the lattice vectors associated to a PDB structure.

        Returns
        -------
        lv : list of ndarrays, shape (3,)

        Raises
        ------
        ParseError
            If the file does not contain a CRYST1 tag.
        """
        self._handle.seek(0)

        for line in filter(lambda l: l.startswith("CRYST1"), self._handle):
            # characters are described in the PDB file content guide.
            a, b, c = float(line[6:15]), float(line[15:24]), float(line[24:33])
            alpha, beta, gamma = (
                float(line[33:40]),
                float(line[40:47]),
                float(line[47:54]),
            )
            break
        else:
            raise ParseError("No CRYST1 line found")

        return Lattice.from_parameters(a, b, c, alpha, beta, gamma).lattice_vectors

    def helices(
        self, ignored: Iterable[str] = ("HOH", "LI1", "SQU")
    ) -> Iterable[Helix]:
        """Returns an iterable of helices present in the protein.

        Parameters
        ----------
        ignored : iterable of str, optional
            3-letter string code for residues to ignore.

        Returns
        -------
        hex : iterable of Helix instances
        """
        residues = self.residues(ignored)
        helices = list()
        self._handle.seek(0)

        # Helices are defined as a list of residues making them up
        # Therefore, we simply filter residues down to the ones in range
        for line in filter(lambda line: line.startswith("HELIX"), self._handle):
            seq_num = int(line[7:10])  # TODO: place in helix
            seq_range = range(int(line[21:25]), int(line[33:37]) + 1)
            helices.append(
                Helix(
                    residues=filter(
                        lambda res: res.sequence_number in seq_range, residues
                    ),
                    sequence_number=seq_num,
                )
            )

        return helices

    def sheets(self, ignored: Iterable[str] = ("HOH", "LI1", "SQU")) -> Iterable[Sheet]:
        """Returns an iterable of sheets present in the protein.

        Parameters
        ----------
        ignored : iterable of str, optional
            3-letter string code for residues to ignore.

        Returns
        -------
        sheets : iterable of Sheet instances
        """
        residues = self.residues(ignored)
        sheets = list()
        self._handle.seek(0)

        # Helices are defined as a list of residues making them up
        # Therefore, we simply filter residues down to the ones in range
        for line in filter(lambda line: line.startswith("SHEET"), self._handle):
            seq_num = int(line[7:10])  # TODO: place in sheet
            seq_range = range(int(line[22:26]), int(line[33:37]) + 1)
            sheets.append(
                Sheet(
                    residues=filter(
                        lambda res: res.sequence_number in seq_range, residues
                    ),
                    sequence_number=seq_num,
                )
            )

        return sheets

    def secondary_structures(
        self, ignored: Iterable[str] = ("HOH", "LI1", "SQU")
    ) -> Iterable[SecondaryStructure]:
        """Iterable of all secondary structures present in this file.

        Parameters
        ----------
        ignored : iterable of str, optional
            3-letter string code for residues to ignore.

        Returns
        -------
        structures : iterable of SecondaryStructure instances
        """
        return list(self.sheets(ignored=ignored)) + list(self.helices(ignored=ignored))

    def residues(
        self, ignored: Iterable[str] = ("HOH", "LI1", "SQU")
    ) -> Iterable[Residue]:
        """
        Iterable of residues present in the structure.

        Parameters
        ----------
        ignored : iterable of str, optional
            3-letter string code for residues to ignore.

        Returns
        -------
        res : iterable of Residue instances
        """
        # Lattice vectors have to be determined first because
        # the file pointer is moved
        lattice_vectors = self.lattice_vectors()

        # Filter lines with start with ATM or HETATM
        is_atom_line = lambda l: l.startswith(("ATOM", "HETATM"))

        # ``residues`` is a dictionary mapping between (sequence numbers, name)
        # and an iterable of ``Atom``. When we have collected all atoms, we then create
        # a ``Residue`` for each sequence number
        residues = dict()

        self._handle.seek(0)
        for line in filter(is_atom_line, self._handle):
            residue_name = str(line[17:20]).replace(" ", "")
            if residue_name in ignored:
                continue

            residue_seq = int(line[22:26])
            if (residue_seq, residue_name) not in residues:
                residues[(residue_seq, residue_name)] = list()

            # TODO: include Atom ID record in Atom objects
            identification = str(line[12:16]).replace(" ", "")
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            coords_fractional = frac_coords([x, y, z], lattice_vectors)
            element = str(line[76:78]).replace(" ", "")

            try:
                occupancy = float(line[54:60])
            except ValueError:
                occupancy = None

            residues[(residue_seq, residue_name)].append(
                Atom(element=element, coords=coords_fractional, occupancy=occupancy)
            )

        if not residues:
            raise ParseError(f"No residues found in {self.filename}")

        structures = list()
        for (seq_number, name), atoms in residues.items():
            structures.append(
                Residue(atoms=atoms, name=name, sequence_number=seq_number)
            )
        return structures

    def atoms(self) -> Iterable[Atom]:
        """
        Returns a list of atoms associated with a PDB structure in fractional coordinates.
        These atoms form the asymmetric unit cell.

        Returns
        -------
        atoms: iterable of Atom
        """
        # Lattice vectors have to be determined first because
        # the file pointer is moved
        lattice_vectors = self.lattice_vectors()

        is_atom_line = lambda l: l.startswith(("ATOM", "HETATM"))

        self._handle.seek(0)
        atoms = list()
        for line in filter(is_atom_line, self._handle):
            element = str(line[76:78]).replace(" ", "").title()

            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            fractional_coordinates = frac_coords(np.array([x, y, z]), lattice_vectors)

            try:
                occupancy = float(line[54:60])
            except ValueError:
                occupancy = None

            atoms.append(
                Atom(
                    element=element, coords=fractional_coordinates, occupancy=occupancy
                )
            )
        return atoms

    def symmetry_operators(self) -> Iterable[np.ndarray]:
        """
        Returns the symmetry operators that map the atomic positions in a
        PDB file to the crystal unit cell.

        Returns
        -------
        sym_ops : iterable of `~numpy.ndarray`, shape (4,4)
            Transformation matrices. Since translations and rotation are combined,
            the transformation matrices are 4x4.
        """
        self._handle.seek(0)

        # This is clunky af
        sym_ops = dict()
        for line in filter(
            lambda l: l.startswith("REMARK 290") and ("SMTRY" in l), self._handle
        ):

            op_num = line[22:23]
            if op_num not in sym_ops:
                sym_ops[op_num] = {"rotation": list(), "translation": list()}

            r1, r2, r3, t = np.fromstring(line[23:], dtype=float, count=4, sep=" ")
            sym_ops[op_num]["rotation"].append([r1, r2, r3])
            sym_ops[op_num]["translation"].append(t)

        if not sym_ops:
            raise ParseError(
                f"No symmetry could be parsed from file {self._handle.filename}"
            )

        operators = list()
        for op in sym_ops.values():
            mat = np.eye(4, dtype=float)
            mat[:3, :3] = np.array(op["rotation"])
            mat[:3, 3] = np.array(op["translation"])
            operators.append(mat)

        return operators


class CIFParser(AbstractStructureParser):
    """
    Collection of methods that parses CIF files based on cif2cell. The preferred method
    of using this object is as a context manager.

    Parameters
    ----------
    filename : str or path-like
        Location of the CIF file.
    """

    def __init__(self, filename: PathLike, **kwargs):
        # ReadCIF would get confused between local files and URLs
        # Therefore, more clear to pass an open file
        self._handle = open(filename, mode="r")
        self.file = ReadCif(self._handle, **kwargs)

    def __exit__(self, *args, **kwargs):
        self._handle.close()

    @property
    def filename(self) -> str:
        return self._handle.name

    @staticmethod
    def sym_ops_from_equiv(equiv_site: Union[str, Iterable[str]]) -> np.ndarray:
        """Parse a symmetry operator from an equivalent-site representation

        Parameters
        ----------
        equiv_site : str or iterable of strings
            Either comma-separated string e.g. "+y, +x, -z + 1/2" or an
            iterable of the comma-separated values, e.g. ["+y", "+x", "-z + 1/2"]

        Returns
        -------
        sym_ops : ndarray, shape (4,4)
            Symmetry operator as a 4x4 affine transformation on the FRACTIONAL
            coordinates.
        """
        symmetry = np.zeros((3, 3))
        translation = np.zeros((3,))

        if isinstance(equiv_site, str):
            equiv_site = equiv_site.split(",")

        equiv_site = tuple(map(lambda s: s.strip().lower(), equiv_site))
        for j in range(3):
            xyz = equiv_site[j].replace("+", " +").replace("-", " -").split()
            for i in xyz:
                if i.strip("+-") == "x":
                    symmetry[0, j] = float(i.strip("x") + "1")
                elif i.strip("+-") == "y":
                    symmetry[1, j] = float(i.strip("y") + "1")
                elif i.strip("+-") == "z":
                    symmetry[2, j] = float(i.strip("z") + "1")

                if i.strip("+-xyz") != "":
                    translation[j] = eval(i)

        symmetry[:] = np.transpose(symmetry)

        # Combination of transform and translation into a single matrix
        # is done in a 4x4 affine transform
        symmetry_operation = affine_map(symmetry)
        symmetry_operation[:3, 3] = translation
        return symmetry_operation

    @property
    def structure_block(self):
        """Retrieve which CIF block has the appropriate structural information"""
        blocks = (self.file[key] for key in self.file.keys())
        for block in blocks:
            try:
                _, _ = get_number_with_esd(block["_cell_length_a"])
            except KeyError:
                continue
            else:
                return block

    @lru_cache(maxsize=1)
    def hall_symbol(self) -> str:
        """Returns the Hall symbol"""
        block = self.structure_block

        hall_symbol = block.get("_symmetry_space_group_name_Hall") or block.get(
            "_space_group_name_Hall"
        )

        # In some rare cases, the given hall symbol in the file isn't standard,
        # otherwise it would be a key in SymOpsHall
        # Then, it is preferable to infer the conventional hall symbol from other info
        if (hall_symbol is None) or (hall_symbol not in SymOpsHall):
            h_m_symbol = block.get("_symmetry_space_group_name_H-M") or block.get(
                "_space_group_name_H-M_alt"
            )

            if h_m_symbol is not None:
                h_m_symbol = re.sub(r"\s+", "", h_m_symbol)
                with suppress(
                    KeyError
                ):  # Symbol could be meaningless, e.g. h_m_symbol = '?' (True story)
                    hall_symbol = HM2Hall[h_m_symbol]

        # Again, if hall_symbol is still missing OR invalid
        if (hall_symbol is None) or (hall_symbol not in SymOpsHall):
            table_number = block.get("_symmetry_Int_Tables_number") or block.get(
                "_space_group_IT_number"
            )

            if table_number is not None:
                hall_symbol = Number2Hall[int(table_number)]

        if hall_symbol is None:
            raise ParseError("Hall symbol could not be inferred")

        if hall_symbol[0] == "-":
            hall_symbol = "-" + hall_symbol[1].upper() + hall_symbol[2:].lower()
        else:
            hall_symbol = hall_symbol[0].upper() + hall_symbol[1:].lower()

        return hall_symbol

    @lru_cache(maxsize=1)
    def lattice_parameters(self) -> Tuple[float, float, float, float, float, float]:
        """
        Returns the lattice parameters associated to a CIF structure.

        Returns
        ----------
        a, b, c : float
            Lengths of lattice vectors [Angstroms]
        alpha, beta, gamma : float
            Angles of lattice vectors [degrees].
        """
        block = self.structure_block

        try:
            a_with_err = block["_cell_length_a"]
        except KeyError:
            raise ParseError(f"No lattice information is present in {self.filename}")

        # In case where b and c are not listed, we use the value of a
        a, _ = get_number_with_esd(a_with_err)
        b, _ = get_number_with_esd(block.get("_cell_length_b", a_with_err))
        c, _ = get_number_with_esd(block.get("_cell_length_c", a_with_err))

        alpha, _ = get_number_with_esd(block["_cell_angle_alpha"])
        beta, _ = get_number_with_esd(block["_cell_angle_beta"])
        gamma, _ = get_number_with_esd(block["_cell_angle_gamma"])

        return a, b, c, alpha, beta, gamma

    @lru_cache(maxsize=1)
    def lattice_vectors(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns the lattice vectors associated to a CIF structure.

        Returns
        -------
        lv : list of ndarrays, shape (3,)
        """
        return Lattice.from_parameters(*self.lattice_parameters()).lattice_vectors

    def symmetry_operators(self) -> Iterable[np.ndarray]:
        """
        Returns the symmetry operators that map the fractional atomic positions in a
        CIF file to the crystal *conventional* unit cell.

        Returns
        -------
        sym_ops : iterable of ndarray, shape (4,4)
            Transformation matrices. Since translations and rotation are combined,
            the transformation matrices are 4x4.
        """
        block = self.structure_block

        equivalent_sites_str = None
        for tag in ["_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz"]:
            with suppress(KeyError):
                equivalent_sites_str = block.GetLoop(tag).get(tag)

        # P1 space group only has a single equivalent site
        if isinstance(equivalent_sites_str, str):
            equivalent_sites_str = [equivalent_sites_str]

        with suppress(ParseError):
            if not equivalent_sites_str:
                equivalent_sites_str = SymOpsHall[self.hall_symbol()]
            elif len(equivalent_sites_str) != len(SymOpsHall[self.hall_symbol()]):
                warnings.warn(
                    "The number of equivalent sites is not in line with the database. The file might be incomplete"
                )

        return list(map(self.sym_ops_from_equiv, equivalent_sites_str))

    def atoms(self) -> Iterable[Atom]:
        """
        Asymmetric unit cell. Combine with CIFParser.symmetry_operators() for a full unit cell.

        Returns
        -------
        atoms : iterable of Atom instance
        """
        block = self.structure_block

        try:
            tmpdata = block.GetLoop("_atom_site_fract_x")
            cartesian = False
        except:
            try:
                tmpdata = block.GetLoop("_atom_site_Cartn_x")
                cartesian = True
            except:
                raise ParseError("Atomic positions could not be found or inferred.")

            t11 = block.get("_atom_sites_Cartn_tran_matrix_11")
            t12 = block.get("_atom_sites_Cartn_tran_matrix_12")
            t13 = block.get("_atom_sites_Cartn_tran_matrix_13")
            t21 = block.get("_atom_sites_Cartn_tran_matrix_21")
            t22 = block.get("_atom_sites_Cartn_tran_matrix_22")
            t23 = block.get("_atom_sites_Cartn_tran_matrix_23")
            t31 = block.get("_atom_sites_Cartn_tran_matrix_13")
            t32 = block.get("_atom_sites_Cartn_tran_matrix_23")
            t33 = block.get("_atom_sites_Cartn_tran_matrix_33")
            cart_trans_matrix_inv = np.array(
                [
                    [float(t11), float(t12), float(t13)],
                    [float(t21), float(t22), float(t23)],
                    [float(t31), float(t32), float(t33)],
                ]
            )
            cart_trans_matrix = inv(cart_trans_matrix_inv)

            if not all([t11, t12, t13, t21, t22, t23, t31, t32, t33]):
                raise ParseError(
                    "Cartesian coordinates in CIF but no transformation matrix given"
                )

        if cartesian:
            xs = tmpdata.get("_atom_site_Cartn_x")
            ys = tmpdata.get("_atom_site_Cartn_y")
            zs = tmpdata.get("_atom_site_Cartn_z")
        else:
            xs = tmpdata.get("_atom_site_fract_x")
            ys = tmpdata.get("_atom_site_fract_y")
            zs = tmpdata.get("_atom_site_fract_z")
        # TODO: handle wildcards like '?', '.' in xs, ys, zs

        elements = tmpdata.get("_atom_site_type_symbol")
        if not elements:
            elements = tmpdata.get("_atom_site_label")
            if not elements:
                raise ParseError("Atom symbols could not be found or inferred.")
        elements = map(lambda s: s.strip(punctuation + digits).title(), elements)

        occupancies = tmpdata.get("_atom_site_occupancy")
        if not occupancies:
            # len(xs) == number of atoms. Therefore, good upper bound for `repeat`, which otherwise
            # might produce an infinitely long iterable
            occupancies = repeat(1.0, len(xs))
        occupancies = map(
            lambda x: get_number_with_esd(x)[0], occupancies
        )  # See issue #7

        atoms = list()
        for e, x, y, z, occ in zip(elements, xs, ys, zs, occupancies):
            coords = np.array(
                [
                    get_number_with_esd(x)[0],
                    get_number_with_esd(y)[0],
                    get_number_with_esd(z)[0],
                ]
            )

            # We normalize atom position to be within the unit cell
            # Therefore we need the fractional coordinates
            if cartesian:
                coords = transform(cart_trans_matrix, coords)
                coords[:] = frac_coords(coords, self.lattice_vectors())

            atoms.append(Atom(element=e, coords=np.mod(coords, 1), occupancy=occ))

        return atoms


class CODParser(CIFParser):
    """
    Collection of methods that parses CIF files retrieved from the Crystallography Open Database.
    The preferred method of using this object is as a context manager.

    Parameters
    ----------
    num : int
        COD identification number.
    revision : int or None, optional
        Revision number. If None (default), the latest revision is used.
    download_dir : path-like object or None, optional
        Directory where to save the CIF file.
    overwrite : bool, optional
        Whether or not to overwrite files in cache if they exist. If no revision
        number is provided, files will always be overwritten.

    Raises
    ------
    RuntimeError : If the file could not be downloaded from any of the mirrors.
    """

    # Database mirrors are made available
    # see http://wiki.crystallography.net/codmirrors/
    mirrors = (
        "http://www.crystallography.net/cod/",
        "http://cod.ibt.lt/cod",
        "http://qiserver.ugr.es/cod/",
    )

    def __init__(
        self,
        num: int,
        revision: Optional[int] = None,
        download_dir: Optional[PathLike] = None,
        overwrite: bool = False,
        **kwargs,
    ):
        if download_dir is None:
            download_dir = STRUCTURE_CACHE

        super().__init__(
            filename=self.download_cif(download_dir, num, revision, overwrite), **kwargs
        )

    @classmethod
    def download_cif(
        cls,
        download_dir: PathLike,
        num: int,
        revision: Optional[int] = None,
        overwrite: bool = False,
    ) -> Path:
        """
        Download a CIF file from the Crystallography Open Database.

        Parameters
        ----------
        download_dir : path-like object
            Directory where to save the CIF file.
        num : int
            COD identification number.
        revision : int or None, optional
            Revision number. If None (default), the latest revision is used.
        overwrite : bool, optional
            Whether or not to overwrite files in cache if they exist. If no revision
            number is provided, files will always be overwritten.

        Returns
        -------
        path : pathlib.Path
            Path to the downloaded file.

        Raises
        ------
        RuntimeError : If the file could not be downloaded from any of the mirrors.

        Notes
        -----
        This function will try three download mirrors. Warnings will be emitted in case
        the file cannot be found in a mirror.
        """
        download_dir = Path(download_dir)
        download_dir.mkdir(exist_ok=True)

        if revision is None:
            # We can never be sure of what the latest revision is
            # Therefore, to be safe, we overwrite
            overwrite = True

            # If latest revision, url should be e.g. (mirror)/1023891.cif
            url_suffix = f"{num}.cif"
            base = f"{num}.cif"
        else:
            # For revision, url should be e.g. (mirror)/9812812.cif@98181234
            url_suffix = str(num) + ".cif@" + str(revision)
            base = f"{num}-{revision}.cif"

        download_path = download_dir / base

        if download_path.exists() and (not overwrite):
            return download_path

        for index, mirror in enumerate(cls.mirrors, start=1):
            url = mirror + url_suffix
            try:
                urlretrieve(url, download_path)
            except URLError as e:
                warn(f"The file {url} could not be downloaded because: {e.reason}")

                # If this is the last mirror, then the file could not
                # be downloaded at all
                if index == len(cls.mirrors):
                    raise RuntimeError(
                        f"Crystallography Open Database ID {num} could not be downloaded from any mirror."
                    ) from None

                continue
            else:
                break

        return download_path


class MPJParser(CIFParser):
    """
    Collection of methods that parses CIF files retrieved from the Materials Project.
    The preferred method of using this object is as a context manager.

    Parameters
    ----------
    query : str
        The query can be a Materials Project material id (e.g., `"mp-1234"`), a
        formula, e.g. (`"Fe2O3"`), or a chemical system ("-" separated list of elemments,
        e.g., `"Li-Fe-O"`).
    api_key : str or None, optional
        An API key for accessing the Materials Project REST interface.
        Please obtain your API key at https://www.materialsproject.org/dashboard.
        If `None` (default), ``crystals`` will look for your API key in the
        `MATERIALS_PROJECT_API_KEY` environment variable.
    download_dir : path-like object or None, optional
        Directory where to save the CIF file. This is used for caching.
    overwrite : bool, optional
        Whether or not to overwrite files in cache if they exist. If True,
        a new file will be downloaded, possibly overwriting previously-downloaded file.

    Raises
    ------
    ValueError : if `api_key` is None and it could not be found in the environment variables.
    """

    def __init__(
        self,
        query: str,
        api_key: Optional[str] = None,
        download_dir: Optional[PathLike] = None,
        overwrite: bool = False,
        **kwargs,
    ):
        if download_dir is None:
            download_dir = STRUCTURE_CACHE

        if api_key is None:
            api_key = environ.get("MATERIALS_PROJECT_API_KEY", None)

            if api_key is None:
                raise ValueError(MISSING_MP_API_KEY_MESSAGE)
        super().__init__(
            filename=self.download_cif(
                api_key=api_key,
                query=query,
                download_dir=download_dir,
                overwrite=overwrite,
            ),
            **kwargs,
        )

    @property
    def material_id(self):
        """Returns the Materials Project material ID from this file."""
        # A comment of the form "Material ID: xxxxxxx" will have been inserted in the first
        # line of the file after download.
        self._handle.seek(0)
        firstline = next(iter(self._handle))
        *_, material_id = firstline.split(":")
        return material_id.strip()

    def download_cif(
        self,
        query: str,
        api_key: Optional[str] = None,
        download_dir: Optional[PathLike] = None,
        overwrite: bool = False,
    ) -> Path:
        """
        Download a CIF file from the Materials Project Database.

        Parameters
        ----------
        api_key : str
            An API key for accessing the Materials Project REST interface.
            Please obtain your API key at https://www.materialsproject.org/dashboard.
        query : str
            The query can be a Materials Project material id (e.g., `"mp-1234"`), a
            formula, e.g. (`"Fe2O3"`), or a chemical system ("-" separated list of elemments,
            e.g., `"Li-Fe-O"`).
        download_dir : path-like object or None, optional
            Directory where to save the CIF file. This is used for caching.
        overwrite : bool, optional
            Whether or not to overwrite files in cache if they exist. If True,
            a new file will be downloaded, possibly overwriting previously-downloaded file.

        Returns
        -------
        path : pathlib.Path
            Path to the downloaded file.

        Raises
        ------
        ConnectionError : If the file could not be downloaded.
        """
        download_dir = Path(download_dir)
        download_dir.mkdir(exist_ok=True)
        target_filename = download_dir / f"{query}.cif"

        if target_filename.exists() and (not overwrite):
            return target_filename

        endpoint = f"https://materialsproject.org/rest/v2/materials/{query}/vasp/cif"
        headers = {
            "x-api-key": api_key,
        }

        with requests.get(endpoint, headers=headers) as response:
            if response.status_code == 403:  # Forbidden access
                raise ConnectionError(
                    "Materials Project API key is not authorized to do this operation."
                    / f"status code {response.status_code}"
                )
            elif response.status_code != 200:
                raise ConnectionError(
                    f"Would not connect: status code {response.status_code}"
                )

            body = response.json()["response"][0]
            material_id = body["material_id"]

            with open(target_filename, "w") as f:
                # We write the material ID as a comment at the top of the CIF file
                # so we can inform the user in the Crystal.source property
                f.write(f"# Material ID: {material_id}\n")
                f.write(body["cif"])

        return target_filename


class PWSCFParser(AbstractStructureParser):
    """
    Collection of methods that parses output files from the Plane-Wave Self-Consistent
    Field (PWSCF) program, part of the Quantum Espresso suite.

    The preferred method of using this object is as a context manager.

    Parameters
    ----------
    filename : str or path-like
        Location of the CIF file.
    """

    # Regular expression pattern for a 3-vector
    # These are numbers separated by whitespace, in parentheses
    # Example:
    #   (  -0.0008701   0.5704561   0.4409210  )
    _vector_pattern = r"[(]\s* (?P<x1>[-]?[0-9]*\.[0-9]+\s*) (?P<x2>[-]?[0-9]*\.[0-9]+\s*) (?P<x3>[-]?[0-9]*\.[0-9]+\s*) [)]"

    # Conversion factor from bohr radius to angstroms
    _bohr_to_angs = 0.529_177_249

    def __init__(self, filename: PathLike, **kwargs):
        self.filename = filename

        with open(filename, mode="r") as f:
            self._filecontent = f.read()

    def __exit__(self, *args, **kwargs):
        pass

    @property
    def alat(self) -> float:
        """Get the lattice parameter [Bohr radius]"""
        match = re.search(
            r"\s*(lattice parameter [(]alat[)])\s*=\s*(?P<alat>\d+[.]\d+)\s*(a.u.)",
            self._filecontent,
        )
        if not match:
            raise ParseError(
                f"Lattice parameter from {self.filename} could not be parsed."
            )

        return float(match.group("alat"))

    @property
    def natoms(self) -> int:
        """Number of atoms defined per cell"""
        match = re.search(
            r"(\s*number of atoms/cell\s*)[=]\s*(?P<natoms>\d+)", self._filecontent
        )
        if not match:
            raise ParseError(
                f"Lattice parameter from {self.filename} could not be parsed."
            )

        return int(match.group("natoms"))

    @property
    def celldm(self) -> Tuple[float, float, float, float, float, float]:
        """Crystallographic constants as defined by INPUT_PW. They are, in order:

        * a,b,c in ANGSTROM

        * cosAB = cosine of the angle between axis a and b (gamma)

        * cosAC = cosine of the angle between axis a and c (beta)

        * cosBC = cosine of the angle between axis b and c (alpha)

        The constants are placed in a dictionary, indexed between 1 and 6.
        """
        params = dict()
        for index in (1, 2, 3, 4, 5, 6):
            pattern = (
                r"(celldm[(]"
                + str(index)
                + r"[)]=\s*)(?P<value>[-]?[0-9]*\.[0-9]+)(\s*)"
            )
            match = re.search(pattern, self._filecontent).group("value")
            params[index] = float(match)

        return params

    def lattice_vectors_alat(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns the lattice vectors associated to a structure. These lattice vectors are in
        units of `lattice parameters` [alat].

        Returns
        -------
        lv : iterable of ndarrays, shape (3,)
        """
        a1 = re.search(
            r"(\s*a[(]1[)]\s*=\s*)" + self._vector_pattern, self._filecontent
        ).group("x1", "x2", "x3")
        a2 = re.search(
            r"(\s*a[(]2[)]\s*=\s*)" + self._vector_pattern, self._filecontent
        ).group("x1", "x2", "x3")
        a3 = re.search(
            r"(\s*a[(]3[)]\s*=\s*)" + self._vector_pattern, self._filecontent
        ).group("x1", "x2", "x3")

        return tuple(np.array(tuple(map(float, a))) for a in (a1, a2, a3))

    def reciprocal_vectors_alat(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns the reciprocal lattice vectors associated to a structure, in units of :math:`2 \\pi / alat`.

        Returns
        -------
        lv : iterable of ndarrays, shape (3,)
        """
        b1 = re.search(
            r"(\s*b[(]1[)]\s*=\s*)" + self._vector_pattern, self._filecontent
        ).group("x1", "x2", "x3")
        b2 = re.search(
            r"(\s*b[(]2[)]\s*=\s*)" + self._vector_pattern, self._filecontent
        ).group("x1", "x2", "x3")
        b3 = re.search(
            r"(\s*b[(]3[)]\s*=\s*)" + self._vector_pattern, self._filecontent
        ).group("x1", "x2", "x3")

        return tuple(np.array(tuple(map(float, b))) for b in (b1, b2, b3))

    def lattice_vectors(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns the lattice vectors associated to a structure [:math:`Å`].

        Returns
        -------
        lv : iterable of ndarrays, shape (3,)
        """
        vectors = np.array(self.lattice_vectors_alat())
        scale = self._bohr_to_angs * self.alat * np.eye(3)
        return tuple(map(np.squeeze, np.vsplit(scale @ vectors, 3)))

    def atoms(self) -> Iterable[Atom]:
        """
        Asymmetric unit cell.

        Returns
        -------
        atoms : iterable of Atom instance
        """
        atoms = list()
        for index in range(1, self.natoms + 1):
            pattern = (
                r"\s*"
                + str(index)
                + r"\s*(?P<element>[A-Z][a-z]?) (\s* tau[(]\s*(?P<atm_index>\d+)[)]\s*=\s*)"
                + self._vector_pattern
            )
            match = re.search(pattern, self._filecontent)
            coords = np.asarray(tuple(map(float, match.group(4, 5, 6))))

            # Convert coordinates in fractional coordinates
            # Note that the atomic coordinates are not quite real-space, but rather
            # real space in units of lattice_parameter (alat) in bohr radius
            atoms.append(
                Atom(
                    element=match.group("element"),
                    coords=frac_coords(coords, self.lattice_vectors_alat()),
                    tag=int(match.group("atm_index")),
                )
            )

        return atoms


class POSCARParser(AbstractStructureParser):
    """
    Collection of methods that parses POSCAR output files from the Vienna Ab initio
    Simulation Package (VASP) suite.

    The preferred method of using this object is as a context manager.

    Parameters
    ----------
    filename : str or path-like
        Location of the POSCAR file.
    """

    def __init__(self, filename: PathLike, **kwargs):
        self.filename = filename

        with open(filename, mode="r") as f:

            next(f)

            scaling_factor = float(next(f))

            self._lattice_vectors = (
                np.array([next(f).split() for _ in range(3)]).astype(float)
                * scaling_factor
            )
            self._atom_types = list(
                zip(
                    next(f).strip().split(),
                    map(int, next(f).strip().split()),
                )
            )
            self._atoms = []
            flag = next(f)

            if flag.startswith("S"):
                raise NotImplementedError(
                    "Selective dynamics tag in POSCAR files are not supported."
                )
            elif flag[0] in ["C", "c", "K", "k"]:
                for element, nat in self._atom_types:
                    for _ in range(nat):
                        coords = np.array(next(f).strip().split()[:3]).astype(float)
                        coords *= scaling_factor
                        coords = coords @ np.linalg.inv(self._lattice_vectors)
                        self._atoms.append(Atom(element, coords))
            else:
                for element, nat in self._atom_types:
                    for _ in range(nat):
                        coords = np.array(next(f).strip().split()[:3]).astype(float)
                        self._atoms.append(Atom(element, coords))

    def __exit__(self, *args, **kwargs):
        pass

    def atoms(self) -> Iterable[Atom]:
        """
        Unit cell.

        Returns
        -------
        atoms : iterable of Atom instance
        """
        return self._atoms

    def lattice_vectors(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns the lattice vectors associated to a POSCAR structure.

        Returns
        -------
        lv : list of ndarrays, shape (3,)
        """
        return self._lattice_vectors
