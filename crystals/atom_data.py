# -*- coding: utf-8 -*-
"""
Atomic data from the `Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase/index.html>`_.

References
----------
  .. [#] Meija, J., Coplen, T., Berglund, M., et al. (2016). Atomic weights of 
    the elements 2013 (IUPAC Technical Report). Pure and Applied Chemistry,
    88(3), pp. 265-291. Retrieved 30 Nov. 2016, doi:10.1515/pac-2015-0305
"""
from math import nan

chemical_symbols = [
    # 1
    "H",
    "He",
    # 2
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    # 3
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    # 4
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    # 5
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    # 6
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    # 7
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]

# IUPAC version dated 28 November 2016
atomic_names = [
    "Hydrogen",
    "Helium",
    "Lithium",
    "Beryllium",
    "Boron",
    "Carbon",
    "Nitrogen",
    "Oxygen",
    "Fluorine",
    "Neon",
    "Sodium",
    "Magnesium",
    "Aluminium",
    "Silicon",
    "Phosphorus",
    "Sulfur",
    "Chlorine",
    "Argon",
    "Potassium",
    "Calcium",
    "Scandium",
    "Titanium",
    "Vanadium",
    "Chromium",
    "Manganese",
    "Iron",
    "Cobalt",
    "Nickel",
    "Copper",
    "Zinc",
    "Gallium",
    "Germanium",
    "Arsenic",
    "Selenium",
    "Bromine",
    "Krypton",
    "Rubidium",
    "Strontium",
    "Yttrium",
    "Zirconium",
    "Niobium",
    "Molybdenum",
    "Technetium",
    "Ruthenium",
    "Rhodium",
    "Palladium",
    "Silver",
    "Cadmium",
    "Indium",
    "Tin",
    "Antimony",
    "Tellurium",
    "Iodine",
    "Xenon",
    "Caesium",
    "Barium",
    "Lanthanum",
    "Cerium",
    "Praseodymium",
    "Neodymium",
    "Promethium",
    "Samarium",
    "Europium",
    "Gadolinium",
    "Terbium",
    "Dysprosium",
    "Holmium",
    "Erbium",
    "Thulium",
    "Ytterbium",
    "Lutetium",
    "Hafnium",
    "Tantalum",
    "Tungsten",
    "Rhenium",
    "Osmium",
    "Iridium",
    "Platinum",
    "Gold",
    "Mercury",
    "Thallium",
    "Lead",
    "Bismuth",
    "Polonium",
    "Astatine",
    "Radon",
    "Francium",
    "Radium",
    "Actinium",
    "Thorium",
    "Protactinium",
    "Uranium",
    "Neptunium",
    "Plutonium",
    "Americium",
    "Curium",
    "Berkelium",
    "Californium",
    "Einsteinium",
    "Fermium",
    "Mendelevium",
    "Nobelium",
    "Lawrencium",
    "Rutherfordium",
    "Dubnium",
    "Seaborgium",
    "Bohrium",
    "Hassium",
    "Meitnerium",
    "Darmastadtium",
    "Roentgenium",
    "Copernicium",
    "Nihonium",
    "Flerovium",
    "Moscovium",
    "Livermorium",
    "Tennessine",
    "Oganesson",
]

NUM_TO_ELEM = dict(enumerate(chemical_symbols, 1))
ELEM_TO_NAME = dict(zip(chemical_symbols, atomic_names))
NAME_TO_ELEM = dict(zip(atomic_names, chemical_symbols))
ELEM_TO_NUM = {v: k for k, v in NUM_TO_ELEM.items()}

# Atomic masses are based on:
#
#   Meija, J., Coplen, T., Berglund, M., et al. (2016). Atomic weights of
#   the elements 2013 (IUPAC Technical Report). Pure and Applied Chemistry,
#   88(3), pp. 265-291. Retrieved 30 Nov. 2016,
#   from doi:10.1515/pac-2015-0305
#
# Standard atomic weights are taken from Table 1: "Standard atomic weights
# 2013", with the uncertainties ignored.
# For hydrogen, helium, boron, carbon, nitrogen, oxygen, magnesium, silicon,
# sulfur, chlorine, bromine and thallium, where the weights are given as a
# range the "conventional" weights are taken from Table 3 and the ranges are
# given in the comments.
# The mass of the most stable isotope (in Table 4) is used for elements
# where there the element has no stable isotopes (to avoid NaNs): Tc, Pm,
# Po, At, Rn, Fr, Ra, Ac, everything after Np
atomic_masses = [
    1.008,  # H [1.00784, 1.00811]
    4.002602,  # He
    6.94,  # Li [6.938, 6.997]
    9.0121831,  # Be
    10.81,  # B [10.806, 10.821]
    12.011,  # C [12.0096, 12.0116]
    14.007,  # N [14.00643, 14.00728]
    15.999,  # O [15.99903, 15.99977]
    18.998403163,  # F
    20.1797,  # Ne
    22.98976928,  # Na
    24.305,  # Mg [24.304, 24.307]
    26.9815385,  # Al
    28.085,  # Si [28.084, 28.086]
    30.973761998,  # P
    32.06,  # S [32.059, 32.076]
    35.45,  # Cl [35.446, 35.457]
    39.948,  # Ar
    39.0983,  # K
    40.078,  # Ca
    44.955908,  # Sc
    47.867,  # Ti
    50.9415,  # V
    51.9961,  # Cr
    54.938044,  # Mn
    55.845,  # Fe
    58.933194,  # Co
    58.6934,  # Ni
    63.546,  # Cu
    65.38,  # Zn
    69.723,  # Ga
    72.630,  # Ge
    74.921595,  # As
    78.971,  # Se
    79.904,  # Br [79.901, 79.907]
    83.798,  # Kr
    85.4678,  # Rb
    87.62,  # Sr
    88.90584,  # Y
    91.224,  # Zr
    92.90637,  # Nb
    95.95,  # Mo
    97.90721,  # 98Tc
    101.07,  # Ru
    102.90550,  # Rh
    106.42,  # Pd
    107.8682,  # Ag
    112.414,  # Cd
    114.818,  # In
    118.710,  # Sn
    121.760,  # Sb
    127.60,  # Te
    126.90447,  # I
    131.293,  # Xe
    132.90545196,  # Cs
    137.327,  # Ba
    138.90547,  # La
    140.116,  # Ce
    140.90766,  # Pr
    144.242,  # Nd
    144.91276,  # 145Pm
    150.36,  # Sm
    151.964,  # Eu
    157.25,  # Gd
    158.92535,  # Tb
    162.500,  # Dy
    164.93033,  # Ho
    167.259,  # Er
    168.93422,  # Tm
    173.054,  # Yb
    174.9668,  # Lu
    178.49,  # Hf
    180.94788,  # Ta
    183.84,  # W
    186.207,  # Re
    190.23,  # Os
    192.217,  # Ir
    195.084,  # Pt
    196.966569,  # Au
    200.592,  # Hg
    204.38,  # Tl [204.382, 204.385]
    207.2,  # Pb
    208.98040,  # Bi
    208.98243,  # 209Po
    209.98715,  # 210At
    222.01758,  # 222Rn
    223.01974,  # 223Fr
    226.02541,  # 226Ra
    227.02775,  # 227Ac
    232.0377,  # Th
    231.03588,  # Pa
    238.02891,  # U
    237.04817,  # 237Np
    244.06421,  # 244Pu
    243.06138,  # 243Am
    247.07035,  # 247Cm
    247.07031,  # 247Bk
    251.07959,  # 251Cf
    252.0830,  # 252Es
    257.09511,  # 257Fm
    258.09843,  # 258Md
    259.1010,  # 259No
    262.110,  # 262Lr
    267.122,  # 267Rf
    268.126,  # 268Db
    271.134,  # 271Sg
    270.133,  # 270Bh
    269.1338,  # 269Hs
    278.156,  # 278Mt
    281.165,  # 281Ds
    281.166,  # 281Rg
    285.177,  # 285Cn
    286.182,  # 286Nh
    289.190,  # 289Fl
    289.194,  # 289Mc
    293.204,  # 293Lv
    293.208,  # 293Ts
    294.214,  # 294Og
]


# set atomic_masses to most recent version
ELEM_TO_MASS = dict(zip(chemical_symbols, atomic_masses))

# http://www.webelements.com
ground_state_magnetic_moments = [
    1.0,  # H
    0.0,  # He
    1.0,  # Li
    0.0,  # Be
    1.0,  # B
    2.0,  # C
    3.0,  # N
    2.0,  # O
    1.0,  # F
    0.0,  # Ne
    1.0,  # Na
    0.0,  # Mg
    1.0,  # Al
    2.0,  # Si
    3.0,  # P
    2.0,  # S
    1.0,  # Cl
    0.0,  # Ar
    1.0,  # K
    0.0,  # Ca
    1.0,  # Sc
    2.0,  # Ti
    3.0,  # V
    6.0,  # Cr
    5.0,  # Mn
    4.0,  # Fe
    3.0,  # Co
    2.0,  # Ni
    1.0,  # Cu
    0.0,  # Zn
    1.0,  # Ga
    2.0,  # Ge
    3.0,  # As
    2.0,  # Se
    1.0,  # Br
    0.0,  # Kr
    1.0,  # Rb
    0.0,  # Sr
    1.0,  # Y
    2.0,  # Zr
    5.0,  # Nb
    6.0,  # Mo
    5.0,  # Tc
    4.0,  # Ru
    3.0,  # Rh
    0.0,  # Pd
    1.0,  # Ag
    0.0,  # Cd
    1.0,  # In
    2.0,  # Sn
    3.0,  # Sb
    2.0,  # Te
    1.0,  # I
    0.0,  # Xe
    1.0,  # Cs
    0.0,  # Ba
    1.0,  # La
    1.0,  # Ce
    3.0,  # Pr
    4.0,  # Nd
    5.0,  # Pm
    6.0,  # Sm
    7.0,  # Eu
    8.0,  # Gd
    5.0,  # Tb
    4.0,  # Dy
    3.0,  # Ho
    2.0,  # Er
    1.0,  # Tm
    0.0,  # Yb
    1.0,  # Lu
    2.0,  # Hf
    3.0,  # Ta
    4.0,  # W
    5.0,  # Re
    4.0,  # Os
    3.0,  # Ir
    2.0,  # Pt
    1.0,  # Au
    0.0,  # Hg
    1.0,  # Tl
    2.0,  # Pb
    3.0,  # Bi
    2.0,  # Po
    1.0,  # At
    0.0,  # Rn
    1.0,  # Fr
    0.0,  # Ra
    1.0,  # Ac
    2.0,  # Th
    3.0,  # Pa
    4.0,  # U
    5.0,  # Np
    6.0,  # Pu
    7.0,  # Am
    8.0,  # Cm
    5.0,  # Bk
    4.0,  # Cf
    4.0,  # Es
    2.0,  # Fm
    1.0,  # Md
    0.0,  # No
    nan,  # Lr
]

ELEM_TO_MAGMOM = dict(zip(chemical_symbols, ground_state_magnetic_moments))

# from Cordero et al, 'Covalent radii revisited' (2008).
# Dalton Transactions. DOI: 10.1039/B801115J
# TODO: implement bond perception based on something like this:
#       https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3557220/
covalent_radii = {
    "H": 0.31,
    "He": 0.28,
    "Li": 1.28,
    "Be": 0.96,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Ne": 0.58,
    "Na": 1.66,
    "Mg": 1.41,
    "Al": 1.21,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Ar": 1.06,
    "K": 2.03,
    "Ca": 1.76,
    "Sc": 1.7,
    "Ti": 1.6,
    "V": 1.53,
    "Cr": 1.39,
    "Mn": 1.61,
    "Fe": 1.52,
    "Co": 1.5,
    "Ni": 1.24,
    "Cu": 1.32,
    "Zn": 1.22,
    "Ga": 1.22,
    "Ge": 1.2,
    "As": 1.19,
    "Se": 1.2,
    "Br": 1.2,
    "Kr": 1.16,
    "Rb": 2.2,
    "Sr": 1.95,
    "Y": 1.9,
    "Zr": 1.75,
    "Nb": 1.64,
    "Mo": 1.54,
    "Tc": 1.47,
    "Ru": 1.46,
    "Rh": 1.42,
    "Pd": 1.39,
    "Ag": 1.45,
    "Cd": 1.44,
    "In": 1.42,
    "Sn": 1.39,
    "Sb": 1.39,
    "Te": 1.38,
    "I": 1.39,
    "Xe": 1.4,
    "Cs": 2.44,
    "Ba": 2.15,
    "La": 2.07,
    "Ce": 2.04,
    "Pr": 2.03,
    "Nd": 2.01,
    "Pm": 1.99,
    "Sm": 1.98,
    "Eu": 1.98,
    "Gd": 1.96,
    "Tb": 1.94,
    "Dy": 1.92,
    "Ho": 1.92,
    "Er": 1.89,
    "Tm": 1.9,
    "Yb": 1.87,
    "Lu": 1.87,
    "Hf": 1.75,
    "Ta": 1.7,
    "W": 1.62,
    "Re": 1.51,
    "Os": 1.44,
    "Ir": 1.41,
    "Pt": 1.36,
    "Au": 1.36,
    "Hg": 1.32,
    "Tl": 1.45,
    "Pb": 1.46,
    "Bi": 1.48,
    "Po": 1.4,
    "At": 1.5,
    "Rn": 1.5,
    "Fr": 2.6,
    "Ra": 2.21,
    "Ac": 2.15,
    "Th": 2.06,
    "Pa": 2.0,
    "U": 1.96,
    "Np": 1.9,
    "Pu": 1.87,
    "Am": 1.8,
    "Cm": 1.69,
}
