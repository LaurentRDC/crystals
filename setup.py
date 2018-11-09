# -*- encoding: utf-8 -*-

import re
from glob import glob
from itertools import chain
from pathlib import Path

from setuptools import find_packages, setup

PACKAGE_NAME = "crystals"
DESCRIPTION = "Data structures for crystallography"
URL = "http://crystals.readthedocs.io"
DOWNLOAD_URL = "http://github.com/LaurentRDC/crystals"
AUTHOR = "Laurent P. René de Cotret"
AUTHOR_EMAIL = "laurent.renedecotret@mail.mcgill.ca"
BASE_PACKAGE = "crystals"

CIF_FILES = chain.from_iterable([glob("crystals\\cifs\\*.cif")])

base_path = Path(__file__).parent
with open(base_path / BASE_PACKAGE / "__init__.py") as f:
    module_content = f.read()
    VERSION = (
        re.compile(r".*__version__ = \"(.*?)\"", re.S).match(module_content).group(1)
    )
    LICENSE = (
        re.compile(r".*__license__ = \"(.*?)\"", re.S).match(module_content).group(1)
    )

with open("README.rst") as f:
    README = f.read()

with open("requirements.txt") as f:
    REQUIREMENTS = [line for line in f.read().split("\n") if len(line.strip())]

if __name__ == "__main__":
    setup(
        name=PACKAGE_NAME,
        description=DESCRIPTION,
        long_description=README,
        license=LICENSE,
        url=URL,
        download_url=DOWNLOAD_URL,
        version=VERSION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        maintainer=AUTHOR,
        maintainer_email=AUTHOR_EMAIL,
        install_requires=REQUIREMENTS,
        keywords=["crystallography"],
        project_urls={
            "Documentation": "https://crystals.readthedocs.io/",
            "Source": "https://github.com/LaurentRDC/crystals",
        },
        python_requires=">=3.6",
        packages=find_packages(PACKAGE_NAME),
        data_files=[("crystals\\cifs", CIF_FILES)],
        include_package_data=True,
        zip_safe=False,
        # list of possible classifiers:
        #  https://pypi.python.org/pypi?%3Aaction=list_classifiers
        classifiers=[
            "Environment :: Console",
            "Development Status :: 4 - Beta",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    )
