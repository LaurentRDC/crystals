# -*- encoding: utf-8 -*-

import platform
import re
from glob import glob
from itertools import chain
from pathlib import Path
import numpy

from setuptools import Extension, find_packages, setup

PACKAGE_NAME = "crystals"
DESCRIPTION = "Data structures for crystallography"
URL = "http://crystals.readthedocs.io"
DOWNLOAD_URL = "http://github.com/LaurentRDC/crystals"
AUTHOR = "Laurent P. René de Cotret"
AUTHOR_EMAIL = "laurent.renedecotret@mail.mcgill.ca"

CIF_FILES = chain.from_iterable([glob("crystals\\cifs\\*.cif")])

base_path = Path(__file__).parent
with open(base_path / PACKAGE_NAME / "__init__.py") as f:
    module_content = f.read()
    VERSION = (
        re.compile(r".*__version__ = \"(.*?)\"", re.S).match(module_content).group(1)
    )
    LICENSE = (
        re.compile(r".*__license__ = \"(.*?)\"", re.S).match(module_content).group(1)
    )

with open("README.md", encoding="utf-8") as f:
    README = f.read()

with open("requirements.txt") as f:
    REQUIREMENTS = [line for line in f.read().split("\n") if len(line.strip())]

ROOT = Path(".") / "pinkindexer"
GCC_COMPILE_ARGS = ["-std=c++11"] if platform.system() != "Windows" else []
pinkindexer_ext = Extension(
    "crystals.indexing._pinkindexer",
    include_dirs=[numpy.get_include()]
    + [ROOT / "src", ROOT / "include", ROOT / "include" / "Eigen"],
    sources=["crystals/indexing/_pinkindexer.cpp"]
    + [str(p) for p in (ROOT / "src").glob("*.cpp")],
    extra_compile_args=[] + GCC_COMPILE_ARGS,
)


if __name__ == "__main__":
    setup(
        name=PACKAGE_NAME,
        description=DESCRIPTION,
        long_description=README,
        long_description_content_type="text/markdown",
        license=LICENSE,
        url=URL,
        download_url=DOWNLOAD_URL,
        version=VERSION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        maintainer=AUTHOR,
        maintainer_email=AUTHOR_EMAIL,
        install_requires=REQUIREMENTS,
        keywords=["crystallography", "material science", "structural biology"],
        project_urls={
            "Documentation": "https://crystals.readthedocs.io/",
            "Source": "https://github.com/LaurentRDC/crystals",
        },
        python_requires=">=3.7",
        packages=find_packages(),
        data_files=[("crystals\\cifs", CIF_FILES)],
        include_package_data=True,
        ext_modules=[pinkindexer_ext],
        zip_safe=False,
        entry_points={"console_scripts": ["crystals = crystals.__main__:main"]},
        # list of possible classifiers:
        #  https://pypi.python.org/pypi?%3Aaction=list_classifiers
        classifiers=[
            "Environment :: Console",
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    )
