# -*- encoding: utf-8 -*-

import os
import platform
import re
import shutil
import tempfile
from distutils.errors import CompileError, LinkError
from glob import glob
from itertools import chain
from pathlib import Path

import numpy
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext

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

# Support for openmp
# This idea is from scikit-image:
# https://github.com/scikit-image/scikit-image/
class BuildExtWithOpenMP(build_ext):
    """
    Try to compile extensions with OpenMP if possible
    """

    def can_compile_link(self, compile_flags, link_flags):
        cc = self.compiler
        fname = "test.c"
        cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp()

        code = "#include <omp.h>" "int main(int argc, char** argv) { return(0); }"

        if self.compiler.compiler_type == "msvc":
            # make sure we build a DLL on Windows
            local_link_flags = link_flags + ["/DLL"]
        else:
            local_link_flags = link_flags

        try:
            os.chdir(tmpdir)
            with open(fname, "wt") as fobj:
                fobj.write(code)
            try:
                objects = cc.compile([fname], extra_postargs=compile_flags)
            except CompileError:
                return False
            try:
                # Link shared lib rather then executable to avoid
                # http://bugs.python.org/issue4431 with MSVC 10+
                cc.link_shared_lib(objects, "testlib", extra_postargs=local_link_flags)
            except (LinkError, TypeError):
                return False
        finally:
            os.chdir(cwd)
            shutil.rmtree(tmpdir)
        return True

    def build_extensions(self):
        """ Hook into extension building to set compiler flags """
        if self.compiler.compiler_type == "msvc":
            compile_flags = ["/openmp"]
            link_flags = []
        else:
            compile_flags = ["-fopenmp"]
            link_flags = ["-fopenmp"]

        if self.can_compile_link(compile_flags, link_flags):
            for ext in self.extensions:
                ext.extra_compile_args += compile_flags
                ext.extra_link_args += link_flags

        return super().build_extensions()


PINKINDEXER = Path(".") / "pinkindexer"
pinkindexer_ext = Extension(
    "crystals.indexing._pinkindexer",
    include_dirs=[
        numpy.get_include(),
        PINKINDEXER / "src",
        PINKINDEXER / "include",
        PINKINDEXER / "include" / "Eigen",
    ],
    sources=["crystals/indexing/_pinkindexer.cpp"]
    + [str(p) for p in (PINKINDEXER / "src").glob("*.cpp")],
    extra_compile_args=["-std=c++11"] if platform.system() != "Windows" else [],
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
        cmdclass=dict(build_ext=BuildExtWithOpenMP),
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
