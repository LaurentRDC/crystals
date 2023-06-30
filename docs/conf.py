# -*- coding: utf-8 -*-

import os
import sys
from datetime import datetime

currentpath = os.path.dirname(__file__)
sys.path.append(os.path.join(currentpath, ".."))

import crystals

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
]

source_suffix = ".rst"
master_doc = "index"
project = "crystals"
year = datetime.now().year
author = crystals.__author__
copyright = f"{year}, {author}"
version = release = crystals.__version__

pygments_style = "sphinx"
templates_path = ["."]
extlinks = {
    "issue": ("https://github.com/LaurentRDC/crystals/issues/%s", "#"),
    "pr": ("https://github.com/LaurentRDC/crystals/pull/%s", "PR #"),
}
# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get("READTHEDOCS", None) == "True"

if not on_rtd:  # only set the theme if we're building docs locally
    html_theme = "sphinx_rtd_theme"

html_use_smartypants = True
html_last_updated_fmt = "%b %d, %Y"
html_split_index = False
html_sidebars = {"**": ["searchbox.html", "globaltoc.html", "sourcelink.html"]}
html_short_title = "%s-%s" % (project, version)

napoleon_use_ivar = True
napoleon_google_docstring = False
napoleon_use_rtype = False
napoleon_use_param = False

autosummary_generate = True

autodoc_default_flags = ["members"]
autoclass_content = "both"
autodoc_member_order = "groupwise"

doctest_global_setup = """
from crystals import Crystal, Atom, ElectronicStructure, Orbital, is_element
"""


def autodoc_skip_member(app, what, name, obj, skip, options):
    exclusions = {"__weakref__", "__doc__", "__module__", "__dict__"}
    exclude = name in exclusions
    return skip or exclude


def setup(app):
    app.connect("autodoc-skip-member", autodoc_skip_member)
