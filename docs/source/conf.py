# Configuration file for the Sphinx documentation builder.

# -- Project information
import os, sys
import subprocess

sys.path.insert(
    0, os.path.abspath("../../bjet_mcmc")
)  # Source code dir relative to this file
sys.path.insert(
    0, os.path.abspath("../../bjet_core")
)  # Source code dir relative to this file

release = "0.2"
version = "0.2.1"

# -- General configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "breathe",
]

# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get("READTHEDOCS", None) == "True"

if read_the_docs_build:
    print("Building on Read the Docs")
    dox_cmd = "doxygen Doxyfile"
    subprocess.run(dox_cmd, shell=True)
    breathe_projects = {"bjet_core": "../build/doxygen/xml"}
    breathe_default_project = "bjet_core"
    breathe_default_members = ("members", "undoc-members")

autosummary_generate = True  # Turn on sphinx.ext.autosummary

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]
pygments_style = "sphinx"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "../../logo/Bjet_MCMC_logo_small_v2.png"
html_theme_options = {
    "logo_only": True,
    "display_version": False,
}

# -- Options for EPUB output
epub_show_urls = "footnote"
