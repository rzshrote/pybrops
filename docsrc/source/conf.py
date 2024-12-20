# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
import pathlib
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())

def setup(app):
    app.add_css_file('wider_borders.css')

# -- Project information -----------------------------------------------------

project = 'pybrops'
copyright = '2024, Robert Z. Shrote'
author = 'Robert Z. Shrote'

# The full version, including alpha/beta/rc tags
release = '1.0.2'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.duration',      # Get build times
    'sphinx.ext.githubpages',   # Add .noj
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints', # Automatically document param types (less noise in class signature)
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# Will change to `root_doc` in Sphinx 4
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for sphinx.ext.autosummary --------------------------------------
autosummary_generate = True
# autosummary_imported_members = True
autoclass_content = "both"  # Add __init__ doc (ie. params) to class summaries
autosummary_generate = True  # Turn on sphinx.ext.autosummary
autoclass_content = "both"  # Add __init__ doc (ie. params) to class summaries
html_show_sourcelink = False  # Remove 'view source code' from top of page (for html, not python)
autodoc_inherit_docstrings = True  # If no docstring, inherit from base class
set_type_checking_flag = True  # Enable 'expensive' imports for sphinx_autodoc_typehints
nbsphinx_allow_errors = True  # Continue through Jupyter errors
# autodoc_typehints = "description" # Sphinx-native method. Not as good as sphinx_autodoc_typehints
# add_module_names = False # Remove namespaces from class/method signatures
add_module_names = True # Remove namespaces from class/method signatures

# -- Options for sphinx.ext.napoleon -----------------------------------------
napoleon_numpy_docstring = True

# -- Options for sphinx.ext.viewcode -----------------------------------------

# -- Options for autoapi.extension -------------------------------------------
# autoapi_dirs = ['../../src']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = "pydata_sphinx_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
