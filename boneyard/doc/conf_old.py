# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

##############################################################################
# -- Path setup --------------------------------------------------------------
##############################################################################

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
import pathlib
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())


##############################################################################
# -- Project information -----------------------------------------------------
##############################################################################

project = 'PyBrOpS'
copyright = '2022, Robert Shrote'
author = 'Robert Shrote'

# The full version, including alpha/beta/rc tags
release = '1.0.0'


##############################################################################
# -- General configuration ---------------------------------------------------
##############################################################################

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',  # Core library for html generation from docstrings
    'sphinx.ext.autosummary',  # Create summary tables
    'sphinx.ext.duration',      # Get build times
    'sphinx.ext.githubpages',   # Add .noj
    'sphinx.ext.napoleon',  # Convert NumPy style comments to Restructured Text
    'sphinx.ext.viewcode',   # Create links to source code
    'sphinx_autodoc_typehints', # Automatically document param types (less noise in class signature)
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# Will change to `root_doc` in Sphinx 4
master_doc = 'index'

# A boolean that decides whether module names are prepended to all object names.
add_module_names = False # Remove namespaces from class/method signatures

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

##############################################################################
# -- Options for sphinx.ext.autodoc ------------------------------------------
##############################################################################

# autoclass_content = "class"  # only use class docstring for class summary
autoclass_content = "both"  # use class + __init__ docstrings for class summary
# autoclass_content = "init"  # only use __init__ docstring for class summary

autodoc_inherit_docstrings = True  # If no docstring, inherit from base class
# autodoc_typehints = "description" # Sphinx-native method. Not as good as sphinx_autodoc_typehints
# autodoc_class_signature = "separated"


##############################################################################
# -- Options for sphinx.ext.autosummary --------------------------------------
##############################################################################
autosummary_generate = True  # Turn on sphinx.ext.autosummary
autosummary_imported_members = True  # document classes and functions imported in modules


##############################################################################
# -- Options for sphinx.ext.napoleon -----------------------------------------
##############################################################################
# napoleon_google_docstring = False
napoleon_numpy_docstring = True
# napoleon_include_init_with_doc = False
# napoleon_include_private_with_doc = False
# napoleon_include_special_with_doc = True
# napoleon_use_admonition_for_examples = False
# napoleon_use_admonition_for_notes = False
# napoleon_use_admonition_for_references = False
# napoleon_use_ivar = False
# napoleon_use_param = True
# napoleon_use_rtype = True
# napoleon_preprocess_types = False
# napoleon_type_aliases = None
# napoleon_attr_annotations = True

##############################################################################
# -- Options for sphinx_autodoc_typehints ------------------------------------
##############################################################################
set_type_checking_flag = True  # Enable 'expensive' imports for sphinx_autodoc_typehints

##############################################################################
# -- Options for HTML output -------------------------------------------------
##############################################################################

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_show_sourcelink = False  # Remove 'view source code' from top of page (for html, not python)
html_theme = "pydata_sphinx_theme"
# html_theme = 'bizstyle'
# html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
