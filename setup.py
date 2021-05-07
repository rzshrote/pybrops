#!/usr/bin/env python3
import pathlib

from setuptools import setup
from setuptools import find_packages

HERE = pathlib.Path(__file__).parent

VERSION = "0.4.0"
PACKAGE_NAME = "PyBrOpt"
AUTHOR = "Robert Z. Shrote"
AUTHOR_EMAIL = "shrotero@msu.edu"
URL = "https://github.com/oceaquaris/PyBrOpt"

LICENSE = "Apache License 2.0"
DESCRIPTION = "Python package for numerical breeding optimizations"
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESCRIPTION_TYPE = "text/markdown"

INSTALL_REQUIRES = [
    "numpy",
    "pandas",
    "cyvcf2",
    "igraph",
    "pytest",
    "pytest-datadir"
]

setup(
    name = PACKAGE_NAME,
    description = DESCRIPTION,
    long_description = LONG_DESCRIPTION,
    long_description_content_type = LONG_DESCRIPTION_TYPE,
    author = AUTHOR,
    license = LICENSE,
    author_email = AUTHOR_EMAIL,
    url = URL,
    install_requires = INSTALL_REQUIRES,
    packages = find_packages()
)
