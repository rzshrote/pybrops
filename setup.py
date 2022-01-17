#!/usr/bin/env python3
import pathlib
import setuptools

# setup.py metadata
setup_location = "." # pathlib.Path(__file__).parents[1]

# package metadata: general descriptors
pybropt_name = "PyBrOpt"
pybropt_version = "1.0.0"
pybropt_author = "Robert Z. Shrote"
pybropt_author_email = "shrotero@msu.edu"
pybropt_description = "Python package for numerical breeding optimizations"
with open("README.md", "r", encoding = "utf-8") as readme_file:
    pybropt_description_long = readme_file.read()
    pybropt_description_long_type = "text/markdown"

# package metadata: project URLs
pybropt_url = "https://github.com/rzshrote/PyBrOpt"
pybropt_project_url = {
    "Bug Tracker": "https://github.com/rzshrote/PyBrOpt/issues",
}

# package metadata: licensing and classifiers
pybropt_license = "Apache License 2.0"
pybropt_classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: Apache Software License",
    "Operating System :: OS Independent",
]

# package metadata: installation requirements
pybropt_requirements_python = ">=3.6"
pybropt_requirements_install = [
    "numpy",
    "pandas",
    "cyvcf2",
    "igraph",
    "rpy2",
    "h5py",
    "deap",
    "scipy",
    "cvxpy",
    "sklearn",
    "matplotlib",
    "pytest",
    "pytest-datadir"
]

# package metadata: package locations
pybropt_package_directory = {"" : setup_location}
pybropt_packages = setuptools.find_packages(where = setup_location)

# setup the package
setuptools.setup(
    name = pybropt_name,
    version = pybropt_version,
    author = pybropt_author,
    author_email = pybropt_author_email,
    description = pybropt_description,
    long_description = pybropt_description_long,
    long_description_content_type = pybropt_description_long_type,
    url = pybropt_url,
    project_urls = pybropt_project_url,
    license = pybropt_license,
    classifiers = pybropt_classifiers,
    package_dir = pybropt_package_directory,
    packages = pybropt_packages,
    python_requires = pybropt_requirements_python,
    install_requires = pybropt_requirements_install
)
