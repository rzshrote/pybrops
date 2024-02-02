#!/usr/bin/env python3
import pathlib
import setuptools

# setup.py metadata
setup_location = "." # pathlib.Path(__file__).parents[1]

# package metadata: general descriptors
pybrops_name = "pybrops"
pybrops_version = "1.0.0"
pybrops_author = "Robert Z. Shrote"
pybrops_author_email = "shrotero@msu.edu"
pybrops_description = "Python package for breeding program numerical optimization and simulation"
with open("README.md", "r", encoding = "utf-8") as readme_file:
    pybrops_description_long = readme_file.read()
    pybrops_description_long_type = "text/markdown"

# package metadata: project URLs
pybrops_url = "https://github.com/rzshrote/pybrops"
pybrops_project_url = {
    "Bug Tracker": "https://github.com/rzshrote/pybrops/issues",
}

# package metadata: licensing and classifiers
pybrops_license = "Apache License 2.0"
pybrops_classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
]

# package metadata: installation requirements
pybrops_requirements_python = ">=3.6"
pybrops_requirements_install = [
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
    "pymoo",
    "cyvcf2",
    "h5py",
    "DEAP",
    "pytest",
]

# package metadata: package locations
pybrops_package_directory = {"" : setup_location}
pybrops_packages = setuptools.find_packages(where = setup_location)

# setup the package
setuptools.setup(
    name = pybrops_name,
    version = pybrops_version,
    author = pybrops_author,
    author_email = pybrops_author_email,
    description = pybrops_description,
    long_description = pybrops_description_long,
    long_description_content_type = pybrops_description_long_type,
    url = pybrops_url,
    project_urls = pybrops_project_url,
    license = pybrops_license,
    classifiers = pybrops_classifiers,
    package_dir = pybrops_package_directory,
    packages = pybrops_packages,
    python_requires = pybrops_requirements_python,
    install_requires = pybrops_requirements_install
)
