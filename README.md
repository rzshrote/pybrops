# PyBrOpS
Python Breeding Optimizer and Simulator: A Python library for simulating
breeding programs and optimizing breeding-related problems.

## General installation requirements
1) `Python 3.8.0+`

### Notes on Linux installation requirements
You will need to install Python developmental headers in addition to
`Python 3.8.0+`. Even though PyBrOpS is written in pure Python, some of its
dependencies require compilation.

Commands for installing Python developmental headers:

| Linux Distro  | Command                           |
| ------------- | --------------------------------- |
| Fedora        | `sudo dnf install python3-devel`  |
| Ubuntu        | `sudo apt install python3-dev`    |

## Python package dependencies
PyBrOpS relies on several dependencies. They are:

1) `cyvcf2-0.30.14+` (reading VCFs)
2) `cvxpy-1.1.18+` (linear programming)
3) `pymoo-0.6.0+` (genetic algorithms)
4) `h5py-3.6.0+` (HDF5 file support)
5) `matplotlib-3.5.1+` (graphics)
6) `numpy-1.22.2+` (matrix storage and algebra)
7) `pandas-1.4.1+` (data frames)
8) `pytest-7.0.1+` (unit tests)
9) `pytest-datadir-1.3.1+` (unit tests)
10) `rpy2-3.4.5+` (interfacing with R)
11) `scipy-1.8.0+` (miscellaneous numerical routines)
12) `setuptools`
13) `wheel`

# Developmental version installation
To install first `clone` the repository:
```
git clone https://github.com/rzshrote/pybrops.git
```

It is a best practice to create a virtual environment where PyBrOpS dependencies
can be installed. To do this, you can use the `Makefile`:
```
make install-devel
```

Alternatively, this may be done manually using the following commands:
```
python3 -m venv env
. env/bin/activate
python3 -m pip install --editable .
```

Next, you must activate the virtual environment using either the `.` command
(for `sh`) or the `source` command (for `bash`):
```
. env/bin/activate
```
or
```
source env/bin/activate
```

Now that the virtual environment is activated, you can access `pybrops`
through the `python3` command-line interface prompt and through script
execution.

# Documentation generation

Non-exhaustive list for building
```
source env/bin/activate

pip3 install sphinx-autodoc-typehints sphinx_rtd_theme pydata_sphinx_theme
```