# PyBrOpS
Python Breeding Optimizer and Simulator: A Python library for simulating
breeding programs and optimizing breeding-related problems.

## Dependencies
1) `cyvcf2-0.30.14+`
2) `cvxpy-1.1.18+`
3) `deap-1.3.1+`
4) `h5py-3.6.0+`
5) `matplotlib-3.5.1+`
6) `numpy-1.22.2+`
7) `pandas-1.4.1+`
8) `pytest-7.0.1+`
9) `pytest-datadir-1.3.1+`
10) `rpy2-3.4.5+`
11) `scipy-1.8.0+`
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
