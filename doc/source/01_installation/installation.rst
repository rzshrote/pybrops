Installation
############

General installation requirements
*********************************

``Python 3.8.0+`` is an installation requirement.

Python package dependencies
===========================

In addition to a base installation of Python, PyBrOpS relies on several dependencies. They are:

#. ``cyvcf2-0.30.14+`` (for reading VCFs)
#. ``cvxpy-1.1.18+`` (for linear programming)
#. ``deap-1.3.1+`` (for genetic algorithms)
#. ``h5py-3.6.0+`` (for HDF5 file support)
#. ``matplotlib-3.5.1+`` (for graphics)
#. ``numpy-1.22.2+`` (for matrix storage and algebra)
#. ``pandas-1.4.1+`` (for data frames)
#. ``pymoo-0.6.0+`` (for genetic algorithms)
#. ``scipy-1.8.0+`` (for miscellaneous numerical routines)

The following additional packages are also required for the developmental installation of PyBrOpS:

#. ``pytest-7.0.1+`` (for unit testing)
#. ``pytest-datadir-1.3.1+`` (for unit testing)
#. ``setuptools`` (for package building)
#. ``wheel`` (for package building)

Additional Linux installation requirements
==========================================

On Linux, you will need to install Python developmental headers in addition to
``Python 3.8.0+``. Even though PyBrOpS is written in pure Python, some of its
dependencies require compilation.

Below are the commands for installing Python developmental headers:

.. list-table:: Title
   :widths: 25 50
   :header-rows: 1

   * - Linux Distro
     - Command
   * - Fedora
     - ``sudo dnf install python3-devel``
   * - Ubuntu
     - ``sudo apt install python3-dev``

Release version installation
****************************

To install PyBrOpS, use ``pip3``:

.. code-block:: console

    $ pip3 install pybrops


Developmental version installation
**********************************

To install the developmental version, first ``clone`` the repository:

.. code-block:: console

    $ git clone https://github.com/rzshrote/pybrops.git

It is a best practice to create a virtual environment where PyBrOpS dependencies
can be installed. To do this, you can use the ``Makefile``:

.. code-block:: console

    $ make install-devel

Alternatively, this may be done manually using the following commands:

.. code-block:: console

    $ python3 -m venv env
    $ . env/bin/activate
    $ python3 -m pip install --editable .

Next, you must activate the virtual environment using either the ``.`` command
(for ``sh``) or the ``source`` command (for ``bash``):

.. code-block:: console

    $ . env/bin/activate

or

.. code-block:: console

    $ source env/bin/activate

Now that the virtual environment is activated, you can access ``pybrops``
through the ``python3`` command-line interface prompt and through script
execution.
