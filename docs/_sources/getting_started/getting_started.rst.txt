Getting Started
###############

Basic Installation
==================

Python installation requirements
--------------------------------

``Python 3.8.0+`` is an installation requirement. To download Python, see the `Python website <https://www.python.org/>`_.

Python package dependencies
---------------------------

In addition to a base installation of Python, PyBrOpS relies on several dependencies. They are:

#. ``numpy-1.22.2+`` (for matrix storage and algebra)
#. ``pandas-1.4.1+`` (for data frames)
#. ``scipy-1.8.0+`` (for miscellaneous numerical routines)
#. ``matplotlib-3.5.1+`` (for graphics)
#. ``pymoo-0.6.0+`` (for genetic algorithms)
#. ``cyvcf2-0.30.14+`` (for reading VCFs)
#. ``h5py-3.6.0+`` (for HDF5 file support)

The following additional packages are also required for the developmental installation of PyBrOpS:

#. ``pytest-7.0.1+`` (for unit testing)
#. ``pytest-datadir-1.3.1+`` (for unit testing)
#. ``setuptools`` (for package building)
#. ``wheel`` (for package building)

Release version installation
----------------------------

To install PyBrOpS, use ``pip3``:

.. code-block:: console

    $ pip3 install pybrops

Basic Programming Knowledge Prerequisites
=========================================

PyBrOpS makes several programming knowledge assumptions of the user. Below is a table of these knowledge prerequisites and a link to a refresher course if you need to understand the basics.

.. list-table:: PyBrOpS Programming Knowledge Prerequisites
    :widths: 25 15 50
    :header-rows: 1

    * - Knowledge Prerequisite
      - Importance
      - Refresher Course
    * - The Python Programming Language
      - Essential
      - See `The Python Tutorial <https://docs.python.org/3/tutorial/>`_
    * - NumPy arrays
      - Essential
      - See `The NumPy User Guide <https://numpy.org/devdocs/user/index.html>`_
    * - Pandas DataFrames
      - Essential
      - See `The Pandas 10 Minute Guide <https://pandas.pydata.org/docs/user_guide/10min.html>`_
    * - PyMOO genetic algorithms
      - Optional
      - See `PyMOO.org <https://pymoo.org/>`_ to learn how to create custom algorithms in PyMOO's framework.

Using PyBrOpS
=============

If you desire to start making your own breeding simulations immediately and want to see complete, working simulation scripts, see :doc:`the examples <../standard_usage/examples/examples>`.

If you desire to gain a more in-depth knowledge about how PyBrOpS works before embarking on creating your own simulation scripts, see :doc:`the fundamentals <../fundamentals/fundamentals>` to understand PyBrOpS's design philosophy, followed by :doc:`the standard usage <../standard_usage/standard_usage>` section.

Advanced Installation
=====================

The following sections detail advanced installation instructions for developing on Linux systems or working with the developmental version of PyBrOpS.

Linux installation requirements
-------------------------------

On Linux, you will need to install Python developmental headers in addition to ``Python 3.8.0+``. Even though PyBrOpS is written in pure Python, some of its dependencies require compilation.

Below are the commands for installing Python developmental headers:

.. list-table:: Linux Developmental Library Installation Commands
   :widths: 25 50
   :header-rows: 1

   * - Linux Distro
     - Command
   * - Fedora
     - ``sudo dnf install python3-devel``
   * - Ubuntu
     - ``sudo apt install python3-dev``

Developmental version installation
----------------------------------

To install the developmental version, first ``clone`` the repository:

.. code-block:: console

    $ git clone https://github.com/rzshrote/pybrops.git

It is a best practice to create a virtual environment where PyBrOpS dependencies can be installed. To do this, you can use the ``Makefile``:

.. code-block:: console

    $ make install-devel

Alternatively, this may be done manually using the following commands:

.. code-block:: console

    $ python3 -m venv env
    $ . env/bin/activate
    $ python3 -m pip install --editable .

Next, you must activate the virtual environment using either the ``.`` command (for ``sh``) or the ``source`` command (for ``bash``):

.. code-block:: console

    $ . env/bin/activate

or

.. code-block:: console

    $ source env/bin/activate

Now that the virtual environment is activated, you can access ``pybrops`` through the ``python3`` command-line interface prompt and through script execution.
