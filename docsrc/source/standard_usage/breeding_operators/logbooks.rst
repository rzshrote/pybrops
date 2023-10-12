Logbooks
########

Class Family Overview
=====================

The purpose of the ``Logbook`` family of classes is to calculate and store simulation metrics. Simulation records may be examined and/or saved after the simulation is complete.

Summary of Logbook Classes
==========================

Logbook classes can be found in the ``pybropos.breed.op.log`` module.

.. list-table:: Summary of classes in the ``pybrops.breed.op.log`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``Logbook``
      - Abstract
      - Interface for all logbook classes.

Logbook Properties
==================

``Logbook`` objects have two properties: a ``data`` property storing simulation records and a ``rep`` property storing the simulation replicate number.

.. list-table:: Summary of ``Logbook`` properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``data``
      - The record data stored in the logbook.
    * - ``rep``
      - The replicate number for the given simulation. Useful for simulations with multiple replicates.

Loading Class Modules
=====================

Logbooks can be imported as demonstrated in the code below.

.. code-block:: python

    # import the Logbook class (an abstract interface class)
    from pybrops.breed.op.log.Logbook import Logbook

Defining Logbooks
=================

Logbooks must be defined by the user since the number of metrics one could collect is extensive. To create a custom logbook class, extend the ``Logbook`` base abstract class and override its abstract properties and methods. Below is an example of how to define a custom logbook class. For the purpose of this example, method code is left empty.

.. code-block:: python

    class MyLogbook(Logbook):
        ################ Special Object Methods ################
        def __init__(self, *args: tuple, **kwargs: dict) -> None:
            """
            Constructor for custom Logbook

            Parameters
            ----------
            args : tuple
                Any user defined arguments.
            
            kwargs : dict
                Any user defined keyword arguments.
            """
            # user defined code
            self.data = {}
            self.rep = 1
        ################## Object Properties ###################
        @property
        def data(self) -> dict:
            """Logbook data."""
            return self._data
        @data.setter
        def data(self, value: dict) -> None:
            """Set logbook data."""
            self._data = value
        @property
        def rep(self) -> int:
            """Replicate number."""
            return self._rep
        @rep.setter
        def rep(self, value: int) -> None:
            """Set replicate number."""
            self._rep = value
        #################### Object Methods ####################
        def log_initialize(
                self, 
                genome: dict, 
                geno: dict, 
                pheno: dict, 
                bval: dict, 
                gmod: dict, 
                t_cur: int, 
                t_max: int, 
                **kwargs: dict
            ) -> None:
            """
            Record information directly after 'InitializationOperator.initialize'
            is called.

            Parameters
            ----------
            genome : dict
                Dictionary of genomes for the breeding program.
            geno : dict
                Dictionary of genotypes for the breeding program.
            pheno : dict
                Dictionary of phenotypes for the breeding program.
            bval : dict
                Dictionary of breeding values for the breeding program.
            gmod : dict
                Dictionary of genomic models for the breeding program.
            t_cur : int
                Current time in the breeding program.
            t_max : int
                Deadline time for the breeding program.
            kwargs : dict
                Additional keyword arguments.
            """
            # user defined code
            pass
        def log_pselect(
                self, 
                mcfg: dict, 
                genome: dict, 
                geno: dict, 
                pheno: dict, 
                bval: dict, 
                gmod: dict, 
                t_cur: int, 
                t_max: int, 
                **kwargs: dict
            ) -> None:
            """
            Record information directly after 'ParentSelectionOperator.pselect'
            is called.

            Parameters
            ----------
            mcfg : dict
                Dictionary of mating configurations for the breeding program.
            genome : dict
                Dictionary of genomes for the breeding program.
            geno : dict
                Dictionary of genotypes for the breeding program.
            pheno : dict
                Dictionary of phenotypes for the breeding program.
            bval : dict
                Dictionary of breeding values for the breeding program.
            gmod : dict
                Dictionary of genomic models for the breeding program.
            t_cur : int
                Current time in the breeding program.
            t_max : int
                Deadline time for the breeding program.
            kwargs : dict
                Additional keyword arguments.
            """
            # user defined code
            pass
        def log_mate(
                self, 
                genome: dict, 
                geno: dict, 
                pheno: dict, 
                bval: dict, 
                gmod: dict, 
                t_cur: int, 
                t_max: int, 
                **kwargs: dict
            ) -> None:
            # user defined code
            """
            Record information directly after 'MatingOperator.mate' is called.

            Parameters
            ----------
            genome : dict
                Dictionary of genomes for the breeding program.
            geno : dict
                Dictionary of genotypes for the breeding program.
            pheno : dict
                Dictionary of phenotypes for the breeding program.
            bval : dict
                Dictionary of breeding values for the breeding program.
            gmod : dict
                Dictionary of genomic models for the breeding program.
            t_cur : int
                Current time in the breeding program.
            t_max : int
                Deadline time for the breeding program.
            kwargs : dict
                Additional keyword arguments.
            """
            pass
        def log_evaluate(
                self, 
                genome: dict, 
                geno: dict, 
                pheno: dict, 
                bval: dict, 
                gmod: dict, 
                t_cur: int, 
                t_max: int, 
                **kwargs: dict
            ) -> None:
            """
            Record information directly after 'EvaluationOperator.evaluate' is
            called.

            Parameters
            ----------
            genome : dict
                Dictionary of genomes for the breeding program.
            geno : dict
                Dictionary of genotypes for the breeding program.
            pheno : dict
                Dictionary of phenotypes for the breeding program.
            bval : dict
                Dictionary of breeding values for the breeding program.
            gmod : dict
                Dictionary of genomic models for the breeding program.
            t_cur : int
                Current time in the breeding program.
            t_max : int
                Deadline time for the breeding program.
            kwargs : dict
                Additional keyword arguments.
            """
            # user defined code
            pass
        def log_sselect(
                self, 
                genome: dict, 
                geno: dict, 
                pheno: dict, 
                bval: dict, 
                gmod: dict, 
                t_cur: int, 
                t_max: int, 
                **kwargs: dict
            ) -> None:
            """
            Record information directly after 'SurvivorSelectionOperator.sselect'
            is called.

            Parameters
            ----------
            genome : dict
                Dictionary of genomes for the breeding program.
            geno : dict
                Dictionary of genotypes for the breeding program.
            pheno : dict
                Dictionary of phenotypes for the breeding program.
            bval : dict
                Dictionary of breeding values for the breeding program.
            gmod : dict
                Dictionary of genomic models for the breeding program.
            t_cur : int
                Current time in the breeding program.
            t_max : int
                Deadline time for the breeding program.
            kwargs : dict
                Additional keyword arguments.
            """
            # user defined code
            pass
        def reset(self) -> None:
            """
            Reset Logbook internals.
            """
            self.data = {}
            self.rep = 1
        def write(self, filename: str) -> None:
            """
            Write Logbook to file

            Parameters
            ----------
            filename : str
                File name to which to write file.
            """
            # user defined code
            pass

Creating Logbooks
=================

To create a ``Logbook`` object, use the corresponding class's constructor. There are no restrictions on how a ``Logbook`` constructor must operate, so object construction is entirely implementation dependent. Below demonstrates the construction of the ``Logbook`` class defined above.

.. code-block:: python

    # create a new logbook
    lbook = MyLogbook()

Logging States in a Breeding Program
====================================

Logbooks have several methods to record data. These recording methods correspond to steps in the universal breeding algorithm.

Logging after breeding program initialization
---------------------------------------------

The ``log_initialize`` method can be used to record simulation statistics directly following the initialization step in the universal breeding algorithm. The code below demonstrates its use.

.. code-block:: python

    # gather data after breeding program initialization
    lbook.log_initialize(
        genome = {},
        geno = {},
        pheno = {},
        bval = {},
        gmod = {},
        t_cur = 0,
        t_max = 0,
    )

Logging after breeding program parent selection
-----------------------------------------------

The ``log_pselect`` method can be used to record simulation statistics directly following the parent selection step in the universal breeding algorithm. The code below demonstrates its use.

.. code-block:: python

    # gather data after breeding program parent selection
    lbook.log_pselect(
        mcfg = {},
        genome = {},
        geno = {},
        pheno = {},
        bval = {},
        gmod = {},
        t_cur = 0,
        t_max = 0,
    )

Logging after breeding program mating
-------------------------------------

The ``log_mate`` method can be used to record simulation statistics directly following the mating step in the universal breeding algorithm. The code below demonstrates its use.

.. code-block:: python

    # gather data after breeding program mating
    lbook.log_mate(
        genome = {},
        geno = {},
        pheno = {},
        bval = {},
        gmod = {},
        t_cur = 0,
        t_max = 0,
    )

Logging after breeding program evaluation
-----------------------------------------

The ``log_evaluate`` method can be used to record simulation statistics directly following the evaluation step in the universal breeding algorithm. The code below demonstrates its use.

.. code-block:: python

    # gather data after breeding program evaluation
    lbook.log_evaluate(
        genome = {},
        geno = {},
        pheno = {},
        bval = {},
        gmod = {},
        t_cur = 0,
        t_max = 0,
    )

Logging after breeding program survivor selection
-------------------------------------------------

The ``log_sselect`` method can be used to record simulation statistics directly following the survivor selection step in the universal breeding algorithm. The code below demonstrates its use.

.. code-block:: python

    # gather data after breeding program survivor selection
    lbook.log_sselect(
        genome = {},
        geno = {},
        pheno = {},
        bval = {},
        gmod = {},
        t_cur = 0,
        t_max = 0,
    )

Resetting a Logbook
===================

Logbooks can be reset and their data erased using the ``reset`` method, demonstrated below.

.. code-block:: python

    # reset logbook internals
    lbook.reset()

Writing a Logbook to a File
===========================

Logbook data can be exported to a file using the ``write`` method, demonstrated below.

.. code-block:: python

    # write logbook to file
    lbook.write("filename.csv")
