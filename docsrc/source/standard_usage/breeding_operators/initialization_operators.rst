Initialization Operators
########################

Class Family Overview
=====================

The ``InitializationOperator`` family of classes is used to initialize a breeding program in the universal breeding algorithm. For breeding programs simulated from scratch, typically the initialization operator is a series of burn-in generations starting from a random population. For breeding program simulations using real-world data, typically the initalization operator involves the reading and processing of breeding program data in preparation for the simulation.

Summary of Initialization Operators
===================================

Initialization operators can be found in the ``pybropos.breed.op.init`` module.

.. list-table:: Summary of classes in the ``pybrops.breed.op.init`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``InitializationOperator``
      - Abstract
      - Interface for all initialization operator classes.

Initialization Operators Properties
===================================

Initialization operators do not have any required properties defined in the ``InitializationOperator`` interface.

Loading Class Modules
=====================

Initialization operators can be imported as demonstrated in the code below.

.. code-block:: python

    # import the InitializationOperator class (an abstract interface class)
    from pybrops.breed.op.init.InitializationOperator import InitializationOperator

Defining Initialization Operators
=================================

Since breeding programs are incredibly diverse, it is the job of the user to define an initialization operator for his or her breeding program simulation. To do this, extend the ``InitializationOperator`` base abstract class and override the ``initialize`` abstract method with the code defining the initialization operation. Below is an example of how to override the ``InitializationOperator`` class. For the purpose of this example, method code is left empty.

.. code-block:: python

    class MyInitializationOperator(InitializationOperator):
        def __init__(self, *args: tuple, **kwargs: dict) -> None:
            """
            Constructor for custom InitializationOperator

            Parameters
            ----------
            args : tuple
                Any user defined arguments.
            
            kwargs : dict
                Any user defined keyword arguments.
            """
            # user defined code
            pass
        def initialize(
                self, 
                miscout: Optional[dict] = None, 
                **kwargs: dict
            ) -> tuple:
            """
            Initialize a breeding program.

            Parameters
            ----------
            miscout : dict, None, default = None
                Pointer to a dictionary for miscellaneous user defined output.
                If dict, write to dict (may overwrite previously defined fields).
                If None, user defined output is not calculated or stored.
            
            kwargs : dict
                Additional keyword arguments.

            Returns
            -------
            out : tuple
                A tuple of length 5: ``(genome, geno, pheno, bval, gmod)``

                Where:

                - ``genome`` is a ``dict`` of genomes for the breeding program.
                - ``geno`` is a ``dict`` of genotypes for the breeding program.
                - ``pheno`` is a ``dict`` of phenotypes for the breeding program.
                - ``bval`` is a ``dict`` of breeding values for the breeding program.
                - ``gmod`` is a ``dict`` of genomic models for the breeding program.
            """
            # user defined code
            return {}, {}, {}, {}, {}

Creating Initialization Operators
=================================

Since ``InitializationOperator`` classes are entirely user defined, object construction is entirely implementation dependent. There are no restrictions on how an ``InitializationOperator`` must be constructed and any number of arguments or keyword arguments may be used in the constructor. Below demonstrates the construction of the ``InitializationOperator`` defined above.

.. code-block:: python

    # create an initalization operator using constructor
    initop = MyInitializationOperator()

Initialization of a Breeding Program
====================================

To initialize a breeding program, use the ``initialize`` method, which returns dictionaries of genomes, genotypes, phenotypes, breeding values, and genomic models for use in a breeding program simulation. The code below demonstrates the use of this method.

.. code-block:: python

    # initialize a breeding program's set of containers
    genome, geno, pheno, bval, gmod = initop.initialize()
