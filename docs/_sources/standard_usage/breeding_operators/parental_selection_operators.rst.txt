Parental Selection Operators
############################

Class Family Overview
=====================

The ``ParentSelectionOperator`` family of classes is used to select parents in a breeding program in the parental selection step of the universal breeding algorithm. Parent selection is extremely diverse. It may involve the selection of individuals based on phenotypes, genotypes, or combinations of the two. There may be one or more traits under selection.

Summary of Parent Selection Operators
=====================================

Parent selection operators can be found in the ``pybropos.breed.op.psel`` module.

.. list-table:: Summary of classes in the ``pybrops.breed.op.psel`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``ParentSelectionOperator``
      - Abstract
      - Interface for all parent selection operator classes.

Parent Selection Operator Properties
====================================

Parent selection operators do not have any required properties defined in the ``ParentSelectionOperator`` interface.

Loading Class Modules
=====================

Parent selection operators can be imported as demonstrated in the code below.

.. code-block:: python

    # import the ParentSelectionOperator class (an abstract interface class)
    from pybrops.breed.op.psel.ParentSelectionOperator import ParentSelectionOperator

Defining Parent Selection Operators
===================================

Since breeding programs are incredibly diverse, it is the job of the user to define a parent selection operator for his or her breeding program simulation. To do this, extend the ``ParentSelectionOperator`` base abstract class and override the ``pselect`` abstract method with the code defining the breeding program. Below is an example of how to override the ``ParentSelectionOperator`` class. For the purpose of this example, method code is left empty.

.. code-block:: python

    class MyParentSelectionOperator(ParentSelectionOperator):
        def __init__(self, *args: tuple, **kwargs: dict) -> None:
            """
            Constructor for custom ParentSelectionOperator

            Parameters
            ----------
            args : tuple
                Any user defined arguments.
            
            kwargs : dict
                Any user defined keyword arguments.
            """
            # user defined code
            pass
        def pselect(
                self, 
                genome: dict, 
                geno: dict, 
                pheno: dict, 
                bval: dict, 
                gmod: dict, 
                t_cur: int, 
                t_max: int, 
                miscout: Optional[dict] = None, 
                **kwargs: dict
            ) -> tuple:
            """
            Select individuals to serve as parents in a breeding program.

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
            miscout : dict, None
                Pointer to a dictionary for miscellaneous user defined output.
                If ``dict``, write to dict (may overwrite previously defined fields).
                If ``None``, user defined output is not calculated or stored.
            kwargs : dict
                Additional keyword arguments.

            Returns
            -------
            out : tuple
                A tuple of length 6: ``(mcfg, genome, geno, pheno, bval, gmod)``

                Where:

                - ``mcfg`` is a ``dict`` of mating configurations for the breeding program.
                - ``genome`` is a ``dict`` of genomes for the breeding program.
                - ``geno`` is a ``dict`` of genotypes for the breeding program.
                - ``pheno`` is a ``dict`` of phenotypes for the breeding program.
                - ``bval`` is a ``dict`` of breeding values for the breeding program.
                - ``gmod`` is a ``dict`` of genomic models for the breeding program.
            """
            # user defined code
            return {}, {}, {}, {}, {}, {}

Creating Parent Selection Operators
===================================

Since ``ParentSelectionOperator`` classes are entirely user defined, object construction is entirely implementation dependent. There are no restrictions on how an ``ParentSelectionOperator`` must be constructed and any number of arguments or keyword arguments may be used in the constructor. Below demonstrates the construction of the ``ParentSelectionOperator`` defined above.

.. code-block:: python

    # create a parent selection operator using constructor
    pselop = MyParentSelectionOperator()

Parental Selection for a Breeding Program
=========================================

To select parents from a breeding program, use the ``pselect`` method, which returns dictionaries of mating configurations, genomes, genotypes, phenotypes, breeding values, and genomic models for use in a breeding program simulation. The code below demonstrates the use of this method.

.. code-block:: python

    # select parents for a breeding program
    mcfg, genome, geno, pheno, bval, gmod = pselop.pselect(
        genome = {}, 
        geno = {}, 
        pheno = {}, 
        bval = {}, 
        gmod = {}, 
        t_cur = 0, 
        t_max = 0, 
    )
