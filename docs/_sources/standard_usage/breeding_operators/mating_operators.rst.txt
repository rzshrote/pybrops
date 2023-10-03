Mating Operators
################

Class Family Overview
=====================

The ``MatingOperator`` family of classes is used to create offspring resulting from parental mating. This operator performs the mating step of the universal breeding algorithm. Since mating and crossing strategies within breeding programs are extremely diverse and species dependent, the user must define this operator in his or her breeding program.

Summary of Mating Operators
=====================================

Mating operators can be found in the ``pybropos.breed.op.mate`` module.

.. list-table:: Summary of classes in the ``pybrops.breed.op.mate`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``MatingOperator``
      - Abstract
      - Interface for all mating operator classes.

Mating Operator Properties
==========================

Mating operators do not have any required properties defined in the ``MatingOperator`` interface.

Loading Class Modules
=====================

Mating operators can be imported as demonstrated in the code below.

.. code-block:: python

    # import the MatingOperator class (an abstract interface class)
    from pybrops.breed.op.mate.MatingOperator import MatingOperator

Defining Mating Operators
=========================

Since breeding programs are incredibly diverse, it is the job of the user to define a mating operator for his or her breeding program simulation. To do this, extend the ``MatingOperator`` base abstract class and override the ``mate`` abstract method with the code defining the mating operation. Below is an example of how to override the ``MatingOperator`` class. For the purpose of this example, method code is left empty.

.. code-block:: python

    class MyMatingOperator(MatingOperator):
        def __init__(self, *args: tuple, **kwargs: dict) -> None:
            """
            Constructor for custom MatingOperator

            Parameters
            ----------
            args : tuple
                Any user defined arguments.
            
            kwargs : dict
                Any user defined keyword arguments.
            """
            # user defined code
            pass
        def mate(
                self, 
                mcfg: dict, 
                genome: dict, 
                geno: dict, 
                pheno: dict, 
                bval: dict, 
                gmod: dict, 
                t_cur: int, 
                t_max: int, 
                miscout: Optional[dict], 
                **kwargs: dict
            ) -> tuple:
            """
            Mate individuals selected as parents in a breeding program.

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
            miscout : dict, None
                Pointer to a dictionary for miscellaneous user defined output.
                If ``dict``, write to dict (may overwrite previously defined fields).
                If ``None``, user defined output is not calculated or stored.
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

Creating Mating Operators
=========================

Since ``MatingOperator`` classes are entirely user defined, object construction is entirely implementation dependent. There are no restrictions on how an ``MatingOperator`` must be constructed and any number of arguments or keyword arguments may be used in the constructor. Below demonstrates the construction of the ``MatingOperator`` defined above.

.. code-block:: python

    # create a mating operator using constructor
    mateop = MyMatingOperator()

Mating Individiuals in a Breeding Program
=========================================

To select parents from a breeding program, use the ``mate`` method, which returns dictionaries of genomes, genotypes, phenotypes, breeding values, and genomic models for use in a breeding program simulation. The code below demonstrates the use of this method.

.. code-block:: python

    # mate individuals in a breeding program
    genome, geno, pheno, bval, gmod = mateop.mate(
        mcfg = {},
        genome = {}, 
        geno = {}, 
        pheno = {}, 
        bval = {}, 
        gmod = {}, 
        t_cur = 0, 
        t_max = 0, 
    )
