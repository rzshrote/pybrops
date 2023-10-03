Survivor Selection Operators
############################

Class Family Overview
=====================

The ``SurvivorSelectionOperator`` family of classes is used to select desired progenies and discard undesired progenies in a breeding program in the survivor selection step of the universal breeding algorithm. Survivor selection is extremely diverse. It may involve the selection of progenies based on phenotypes, genotypes, or combinations of the two. Additionally, there may be one or more traits under selection.

Summary of Survivor Selection Operators
=======================================

Survivor selection operators can be found in the ``pybropos.breed.op.ssel`` module.

.. list-table:: Summary of classes in the ``pybrops.breed.op.ssel`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``SurvivorSelectionOperator``
      - Abstract
      - Interface for all survivor selection operator classes.

Survivor Selection Operator Properties
======================================

Survivor selection operators do not have any required properties defined in the ``SurvivorSelectionOperator`` interface.

Loading Class Modules
=====================

Survivor selection operators can be imported as demonstrated in the code below.

.. code-block:: python

    # import the SurvivorSelectionOperator class (an abstract interface class)
    from pybrops.breed.op.ssel.SurvivorSelectionOperator import SurvivorSelectionOperator

Defining Survivor Selection Operators
=====================================

Since breeding programs are incredibly diverse, it is the job of the user to define a survivor selection operator for his or her breeding program simulation. To do this, extend the ``SurvivorSelectionOperator`` base abstract class and override the ``sselect`` abstract method with the code defining the surivor selection operation. Below is an example of how to override the ``SurvivorSelectionOperator`` class. For the purpose of this example, method code is left empty.

.. code-block:: python

    class MySurvivorSelectionOperator(SurvivorSelectionOperator):
        def __init__(self, *args: tuple, **kwargs: dict) -> None:
            """
            Constructor for custom SurvivorSelectionOperator

            Parameters
            ----------
            args : tuple
                Any user defined arguments.
            
            kwargs : dict
                Any user defined keyword arguments.
            """
            # user defined code
            pass
        def sselect(
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
            Select progeny survivors in a breeding program.

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

Creating Survivor Selection Operators
=====================================

Since ``SurvivorSelectionOperator`` classes are entirely user defined, object construction is entirely implementation dependent. There are no restrictions on how an ``SurvivorSelectionOperator`` must be constructed and any number of arguments or keyword arguments may be used in the constructor. Below demonstrates the construction of the ``SurvivorSelectionOperator`` defined above.

.. code-block:: python

    # create a survivor selection operator using constructor
    sselop = MySurvivorSelectionOperator()

Survivor Selection for a Breeding Program
=========================================

To select survivors from a breeding program, use the ``sselect`` method, which returns dictionaries of genomes, genotypes, phenotypes, breeding values, and genomic models for use in a breeding program simulation. The code below demonstrates the use of this method.

.. code-block:: python

    # select survivors for a breeding program
    genome, geno, pheno, bval, gmod = sselop.sselect(
        genome = {}, 
        geno = {}, 
        pheno = {}, 
        bval = {}, 
        gmod = {}, 
        t_cur = 0, 
        t_max = 0, 
    )
