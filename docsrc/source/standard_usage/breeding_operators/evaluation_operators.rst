Evaluation Operators
####################

Class Family Overview
=====================

The ``EvaluationOperator`` family of classes is used to evaluate individuals in a breeding program. This family of classes performs the evaluation step of the universal breeding algorithm. Evaluation methodology is extremely diverse across breeding programs and species. It may involve the evaluation individual phenotypes, genotypes, or combinations of the two. It is the job of the user to define these evalution activities within an ``EvalutationOperator``.

Summary of Evaluation Operators
===============================

Evaluation operators can be found in the ``pybropos.breed.op.eval`` module.

.. list-table:: Summary of classes in the ``pybrops.breed.op.eval`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``EvaluationOperator``
      - Abstract
      - Interface for all evaluation operator classes.

Evaluation Operator Properties
==============================

Evaluation operators do not have any required properties defined in the ``EvaluationOperator`` interface.

Loading Class Modules
=====================

Evaluation operators can be imported as demonstrated in the code below.

.. code-block:: python

    # import the EvaluationOperator class (an abstract interface class)
    from pybrops.breed.op.eval.EvaluationOperator import EvaluationOperator

Defining Evaluation Operators
=============================

Since breeding programs are incredibly diverse, it is the job of the user to define a evaluation operator for his or her breeding program simulation. To do this, extend the ``EvaluationOperator`` base abstract class and override the ``evaluate`` abstract method with the code defining the evaluation operation. Below is an example of how to override the ``EvaluationOperator`` class. For the purpose of this example, method code is left empty.

.. code-block:: python

    class MyEvaluationOperator(EvaluationOperator):
        def __init__(self, *args: tuple, **kwargs: dict) -> None:
            """
            Constructor for custom EvaluationOperator

            Parameters
            ----------
            args : tuple
                Any user defined arguments.
            
            kwargs : dict
                Any user defined keyword arguments.
            """
            # user defined code
            pass
        def evaluate(
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
            Evaluate individuals in a breeding program.

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

Creating Evaluation Operators
=============================

Since ``EvaluationOperator`` classes are entirely user defined, object construction is entirely implementation dependent. There are no restrictions on how an ``EvaluationOperator`` must be constructed and any number of arguments or keyword arguments may be used in the constructor. Below demonstrates the construction of the ``EvaluationOperator`` defined above.

.. code-block:: python

    # create an evaluation operator using constructor
    evalop = MyEvaluationOperator()

Evaluation of Individuals in a Breeding Program
===============================================

To select parents from a breeding program, use the ``evaluate`` method, which returns dictionaries of genomes, genotypes, phenotypes, breeding values, and genomic models for use in a breeding program simulation. The code below demonstrates the use of this method.

.. code-block:: python

    # evaluate a breeding program
    genome, geno, pheno, bval, gmod = evalop.evaluate(
        genome = {}, 
        geno = {}, 
        pheno = {}, 
        bval = {}, 
        gmod = {}, 
        t_cur = 0, 
        t_max = 0, 
    )
