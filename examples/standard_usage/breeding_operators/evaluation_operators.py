#!/usr/bin/env python3

from typing import Optional

###
### Loading Class Modules
### =====================

# import the EvaluationOperator class (an abstract interface class)
from pybrops.breed.op.eval.EvaluationOperator import EvaluationOperator

###
### Defining Evaluation Operators
### =============================

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

###
### Creating Evaluation Operators
### =============================

# create an evaluation operator using constructor
evalop = MyEvaluationOperator()

###
### Evaluation of Individuals in a Breeding Program
### ===============================================

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
