#!/usr/bin/env python3

from typing import Optional

###
### Loading Class Modules
### =====================

# import the ParentSelectionOperator class (an abstract interface class)
from pybrops.breed.op.psel.ParentSelectionOperator import ParentSelectionOperator

###
### Defining Parent Selection Operators
### ===================================

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

###
### Creating Parent Selection Operators
### =================================

# create a parent selection operator using constructor
pselop = MyParentSelectionOperator()

###
### Parental Selection for a Breeding Program
### =========================================

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