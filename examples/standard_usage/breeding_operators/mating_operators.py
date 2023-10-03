#!/usr/bin/env python3

from typing import Optional

###
### Loading Class Modules
### =====================

# import the MatingOperator class (an abstract interface class)
from pybrops.breed.op.mate.MatingOperator import MatingOperator

###
### Defining Mating Operators
### =========================

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

###
### Creating Mating Operators
### =========================

# create a mating operator using constructor
mateop = MyMatingOperator()

###
### Mating Individiuals in a Breeding Program
### =========================================

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
