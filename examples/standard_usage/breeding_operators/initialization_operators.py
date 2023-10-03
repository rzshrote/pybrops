#!/usr/bin/env python3

from typing import Optional

###
### Loading Class Modules
### =====================

# import the InitializationOperator class (an abstract interface class)
from pybrops.breed.op.init.InitializationOperator import InitializationOperator

###
### Defining Initialization Operators
### =================================

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

###
### Creating Initialization Operators
### =================================

# create an initalization operator using constructor
initop = MyInitializationOperator()

###
### Initialization of a Breeding Program
### ====================================

# initialize a breeding program's set of containers
genome, geno, pheno, bval, gmod = initop.initialize()
