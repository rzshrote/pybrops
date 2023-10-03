#!/usr/bin/env python3

###
### Loading Class Modules
### =====================

# import the Logbook class (an abstract interface class)
from pybrops.breed.op.log.Logbook import Logbook

###
### Defining Logbooks
### =================

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
    
###
### Creating Logbooks
### =================

# create a new logbook
lbook = MyLogbook()

###
### Logging States in a Breeding Program
### ====================================

##
## Logging after breeding program initialization
## ---------------------------------------------

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

##
## Logging after breeding program parent selection
## -----------------------------------------------

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

##
## Logging after breeding program mating
## -------------------------------------

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

##
## Logging after breeding program evaluation
## -----------------------------------------

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

##
## Logging after breeding program survivor selection
## -------------------------------------------------

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

## 
## Resetting a logbook
## -------------------

# reset logbook internals
lbook.reset()

## 
## Writing a logbook to file
## -------------------------

# write logbook to file
lbook.write("filename.csv")