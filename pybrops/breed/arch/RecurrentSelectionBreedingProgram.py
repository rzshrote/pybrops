"""
Module implementing recurrent selection.
"""

import copy
from typing import Any
import numpy

from pybrops.breed.arch.BreedingProgram import BreedingProgram
from pybrops.breed.op.init.InitializationOperator import check_is_InitializationOperator
from pybrops.breed.op.eval.EvaluationOperator import check_is_EvaluationOperator
from pybrops.breed.op.mate.MatingOperator import check_is_MatingOperator
from pybrops.breed.op.psel.ParentSelectionOperator import check_is_ParentSelectionOperator
from pybrops.breed.op.ssel.SurvivorSelectionOperator import check_is_SurvivorSelectionOperator
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_int
from pybrops.core.error import check_keys_in_dict

class RecurrentSelectionBreedingProgram(BreedingProgram):
    """
    Class implementing recurrent selection. This class is very generic and
    highly modular to facilitate rapid prototyping.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, initop, pselop, mateop, evalop, sselop, t_max, start_genome = None, start_geno = None, start_pheno = None, start_bval = None, start_gmod = None, **kwargs):
        """
        Constructor for the concrete class RecurrentSelectionBreedingProgram.

        This class simulates a recurrent selection breeding program.

        Parameters
        ----------
        initop : InitializationOperator
        pselop : ParentSelectionOperator
        mateop : MatingOperator
        evalop : EvaluationOperator
        sselop : SurvivorSelectionOperator
        t_max : int
        """
        super(RecurrentSelectionBreedingProgram, self).__init__(**kwargs)

        # save operators
        self.initop = initop
        self.pselop = pselop
        self.mateop = mateop
        self.evalop = evalop
        self.sselop = sselop

        # set time variables
        self.t_cur = 0
        self.t_max = t_max

        # TODO: go through set methods properly
        self.start_genome = start_genome
        self.start_geno = start_geno
        self.start_pheno = start_pheno
        self.start_bval = start_bval
        self.start_gmod = start_gmod

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############ Starting condition containers #############
    def start_genome():
        doc = "Starting genomes for individuals in the breeding program."
        def fget(self):
            """Get starting genomes for individuals in the breeding program"""
            return self._start_genome
        def fset(self, value):
            """Set starting genomes for individuals in the breeding program"""
            if value is not None:
                check_is_dict(value, "start_genome")
            self._start_genome = value
        def fdel(self):
            """Delete starting genomes for individuals in the breeding program"""
            del self._start_genome
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    start_genome = property(**start_genome())

    def start_geno():
        doc = "Starting genotypes for individuals in the breeding program."
        def fget(self):
            """Get starting genotypes for individuals in the breeding program"""
            return self._start_geno
        def fset(self, value):
            """Set starting genotypes for individuals in the breeding program"""
            if value is not None:
                check_is_dict(value, "start_geno")
            self._start_geno = value
        def fdel(self):
            """Delete starting genotypes for individuals in the breeding program"""
            del self._start_geno
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    start_geno = property(**start_geno())

    def start_pheno():
        doc = "Starting phenotypes for individuals in the breeding program."
        def fget(self):
            """Get starting phenotypes for individuals in the breeding program"""
            return self._start_pheno
        def fset(self, value):
            """Set starting phenotypes for individuals in the breeding program"""
            if value is not None:
                check_is_dict(value, "start_pheno")
            self._start_pheno = value
        def fdel(self):
            """Delete starting phenotypes for individuals in the breeding program"""
            del self._start_pheno
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    start_pheno = property(**start_pheno())

    def start_bval():
        doc = "Starting breeding values for individuals in the breeding program."
        def fget(self):
            """Get starting breeding values for individuals in the breeding program"""
            return self._start_bval
        def fset(self, value):
            """Set starting breeding values for individuals in the breeding program"""
            if value is not None:
                check_is_dict(value, "start_bval")
            self._start_bval = value
        def fdel(self):
            """Delete starting breeding values for individuals in the breeding program"""
            del self._start_bval
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    start_bval = property(**start_bval())

    def start_gmod():
        doc = "Starting genomic models for individuals in the breeding program."
        def fget(self):
            """Get starting genomic models for individuals in the breeding program"""
            return self._start_gmod
        def fset(self, value):
            """Set starting genomic models for individuals in the breeding program"""
            if value is not None:
                check_is_dict(value, "start_gmod")
            self._start_gmod = value
        def fdel(self):
            """Delete starting genomic models for individuals in the breeding program"""
            del self._start_gmod
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    start_gmod = property(**start_gmod())

    ############ Program information containers ############
    def genome():
        doc = "Genomes for individuals in the breeding program."
        def fget(self):
            """Get genomes for individuals in the breeding program"""
            return self._genome
        def fset(self, value):
            """Set genomes for individuals in the breeding program"""
            check_is_dict(value, "genome")
            self._genome = value
        def fdel(self):
            """Delete genomes for individuals in the breeding program"""
            del self._genome
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    genome = property(**genome())

    def geno():
        doc = "Genotypes for individuals in the breeding program."
        def fget(self):
            """Get genotypes for individuals in the breeding program"""
            return self._geno
        def fset(self, value):
            """Set genotypes for individuals in the breeding program"""
            check_is_dict(value, "geno")
            self._geno = value
        def fdel(self):
            """Delete genotypes for individuals in the breeding program"""
            del self._geno
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    geno = property(**geno())

    def pheno():
        doc = "Phenotypes for individuals in the breeding program."
        def fget(self):
            """Get phenotypes for individuals in the breeding program"""
            return self._pheno
        def fset(self, value):
            """Set phenotypes for individuals in the breeding program"""
            check_is_dict(value, "pheno")
            self._pheno = value
        def fdel(self):
            """Delete phenotypes for individuals in the breeding program"""
            del self._pheno
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    pheno = property(**pheno())

    def bval():
        doc = "Breeding values for individuals in the breeding program."
        def fget(self):
            """Get breeding values for individuals in the breeding program"""
            return self._bval
        def fset(self, value):
            """Set breeding values for individuals in the breeding program"""
            check_is_dict(value, "bval")
            self._bval = value
        def fdel(self):
            """Delete breeding values for individuals in the breeding program"""
            del self._bval
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    bval = property(**bval())

    def gmod():
        doc = "Genomic models for individuals in the breeding program."
        def fget(self):
            """Get genomic models for individuals in the breeding program"""
            return self._gmod
        def fset(self, value):
            """Set genomic models for individuals in the breeding program"""
            check_is_dict(value, "gmod")
            self._gmod = value
        def fdel(self):
            """Delete genomic models for individuals in the breeding program"""
            del self._gmod
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    gmod = property(**gmod())

    ######### Breeding program operator properties #########
    def initop():
        doc = "Initialization operator"
        def fget(self):
            """Get the initialization operator"""
            return self._initop
        def fset(self, value):
            """Set the initialization operator"""
            check_is_InitializationOperator(value, "initop")
            self._initop = value
        def fdel(self):
            """Delete the initialization operator"""
            del self._initop
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    initop = property(**initop())

    def pselop():
        doc = "Parent selection operator"
        def fget(self):
            """Get the parent selection operator"""
            return self._pselop
        def fset(self, value):
            """Set the parent selection operator"""
            check_is_ParentSelectionOperator(value, "pselop")
            self._pselop = value
        def fdel(self):
            """Delete the parent selection operator"""
            del self._pselop
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    pselop = property(**pselop())

    def mateop():
        doc = "Mating operator"
        def fget(self):
            """Get the mating operator"""
            return self._mateop
        def fset(self, value):
            """Set the mating operator"""
            check_is_MatingOperator(value, "mateop")
            self._mateop = value
        def fdel(self):
            """Delete the mating operator"""
            del self._mateop
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    mateop = property(**mateop())

    def evalop():
        doc = "Evaluation operator"
        def fget(self):
            """Get the evaluation operator"""
            return self._evalop
        def fset(self, value):
            """Set the evaluation operator"""
            check_is_EvaluationOperator(value, "evalop")
            self._evalop = value
        def fdel(self):
            """Delete the evaluation operator"""
            del self._evalop
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    evalop = property(**evalop())

    def sselop():
        doc = "Survivor selection operator"
        def fget(self):
            """Get the survivor selection operator"""
            return self._sselop
        def fset(self, value):
            """Set the survivor selection operator"""
            check_is_SurvivorSelectionOperator(value, "sselop")
            self._sselop = value
        def fdel(self):
            """Delete the survivor selection operator"""
            del self._sselop
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    sselop = property(**sselop())

    ############# Generation number properties #############
    def t_cur():
        doc = "Current time of the BreedingNode."
        def fget(self):
            """Get the current time for the breeding program"""
            return self._t_cur
        def fset(self, value):
            """Set the current time for the breeding program"""
            check_is_int(value, "t_cur")
            self._t_cur = value
        def fdel(self):
            """Delete the current time for the breeding program"""
            del self._t_cur
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    t_cur = property(**t_cur())

    def t_max():
        doc = "Maximum time of the BreedingNode."
        def fget(self):
            """Get the maximum time for the breeding program"""
            return self._t_max
        def fset(self, value):
            """Set the maximum time for the breeding program"""
            check_is_int(value, "t_max")
            self._t_max = value
        def fdel(self):
            """Delete the maximum time for the breeding program"""
            del self._t_max
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    t_max = property(**t_max())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############# Initialize breeding program ##############
    def initialize(self, **kwargs):
        """
        Initialize the breeding program with genotypes, phenotypes, and genomic
        models.
        """
        self.start_genome, self.start_geno, self.start_pheno, self.start_bval, self.start_gmod = self._initop.initialize(**kwargs)

    def is_initialized(self, **kwargs):
        """
        Return whether or not the BreedingProgram has been initialized with a
        starting set of conditions.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : boolean
            True if the BreedingProgram has been initialized.
            False if the BreedingProgram has not been initialized.
        """
        return (
            (self._start_genome is not None) and
            (self._start_geno is not None) and
            (self._start_pheno is not None) and
            (self._start_bval is not None) and
            (self._start_gmod is not None)
        )

    ################ Whole breeding program ################
    def reset(self, **kwargs):
        """
        Reset the evolution of the breeding program back to starting conditions.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        self.genome = copy.deepcopy(self.start_genome)  # reset genomes container
        self.geno = copy.deepcopy(self.start_geno)      # reset genotypes container
        self.pheno = copy.deepcopy(self.start_pheno)    # reset phenotypes container
        self.bval = copy.deepcopy(self.start_bval)      # reset breeding values container
        self.gmod = copy.deepcopy(self.start_gmod)      # reset genomic model container
        self.t_cur = 0                                  # reset time

    def advance(self, ngen, lbook, verbose = False, **kwargs):
        """
        Advance the breeding program by a specified number of generations.

        Parameters
        ----------
        ngen : int
            Number of generations to advance the BreedingProgram.
        lbook : Logbook
            Logbook into which to write statistics.
        kwargs : dict
            Additional keyword arguments.
        """
        # iterate through main breeding loop for ngen generations
        for _ in range(ngen):
            # display verobse messages if needed
            if verbose:
                print("Simulating rep {0}: gen {1} of {2}".format(lbook.rep, _+1, ngen,))

            ####################################################################
            ########################## select parents ##########################
            ####################################################################
            misc = {}
            mcfg, self.genome, self.geno, self.pheno, self.bval, self.gmod = self._pselop.pselect(
                genome = self._genome,
                geno = self._geno,
                pheno = self._pheno,
                bval = self._bval,
                gmod = self._gmod,
                t_cur = self._t_cur,
                t_max = self._t_max,
                miscout = misc
            )
            lbook.log_pselect(
                mcfg = mcfg,
                genome = self._genome,
                geno = self._geno,
                pheno = self._pheno,
                bval = self._bval,
                gmod = self._gmod,
                t_cur = self._t_cur,
                t_max = self._t_max,
                **misc
            )

            ####################################################################
            ########################### mate parents ###########################
            ####################################################################
            misc = {}
            self.genome, self.geno, self.pheno, self.bval, self.gmod = self._mateop.mate(
                mcfg = mcfg,
                genome = self._genome,
                geno = self._geno,
                pheno = self._pheno,
                bval = self._bval,
                gmod = self._gmod,
                t_cur = self._t_cur,
                t_max = self._t_max,
                miscout = misc
            )
            lbook.log_mate(
                mcfg = mcfg,
                genome = self._genome,
                geno = self._geno,
                pheno = self._pheno,
                bval = self._bval,
                gmod = self._gmod,
                t_cur = self._t_cur,
                t_max = self._t_max,
                **misc
            )

            ####################################################################
            ######################## evaluate genotypes ########################
            ####################################################################
            misc = {}
            self.genome, self.geno, self.pheno, self.bval, self.gmod = self._evalop.evaluate(
                genome = self._genome,
                geno = self._geno,
                pheno = self._pheno,
                bval = self._bval,
                gmod = self._gmod,
                t_cur = self._t_cur,
                t_max = self._t_max,
                miscout = misc
            )
            lbook.log_evaluate(
                genome = self._genome,
                geno = self._geno,
                pheno = self._pheno,
                bval = self._bval,
                gmod = self._gmod,
                t_cur = self._t_cur,
                t_max = self._t_max,
                **misc
            )

            ####################################################################
            ######################### select survivors #########################
            ####################################################################
            misc = {}
            self.genome, self.geno, self.pheno, self.bval, self.gmod = self._sselop.sselect(
                genome = self._genome,
                geno = self._geno,
                pheno = self._pheno,
                bval = self._bval,
                gmod = self._gmod,
                t_cur = self._t_cur,
                t_max = self._t_max,
                miscout = misc
            )
            lbook.log_sselect(
                genome = self._genome,
                geno = self._geno,
                pheno = self._pheno,
                bval = self._bval,
                gmod = self._gmod,
                t_cur = self._t_cur,
                t_max = self._t_max,
                **misc
            )

            ####################################################################
            ######################### variable updates #########################
            ####################################################################
            # increment time variables
            self._t_cur += 1

    def evolve(self, nrep, ngen, lbook, loginit = True, verbose = False, **kwargs):
        """
        Evolve the breeding program for a set number of replications and
        generations. The BreedingProgram is restarted using the starting geno,
        bval, gmod containers.

        Parameters
        ----------
        nrep : int
            Number of evolution replicates.
        ngen : int, None
            Number of generations to evolve the population for each replicate.
            If None, use 't_max'.
            Note: if specified this does not modify 't_max' which may affect
            operators that utilize 't_max'.
        lbook : Logbook
            Logbook into which to write statistics.
        loginit : bool
            Whether to log the initial state before main loop.
        verbose : bool
            Whether to print the rep number.
        """
        # initialize if needed
        if not self.is_initialized():
            self.initialize()

        # main replication loop
        for r in range(nrep):
            # increment rep in Logbook
            lbook.rep += 1

            # verbose messages
            if verbose:
                print("Simulating rep {0} of {1} (Logbook entry {2})".format(r+1, nrep, lbook.rep))

            # reset simulation to starting conditions
            self.reset()

            # evaluate breeding program at starting conditions
            misc = {}
            self.genome, self.geno, self.pheno, self.bval, self.gmod = self._evalop.evaluate(
                genome = self._genome,
                geno = self._geno,
                pheno = self._pheno,
                bval = self._bval,
                gmod = self._gmod,
                t_cur = self._t_cur,
                t_max = self._t_max,
                miscout = misc
            )
            if loginit:
                lbook.log_initialize(
                    genome = self._genome,
                    geno = self._geno,
                    pheno = self._pheno,
                    bval = self._bval,
                    gmod = self._gmod,
                    t_cur = self._t_cur,
                    t_max = self._t_max,
                    **misc
                )

            # increment t_cur from 0 to 1 (first generation)
            self.t_cur += 1

            # evolve the population
            self.advance(
                ngen = ngen,
                lbook = lbook,
                verbose = verbose,
                **kwargs
            )



################################################################################
################################## Utilities ###################################
################################################################################
def is_RecurrentSelectionBreedingProgram(v: Any) -> bool:
    """
    Determine whether an object is a RecurrentSelectionBreedingProgram.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a RecurrentSelectionBreedingProgram object instance.
    """
    return isinstance(v, RecurrentSelectionBreedingProgram)

def check_is_RecurrentSelectionBreedingProgram(v: Any, varname: str) -> None:
    """
    Check if object is of type RecurrentSelectionBreedingProgram. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RecurrentSelectionBreedingProgram):
        raise TypeError("'%s' must be a RecurrentSelectionBreedingProgram." % varname)
