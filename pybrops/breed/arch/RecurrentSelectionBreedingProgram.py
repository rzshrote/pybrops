"""
Module implementing recurrent selection.
"""

__all__ = [
    "RecurrentSelectionBreedingProgram",
    "check_is_RecurrentSelectionBreedingProgram",
]

import copy
from typing import Union

from pybrops.breed.arch.BreedingProgram import BreedingProgram
from pybrops.breed.op.init.InitializationOperator import InitializationOperator, check_is_InitializationOperator
from pybrops.breed.op.eval.EvaluationOperator import EvaluationOperator, check_is_EvaluationOperator
from pybrops.breed.op.mate.MatingOperator import MatingOperator, check_is_MatingOperator
from pybrops.breed.op.psel.ParentSelectionOperator import ParentSelectionOperator, check_is_ParentSelectionOperator
from pybrops.breed.op.ssel.SurvivorSelectionOperator import SurvivorSelectionOperator, check_is_SurvivorSelectionOperator
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.core.error.error_type_python import check_is_int
from pybrops.core.error.error_value_python import check_dict_has_keys

class RecurrentSelectionBreedingProgram(BreedingProgram):
    """
    Class implementing recurrent selection. This class is very generic and
    highly modular to facilitate rapid prototyping.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            initop: InitializationOperator, 
            pselop: ParentSelectionOperator, 
            mateop: MatingOperator, 
            evalop: EvaluationOperator, 
            sselop: SurvivorSelectionOperator, 
            t_max, 
            start_genome = None, 
            start_geno = None, 
            start_pheno = None, 
            start_bval = None, 
            start_gmod = None, 
            **kwargs: dict
        ) -> None:
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

    ############################ Object Properties #############################

    ############ Starting condition containers #############
    @property
    def start_genome(self) -> Union[dict,None]:
        """Starting genomes for individuals in the breeding program."""
        return self._start_genome
    @start_genome.setter
    def start_genome(self, value: Union[dict,None]) -> None:
        """Set starting genomes for individuals in the breeding program"""
        if value is not None:
            check_is_dict(value, "start_genome")
        self._start_genome = value

    @property
    def start_geno(self) -> Union[dict,None]:
        """Starting genotypes for individuals in the breeding program."""
        return self._start_geno
    @start_geno.setter
    def start_geno(self, value: Union[dict,None]) -> None:
        """Set starting genotypes for individuals in the breeding program"""
        if value is not None:
            check_is_dict(value, "start_geno")
        self._start_geno = value

    @property
    def start_pheno(self) -> Union[dict,None]:
        """Starting phenotypes for individuals in the breeding program."""
        return self._start_pheno
    @start_pheno.setter
    def start_pheno(self, value: Union[dict,None]) -> None:
        """Set starting phenotypes for individuals in the breeding program"""
        if value is not None:
            check_is_dict(value, "start_pheno")
        self._start_pheno = value

    @property
    def start_bval(self) -> Union[dict,None]:
        """Starting breeding values for individuals in the breeding program."""
        return self._start_bval
    @start_bval.setter
    def start_bval(self, value: Union[dict,None]) -> None:
        """Set starting breeding values for individuals in the breeding program"""
        if value is not None:
            check_is_dict(value, "start_bval")
        self._start_bval = value

    @property
    def start_gmod(self) -> Union[dict,None]:
        """Starting genomic models for individuals in the breeding program."""
        return self._start_gmod
    @start_gmod.setter
    def start_gmod(self, value: Union[dict,None]) -> None:
        """Set starting genomic models for individuals in the breeding program"""
        if value is not None:
            check_is_dict(value, "start_gmod")
        self._start_gmod = value


    ############ Program information containers ############
    @property
    def genome(self) -> object:
        """Genomes for individuals in the breeding program."""
        return self._genome
    @genome.setter
    def genome(self, value: object) -> None:
        """Set genomes for individuals in the breeding program"""
        check_is_dict(value, "genome")
        self._genome = value

    @property
    def geno(self) -> object:
        """Genotypes for individuals in the breeding program."""
        return self._geno
    @geno.setter
    def geno(self, value: object) -> None:
        """Set genotypes for individuals in the breeding program"""
        check_is_dict(value, "geno")
        self._geno = value

    @property
    def pheno(self) -> object:
        """Phenotypes for individuals in the breeding program."""
        return self._pheno
    @pheno.setter
    def pheno(self, value: object) -> None:
        """Set phenotypes for individuals in the breeding program"""
        check_is_dict(value, "pheno")
        self._pheno = value

    @property
    def bval(self) -> object:
        """Breeding values for individuals in the breeding program."""
        return self._bval
    @bval.setter
    def bval(self, value: object) -> None:
        """Set breeding values for individuals in the breeding program"""
        check_is_dict(value, "bval")
        self._bval = value

    @property
    def gmod(self) -> object:
        """Genomic models for individuals in the breeding program."""
        return self._gmod
    @gmod.setter
    def gmod(self, value: object) -> None:
        """Set genomic models for individuals in the breeding program"""
        check_is_dict(value, "gmod")
        self._gmod = value

    ######### Breeding program operator properties #########
    @property
    def initop(self) -> InitializationOperator:
        """Initialization operator."""
        return self._initop
    @initop.setter
    def initop(self, value: InitializationOperator) -> None:
        """Set the initialization operator"""
        check_is_InitializationOperator(value, "initop")
        self._initop = value

    @property
    def pselop(self) -> ParentSelectionOperator:
        """Parent selection operator."""
        return self._pselop
    @pselop.setter
    def pselop(self, value: ParentSelectionOperator) -> None:
        """Set the parent selection operator"""
        check_is_ParentSelectionOperator(value, "pselop")
        self._pselop = value

    @property
    def mateop(self) -> MatingOperator:
        """Mating operator."""
        return self._mateop
    @mateop.setter
    def mateop(self, value: MatingOperator) -> None:
        """Set the mating operator"""
        check_is_MatingOperator(value, "mateop")
        self._mateop = value

    @property
    def evalop(self) -> EvaluationOperator:
        """Evaluation operator."""
        return self._evalop
    @evalop.setter
    def evalop(self, value: EvaluationOperator) -> None:
        """Set the evaluation operator"""
        check_is_EvaluationOperator(value, "evalop")
        self._evalop = value

    @property
    def sselop(self) -> SurvivorSelectionOperator:
        """Survivor selection operator."""
        return self._sselop
    @sselop.setter
    def sselop(self, value: SurvivorSelectionOperator) -> None:
        """Set the survivor selection operator"""
        check_is_SurvivorSelectionOperator(value, "sselop")
        self._sselop = value

    ############# Generation number properties #############
    @property
    def t_cur(self) -> int:
        """Current time of the BreedingNode."""
        return self._t_cur
    @t_cur.setter
    def t_cur(self, value: int) -> None:
        """Set the current time for the BreedingNode"""
        check_is_int(value, "t_cur")
        self._t_cur = value

    @property
    def t_max(self) -> int:
        """Maximum time of the BreedingNode."""
        return self._t_max
    @t_max.setter
    def t_max(self, value: int) -> None:
        """Set the maximum time for the BreedingNode"""
        check_is_int(value, "t_max")
        self._t_max = value

    ############################## Object Methods ##############################

    ############# Initialize breeding program ##############
    def initialize(self, **kwargs: dict):
        """
        Initialize the breeding program with genotypes, phenotypes, and genomic
        models.
        """
        self.start_genome, self.start_geno, self.start_pheno, self.start_bval, self.start_gmod = self._initop.initialize(**kwargs)

    def is_initialized(self, **kwargs: dict):
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
    def reset(self, **kwargs: dict):
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

    def advance(self, ngen, lbook, verbose = False, **kwargs: dict):
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

    def evolve(self, nrep, ngen, lbook, loginit = True, verbose = False, **kwargs: dict):
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



################################## Utilities ###################################
def check_is_RecurrentSelectionBreedingProgram(v: object, vname: str) -> None:
    """
    Check if object is of type RecurrentSelectionBreedingProgram. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RecurrentSelectionBreedingProgram):
        raise TypeError("'%s' must be a RecurrentSelectionBreedingProgram." % vname)
