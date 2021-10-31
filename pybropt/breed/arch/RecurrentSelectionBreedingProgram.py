import copy
import numpy

from . import BreedingProgram

from pybropt.core.error import check_is_int
from pybropt.core.error import check_is_dict
from pybropt.core.error import check_keys_in_dict
# from pybropt.breed.op.init import check_is_InitializationOperator
# from pybropt.breed.op.psel import check_is_ParentSelectionOperator
# from pybropt.breed.op.mate import check_is_MatingOperator
# from pybropt.breed.op.intg import check_is_GenotypeIntegrationOperator
# from pybropt.breed.op.eval import check_is_EvaluationOperator
# from pybropt.breed.op.intg import check_is_BreedingValueIntegrationOperator
# from pybropt.breed.op.calibr import check_is_GenomicModelCalibrationOperator
# from pybropt.breed.op.ssel import check_is_SurvivorSelectionOperator

class RecurrentSelectionBreedingProgram(BreedingProgram):
    """docstring for RecurrentSelectionBreedingProgram."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, initop, pselop, mateop, evalop, sselop, t_max, start_genome = None, start_geno = None, start_pheno = None, start_bval = None, start_gmod = None, **kwargs):
        """
        Constructor for the concrete class RecurrentSelectionBreedingProgram.

        This class simulates a recurrent selection breeding program.

        Parameters
        ----------
        initop :
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
        return locals()
    start_genome = property(**start_genome())

    def start_geno():
        doc = "Starting genotypes for individuals in the breeding program."
        def fget(self):
            """Get starting genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set starting genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete starting genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    start_geno = property(**start_geno())

    def start_pheno():
        doc = "Starting phenotypes for individuals in the breeding program."
        def fget(self):
            """Get starting phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set starting phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete starting phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    start_pheno = property(**start_pheno())

    def start_bval():
        doc = "Starting breeding values for individuals in the breeding program."
        def fget(self):
            """Get starting breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set starting breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete starting breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    start_bval = property(**start_bval())

    def start_gmod():
        doc = "Starting genomic models for individuals in the breeding program."
        def fget(self):
            """Get starting genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set starting genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete starting genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    start_gmod = property(**start_gmod())

    ############ Program information containers ############
    def genome():
        doc = "Genomes for individuals in the breeding program."
        def fget(self):
            """Get genomes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set genomes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete genomes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    genome = property(**genome())

    def geno():
        doc = "Genotypes for individuals in the breeding program."
        def fget(self):
            """Get genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    geno = property(**geno())

    def pheno():
        doc = "Phenotypes for individuals in the breeding program."
        def fget(self):
            """Get phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    pheno = property(**pheno())

    def bval():
        doc = "Breeding values for individuals in the breeding program."
        def fget(self):
            """Get breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    bval = property(**bval())

    def gmod():
        doc = "Genomic models for individuals in the breeding program."
        def fget(self):
            """Get genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    gmod = property(**gmod())

    ######### Breeding program operator properties #########
    def initop():
        doc = "Initialization operator."
        def fget(self):
            return self._initop
        def fset(self, value):
            # check_is_InitializationOperator(value, "initop")
            self._initop = value
        def fdel(self):
            del self._initop
        return locals()
    initop = property(**initop())

    def pselop():
        doc = "Parental selection operator."
        def fget(self):
            return self._pselop
        def fset(self, value):
            check_is_ParentSelectionOperator(value, "pselop")
            self._pselop = value
        def fdel(self):
            del self._pselop
        return locals()
    pselop = property(**pselop())

    def mateop():
        doc = "Mating operator."
        def fget(self):
            return self._mateop
        def fset(self, value):
            check_is_MatingOperator(value, "mateop")
            self._mateop = value
        def fdel(self):
            del self._mateop
        return locals()
    mateop = property(**mateop())

    def evalop():
        doc = "Evaluation operator."
        def fget(self):
            return self._evalop
        def fset(self, value):
            check_is_EvaluationOperator(value, "evalop")
            self._evalop = value
        def fdel(self):
            del self._evalop
        return locals()
    evalop = property(**evalop())

    def sselop():
        doc = "Survivor selection operator."
        def fget(self):
            return self._sselop
        def fset(self, value):
            check_is_SurvivorSelectionOperator(value, "sselop")
            self._sselop = value
        def fdel(self):
            del self._sselop
        return locals()
    sselop = property(**sselop())

    ############# Generation number properties #############
    def t_cur():
        doc = "Current generation number of the BreedingNode."
        def fget(self):
            return self._t_cur
        def fset(self, value):
            check_is_int(value, "t_cur")
            self._t_cur = value
        def fdel(self):
            del self._t_cur
        return locals()
    t_cur = property(**t_cur())

    def t_max():
        doc = "Maximum generation number of the BreedingNode."
        def fget(self):
            return self._t_max
        def fset(self, value):
            check_is_int(value, "t_max")
            self._t_max = value
        def fdel(self):
            del self._t_max
        return locals()
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
        self._start_geno, self._start_bval, self._start_gmod = self._initop.initialize(
            **kwargs
        )

    def is_initialized(self, **kwargs):
        """
        Return whether or not the BreedingProgram has been initialized with a
        starting set of conditions.

        Parameters
        ----------
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : boolean
            True if the BreedingProgram has been initialized.
            False if the BreedingProgram has not been initialized.
        """
        raise NotImplementedError("method is abstract")

    ################ Whole breeding program ################
    def advance(self, ngen, lbook, **kwargs):
        """
        Advance the breeding program by a specified number of generations.

        Parameters
        ----------
        ngen : int
            Number of generations to advance the BreedingProgram.
        lbook : Logbook
            Logbook into which to write statistics.
        **kwargs
            Additional keyword arguments.
        """
        # iterate through main breeding loop for ngen generations
        for t in range(self._t_cur, self._t_cur + ngen):
            ####################################################################
            ########################## select parents ##########################
            ####################################################################
            parent_gmat, sel, ncross, nprogeny, misc = self._pselop.pselect(
                t_cur = t,
                t_max = self._t_max,
                geno = self._geno,
                bval = self._bval,
                gmod = self._gmod
            )
            lbook.log_pselect(
                t_cur = t,
                t_max = self._t_max,
                pgvmat = parent_gmat,
                sel = sel,
                ncross = ncross,
                nprogeny = nprogeny,
                misc = misc
            )

            ####################################################################
            ########################### mate parents ###########################
            ####################################################################
            progeny_gmat, misc = self._mateop.mate(
                t_cur = t,
                t_max = self._t_max,
                pgvmat = parent_gmat,
                sel = sel,
                ncross = ncross,
                nprogeny = nprogeny
            )
            lbook.log_mate(
                t_cur = t,
                t_max = self._t_max,
                pgvmat = progeny_gmat,
                misc = misc
            )

            ####################################################################
            ####################### integrate genotypes ########################
            ####################################################################
            self._geno, misc = self._gintgop.gintegrate(
                t_cur = t,
                t_max = self._t_max,
                pgvmat = progeny_gmat,
                geno = self._geno,
            )
            lbook.log_gintegrate(
                t_cur = t,
                t_max = self._t_max,
                geno = self._geno,
                misc = misc
            )

            ####################################################################
            ######################## evaluate genotypes ########################
            ####################################################################
            bvmat, bvmat_true, misc = self._evalop.evaluate(
                t_cur = t,
                t_max = self._t_max,
                pgvmat = self._geno["main"],
                gmod_true = self._gmod["true"]
            )
            lbook.log_evaluate(
                t_cur = t,
                t_max = self._t_max,
                bvmat = bvmat,
                bvmat_true = bvmat_true,
                misc = misc
            )

            ####################################################################
            #################### integrate breeding values #####################
            ####################################################################
            self._bval, misc = self._bvintgop.bvintegrate(
                t_cur = t,
                t_max = self._t_max,
                bvmat = bvmat,
                bvmat_true = bvmat_true,
                bval = self._bval,
            )
            lbook.log_bvintegrate(
                t_cur = t,
                t_max = self._t_max,
                bval = self._bval,
                misc = misc
            )

            ####################################################################
            ######################### calibrate models #########################
            ####################################################################
            self._gmod, misc = self._calop.calibrate(
                t_cur = t,
                t_max = self._t_max,
                geno = self._geno,
                bval = self._bval,
                gmod = self._gmod
            )
            lbook.log_calibrate(
                t_cur = t,
                t_max = self._t_max,
                gmod = self._gmod,
                misc = misc
            )

            ####################################################################
            ######################### select survivors #########################
            ####################################################################
            self._geno, self._bval, self._gmod, misc = self._sselop.sselect(
                t_cur = t,
                t_max = self._t_max,
                geno = self._geno,
                bval = self._bval,
                gmod = self._gmod
            )
            lbook.log_sselect(
                t_cur = t,
                t_max = self._t_max,
                geno = self._geno,
                bval = self._bval,
                gmod = self._gmod,
                misc = misc
            )

            ####################################################################
            ######################### variable updates #########################
            ####################################################################
            # increment time variables
            self._t_cur += 1

    def evolve(self, nrep, ngen, lbook, loginit, verbose = False, **kwargs):
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
        if any(e is None for e in (self._start_geno, self._start_bval, self._start_gmod)):
            self.initialize()

        # main replication loop
        for r in range(1, nrep+1):
            # verbose messages
            if verbose:
                print("Simulating rep {0}".format(r))

            # copy starting pointers (okay since operators replace, not overwrite)
            self._geno = self._start_geno
            self._bval = self._start_bval
            self._gmod = self._start_gmod

            # set t_cur to zero (starting population)
            self._t_cur = 0

            # set rep in Logbook
            lbook.rep = r

            # cand_bvmat, cand_bvmat_true, misc = self._evalop.evaluate(
            #     t_cur = self._t_cur,
            #     t_max = self._t_max,
            #     pgvmat = self._geno["cand"],
            #     gmod_true = self._gmod["true"]
            # )
            #
            # self._bval["cand"] = cand_bvmat
            # self._bval["cand_true"] = cand_bvmat_true

            # evaluate main population starting genotypes using evalop
            main_bvmat, main_bvmat_true, misc = self._evalop.evaluate(
                t_cur = self._t_cur,
                t_max = self._t_max,
                pgvmat = self._geno["main"],
                gmod_true = self._gmod["true"]
            )

            self._bval["main"] = main_bvmat
            self._bval["main_true"] = main_bvmat_true

            # log initial conditions if needed
            if loginit:
                lbook.log_initialize(
                    t_cur = self._t_cur,
                    t_max = self._t_max,
                    geno = self._geno,
                    bval = self._bval,
                    gmod = self._gmod
                )

            # set t_cur to 1 (first generation)
            self._t_cur = 1

            # evolve the population
            self.advance(
                ngen = ngen,
                lbook = lbook,
                **kwargs
            )



################################################################################
################################## Utilities ###################################
################################################################################
def is_RecurrentSelectionBreedingProgram(v):
    """
    Determine whether an object is a RecurrentSelectionBreedingProgram.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a RecurrentSelectionBreedingProgram object instance.
    """
    return isinstance(v, RecurrentSelectionBreedingProgram)

def check_is_RecurrentSelectionBreedingProgram(v, varname):
    """
    Check if object is of type RecurrentSelectionBreedingProgram. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RecurrentSelectionBreedingProgram):
        raise TypeError("'%s' must be a RecurrentSelectionBreedingProgram." % varname)

def cond_check_is_RecurrentSelectionBreedingProgram(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type RecurrentSelectionBreedingProgram. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a RecurrentSelectionBreedingProgram.
    """
    if cond(v):
        check_is_RecurrentSelectionBreedingProgram(v, varname)
