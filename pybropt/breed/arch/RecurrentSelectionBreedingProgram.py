from . import BreedingProgram

from pybropt.core.error import check_is_int
from pybropt.core.error import check_is_dict
from pybropt.core.error import check_keys_in_dict
from pybropt.breed.init import check_is_InitializationOperator
from pybropt.breed.sel import check_is_ParentSelectionOperator
from pybropt.breed.mate import check_is_MatingOperator
from pybropt.breed.intg import check_is_GenotypeIntegrationOperator
from pybropt.breed.eval import check_is_EvaluationOperator
from pybropt.breed.intg import check_is_BreedingValueIntegrationOperator
from pybropt.breed.calibr import check_is_GenomicModelCalibrationOperator
from pybropt.breed.sel import check_is_SurvivorSelectionOperator

class RecurrentSelectionBreedingProgram(BreedingProgram):
    """docstring for RecurrentSelectionBreedingProgram."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, t_max, initop, pselop, mateop, gintgop, evalop, bvintgop, calop, sselop, **kwargs):
        super(RecurrentSelectionBreedingProgram, self).__init__(
            t_max = t_max,
            initop = initop,
            pselop = pselop,
            mateop = mateop,
            evalop = evalop,
            intgop = intgop,
            calop = calop,
            sselop = sselop
            **kwargs
        )

        # set time variables
        self.t_cur = 0
        self.t_max = t_max

        # save operators
        self.initop = initop
        self.pselop = pselop
        self.mateop = mateop
        self.gintgop = gintgop
        self.evalop = evalop
        self.bvintgop = bvintgop
        self.calop = calop
        self.sselop = sselop

        # TODO: go through set methods properly
        self._geno = None
        self._bval = None
        self._gmod = None

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

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

    ################ Population properties #################
    def geno():
        doc = "Main breeding population of the BreedingNode."
        def fget(self):
            return self._pop
        def fset(self, value):
            check_is_dict(value, "geno")
            check_keys_in_dict(value, "geno", "cand", "main", "queue")
            self._pop = value
        def fdel(self):
            del self._pop
        return locals()
    geno = property(**geno())

    ############## Breeding value properties ###############
    def bval():
        doc = "Estimated breeding values for the main breeding population of the BreedingNode."
        def fget(self):
            return self._bval
        def fset(self, value):
            check_is_dict(value, "bval")
            check_keys_in_dict(value, "bval", "cand", "cand_true", "main", "main_true", "queue", "queue_true")
            self._bval = value
        def fdel(self):
            del self._bval
        return locals()
    bval = property(**bval())

    ############### Genomic model properties ###############
    def gmod():
        doc = "Estimated genomic model for the main breeding population of the BreedingNode."
        def fget(self):
            return self._gmod
        def fset(self, value):
            check_is_dict(value, "gmod")
            check_keys_in_dict(value, "gmod", "cand", "main", "queue", "true")
            self._gmod = value
        def fdel(self):
            del self._gmod
        return locals()
    gmod = property(**gmod())

    ######### Breeding program operator properties #########
    def initop():
        doc = "Initialization operator."
        def fget(self):
            return self._initop
        def fset(self, value):
            check_is_InitializationOperator(value, "initop")
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

    def gintgop():
        doc = "Genotype integration operator."
        def fget(self):
            return self._gintgop
        def fset(self, value):
            check_is_GenotypeIntegrationOperator(value, "gintgop")
            self._gintgop = value
        def fdel(self):
            del self._gintgop
        return locals()
    gintgop = property(**gintgop())

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

    def bvintgop():
        doc = "Breeding value integration operator."
        def fget(self):
            return self._bvintgop
        def fset(self, value):
            check_is_BreedingValueIntegrationOperator(value, "bvintgop")
            self._bvintgop = value
        def fdel(self):
            del self._bvintgop
        return locals()
    bvintgop = property(**bvintgop())

    def calop():
        doc = "Genomic model calibration operator."
        def fget(self):
            return self._calop
        def fset(self, value):
            check_is_GenomicModelCalibrationOperator(value, "calop")
            self._calop = value
        def fdel(self):
            del self._calop
        return locals()
    calop = property(**calop())

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

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############# Initialize breeding program ##############
    def initialize(self, **kwargs):
        """
        Initialize the breeding program with genotypes, phenotypes, and genomic
        models.
        """
        self.geno, self.bval, self.gmod = self.initop.initialize(
            **kwargs
        )

    ################ Whole breeding program ################
    def evolve(self, ngen, lbook, loginit, **kwargs):
        """
        Evolve the breeding program for a number of generations.

        Parameters
        ----------
        ngen : int
            Number of generations to evolve the population.
        lbook : LogBook
            LogBook into which to write statistics.
        loginit : bool
            Whether to log the initial state before main loop.
        """
        # initialize if needed
        if any(e is None for e in (self._geno, self._bval, self._gmod)):
            self.initialize()

        # log initial conditions if needed
        if loginit:
            lbook.log_initialize(
                t_cur = self._t_cur,
                t_max = self._t_max,
                geno = self._geno,
                bval = self._bval,
                gmod = self._gmod
            )

        # main breeding loop
        for _ in range(ngen):
            ####################################################################
            ########################## select parents ##########################
            ####################################################################
            pgvmat, sel, ncross, nprogeny, misc = self.pselop.pselect(
                t_cur = self._t_cur,
                t_max = self._t_max,
                geno = self._geno,
                bval = self._bval,
                gmod = self._gmod
            )
            lbook.log_pselect(
                pgvmat = pgvmat,
                sel = sel,
                ncross = ncross,
                nprogeny = nprogeny,
                misc = misc
            )

            ####################################################################
            ########################### mate parents ###########################
            ####################################################################
            pgvmat, misc = self.mateop.mate(
                t_cur = self._t_cur,
                t_max = self._t_max,
                pgvmat = pgvmat,
                sel = sel,
                ncross = ncross,
                nprogeny = nprogeny
            )
            lbook.log_mate(
                t_cur = self._t_cur,
                t_max = self._t_max,
                pgvmat = pgvmat,
                misc = misc
            )

            ####################################################################
            ####################### integrate genotypes ########################
            ####################################################################
            geno_tmp, misc = self.gintgop.gintegrate(
                t_cur = self._t_cur,
                t_max = self._t_max,
                pgvmat = pgvmat,
                geno = self._geno,
            )
            lbook.log_gintegrate(
                t_cur = self._t_cur,
                t_max = self._t_max,
                geno = geno_tmp,
                bval = bval_tmp,
                misc = misc
            )

            ####################################################################
            ######################## evaluate genotypes ########################
            ####################################################################
            bvmat, bvmat_true, misc = self.evalop.evaluate(
                t_cur = self._t_cur,
                t_max = self._t_max,
                pgvmat = geno_tmp["main"],
                gmod_true = self._gmod["true"]
            )
            lbook.log_evaluate(
                t_cur = self._t_cur,
                t_max = self._t_max,
                bvmat = bvmat,
                bvmat_true = bvmat_true,
                misc = misc
            )

            ####################################################################
            #################### integrate breeding values #####################
            ####################################################################
            bval_tmp, misc = self.bvintgop.bvintegrate(
                t_cur = self._t_cur,
                t_max = self._t_max,
                bvmat = bvmat,
                bvmat_true = bvmat_true,
                bval = self._geno,
            )
            lbook.log_bvintegrate(
                t_cur = self._t_cur,
                t_max = self._t_max,
                bval = bval_tmp,
                misc = misc
            )

            ####################################################################
            ######################### calibrate models #########################
            ####################################################################
            gmod_tmp, misc = self.calop.calibrate(
                t_cur = self._t_cur,
                t_max = self._t_max,
                geno = geno_tmp,
                bval = bval_tmp,
                gmod = self._gmod
            )
            lbook.log_calibrate(
                t_cur = self._t_cur,
                t_max = self._t_max,
                gmod = gmod_tmp,
                misc = misc
            )

            ####################################################################
            ######################### select survivors #########################
            ####################################################################
            geno_new, bval_new, gmod_new, misc = self.sselop.sselect(
                t_cur = self._t_cur,
                t_max = self._t_max,
                geno = geno_tmp,
                bval = bval_tmp,
                gmod = gmod_tmp
            )
            lbook.log_sselect(
                t_cur = self._t_cur,
                t_max = self._t_max,
                geno = geno_new,
                bval = bval_new,
                gmod = gmod_new,
                misc = misc
            )

            ####################################################################
            ######################### variable updates #########################
            ####################################################################
            # update population variables
            self.geno = geno_new
            self.bval = bval_new
            self.gmod = gmod_new

            # increment time variables
            self._t_cur += 1



################################################################################
################################## Utilities ###################################
################################################################################
def is_RecurrentSelectionBreedingProgram(v):
    return isinstance(v, RecurrentSelectionBreedingProgram)

def check_is_RecurrentSelectionBreedingProgram(v, vname):
    if not isinstance(v, RecurrentSelectionBreedingProgram):
        raise TypeError("variable '{0}' must be a RecurrentSelectionBreedingProgram".format(vname))

def cond_check_is_BreedingProgram(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_RecurrentSelectionBreedingProgram(v, vname)
