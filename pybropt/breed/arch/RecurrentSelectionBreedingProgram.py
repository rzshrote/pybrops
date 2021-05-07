import copy
import numpy

from . import BreedingProgram

from pybropt.core.error import check_is_int
from pybropt.core.error import check_is_dict
from pybropt.core.error import check_keys_in_dict
from pybropt.breed.init import check_is_InitializationOperator
from pybropt.breed.psel import check_is_ParentSelectionOperator
from pybropt.breed.mate import check_is_MatingOperator
from pybropt.breed.intg import check_is_GenotypeIntegrationOperator
from pybropt.breed.eval import check_is_EvaluationOperator
from pybropt.breed.intg import check_is_BreedingValueIntegrationOperator
from pybropt.breed.calibr import check_is_GenomicModelCalibrationOperator
from pybropt.breed.ssel import check_is_SurvivorSelectionOperator

# def dprint(t, lab, d, mod, ID = True, USL = True):
#     if ID:
#         x = "cand:" + str(id(d["cand"]))
#         y = "main:" + str(id(d["main"]))
#         z = "queue:" + str([id(a) for a in d["queue"]])
#         print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(t, "ID", lab, x, y, z))
#     if USL:
#         x = "cand:" + str(mod.usl(d["cand"]))
#         y = "main:" + str(mod.usl(d["main"]))
#         z = "queue:" + str([mod.usl(a) for a in d["queue"]])
#         print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(t, "USL", lab, x, y, z))

class RecurrentSelectionBreedingProgram(BreedingProgram):
    """docstring for RecurrentSelectionBreedingProgram."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, t_max, initop, pselop, mateop, gintgop, evalop, bvintgop, calop, sselop, **kwargs):
        super(RecurrentSelectionBreedingProgram, self).__init__(**kwargs)

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
        self._start_geno = None
        self._start_bval = None
        self._start_gmod = None
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
    def start_geno():
        doc = "The start_geno property."
        def fget(self):
            return self._start_geno
        def fset(self, value):
            self._start_geno = value
        def fdel(self):
            del self._start_geno
        return locals()
    start_geno = property(**start_geno())

    def geno():
        doc = "Main breeding population of the BreedingNode."
        def fget(self):
            return self._geno
        def fset(self, value):
            check_is_dict(value, "geno")
            # TODO:
            # check_keys_in_dict(value, "geno", "cand", "main", "queue")
            self._geno = value
        def fdel(self):
            del self._geno
        return locals()
    geno = property(**geno())

    ############## Breeding value properties ###############
    def start_bval():
        doc = "The start_bval property."
        def fget(self):
            return self._start_bval
        def fset(self, value):
            self._start_bval = value
        def fdel(self):
            del self._start_bval
        return locals()
    start_bval = property(**start_bval())

    def bval():
        doc = "Estimated breeding values for the main breeding population of the BreedingNode."
        def fget(self):
            return self._bval
        def fset(self, value):
            check_is_dict(value, "bval")
            # TODO:
            # check_keys_in_dict(value, "bval", "cand", "cand_true", "main", "main_true", "queue", "queue_true")
            self._bval = value
        def fdel(self):
            del self._bval
        return locals()
    bval = property(**bval())

    ############### Genomic model properties ###############
    def start_gmod():
        doc = "The start_gmod property."
        def fget(self):
            return self._start_gmod
        def fset(self, value):
            self._start_gmod = value
        def fdel(self):
            del self._start_gmod
        return locals()
    start_gmod = property(**start_gmod())

    def gmod():
        doc = "Estimated genomic model for the main breeding population of the BreedingNode."
        def fget(self):
            return self._gmod
        def fset(self, value):
            check_is_dict(value, "gmod")
            # TODO:
            # check_keys_in_dict(value, "gmod", "cand", "main", "queue", "true")
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
        self._start_geno, self._start_bval, self._start_gmod = self._initop.initialize(
            **kwargs
        )

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
    return isinstance(v, RecurrentSelectionBreedingProgram)

def check_is_RecurrentSelectionBreedingProgram(v, vname):
    if not isinstance(v, RecurrentSelectionBreedingProgram):
        raise TypeError("variable '{0}' must be a RecurrentSelectionBreedingProgram".format(vname))

def cond_check_is_BreedingProgram(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_RecurrentSelectionBreedingProgram(v, vname)
