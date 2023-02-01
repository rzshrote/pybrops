import numpy

from . import InitializationOperator

from pybrops.core import random as pbo_rng
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_keys_in_dict
from pybrops.core.error import cond_check_is_Generator
from pybrops.breed.calibr import check_is_GenomicModelCalibrationOperator
from pybrops.breed.eval import check_is_EvaluationOperator
from pybrops.breed.intg import check_is_GenotypeIntegrationOperator
from pybrops.breed.intg import check_is_BreedingValueIntegrationOperator
from pybrops.breed.mate import check_is_MatingOperator
from pybrops.breed.prot.sel import check_is_ParentSelectionOperator
from pybrops.breed.ssel import check_is_SurvivorSelectionOperator
from pybrops.model.gmod import check_is_GenomicModel
from pybrops.popgen.gmat import DensePhasedGenotypeMatrix
from pybrops.popgen.gmat import check_is_DensePhasedGenotypeMatrix

class SimpleInitializationOperator(InitializationOperator):
    """docstring for SimpleInitializationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, burnin, founder_geno, founder_bval, founder_gmod, pselop, mateop, gintgop, evalop, bvintgop, calop, sselop, **kwargs: dict):
        """
        Parameters
        ----------
        burnin : int
            Number of generations to burnin.
        founder_geno : dict
        founder_bval : dict
        founder_gmod : dict

        """
        super(SimpleInitializationOperator, self).__init__(**kwargs)

        # error checks
        check_is_int(burnin, "burnin")

        check_is_dict(founder_geno, "founder_geno")
        check_keys_in_dict(founder_geno, "founder_geno", "cand", "main", "queue")
        check_is_dict(founder_bval, "founder_bval")
        check_keys_in_dict(founder_bval, "founder_bval", "cand", "cand_true", "main", "main_true")
        check_is_dict(founder_gmod, "founder_gmod")
        check_keys_in_dict(founder_gmod, "founder_gmod", "cand", "main", "true")

        check_is_ParentSelectionOperator(pselop, "pselop")
        check_is_MatingOperator(mateop, "mateop")
        check_is_GenotypeIntegrationOperator(gintgop, "gintgop")
        check_is_EvaluationOperator(evalop, "evalop")
        check_is_BreedingValueIntegrationOperator(bvintgop, "bvintgop")
        check_is_GenomicModelCalibrationOperator(calop, "calop")
        check_is_SurvivorSelectionOperator(sselop, "sselop")

        # variable assign
        self.burnin = burnin

        self.founder_geno = founder_geno
        self.founder_bval = founder_bval
        self.founder_gmod = founder_gmod

        self.pselop = pselop
        self.mateop = mateop
        self.gintgop = gintgop
        self.evalop = evalop
        self.bvintgop = bvintgop
        self.calop = calop
        self.sselop = sselop

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def initialize(self, burnin = None, **kwargs: dict):
        """
        Initialize a breeding program.

        Parameters
        ----------
        burnin : int, None
            Number of generations to burnin. If None, use default burnin
            generations.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing three dict objects: (geno, bval, gmod)
            Elements of this tuple are as follows:
            Element | Description
            --------+-----------------------------------
            geno    | A dict of genotypic data.
            bval    | A dict of breeding value data.
            gmod    | A dict of genomic models.

            Dictionaries must have the following fields:
            ============================================
            geno : dict
                Field | Type                         | Description
                ------+------------------------------+--------------------------
                cand  | PhasedGenotypeMatrix         | Parental candidate breeding population
                main  | PhasedGenotypeMatrix         | Main breeding population
                queue | List of PhasedGenotypeMatrix | Breeding populations on queue
                ""
            bval : dict
                Field      | Type                        | Description
                -----------+-----------------------------+----------------------
                cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                main       | BreedingValueMatrix         | Main breeding population breeding values
                main_true  | BreedingValueMatrix         | Main breeding population true breeding values
            gmod : dict
                Field | Type                 | Description
                ------+----------------------+----------------------------------
                cand  | GenomicModel         | Parental candidate breeding population genomic model
                main  | GenomicModel         | Main breeding population genomic model
                queue | List of GenomicModel | Genomic models for populations on queue
                true  | GenomicModel         | True genomic model for trait(s)
        """
        # get number of burnin generations
        if burnin is None:
            burnin = self.burnin

        # copy dict's and any list's within dict
        geno = dict((p[0],list(p[1])) if isinstance(p[1],list) else p for p in self.founder_geno.items())
        bval = dict((p[0],list(p[1])) if isinstance(p[1],list) else p for p in self.founder_bval.items())
        gmod = dict((p[0],list(p[1])) if isinstance(p[1],list) else p for p in self.founder_gmod.items())

        for t in range(-(self.burnin-1), 1):
            ####################################################################
            ########################## select parents ##########################
            ####################################################################
            pgvmat, sel, ncross, nprogeny, misc = self.pselop.pselect(
                t_cur = t,
                t_max = 0,
                geno = geno,
                bval = bval,
                gmod = gmod
            )

            ####################################################################
            ########################### mate parents ###########################
            ####################################################################
            pgvmat, misc = self.mateop.mate(
                t_cur = t,
                t_max = 0,
                pgvmat = pgvmat,
                sel = sel,
                ncross = ncross,
                nprogeny = nprogeny
            )

            ####################################################################
            ####################### integrate genotypes ########################
            ####################################################################
            geno, misc = self.gintgop.gintegrate(
                t_cur = t,
                t_max = 0,
                pgvmat = pgvmat,
                geno = geno,
            )

            ####################################################################
            ######################## evaluate genotypes ########################
            ####################################################################
            bvmat, bvmat_true, misc = self.evalop.evaluate(
                t_cur = t,
                t_max = 0,
                pgvmat = geno["main"],
                gmod_true = gmod["true"]
            )

            ####################################################################
            #################### integrate breeding values #####################
            ####################################################################
            bval, misc = self.bvintgop.bvintegrate(
                t_cur = t,
                t_max = 0,
                bvmat = bvmat,
                bvmat_true = bvmat_true,
                bval = bval,
            )

            ####################################################################
            ######################### calibrate models #########################
            ####################################################################
            gmod, misc = self.calop.calibrate(
                t_cur = t,
                t_max = 0,
                geno = geno,
                bval = bval,
                gmod = gmod
            )

            ####################################################################
            ######################### select survivors #########################
            ####################################################################
            geno, bval, gmod, misc = self.sselop.sselect(
                t_cur = t,
                t_max = 0,
                geno = geno,
                bval = bval,
                gmod = gmod
            )

            ####################################################################
            ######################### variable updates #########################
            ####################################################################
            # nothing to do!

        return geno, bval, gmod

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    @staticmethod
    def from_dpgvmat(dpgvmat, nfounder, founder_ncross, founder_nprogeny, gmod_true, burnin, pselop, mateop, gintgop, evalop, bvintgop, calop, sselop, rng = None):
        """
        Create an initialization operator from a DensePhasedGenotypeMatrix.

        Initialization is conducted as follows:
            1) 'nfounder' lines are selected as founders.
            2) founders are randomly intermated to create base population.

        Parameters
        ----------
        dpgvmat : DensePhasedGenotypeMatrix
        rng : numpy.Generator
        nfounder : int
            Number of founders to select. Sampling is done without replacement.
        """
        # perform error checks
        check_is_DensePhasedGenotypeMatrix(dpgvmat, "dpgvmat")
        check_is_int(nfounder, "nfounder")
        check_is_int(founder_ncross, "founder_ncross")
        check_is_int(founder_nprogeny, "founder_nprogeny")
        check_is_GenomicModel(gmod_true, "gmod_true")
        check_is_int(burnin, "burnin")
        check_is_ParentSelectionOperator(pselop, "pselop")
        check_is_MatingOperator(mateop, "mateop")
        check_is_GenotypeIntegrationOperator(gintgop, "gintgop")
        check_is_EvaluationOperator(evalop, "evalop")
        check_is_BreedingValueIntegrationOperator(bvintgop, "bvintgop")
        check_is_GenomicModelCalibrationOperator(calop, "calop")
        check_is_SurvivorSelectionOperator(sselop, "sselop")
        cond_check_is_Generator(rng, "rng")

        # assign default random number generator
        if rng is None:
            rng = pbo_rng

        ####################################################
        ### step 1: count available taxa ###
        ####################################################
        ntaxa = dpgvmat.ntaxa

        ###################################################################
        ### step 2: define and populate founder_geno, founder_bval, founder_gmod ###
        ###################################################################
        # define founder_geno
        founder_geno = {
            "cand" : None,
            "main" : None,
            "queue" : []
        }

        # define founder_bval
        founder_bval = {
            "cand" : None,
            "cand_true" : None,
            "main" : None,
            "main_true" : None
        }

        # define founder_gmod
        founder_gmod = {
            "cand" : gmod_true,
            "main" : gmod_true,
            "true" : gmod_true
        }

        # pselect
        # mate
        # gintegrate
        # evaluate
        # bvintegrate
        # calibrate
        # sselect
        #######################################
        ### Main populataion initialization ###
        #######################################
        # get random selections (pselect)
        sel = rng.choice(ntaxa, nfounder, replace = False)

        # mate random selections (mate)
        pgvmat, misc = mateop.mate(
            t_cur = -burnin,
            t_max = 0,
            pgvmat = dpgvmat,
            sel = sel,
            ncross = founder_ncross,
            nprogeny = founder_nprogeny
        )

        # assign genotypes (gintegrate)
        founder_geno["main"] = pgvmat

        #################################
        ### Evaluate main populataion ###
        #################################
        # phenotype main population (evaluate)
        bvmat, bvmat_true, misc = evalop.evaluate(
            t_cur = -burnin,
            t_max = 0,
            pgvmat = founder_geno["main"],
            gmod_true = founder_gmod["true"]
        )

        # assign breeding values (bvintegrate)
        founder_bval, misc = bvintgop.bvintegrate(
            t_cur = -burnin,
            t_max = 0,
            bvmat = bvmat,
            bvmat_true = bvmat_true,
            bval = founder_bval,
        )

        # calibrate genomic model (calibrate)
        founder_gmod, misc = calop.calibrate(
            t_cur = -burnin,
            t_max = 0,
            geno = founder_geno,
            bval = founder_bval,
            gmod = founder_gmod
        )

        # select survivors to fill breeding candidates (sselect)
        founder_geno, founder_bval, founder_gmod, misc = sselop.sselect(
            t_cur = -burnin,
            t_max = 0,
            geno = founder_geno,
            bval = founder_bval,
            gmod = founder_gmod
        )

        geninitop = SimpleInitializationOperator(
            founder_geno = founder_geno,
            founder_bval = founder_bval,
            founder_gmod = founder_gmod,
            burnin = burnin,
            pselop = pselop,
            mateop = mateop,
            gintgop = gintgop,
            evalop = evalop,
            bvintgop = bvintgop,
            calop = calop,
            sselop = sselop,
        )

        return geninitop

    @staticmethod
    def from_vcf(fname, nfounder, founder_ncross, founder_nprogeny, gmod_true, burnin, pselop, mateop, gintgop, evalop, bvintgop, calop, sselop, rng = None):
        """
        Create a SimpleInitializationOperator from a VCF file.

        Initializes a "main" population with genotypes, and queue populations
        of various lengths. Uses the individuals in "main" to select a set of
        parental candidates using a provided SurvivorSelectionOperator. Then,
        a provided ParentSelectionOperator, and Mating operator is used to
        select and mate parents to create one additional generation that is
        added to the queue.

        Parameters
        ----------
        fname : str
            VCF file name.
        size : dict
            Field | Type        | Description
            ------+-------------+-----------------------------------------------
            cand  | None        | None
            main  | int         | Number of taxa in main breeding population.
            queue | list of int | Number of taxa in breeding populations on queue.
        rng : numpy.random.Generator
            A random number generator object.
        replace : bool, default = True
            Whether genotype sampling is with or without replacement.
            If replace == False:
                If the number of genotypes in the provided file is less than
                the sum of the required initial number of genotypes
                (main + sum(queue)), then sample all genotypes and fill
                remaining genotypes with genotypes sampled without replacement.
        """
        # step 1: load genotype matrix
        dpgvmat = DensePhasedGenotypeMatrix.from_vcf(fname)

        # step 2: create from genotype matrix
        geninitop = SimpleInitializationOperator.from_dpgvmat(
            dpgvmat = dpgvmat,
            nfounder = nfounder,
            founder_ncross = founder_ncross,
            founder_nprogeny = founder_nprogeny,
            gmod_true = gmod_true,
            burnin = burnin,
            pselop = pselop,
            mateop = mateop,
            gintgop = gintgop,
            evalop = evalop,
            bvintgop = bvintgop,
            calop = calop,
            sselop = sselop,
            rng = rng
        )

        return geninitop



################################################################################
################################## Utilities ###################################
################################################################################
def is_SimpleInitializationOperator(v):
    return isinstance(v, SimpleInitializationOperator)

def check_is_SimpleInitializationOperator(v, varname):
    if not isinstance(v, SimpleInitializationOperator):
        raise TypeError("'%s' must be a SimpleInitializationOperator." % varname)

def cond_check_is_SimpleInitializationOperator(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_SimpleInitializationOperator(v, varname)
