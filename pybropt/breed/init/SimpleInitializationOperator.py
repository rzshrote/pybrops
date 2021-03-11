import numpy

from . import InitializationOperator

from pybropt.core.error import check_is_int
from pybropt.core.error import check_is_Generator
from pybropt.core.error import check_is_dict
from pybropt.core.error import check_keys_in_dict

from pybropt.breed.sel import check_is_ParentSelectionOperator
from pybropt.breed.mate import check_is_MatingOperator
from pybropt.breed.intg import check_is_GenotypeIntegrationOperator
from pybropt.breed.eval import check_is_EvaluationOperator
from pybropt.breed.intg import check_is_BreedingValueIntegrationOperator
from pybropt.breed.calibr import check_is_GenomicModelCalibrationOperator
from pybropt.breed.sel import check_is_SurvivorSelectionOperator

from pybropt.model.gmod import check_is_GenomicModel

from pybropt.popgen.gmat import DensePhasedGenotypeVariantMatrix
from pybropt.breed.eval import NoGxEEvaluationOperator

class SimpleInitializationOperator(InitializationOperator):
    """docstring for SimpleInitializationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, burnin, founder_geno, founder_bval, founder_gmod, gmult, pselop, mateop, gintgop, evalop, bvintgop, calop, sselop, **kwargs):
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

        self.gmult = gmult
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
    def initialize(self, **kwargs):
        """
        Initialize a breeding program.

        Parameters
        ----------
        **kwargs

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
                queue      | List of BreedingValueMatrix | Breeding values for populations on queue
                queue_true | List of BreedingValueMatrix | True breeding values for populations on queue
            gmod : dict
                Field | Type                 | Description
                ------+----------------------+----------------------------------
                cand  | GenomicModel         | Parental candidate breeding population genomic model
                main  | GenomicModel         | Main breeding population genomic model
                queue | List of GenomicModel | Genomic models for populations on queue
                true  | GenomicModel         | True genomic model for trait(s)
        """
        geno = self.founder_geno
        bval = self.founder_bval
        gmod = self.founder_gmod
        # print("cand:", geno["cand"].taxa_grp)

        for t in range(-(self.burnin-1), 1):
            # print("################################################################################")
            # print("iteration:", t)
            ####################################################################
            ########################## select parents ##########################
            ####################################################################
            # print("main:", geno["main"].taxa_grp)
            # print("cand:", geno["cand"].taxa_grp)
            pgvmat, sel, ncross, nprogeny, misc = self.pselop.pselect(
                t_cur = t,
                t_max = 0,
                geno = geno,
                bval = bval,
                gmod = gmod
            )
            # print("pgvmat:", pgvmat.taxa_grp)
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
            # print("gintegrate:",geno["main"].mat.shape)
            ####################################################################
            ######################## evaluate genotypes ########################
            ####################################################################
            bvmat, bvmat_true, misc = self.evalop.evaluate(
                t_cur = t,
                t_max = 0,
                pgvmat = geno["main"],
                gmod_true = gmod["true"]
            )
            # print("evaluate:",geno["main"].mat.shape)
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
            # print("bvintegrate:",bval["main"].mat.shape)
            ####################################################################
            ######################### calibrate models #########################gmod
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
            geno_new, bval_new, gmod_new, misc = self.sselop.sselect(
                t_cur = t,
                t_max = 0,
                geno = geno,
                bval = bval,
                gmod = gmod
            )

            ####################################################################
            ######################### variable updates #########################
            ####################################################################
            # update population variables
            geno = geno_new
            bval = bval_new
            gmod = gmod_new
            # print("cand:", geno["cand"].taxa_grp)
            # print("cand:", geno["cand"].mat.shape)
            # cand_mean = bval["cand"].mat.mean(0)
            # print("cand mean:", cand_mean)
            # print("cand mean sum:", cand_mean.sum())

        # HACK: # FIXME:
        offset = ((geno["main"].taxa_grp.max() // self.gmult) ) * (self.gmult)

        # adjust generations to zero
        geno["cand"].taxa_grp = geno["cand"].taxa_grp - offset
        geno["main"].taxa_grp = geno["main"].taxa_grp - offset
        for i in range(len(geno["queue"])):
            geno["queue"][i].taxa_grp = geno["queue"][i].taxa_grp - offset
        bval["cand"].taxa_grp = bval["cand"].taxa_grp - offset
        bval["main"].taxa_grp = bval["main"].taxa_grp - offset

        return geno, bval, gmod

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    @staticmethod
    def from_dpgvmat(dpgvmat, rng, nfounder, founder_ncross, founder_nprogeny, gmult, gmod_true, burnin, pselop, mateop, gintgop, evalop, bvintgop, calop, sselop):
        """
        Create an initialization operator from a DensePhasedGenotypeVariantMatrix.

        Initialization is conducted as follows:
            1) 'nfounder' lines are selected as founders.
            2) founders are randomly intermated to create base population.

        Parameters
        ----------
        dpgvmat : DensePhasedGenotypeVariantMatrix
        rng : numpy.Generator
        nfounder : int
            Number of founders to select. Sampling is done without replacement.
        """
        # perform error checks
        check_is_Generator(rng, "rng")
        check_is_GenomicModel(gmod_true, "gmod_true")

        check_is_int(burnin, "burnin")

        check_is_ParentSelectionOperator(pselop, "pselop")
        check_is_MatingOperator(mateop, "mateop")
        check_is_GenotypeIntegrationOperator(gintgop, "gintgop")
        check_is_EvaluationOperator(evalop, "evalop")
        check_is_BreedingValueIntegrationOperator(bvintgop, "bvintgop")
        check_is_GenomicModelCalibrationOperator(calop, "calop")
        check_is_SurvivorSelectionOperator(sselop, "sselop")

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
        founder_bval["main"] = bvmat
        founder_bval["main_true"] = bvmat_true

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
            gmult = gmult,
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
    def from_vcf(fname, rng, nfounder, founder_ncross, founder_nprogeny, gmult, gmod_true, burnin, pselop, mateop, gintgop, evalop, bvintgop, calop, sselop):
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
        dpgvmat = DensePhasedGenotypeVariantMatrix.from_vcf(fname)

        # step 2: create from genotype matrix
        geninitop = SimpleInitializationOperator.from_dpgvmat(
            dpgvmat = dpgvmat,
            rng = rng,
            nfounder = nfounder,
            founder_ncross = founder_ncross,
            founder_nprogeny = founder_nprogeny,
            gmult = gmult,
            gmod_true = gmod_true,
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
