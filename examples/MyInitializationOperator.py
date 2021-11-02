import pybropt.core.random
from pybropt.breed.op.init import InitializationOperator
from pybropt.breed.op.psel import ParentSelectionOperator
from pybropt.breed.op.mate import MatingOperator
from pybropt.breed.op.eval import EvaluationOperator
from pybropt.breed.op.ssel import SurvivorSelectionOperator

from pybropt.breed.prot.sel import ConventionalPhenotypicSelection
from pybropt.breed.prot.sel.transfn import trans_sum
from pybropt.breed.prot.mate import GenerationalTwoWayDHCross

class MyInitParentSelectionOperator(ParentSelectionOperator):
    def __init__(self, nparent, ncross, nprogeny, **kwargs):
        super(InitParentSelectionOperator, self).__init__(**kwargs)
        self.cps = ConventionalPhenotypicSelection( # create the selection protocol
            nparent = nparent,                      # number of parents
            ncross = ncross,                        # number of times to perform cross
            nprogeny = nprogeny,                    # number of progeny per cross
            objfn_trans = trans_sum                 # sum across all traits (equal weight)
        )
    def pselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        pgmat, sel, ncross, nprogeny, misc = self.cps.select(   # make selections from candidates
            pgmat = genome["cand"],                             # not used by CPS
            gmat = geno["cand"],                                # not used by CPS
            ptdf = pheno["cand"],                               # not used by CPS
            bvmat = bval["cand"],                               # breeding values used by CPS
            gpmod = gmod["cand"],                               # not used by CPS
            t_cur = t_cur,
            t_max = t_max,
            method = "single"                                   # just use equal weights for each trait
        )
        mcfg = {                    # store parent selections configurations in dict
            "pgmat" = pgmat,
            "sel" = sel,
            "ncross" = ncross,
            "nprogeny" = nprogeny
        }
        return mcfg, genome, geno, pheno, bval, gmod, misc

class MyInitMatingOperator(MatingOperator):
    def __init__(self, **kwargs):
        super(MyInitMatingOperator, self).__init__(**kwargs)
        self.mprot = GenerationalTwoWayDHCross()
    def mate(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        progeny = self.mprot.mate(**mcfg, s = 0)    # mate parents
        genome["progeny"] = progeny                 # add progeny to genome dict
        misc = {}
        return genome, geno, pheno, bval, gmod, misc

class MyInitEvaluationOperator(EvaluationOperator):
    def __init__(self, **kwargs):
        super(MyInitEvaluationOperator, self).__init__(**kwargs)
        self.pt = G_E_Phenotyping(
            gpmod = gmod["true"],
            nenv = 4,

        )

    def evaluate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        progeny_bv = gmod["true"].gebv(genome["progeny"])
        pass


class MyInitializationOperator(InitializationOperator):
    def __init__(self, arg):
        super(MyInitializationOperator, self).__init__()
        self.pselop = MyInitParentSelectionOperator(20, 1, 8)
        self.mateop = MyInitMatingOperator()
        self.evalop = MyInitEvaluationOperator()
        self.arg = arg
    @staticmethod
    def from_dpgmat(dpgmat, nfounder, founder_ncross, founder_nprogeny, gqlen, gmod_true, burnin):
        # perform error checks
        check_is_DensePhasedGenotypeMatrix(dpgmat, "dpgmat")
        check_is_int(nfounder, "nfounder")
        check_is_int(founder_ncross, "founder_ncross")
        check_is_int(founder_nprogeny, "founder_nprogeny")
        check_is_int(gqlen, "gqlen")
        check_is_GenomicModel(gmod_true, "gmod_true")
        check_is_int(burnin, "burnin")
        rng = pybropt.core.random

        ####################################################
        ### step 1: count available taxa ###
        ####################################################
        ntaxa = dpgmat.ntaxa

        ###################################################################
        ### step 2: define and populate founder_geno, founder_bval, founder_gmod ###
        ###################################################################
        # define founding dictionaries
        founder_genome = {"cand":None,      "main":None,      "queue":[]}
        founder_geno =   {"cand":None,      "main":None,      "queue":[]}
        founder_pheno =  {"cand":None,      "main":None}
        founder_bval =   {"cand":None,      "main":None}
        founder_gmod =   {"cand":gmod_true, "main":gmod_true, "true":gmod_true}

        #######################################
        ### Main populataion initialization ###
        #######################################
        # get random selections (pselect)
        sel = rng.choice(ntaxa, nfounder, replace = False)

        # fill queue with random matings of founders
        for t in range(-(burnin+gqlen), -burnin):
            # mate random selections (mate)
            pgvmat, misc = mateop.mate(
                t_cur = t,
                t_max = 0,
                pgvmat = dpgmat,
                sel = sel,
                ncross = founder_ncross,
                nprogeny = founder_nprogeny
            )

            # add progeny to queue
            founder_geno["queue"].append(pgvmat)

            # randomly shuffle founders again to get new crosses
            rng.shuffle(sel)

        # mate random selections (mate)
        pgvmat, misc = mateop.mate(
            t_cur = -burnin,
            t_max = 0,
            pgvmat = dpgmat,
            sel = sel,
            ncross = founder_ncross,
            nprogeny = founder_nprogeny
        )

        # integrate genotypes (gintegrate)
        founder_geno, misc = gintgop.gintegrate(
            t_cur = -burnin,
            t_max = 0,
            pgvmat = pgvmat,
            geno = founder_geno,
        )

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

        geninitop = GenerationalInitializationOperator(
            burnin = burnin,
            founder_geno = founder_geno,
            founder_bval = founder_bval,
            founder_gmod = founder_gmod,
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
    def from_vcf(fname, nfounder, founder_ncross, founder_nprogeny, gqlen, gmod_true, burnin, pselop, mateop, gintgop, evalop, bvintgop, calop, sselop, rng = None):
        """
        Create a GenerationalInitializationOperator from a VCF file.

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
        """
        # step 1: load genotype matrix
        dpgmat = DensePhasedGenotypeMatrix.from_vcf(fname)

        # step 2: create from genotype matrix
        geninitop = GenerationalInitializationOperator.from_dpgmat(
            dpgmat = dpgmat,
            nfounder = nfounder,
            founder_ncross = founder_ncross,
            founder_nprogeny = founder_nprogeny,
            gqlen = gqlen,
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
