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
    def __init__(self, ptprot, **kwargs):
        super(MyInitEvaluationOperator, self).__init__(**kwargs)
        self.ptprot = ptprot
    def evaluate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        pass
    def set_h2(self, ):
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

        # phenotype main population (evaluate)
        bvmat, bvmat_true, misc = evalop.evaluate(
            t_cur = -burnin,
            t_max = 0,
            pgvmat = founder_geno["main"],
            gmod_true = founder_gmod["true"]
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
