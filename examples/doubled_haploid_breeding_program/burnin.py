#!/usr/bin/env python3

import numpy
import pybropt.core.random
from pybropt.popgen.gmat import DensePhasedGenotypeMatrix
from pybropt.popgen.gmap import ExtendedGeneticMap
from pybropt.popgen.gmap import HaldaneMapFunction
from pybropt.model.gmod import GenericLinearGenomicModel
from pybropt.breed.op.init import InitializationOperator
from pybropt.breed.op.mate import MatingOperator
from pybropt.breed.prot.mate import FamilyGroupTwoWayDHCross
from pybropt.breed.op.eval import EvaluationOperator

# seed python random and numpy random
pybropt.core.random.seed(941)

def build_gmap(path):
    """Load genetic map."""
    genetic_map = ExtendedGeneticMap.from_egmap(path)   # read from file
    genetic_map.group()                                 # group markers
    genetic_map.build_spline()                          # construct spline
    assert genetic_map.is_congruent()                   # make sure we don't have map inversions
    return genetic_map

def build_gmapfn():
    """Create recombination mapping function."""
    gmapfn = HaldaneMapFunction()                   # use generic Haldane function
    return gmapfn

def build_dpgmat(path, gmap, gmapfn):
    """Load genotype matrix."""
    gmat = DensePhasedGenotypeMatrix.from_vcf(path) # read from file
    gmat.group()                                    # group markers
    gmat.interp_xoprob(gmap, gmapfn)                # interpolate genetic positions
    return gmat

def build_gmod_true(gmat):
    nvrnt = gmat.nvrnt                          # number of markers
    beta = numpy.float64([[100.0]])             # (1,1); model intercept
    u = numpy.random.normal(0, 0.05, (nvrnt,1)) # (nvrnt,1); randomly assign marker weights
    trait = numpy.object_(["protein"])          # (1,); trait name
    model_name = "protein_test"
    glgmod = GenericLinearGenomicModel(
        beta = beta,
        u = u,
        trait = trait,
        model_name = model_name,
        params = None
    )
    return glgmod

gmap = build_gmap("McMullen_2009_US_NAM_corrected.M.egmap")
gmapfn = build_gmapfn()
dpgmat = build_dpgmat("widiv_2000SNPs_imputed_chr1-10_APGv4_noNA_noHet_q0.2_Q0.8.vcf.gz", gmap, gmapfn)
gmod_true = build_gmod_true(dpgmat)

class MyInitializationOperator(InitializationOperator):
    def __init__(self, founder_genome, founder_geno, founder_pheno, founder_bval, founder_gmod, **kwargs):
        super(MyInitializationOperator, self).__init__(**kwargs)
        self.founder_genome = founder_genome
        self.founder_geno = founder_geno
        self.founder_pheno = founder_pheno
        self.founder_bval = founder_bval
        self.founder_gmod = founder_gmod
    @classmethod
    def from_dpgmat(cls, dpgmat, mateprot, nfounder, founder_ncross, founder_nprogeny, gqlen, gmod_true, burnin, mateop):
        rng = pybropt.core.random   # get default random number generator
        ####################################################
        ### step 1: count available taxa
        ####################################################
        ntaxa = dpgmat.ntaxa
        ###################################################################
        ### step 2: define and populate dictionaries containing founder:
        ###             1) genomes
        ###             2) genotypes
        ###             3) phenotypes
        ###             4) breeding values
        ###             5) genomic models
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
        sel = rng.choice(ntaxa, nfounder, replace = False)  # get random selections (pselect)
        sel.sort()                                          # sort indices for random founder selections
        dpgmat = dpgmat.select_taxa(sel)                    # select founder individuals
        sel = numpy.arange(dpgmat.ntaxa)                    # create indices [0,1,...,nfounder-1]
        for _ in range(gqlen):                              # fill queue with random matings of founders
            rng.shuffle(sel)                                # randomly shuffle indices
            pgmat, misc = mateprot.mate(                    # mate random selections (mate)
                pgmat = dpgmat,
                sel = sel,
                ncross = founder_ncross,
                nprogeny = founder_nprogeny
            )
            founder_genome["queue"].append(pgmat)           # add progeny to queue
        rng.shuffle(sel)                                    # randomly shuffle indices
        mcfg = {                                            # random mating initially
            "pgmat" : dpgmat,
            "sel" : sel,
            "ncross" : founder_ncross,
            "nprogeny" : founder_nprogeny
        }
        genome, geno, pheno, bval, gmod, misc = mateop.mate(    # initial bootstrap mate
            mcfg = mcfg,
            genome = founder_genome,
            geno = founder_geno,
            pheno = founder_pheno,
            bval = founder_bval,
            gmod = founder_gmod,
            t_cur = -burnin,
            t_max = 0
        )
        out = cls(                                          # construct object
            founder_genome = founder_genome,
            founder_geno = founder_geno,
            founder_pheno = founder_pheno,
            founder_bval = founder_bval,
            founder_gmod = founder_gmod
        )
        return out

class MyInitMatingOperator(MatingOperator):
    def __init__(self, mprot, **kwargs):
        super(MyInitMatingOperator, self).__init__(**kwargs)
        self.mprot = mprot
    def mate(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        progeny, misc = self.mprot.mate(**mcfg, s = 0)  # mate parents
        genome["queue"].append(progeny)                 # add progeny to queue in genome dict
        misc = {}
        return genome, geno, pheno, bval, gmod, misc

class MyInitEvaluationOperator(EvaluationOperator):
    def __init__(self, ptprot, **kwargs):
        super(MyInitEvaluationOperator, self).__init__(**kwargs)
        self.ptprot = ptprot
    def evaluate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        pass
    def set_h2(self, h2, pgmat, gpmod = None, **kwargs):
        self.ptprot.set_h2(h2, pgmat, gpmod, **kwargs)

mateprot = FamilyGroupTwoWayDHCross()       # make mating protocol
imateop = MyInitMatingOperator(mateprot)    # make init mating operator

initop = MyInitializationOperator.from_dpgmat(
    dpgmat = dpgmat,
    mateprot = mateprot,
    nfounder = 40,
    founder_ncross = 1,
    founder_nprogeny = 80,
    gqlen = 6,
    gmod_true = gmod_true,
    burnin = 20,
    mateop = imateop,
)
