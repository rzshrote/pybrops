#!/usr/bin/env python3

import numpy
import pybropt.core.random
from pybropt.breed.op.init import InitializationOperator
from pybropt.breed.op.mate import MatingOperator
from pybropt.breed.op.eval import EvaluationOperator
from pybropt.breed.prot.bv import MeanPhenotypicBreedingValue
from pybropt.breed.prot.gt import DenseUnphasedGenotyping
from pybropt.breed.prot.mate import TwoWayDHCross
from pybropt.breed.prot.pt import G_E_Phenotyping
from pybropt.breed.prot.sel import FamilyPhenotypicSelection
from pybropt.breed.prot.sel.transfn import trans_sum
from pybropt.model.gmod import GenericLinearGenomicModel
from pybropt.popgen.gmat import DensePhasedGenotypeMatrix
from pybropt.popgen.gmap import ExtendedGeneticMap
from pybropt.popgen.gmap import HaldaneMapFunction

# seed python random and numpy random
pybropt.core.random.seed(941)

class MyInitializationOperator(InitializationOperator):
    def __init__(self, founder_genome, founder_geno, founder_pheno, founder_bval, founder_gmod, mateop, evalop, **kwargs):
        super(MyInitializationOperator, self).__init__(**kwargs)
        self.founder_genome = founder_genome
        self.founder_geno = founder_geno
        self.founder_pheno = founder_pheno
        self.founder_bval = founder_bval
        self.founder_gmod = founder_gmod
        self.mateop = mateop
        self.evalop = evalop
    @classmethod
    def from_dpgmat(cls, dpgmat, mateprot, nfounder, founder_ncross, founder_nprogeny, gqlen, gmod_true, burnin, mateop, evalop):
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
            pgmat = mateprot.mate(                    # mate random selections (mate)
                pgmat = dpgmat,
                sel = sel,
                ncross = founder_ncross,
                nprogeny = founder_nprogeny
            )
            founder_genome["queue"].append(pgmat)           # add progeny to queue
        evalop.set_h2(0.4, dpgmat)                          # set heritability
        rng.shuffle(sel)                                    # randomly shuffle indices
        mcfg = {                                            # random mating initially
            "pgmat" : dpgmat,
            "sel" : sel,
            "ncross" : founder_ncross,
            "nprogeny" : founder_nprogeny
        }
        t_cur = -burnin                                     # set t_cur
        t_max = 0                                           # set t_max
        # initial bootstrap mate
        genome, geno, pheno, bval, gmod = mateop.mate(
            mcfg = mcfg,
            genome = founder_genome,
            geno = founder_geno,
            pheno = founder_pheno,
            bval = founder_bval,
            gmod = founder_gmod,
            t_cur = t_cur,
            t_max = t_max
        )
        # initial bootstrap evaluate
        genome, geno, pheno, bval, gmod = evalop.evaluate(
            genome = genome,
            geno = geno,
            pheno = pheno,
            bval = bval,
            gmod = gmod,
            t_cur = t_cur,
            t_max = t_max
        )
        # construct object
        out = cls(
            founder_genome = founder_genome,
            founder_geno = founder_geno,
            founder_pheno = founder_pheno,
            founder_bval = founder_bval,
            founder_gmod = founder_gmod,
            mateop = mateop,
            evalop = evalop
        )
        return out

class MyInitMatingOperator(MatingOperator):
    def __init__(self, mprot, **kwargs):
        super(MyInitMatingOperator, self).__init__(**kwargs)
        self.mprot = mprot
    def mate(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout = None, **kwargs):
        progeny = self.mprot.mate(**mcfg, s = 0)  # mate parents
        genome["queue"].append(progeny)                 # add progeny to queue in genome dict
        return genome, geno, pheno, bval, gmod

class MyInitEvaluationOperator(EvaluationOperator):
    def __init__(self, ptprot, **kwargs):
        super(MyInitEvaluationOperator, self).__init__(**kwargs)
        self.ptprot = ptprot
    def evaluate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        return genome, geno, pheno, bval, gmod, {}
    def set_h2(self, h2, pgmat, gpmod = None, **kwargs):
        self.ptprot.set_h2(h2, pgmat, gpmod, **kwargs)

################################################################################

##################### Read genetic map #####################
gmap = ExtendedGeneticMap.from_egmap("McMullen_2009_US_NAM_corrected.M.egmap")  # read from file
gmap.group()                                                                    # group markers
gmap.build_spline()                                                             # construct spline

############### Create genetic map function ################
gmapfn = HaldaneMapFunction()                                                   # use generic Haldane function

#################### Load genetic data #####################
dpgmat = DensePhasedGenotypeMatrix.from_vcf("widiv_2000SNPs_imputed_chr1-10_APGv4_noNA_noHet_q0.2_Q0.8.vcf.gz") # read from file
dpgmat.group()                                                                  # group markers
dpgmat.interp_xoprob(gmap, gmapfn)                                              # interpolate crossover probabilies

################# Construct genomic model ##################
gmod_true = GenericLinearGenomicModel(                                          # create model
    beta = numpy.float64([[10.0]]),                                            # model intercepts
    u = pybropt.core.random.normal(0, 0.01, (dpgmat.nvrnt,1)),                  # random marker weights
    trait = numpy.object_(["yield"]),                                           # trait names
    model_name = "yield_model",
    params = None
)

################################################################################

################### Founding parameters ####################
founder_heritability = 0.4                                                      # heritability of founder lines
burnin = 20
t_cur = -burnin                                                                 # set t_cur
t_max = 0                                                                       # set t_max
gqlen = 6
nfounder = 40
founder_ncross = 1
founder_nprogeny = 80
fps_nparent = 4
fps_ncross = 1
fps_nprogeny = 80

#################### Determine founders ####################
sel = pybropt.core.random.choice(dpgmat.ntaxa, nfounder, replace = False)       # get random selections (pselect)
sel.sort()                                                                      # sort indices for random founder selections
dpgmat = dpgmat.select_taxa(sel)                                                # select founder individuals

################ Build founder populations #################
mateprot = TwoWayDHCross()                                                      # make mating protocol
gtprot = DenseUnphasedGenotyping()                                              # genotyping protocols
ptprot = G_E_Phenotyping(gmod_true, nenv = 4)                                   # make phenotyping protocol
ptprot.set_h2(founder_heritability, dpgmat)                                     # set heritability
bvprot = MeanPhenotypicBreedingValue("taxa", ["yield"])                         # make breeding value protocol
sselprot = FamilyPhenotypicSelection(
    nparent = fps_nparent,
    ncross = fps_ncross,
    nprogeny = fps_nprogeny,
    objfn_trans = trans_sum
)

################ Build founder populations #################
founder_genome = {"cand":None,      "main":None,      "queue":[]}
founder_geno =   {"cand":None,      "main":None,      "queue":[]}
founder_pheno =  {"cand":None,      "main":None}
founder_bval =   {"cand":None,      "main":None}
founder_gmod =   {"cand":gmod_true, "main":gmod_true, "true":gmod_true}

################ Build founder populations #################
sel = numpy.arange(dpgmat.ntaxa)                                                # create indices [0,1,...,nfounder-1]
for _ in range(gqlen):                                                          # fill queue with random matings of founders
    pybropt.core.random.shuffle(sel)                                            # randomly shuffle indices
    pgmat = mateprot.mate(                                                      # mate random selections (mate)
        pgmat = dpgmat,
        sel = sel,
        ncross = founder_ncross,
        nprogeny = founder_nprogeny
    )
    founder_genome["queue"].append(pgmat)                                       # add progeny to queue
    gmat = gtprot.genotype(                                                     # genotype progeny
        pgmat = pgmat
    )
    founder_geno["queue"].append(gmat)                                          # add progeny genotypes to queue

founder_genome["main"] = founder_genome["queue"][0].concat_taxa(                # construct main population
    founder_genome["queue"][0:3]
)
founder_geno["main"] = founder_geno["queue"][0].concat_taxa(                    # construct main population genotypes
    founder_geno["queue"][0:3]
)

################ Phenotype founder populations #################
founder_pheno["main"] = ptprot.phenotype(                                       # phenotype main population
    founder_genome["main"]
)

################ Calculate breeding values for founder populations #################
founder_bval["main"] = bvprot.estimate(                                         # estimate breeding values
    founder_pheno["main"],
    founder_geno["main"]
)

#### Develop first set of parental candidates ####
pgmat, sel, ncross, nprogeny = sselprot.select(
    pgmat = founder_genome["main"],
    gmat = founder_geno["main"],
    ptdf = founder_pheno["main"],
    bvmat = founder_bval["main"],
    gpmod = founder_gmod["main"],
    t_cur = t_cur,
    t_max = t_max,
    miscout = None,
    method = "single"
)

### select parental candidates ###
sel.sort()
founder_genome["cand"] = founder_genome["main"].select_taxa(sel)
founder_geno["cand"] = founder_geno["main"].select_taxa(sel)
founder_bval["cand"] = founder_bval["main"].select_taxa(sel)
# print(founder_genome["cand"].taxa)
# print(founder_bval["cand"].taxa)
# print(founder_genome["cand"].taxa == founder_bval["cand"].taxa)
