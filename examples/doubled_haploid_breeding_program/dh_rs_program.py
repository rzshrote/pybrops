#!/usr/bin/env python3

import numpy
import pandas
import pybrops.core.random
from pybrops.breed.arch.RecurrentSelectionBreedingProgram import RecurrentSelectionBreedingProgram
from pybrops.breed.op.init.InitializationOperator import InitializationOperator
from pybrops.breed.op.psel.ParentSelectionOperator import ParentSelectionOperator
from pybrops.breed.op.mate.MatingOperator import MatingOperator
from pybrops.breed.op.eval.EvaluationOperator import EvaluationOperator
from pybrops.breed.op.ssel.SurvivorSelectionOperator import SurvivorSelectionOperator
from pybrops.breed.prot.bv.MeanPhenotypicBreedingValue import MeanPhenotypicBreedingValue
from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.breed.prot.mate.TwoWayDHCross import TwoWayDHCross
from pybrops.breed.prot.pt.G_E_Phenotyping import G_E_Phenotyping
from pybrops.breed.prot.sel.FamilyPhenotypicSelection import FamilyPhenotypicSelection
from pybrops.breed.prot.sel.ConventionalPhenotypicSelection import ConventionalPhenotypicSelection
from pybrops.breed.prot.sel.transfn import trans_sum
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.breed.op.log.Logbook import Logbook

# seed python random and numpy random
pybrops.core.random.seed(941)

class MyInitParentSelectionOperator(ParentSelectionOperator):
    def __init__(self):
        super(MyInitParentSelectionOperator, self).__init__()
        self.pselprot = ConventionalPhenotypicSelection(
            nparent = 40,
            ncross = 1,
            nprogeny = 80,
            objfn_trans = trans_sum
        )
    def pselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        pgmat, sel, ncross, nprogeny = self.pselprot.select(
            pgmat = genome["cand"],
            gmat = geno["cand"],
            ptdf = pheno["cand"],
            bvmat = bval["cand"],
            gpmod = gmod["cand"],
            t_cur = t_cur,
            t_max = t_max,
            miscout = miscout,
            method = "single"
        )
        mcfg = {
            "pgmat": pgmat,
            "sel": sel,
            "ncross": ncross,
            "nprogeny": nprogeny
        }
        return mcfg, genome, geno, pheno, bval, gmod

class MyInitMatingOperator(MatingOperator):
    def __init__(self, pcnt, fcnt, **kwargs):
        super(MyInitMatingOperator, self).__init__(**kwargs)
        self.mprot = TwoWayDHCross(
            progeny_counter = pcnt,
            family_counter = fcnt
        )
    def mate(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout = None, **kwargs):
        progeny = self.mprot.mate(**mcfg, s = 0, miscout = miscout)  # mate parents
        genome["queue"].append(progeny)                 # add progeny to queue in genome dict
        return genome, geno, pheno, bval, gmod

class MyInitEvaluationOperator(EvaluationOperator):
    def __init__(self, gpmod, var_err, **kwargs):
        super(MyInitEvaluationOperator, self).__init__(**kwargs)
        self.gtprot = DenseUnphasedGenotyping()
        self.ptprot = G_E_Phenotyping(
            gpmod = gpmod,
            nenv = 4,
            var_err = var_err
        )
        self.bvprot = MeanPhenotypicBreedingValue("taxa", ["yield"])
    def evaluate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        geno["queue"].append(self.gtprot.genotype(genome["queue"][-1]))         # genotype incoming inbreds
        genome["queue"].pop(0)                                                  # remove oldest inbreds
        geno["queue"].pop(0)                                                    # remove oldest inbreds
        genome["main"] = genome["queue"][0].concat_taxa(genome["queue"][0:3])   # construct main population using 3 oldest cohorts
        geno["main"] = geno["queue"][0].concat_taxa(geno["queue"][0:3])         # construct main population using 3 oldest cohorts
        pheno["main"] = self.ptprot.phenotype(genome["main"])                   # phenotype main population
        bval["main"] = self.bvprot.estimate(pheno["main"], geno["main"])        # estimate breeding values after phenotyping
        return genome, geno, pheno, bval, gmod

class MyInitSurvivorSelectionOperator(SurvivorSelectionOperator):
    def __init__(self):
        super(MyInitSurvivorSelectionOperator, self).__init__()
        self.sselprot = FamilyPhenotypicSelection(
            nparent = 4,
            ncross = fps_ncross,
            nprogeny = fps_nprogeny,
            objfn_trans = trans_sum
        )
    def sselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        pgmat, sel, ncross, nprogeny = self.sselprot.select(
            pgmat = genome["main"],
            gmat = geno["main"],
            ptdf = pheno["main"],
            bvmat = bval["main"],
            gpmod = gmod["main"],
            t_cur = t_cur,
            t_max = t_max,
            miscout = miscout,
            method = "single"
        )
        sel.sort()                                          # sort selection indices
        genome["cand"] = genome["main"].select_taxa(sel)    # select parent candidates
        geno["cand"] = geno["main"].select_taxa(sel)        # select parent candidates
        bval["cand"] = bval["main"].select_taxa(sel)        # select parent candidates
        return genome, geno, pheno, bval, gmod

class MyInitializationOperator(InitializationOperator):
    def __init__(self, fndr_genome, fndr_geno, fndr_pheno, fndr_bval, fndr_gmod, pselop, mateop, evalop, sselop, burnin, **kwargs):
        super(MyInitializationOperator, self).__init__(**kwargs)
        self.genome = fndr_genome
        self.geno = fndr_geno
        self.pheno = fndr_pheno
        self.bval = fndr_bval
        self.gmod = fndr_gmod
        self.pselop = pselop
        self.mateop = mateop
        self.evalop = evalop
        self.sselop = sselop
        self.burnin = burnin
        self.t_cur = -burnin
    def initialize(self, miscout = None, verbose = True, **kwargs):
        for _ in range(self.burnin): # iterate through main breeding loop for burnin generations
            mcfg, self.genome, self.geno, self.pheno, self.bval, self.gmod = self.pselop.pselect(
                genome = self.genome,
                geno = self.geno,
                pheno = self.pheno,
                bval = self.bval,
                gmod = self.gmod,
                t_cur = self.t_cur,
                t_max = 0,
                miscout = None
            )
            self.genome, self.geno, self.pheno, self.bval, self.gmod = self.mateop.mate(
                mcfg = mcfg,
                genome = self.genome,
                geno = self.geno,
                pheno = self.pheno,
                bval = self.bval,
                gmod = self.gmod,
                t_cur = self.t_cur,
                t_max = 0,
                miscout = None
            )
            self.genome, self.geno, self.pheno, self.bval, self.gmod = self.evalop.evaluate(
                genome = self.genome,
                geno = self.geno,
                pheno = self.pheno,
                bval = self.bval,
                gmod = self.gmod,
                t_cur = self.t_cur,
                t_max = 0,
                miscout = None
            )
            self.genome, self.geno, self.pheno, self.bval, self.gmod = self.sselop.sselect(
                genome = self.genome,
                geno = self.geno,
                pheno = self.pheno,
                bval = self.bval,
                gmod = self.gmod,
                t_cur = self.t_cur,
                t_max = 0,
                miscout = None
            )
            self.t_cur += 1     # increment time variables
            if verbose:
                print("Burn-in generation {0} of {1}".format(_+1, self.burnin))
        return self.genome, self.geno, self.pheno, self.bval, self.gmod

class MyParentSelectionOperator(ParentSelectionOperator):
    def __init__(self, pselprot = None):
        super(MyParentSelectionOperator, self).__init__()
        self.pselprot = pselprot
        if self.pselprot is None:
            self.pselprot = ConventionalPhenotypicSelection(
                nparent = 40,
                ncross = 1,
                nprogeny = 80,
                objfn_trans = trans_sum
            )
    def pselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        pgmat, sel, ncross, nprogeny = self.pselprot.select(
            pgmat = genome["cand"],
            gmat = geno["cand"],
            ptdf = pheno["cand"],
            bvmat = bval["cand"],
            gpmod = gmod["cand"],
            t_cur = t_cur,
            t_max = t_max,
            miscout = miscout,
            method = "single"
        )
        mcfg = {
            "pgmat": pgmat,
            "sel": sel,
            "ncross": ncross,
            "nprogeny": nprogeny
        }
        return mcfg, genome, geno, pheno, bval, gmod

class MyMatingOperator(MatingOperator):
    def __init__(self, pcnt, fcnt, **kwargs):
        super(MyMatingOperator, self).__init__(**kwargs)
        self.mprot = TwoWayDHCross(
            progeny_counter = pcnt,
            family_counter = fcnt
        )
    def mate(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout = None, **kwargs):
        progeny = self.mprot.mate(**mcfg, s = 0, miscout = miscout)  # mate parents
        genome["queue"].append(progeny)                 # add progeny to queue in genome dict
        return genome, geno, pheno, bval, gmod

class MyEvaluationOperator(EvaluationOperator):
    def __init__(self, gpmod, var_err, **kwargs):
        super(MyEvaluationOperator, self).__init__(**kwargs)
        self.gtprot = DenseUnphasedGenotyping()
        self.ptprot = G_E_Phenotyping(
            gpmod = gpmod,
            nenv = 4,
            var_err = var_err
        )
        self.bvprot = MeanPhenotypicBreedingValue("taxa", ["yield"])
    def evaluate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        geno["queue"].append(self.gtprot.genotype(genome["queue"][-1]))         # genotype incoming inbreds
        genome["queue"].pop(0)                                                  # remove oldest inbreds
        geno["queue"].pop(0)                                                    # remove oldest inbreds
        genome["main"] = genome["queue"][0].concat_taxa(genome["queue"][0:3])   # construct main population using 3 oldest cohorts
        geno["main"] = geno["queue"][0].concat_taxa(geno["queue"][0:3])         # construct main population using 3 oldest cohorts
        pheno["main"] = self.ptprot.phenotype(genome["main"])                   # phenotype main population
        bval["main"] = self.bvprot.estimate(pheno["main"], genome["main"])      # estimate breeding values after phenotyping
        return genome, geno, pheno, bval, gmod

class MySurvivorSelectionOperator(SurvivorSelectionOperator):
    def __init__(self):
        super(MySurvivorSelectionOperator, self).__init__()
        self.sselprot = FamilyPhenotypicSelection(
            nparent = 4,
            ncross = fps_ncross,
            nprogeny = fps_nprogeny,
            objfn_trans = trans_sum
        )
    def sselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        pgmat, sel, ncross, nprogeny = self.sselprot.select(
            pgmat = genome["main"],
            gmat = geno["main"],
            ptdf = pheno["main"],
            bvmat = bval["main"],
            gpmod = gmod["main"],
            t_cur = t_cur,
            t_max = t_max,
            miscout = miscout,
            method = "single"
        )
        sel.sort()                                          # sort selection indices
        genome["cand"] = genome["main"].select_taxa(sel)    # select parent candidates
        geno["cand"] = geno["main"].select_taxa(sel)        # select parent candidates
        bval["cand"] = bval["main"].select_taxa(sel)        # select parent candidates
        return genome, geno, pheno, bval, gmod

class MyLogbook(Logbook):
    def __init__(self):
        super(MyLogbook, self).__init__()
        self.reset()
    def data():
        doc = "The data property."
        def fget(self):
            return self._data
        def fset(self, value):
            self._data = value
        def fdel(self):
            del self._data
        return locals()
    data = property(**data())
    def rep():
        doc = "The rep property."
        def fget(self):
            return self._rep
        def fset(self, value):
            self._rep = value
        def fdel(self):
            del self._rep
        return locals()
    rep = property(**rep())
    def log_initialize(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        gpmod = gmod["true"]
        cand_bval_true = gpmod.gebv(genome["cand"])
        main_bval_true = gpmod.gebv(genome["main"])
        self.data["rep"].append(self.rep)
        self.data["t_cur"].append(t_cur)
        self.data["cand_mehe"].append(genome["cand"].mehe())
        self.data["cand_mean"].append(bval["cand"].tmean(descale = True)[0])
        self.data["cand_true_mean"].append(cand_bval_true.tmean(descale = True)[0])
        self.data["cand_std"].append(bval["cand"].tstd(descale = True)[0])
        self.data["cand_true_std"].append(cand_bval_true.tstd(descale = True)[0])
        self.data["cand_true_var_A"].append(gpmod.var_A(genome["cand"])[0])
        self.data["cand_true_var_a"].append(gpmod.var_a(genome["cand"])[0])
        self.data["cand_true_bulmer"].append(gpmod.bulmer(genome["cand"])[0])
        self.data["cand_true_usl"].append(gpmod.usl(genome["cand"], descale = True)[0])
        self.data["cand_true_lsl"].append(gpmod.lsl(genome["cand"], descale = True)[0])
        self.data["main_mehe"].append(genome["main"].mehe())
        self.data["main_mean"].append(bval["main"].tmean(descale = True)[0])
        self.data["main_true_mean"].append(main_bval_true.tmean(descale = True)[0])
        self.data["main_std"].append(bval["main"].tstd(descale = True)[0])
        self.data["main_true_std"].append(main_bval_true.tstd(descale = True)[0])
        self.data["main_true_var_A"].append(gpmod.var_A(genome["main"])[0])
        self.data["main_true_var_a"].append(gpmod.var_a(genome["main"])[0])
        self.data["main_true_bulmer"].append(gpmod.bulmer(genome["main"])[0])
        self.data["main_true_usl"].append(gpmod.usl(genome["main"], descale = True)[0])
        self.data["main_true_lsl"].append(gpmod.lsl(genome["main"], descale = True)[0])
    def log_pselect(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        pass
    def log_mate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        pass
    def log_evaluate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        pass
    def log_sselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        gpmod = gmod["true"]
        cand_bval_true = gpmod.gebv(genome["cand"])
        main_bval_true = gpmod.gebv(genome["main"])
        self.data["rep"].append(self.rep)
        self.data["t_cur"].append(t_cur)
        self.data["cand_mehe"].append(genome["cand"].mehe())
        self.data["cand_mean"].append(bval["cand"].tmean(descale = True)[0])
        self.data["cand_true_mean"].append(cand_bval_true.tmean(descale = True)[0])
        self.data["cand_std"].append(bval["cand"].tstd(descale = True)[0])
        self.data["cand_true_std"].append(cand_bval_true.tstd(descale = True)[0])
        self.data["cand_true_var_A"].append(gpmod.var_A(genome["cand"])[0])
        self.data["cand_true_var_a"].append(gpmod.var_a(genome["cand"])[0])
        self.data["cand_true_bulmer"].append(gpmod.bulmer(genome["cand"])[0])
        self.data["cand_true_usl"].append(gpmod.usl(genome["cand"], descale = True)[0])
        self.data["cand_true_lsl"].append(gpmod.lsl(genome["cand"], descale = True)[0])
        self.data["main_mehe"].append(genome["main"].mehe())
        self.data["main_mean"].append(bval["main"].tmean(descale = True)[0])
        self.data["main_true_mean"].append(main_bval_true.tmean(descale = True)[0])
        self.data["main_std"].append(bval["main"].tstd(descale = True)[0])
        self.data["main_true_std"].append(main_bval_true.tstd(descale = True)[0])
        self.data["main_true_var_A"].append(gpmod.var_A(genome["main"])[0])
        self.data["main_true_var_a"].append(gpmod.var_a(genome["main"])[0])
        self.data["main_true_bulmer"].append(gpmod.bulmer(genome["main"])[0])
        self.data["main_true_usl"].append(gpmod.usl(genome["main"], descale = True)[0])
        self.data["main_true_lsl"].append(gpmod.lsl(genome["main"], descale = True)[0])
    def reset(self):
        self.data = {
            "rep": [],
            "t_cur": [],
            "cand_mehe": [],
            "cand_mean": [],
            "cand_true_mean": [],
            "cand_std": [],
            "cand_true_std": [],
            "cand_true_var_A": [],
            "cand_true_var_a": [],
            "cand_true_bulmer": [],
            "cand_true_usl": [],
            "cand_true_lsl": [],
            "main_mehe": [],
            "main_mean": [],
            "main_true_mean": [],
            "main_std": [],
            "main_true_std": [],
            "main_true_var_A": [],
            "main_true_var_a": [],
            "main_true_bulmer": [],
            "main_true_usl": [],
            "main_true_lsl": [],
        }
        self.rep = 0
    def write(self, fname):
        pandas_df = pandas.DataFrame(self.data)
        pandas_df.to_csv(fname, index = False)

################################################################################
################################################################################
################################################################################

##################### Read genetic map #####################
gmap = ExtendedGeneticMap.from_egmap(                           # read from file
    "McMullen_2009_US_NAM_corrected.M.egmap"
)
gmap.group()                                                    # group markers
gmap.build_spline()                                             # construct spline

############### Create genetic map function ################
gmapfn = HaldaneMapFunction()                                   # use generic Haldane function

#################### Load genetic data #####################
dpgmat = DensePhasedGenotypeMatrix.from_vcf(                    # read from file
    "widiv_2000SNPs_imputed_chr1-10_APGv4_noNA_noHet_q0.2_Q0.8.vcf.gz"
)
dpgmat.group_vrnt()                                             # group markers
dpgmat.interp_xoprob(gmap, gmapfn)                              # interpolate crossover probabilies

################# Construct genomic model ##################
gmod_true = DenseAdditiveLinearGenomicModel(                    # create model
    beta = numpy.float64([[10.0]]),                             # model intercepts
    u_misc = None,                                              # miscellaneous random effects
    u_a = pybrops.core.random.normal(0, 0.01, (dpgmat.nvrnt,1)),# random marker effects
    trait = numpy.object_(["yield"]),                           # trait names
    model_name = "yield_model",                                 # name of the model
    params = None
)

################### Founding parameters ####################
fndr_heritability = 0.4                                         # heritability of founder lines
burnin = 20                                                     # number of burnin generations
t_cur = -burnin                                                 # set t_cur
t_max = 0                                                       # set t_max
gqlen = 6                                                       # breeding pipeline queue length
nfounder = 40                                                   # number of random founders to select
fndr_ncross = 1                                                 # number of founder crosses
fndr_nprogeny = 80
fps_nparent = 4
fps_ncross = 1
fps_nprogeny = 80

#################### Determine founders ####################
sel = pybrops.core.random.choice(                               # get random selections (pselect)
    dpgmat.ntaxa,
    nfounder,
    replace = False
)
sel.sort()                                                      # sort indices for random founder selections
dpgmat = dpgmat.select_taxa(sel)                                # select founder individuals

################ Build founder populations #################
mateprot = TwoWayDHCross()                                      # make mating protocol
gtprot = DenseUnphasedGenotyping()                              # genotyping protocols
ptprot = G_E_Phenotyping(gmod_true, nenv = 4)                   # make phenotyping protocol
ptprot.set_h2(fndr_heritability, dpgmat)                        # set heritability
bvprot = MeanPhenotypicBreedingValue("taxa", ["yield"])         # make breeding value protocol
sselprot = FamilyPhenotypicSelection(                           # family-based survivor selection
    nparent = fps_nparent,
    ncross = fps_ncross,
    nprogeny = fps_nprogeny,
    objfn_trans = trans_sum
)

################ Build founder populations #################
fndr_genome = {"cand":None,      "main":None,      "queue":[]}
fndr_geno =   {"cand":None,      "main":None,      "queue":[]}
fndr_pheno =  {"cand":None,      "main":None}
fndr_bval =   {"cand":None,      "main":None}
fndr_gmod =   {"cand":gmod_true, "main":gmod_true, "true":gmod_true}

################ Build founder populations #################
sel = numpy.arange(dpgmat.ntaxa)                                                # create indices [0,1,...,nfounder-1]
for _ in range(gqlen):                                                          # fill queue with random matings of founders
    pybrops.core.random.shuffle(sel)                                            # randomly shuffle indices
    pgmat = mateprot.mate(dpgmat, sel, fndr_ncross, fndr_nprogeny)              # mate random selections (mate)
    fndr_genome["queue"].append(pgmat)                                          # add progeny to queue
    gmat = gtprot.genotype(pgmat = pgmat)                                       # genotype progeny
    fndr_geno["queue"].append(gmat)                                             # add progeny genotypes to queue

fndr_genome["main"] = fndr_genome["queue"][0].concat_taxa(fndr_genome["queue"][0:3])    # construct main population
fndr_geno["main"] = fndr_geno["queue"][0].concat_taxa(fndr_geno["queue"][0:3])  # construct main population genotypes

################ Phenotype founder populations #################
fndr_pheno["main"] = ptprot.phenotype(fndr_genome["main"])                      # phenotype main population

################ Calculate breeding values for founder populations #################
fndr_bval["main"] = bvprot.estimate(fndr_pheno["main"], fndr_geno["main"])      # estimate breeding values

#### Develop first set of parental candidates ####
pgmat, sel, ncross, nprogeny = sselprot.select(
    pgmat = fndr_genome["main"],
    gmat = fndr_geno["main"],
    ptdf = fndr_pheno["main"],
    bvmat = fndr_bval["main"],
    gpmod = fndr_gmod["main"],
    t_cur = t_cur,
    t_max = t_max,
    miscout = None,
    method = "single"
)

### select parental candidates ###
sel.sort()
fndr_genome["cand"] = fndr_genome["main"].select_taxa(sel)
fndr_geno["cand"] = fndr_geno["main"].select_taxa(sel)
fndr_bval["cand"] = fndr_bval["main"].select_taxa(sel)


################################################################################
################################################################################
################################################################################

### create init operators ###
init_pselop = MyInitParentSelectionOperator()
init_mateop = MyInitMatingOperator(mateprot.progeny_counter, mateprot.family_counter)
init_evalop = MyInitEvaluationOperator(gmod_true, ptprot.var_err)
init_sselop = MyInitSurvivorSelectionOperator()

### create main operators ###
initop = MyInitializationOperator(
    fndr_genome = fndr_genome,
    fndr_geno = fndr_geno,
    fndr_pheno = fndr_pheno,
    fndr_bval = fndr_bval,
    fndr_gmod = fndr_gmod,
    pselop = init_pselop,
    mateop = init_mateop,
    evalop = init_evalop,
    sselop = init_sselop,
    burnin = burnin
)
pselop = MyParentSelectionOperator()
mateop = MyMatingOperator(mateprot.progeny_counter, mateprot.family_counter)
evalop = MyEvaluationOperator(gmod_true, ptprot.var_err)
sselop = MySurvivorSelectionOperator()
lbook = MyLogbook()

rsprog = RecurrentSelectionBreedingProgram(
    initop = initop,
    pselop = pselop,
    mateop = mateop,
    evalop = evalop,
    sselop = sselop,
    t_max = 20
)

# evolve the population
rsprog.evolve(nrep = 4, ngen = 20, lbook = lbook, verbose = True)

lbook.write("dh_rs_program.csv")
