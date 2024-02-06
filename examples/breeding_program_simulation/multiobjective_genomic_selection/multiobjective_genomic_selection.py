#!/usr/bin/env python3

###
### Multi-objective Genomic Selection
### #################################

##
## Simulation Preliminaries
## ========================

#
# Loading Required Modules and Seeding the global PRNG
# ----------------------------------------------------

import numpy
import pandas
import pybrops
from pybrops.breed.arch.RecurrentSelectionBreedingProgram import RecurrentSelectionBreedingProgram
from pybrops.breed.op.eval.EvaluationOperator import EvaluationOperator
from pybrops.breed.op.init.InitializationOperator import InitializationOperator
from pybrops.breed.op.log.Logbook import Logbook
from pybrops.breed.op.mate.MatingOperator import MatingOperator
from pybrops.breed.op.psel.ParentSelectionOperator import ParentSelectionOperator
from pybrops.breed.op.ssel.SurvivorSelectionOperator import SurvivorSelectionOperator
from pybrops.breed.prot.bv.MeanPhenotypicBreedingValue import MeanPhenotypicBreedingValue
from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol
from pybrops.breed.prot.mate.TwoWayDHCross import TwoWayDHCross
from pybrops.breed.prot.pt.G_E_Phenotyping import G_E_Phenotyping
from pybrops.breed.prot.sel.EstimatedBreedingValueSelection import EstimatedBreedingValueSubsetSelection
from pybrops.breed.prot.sel.GenomicEstimatedBreedingValueSelection import GenomicEstimatedBreedingValueSubsetSelection
from pybrops.breed.prot.sel.transfn import trans_ndpt_to_vec_dist
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.gmap.StandardGeneticMap import StandardGeneticMap
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

# seed python ``random`` and ``numpy.random`` simultaneously for reproducability
pybrops.core.random.prng.seed(17320132226)

#
# Simulation Parameters
# ---------------------

nfndr = 40          # number of random founders to select
nqtl = 1000         # number of QTL
qlen = 6            # breeding pipeline queue length
ncross = 20         # number of cross configurations (20, 2-way crosses)
nparent = 2         # number of parents per cross configuration (2-way crosses)
nmating = 1         # number of times to perform cross configuration
nprogeny = 80       # number of progenies per cross attempt
nrandmate = 20      # number of random intermatings
nburnin = 10        # number of burnin generations
ngen = 10           # number of simulation generations
nrep = 1            # number of simulation replications
fndr_heritability = numpy.array([0.4, 0.6]) # heritability of founder lines
t_cur = -nburnin    # current time for burn-in
t_max = 0           # maximum time for burn-in
nindiv_fam = 4      # number of individuals to select per family for within-family selection
nmating_fam = 1     # number of times to perform cross configuration for within-family selection 
nprogeny_fam = 80   # number of progenies to generate for within-family selection

##
## Loading Genetic Map Data from a Text File
## =========================================

# read genetic map from a CSV-like file
gmap = StandardGeneticMap.from_csv(
    "McMullen_2009_US_NAM.gmap",
    vrnt_chrgrp_col     = "chr", # name of marker chromosome number column
    vrnt_phypos_col     = "pos", # name of marker physical position column
    vrnt_genpos_col     = "cM",  # name of marker genetic position column
    vrnt_genpos_units   = "cM",  # units of the genetic position column (cM, M)
    auto_group          = True,  # whether to auto group chromosomes on load
    auto_build_spline   = True,  # whether to auto build interpolation spline on load
    sep                 = "\t",  # separator for file format
    header              = 0,     # index of the header row
)

##
## Creating a Genetic Map Function
## ===============================

# use Haldane map function to calculate crossover probabilities
gmapfn = HaldaneMapFunction()

##
## Loading Genome Data from a VCF File
## ===================================

# read phased genetic markers from a vcf file
fndr_pgmat = DensePhasedGenotypeMatrix.from_vcf(
    "widiv_2000SNPs_q0.2_Q0.8.vcf.gz", # file name to load
    auto_group_vrnt = True,            # automatically sort and group variants
)

# interpolate genetic map positions
fndr_pgmat.interp_xoprob(gmap, gmapfn)

##
## Constructing a Bi-Trait Genomic Model
## =====================================

# create trait means (intercepts) for model
beta = numpy.array([[10.0, 25.0]], dtype = float)

# create marker effects from MVN distribution for model:
# 1) marker effects have mean zero
# 2) marker effects have a covariance structure with negative marker covariance
#    meaning that markers are pleiotropic and competing in nature.
# marker effect matrix is shape (nvrnt,2)
mkreffect = numpy.random.multivariate_normal(
    mean = numpy.array([0.0, 0.0]),
    cov = numpy.array([
        [ 1.0, -0.4],
        [-0.4,  1.0]
    ]),
    size = fndr_pgmat.nvrnt
)

# create trait names for model:
trait = numpy.array(["syn1","syn2"], dtype = object)

# create the true genomic model
algmod_true = DenseAdditiveLinearGenomicModel(
    beta   = beta,      # model intercepts
    u_misc = None,      # miscellaneous random effects
    u_a    = mkreffect, # random marker effects
    trait  = trait,     # trait names
)

##
## Determine the Founder Population
## --------------------------------

# get indices of random selections
# randomly choose ``nfounder`` indices from ``ntaxa``
sel = numpy.random.choice(fndr_pgmat.ntaxa, nfndr, replace=False)

# sort indices for random founder selections
sel.sort()

# select founder individuals
fndr_pgmat = fndr_pgmat.select_taxa(sel)

#
# Create Mating Protocols for Burn-In
# -----------------------------------

# create a 2-way DH cross object
mate2waydh = TwoWayDHCross()

#
# Create Genotyping Protocols for Burn-In
# ---------------------------------------

# create a genotyping protocol
gtprot = DenseUnphasedGenotyping()

#
# Create Phenotyping Protocols for Burn-In
# ----------------------------------------

# create a phenotyping protocol using the true genomic model
ptprot = G_E_Phenotyping(algmod_true, 4, 1)

# set the heritability using the founder population
ptprot.set_h2(0.4, fndr_pgmat)

#
# Create Breeing Value Estimation Protocols for Burn-In
# -----------------------------------------------------

# create a breeding value estimation protocol
bvprot = MeanPhenotypicBreedingValue("taxa", "taxa_grp", trait)

#
# Create a Within-Family Selection Helper Function
# ------------------------------------------------

# define function to do within family selection based on yield
def within_family_selection(bvmat: DenseBreedingValueMatrix, nindiv: int) -> numpy.ndarray:
    """
    Select individuals within a family based on the sum of their breeding values.

    Parameters
    ----------
    bvmat : DenseBreedingValueMatrix
        Input breeding value matrix from which to calculate values.
    nindiv : int
        Number of individuals to select from each family.
    
    Returns
    -------
    indices : numpy.ndarray
        An array of array indices indicating which individuals are to be selected.
    """
    order = numpy.arange(bvmat.ntaxa)
    value = bvmat.mat.sum(1)
    indices = []
    groups = numpy.unique(bvmat.taxa_grp)
    for group in groups:
        mask = bvmat.taxa_grp == group
        tmp_order = order[mask]
        tmp_value = value[mask]
        value_argsort = tmp_value.argsort()
        ix = value_argsort[::-1][:nindiv]
        indices.append(tmp_order[ix])
    indices = numpy.concatenate(indices)
    return indices

#
# Create a Cohort Construction Helper function
# --------------------------------------------

# define a helper function to help make cohorts of individuals
def cohort(
        mateprot: MatingProtocol, 
        fndr_pgmat: DensePhasedGenotypeMatrix, 
        ncross: int, 
        nparent: int,
        nmating: int, 
        nprogeny: int
    ) -> DensePhasedGenotypeMatrix:
    """
    Randomly sample individuals from a founder population
    """
    # sample indicies of individuals and reshape for input into mating protocol
    xconfix = numpy.random.choice(fndr_pgmat.ntaxa, ncross * nparent, replace = False)
    xconfig = xconfix.reshape(ncross, nparent)
    # mate individuals
    out = mateprot.mate(fndr_pgmat, xconfig, nmating, nprogeny)
    return out

#
# Create Cohort Structure
# -----------------------

################ Build founder populations #################
fndr_genome = {"cand":None,      "main":None,      "queue":[]}
fndr_geno =   {"cand":None,      "main":None,      "queue":[]}
fndr_pheno =  {"cand":None,      "main":None}
fndr_bval =   {"cand":None,      "main":None}
fndr_gmod =   {"cand":algmod_true, "main":algmod_true, "true":algmod_true}

# fill queue with cohort genomes derived from randomly mating the founders
fndr_genome["queue"] = [cohort(mate2waydh, fndr_pgmat, ncross, nparent, nmating, nprogeny) for _ in range(qlen)]

# construct the main population genomes from the first three cohorts in the queue
fndr_genome["main"] = DensePhasedGenotypeMatrix.concat_taxa(fndr_genome["queue"][0:3])

# genotype individuals to fill the genotyping queue
fndr_geno["queue"] = [gtprot.genotype(genome) for genome in fndr_genome["queue"]]

# construct the main population genotypes from the first three cohorts in the queue
fndr_geno["main"] = DenseGenotypeMatrix.concat_taxa(fndr_geno["queue"][0:3])

# phenotype the main population
fndr_pheno["main"] = ptprot.phenotype(fndr_genome["main"])

# calculate breeding values for the main population
fndr_bval["main"] = bvprot.estimate(fndr_pheno["main"], fndr_geno["main"])

# # calculate indices for within family selection to get parental candidates
ix = within_family_selection(fndr_bval["main"], nindiv_fam) # select top 5%

# # select parental candidates
fndr_genome["cand"] = fndr_genome["main"].select_taxa(ix)
fndr_geno["cand"]   = fndr_geno["main"].select_taxa(ix)
fndr_bval["cand"]   = fndr_bval["main"].select_taxa(ix) # breeding values have been recentered and rescaled

###
### Define breeding program operators
### =================================

class MyInitParentSelectionOperator(ParentSelectionOperator):
    """
    Custom Parent Selection Operator class for our custom Initialization Operator
    """
    def __init__(self):
        # create parental selection protocol that selects individuals based on
        # their estimated breeding values; this is phenotypic selection
        self.pselprot = EstimatedBreedingValueSubsetSelection(
            ntrait = 2,
            unscale = True,
            ncross = 20,
            nparent = 2,
            nmating = 1,
            nprogeny = 80,
            nobj = 2,
            ndset_wt = 1.0,
            ndset_trans = trans_ndpt_to_vec_dist,
            ndset_trans_kwargs = {
                "objfn_wt": numpy.array([1.0, 1.0]),    # all objectives maximizing
                "wt": numpy.array([0.5, 0.5])           # 1/2; equal weight to all
            },
        )
    def pselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        mcfg = {}
        mcfg["cand"] = self.pselprot.select(
            pgmat = genome["cand"],
            gmat = geno["cand"],
            ptdf = pheno["cand"],
            bvmat = bval["cand"],
            gpmod = gmod["cand"],
            t_cur = t_cur,
            t_max = t_max,
            miscout = miscout
        )
        return mcfg, genome, geno, pheno, bval, gmod

class MyInitMatingOperator(MatingOperator):
    def __init__(self, pcnt, fcnt, **kwargs):
        super(MyInitMatingOperator, self).__init__(**kwargs)
        self.mprot = TwoWayDHCross(
            progeny_counter = pcnt,
            family_counter = fcnt
        )
    def mate(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout = None, **kwargs):
        progeny = self.mprot.mate(
            pgmat = mcfg["cand"].pgmat,
            xconfig = mcfg["cand"].xconfig,
            nmating = mcfg["cand"].nmating,
            nprogeny = mcfg["cand"].nprogeny,
            miscout = miscout,
            nself = 0,
        )
        genome["queue"].append(progeny)                 # add progeny to queue in genome dict
        return genome, geno, pheno, bval, gmod

class MyInitEvaluationOperator(EvaluationOperator):
    def __init__(self, gpmod, var_err, **kwargs):
        self.gtprot = DenseUnphasedGenotyping()
        self.ptprot = G_E_Phenotyping(gpmod = gpmod, nenv = 4, var_err = var_err)
        self.bvprot = MeanPhenotypicBreedingValue("taxa", "taxa_grp", trait)
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
    def __init__(self, nindiv_fam):
        self.nindiv_fam = nindiv_fam
    def sselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        # calculate indices for within family selection to get parental candidates
        ix = within_family_selection(bval["main"], self.nindiv_fam) # select top 5%
        # select parental candidates
        genome["cand"] = genome["main"].select_taxa(ix)
        geno["cand"]   = geno["main"].select_taxa(ix)
        bval["cand"]   = bval["main"].select_taxa(ix) # breeding values have been recentered and rescaled
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
    def __init__(self):
        self.pselprot = GenomicEstimatedBreedingValueSubsetSelection(
            ntrait = 2,
            unscale = True,
            ncross = 20,
            nparent = 2,
            nmating = 1,
            nprogeny = 80,
            nobj = 2,
            ndset_wt = 1.0,
            ndset_trans = trans_ndpt_to_vec_dist,
            ndset_trans_kwargs = {
                "objfn_wt": numpy.array([1.0, 1.0]),    # all objectives maximizing
                "wt": numpy.array([0.5, 0.5])           # 1/2; equal weight to all
            },
        )
    def pselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        mcfg = {}
        mcfg["cand"] = self.pselprot.select(
            pgmat = genome["cand"],
            gmat = geno["cand"],
            ptdf = pheno["cand"],
            bvmat = bval["cand"],
            gpmod = gmod["cand"],
            t_cur = t_cur,
            t_max = t_max,
            miscout = miscout
        )
        return mcfg, genome, geno, pheno, bval, gmod

class MyMatingOperator(MatingOperator):
    def __init__(self, pcnt, fcnt, **kwargs):
        self.mprot = TwoWayDHCross(
            progeny_counter = pcnt,
            family_counter = fcnt
        )
    def mate(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout = None, **kwargs):
        # mate parents
        progeny = self.mprot.mate(
            pgmat = mcfg["cand"].pgmat,
            xconfig = mcfg["cand"].xconfig,
            nmating = mcfg["cand"].nmating,
            nprogeny = mcfg["cand"].nprogeny,
            miscout = miscout,
            nself = 0,
        )
        genome["queue"].append(progeny)                 # add progeny to queue in genome dict
        return genome, geno, pheno, bval, gmod

class MyEvaluationOperator(EvaluationOperator):
    def __init__(self, gpmod, var_err, **kwargs):
        self.gtprot = DenseUnphasedGenotyping()
        self.ptprot = G_E_Phenotyping(gpmod = gpmod, nenv = 4, var_err = var_err)
        self.bvprot = MeanPhenotypicBreedingValue("taxa", "taxa_grp", trait)
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
    def __init__(self, nindiv_fam):
        self.nindiv_fam = nindiv_fam
    def sselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        # calculate indices for within family selection to get parental candidates
        ix = within_family_selection(bval["main"], self.nindiv_fam) # select top 5%
        # select parental candidates
        genome["cand"] = genome["main"].select_taxa(ix)
        geno["cand"]   = geno["main"].select_taxa(ix)
        bval["cand"]   = bval["main"].select_taxa(ix) # breeding values have been recentered and rescaled
        return genome, geno, pheno, bval, gmod

class MyLogbook(Logbook):
    def __init__(self):
        super(MyLogbook, self).__init__()
        self.reset()
    @property
    def data(self) -> dict:
        """Logbook data."""
        return self._data
    @data.setter
    def data(self, value: dict) -> None:
        """Set logbook data."""
        self._data = value
    @property
    def data_frontier(self) -> dict:
        """Logbook frontier data."""
        return self._data_frontier
    @data_frontier.setter
    def data_frontier(self, value: dict) -> None:
        """Set logbook frontier data."""
        self._data_frontier = value
    @property
    def rep(self) -> int:
        """Replication number."""
        return self._rep
    @rep.setter
    def rep(self, value: int) -> None:
        """Set replication number."""
        self._rep = value
    def log_initialize(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        gpmod = gmod["true"]
        cand_bval_true = gpmod.gebv(genome["cand"])
        main_bval_true = gpmod.gebv(genome["main"])
        self.data["rep"].append(self.rep)
        self.data["t_cur"].append(t_cur)
        ################ candidate mean expected heterozygosity ################
        self.data["cand_meh"].append(genome["cand"].meh())
        ########################### candidate means ############################
        tmp = bval["cand"].tmean(unscale = True)
        self.data["cand_mean_syn1"].append(tmp[0])
        self.data["cand_mean_syn2"].append(tmp[1])
        ######################### candidate true means #########################
        tmp = cand_bval_true.tmean(unscale = True)
        self.data["cand_true_mean_syn1"].append(tmp[0])
        self.data["cand_true_mean_syn2"].append(tmp[1])
        #################### candidate standard deviations #####################
        tmp = bval["cand"].tstd(unscale = True)
        self.data["cand_std_syn1"].append(tmp[0])
        self.data["cand_std_syn2"].append(tmp[1])
        ################## candidate true standard deviations ##################
        tmp = cand_bval_true.tstd(unscale = True)
        self.data["cand_true_std_syn1"].append(tmp[0])
        self.data["cand_true_std_syn2"].append(tmp[1])
        ############### candidate true additive genetic variance ###############
        tmp = gpmod.var_A(genome["cand"])
        self.data["cand_true_var_A_syn1"].append(tmp[0])
        self.data["cand_true_var_A_syn2"].append(tmp[1])
        ################ candidate true additive genic variance ################
        tmp = gpmod.var_a(genome["cand"])
        self.data["cand_true_var_a_syn1"].append(tmp[0])
        self.data["cand_true_var_a_syn2"].append(tmp[1])
        ##################### candidate true bulmer ratio ######################
        tmp = gpmod.bulmer(genome["cand"])
        self.data["cand_true_bulmer_syn1"].append(tmp[0])
        self.data["cand_true_bulmer_syn2"].append(tmp[1])
        ################# candidate true upper selection limit #################
        tmp = gpmod.usl(genome["cand"], unscale = True)
        self.data["cand_true_usl_syn1"].append(tmp[0])
        self.data["cand_true_usl_syn2"].append(tmp[1])
        ################# candidate true lower selection limit #################
        tmp = gpmod.lsl(genome["cand"], unscale = True)
        self.data["cand_true_lsl_syn1"].append(tmp[0])
        self.data["cand_true_lsl_syn2"].append(tmp[1])
        ########################################################################
        ################## main mean expected heterozygosity ###################
        self.data["main_meh"].append(genome["main"].meh())
        ############################## main means ##############################
        tmp = bval["main"].tmean(unscale = True)
        self.data["main_mean_syn1"].append(tmp[0])
        self.data["main_mean_syn2"].append(tmp[1])
        ########################### main true means ############################
        tmp = main_bval_true.tmean(unscale = True)
        self.data["main_true_mean_syn1"].append(tmp[0])
        self.data["main_true_mean_syn2"].append(tmp[1])
        ####################### main standard deviations #######################
        tmp = bval["main"].tstd(unscale = True)
        self.data["main_std_syn1"].append(tmp[0])
        self.data["main_std_syn2"].append(tmp[1])
        #################### main true standard deviations #####################
        tmp = main_bval_true.tstd(unscale = True)
        self.data["main_true_std_syn1"].append(tmp[0])
        self.data["main_true_std_syn2"].append(tmp[1])
        ##################### main true genetic variances ######################
        tmp = gpmod.var_A(genome["main"])
        self.data["main_true_var_A_syn1"].append(tmp[0])
        self.data["main_true_var_A_syn2"].append(tmp[1])
        ###################### main true genic variances #######################
        tmp = gpmod.var_a(genome["main"])
        self.data["main_true_var_a_syn1"].append(tmp[0])
        self.data["main_true_var_a_syn2"].append(tmp[1])
        ####################### main true bulmer ratios ########################
        tmp = gpmod.bulmer(genome["main"])
        self.data["main_true_bulmer_syn1"].append(tmp[0])
        self.data["main_true_bulmer_syn2"].append(tmp[1])
        ################### main true lower selection limits ###################
        tmp = gpmod.usl(genome["main"], unscale = True)
        self.data["main_true_usl_syn1"].append(tmp[0])
        self.data["main_true_usl_syn2"].append(tmp[1])
        ################### main true lower selection limits ###################
        tmp = gpmod.lsl(genome["main"], unscale = True)
        self.data["main_true_lsl_syn1"].append(tmp[0])
        self.data["main_true_lsl_syn2"].append(tmp[1])
    def log_pselect(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        self.data_frontier["rep"].append(self.rep)
        self.data_frontier["t_cur"].append(t_cur)
        if "mosoln" in kwargs:
            self.data_frontier["frontier"].append(kwargs["mosoln"])
        else:
            self.data_frontier["frontier"].append(None)
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
        ################ candidate mean expected heterozygosity ################
        self.data["cand_meh"].append(genome["cand"].meh())
        ########################### candidate means ############################
        tmp = bval["cand"].tmean(unscale = True)
        self.data["cand_mean_syn1"].append(tmp[0])
        self.data["cand_mean_syn2"].append(tmp[1])
        ######################### candidate true means #########################
        tmp = cand_bval_true.tmean(unscale = True)
        self.data["cand_true_mean_syn1"].append(tmp[0])
        self.data["cand_true_mean_syn2"].append(tmp[1])
        #################### candidate standard deviations #####################
        tmp = bval["cand"].tstd(unscale = True)
        self.data["cand_std_syn1"].append(tmp[0])
        self.data["cand_std_syn2"].append(tmp[1])
        ################## candidate true standard deviations ##################
        tmp = cand_bval_true.tstd(unscale = True)
        self.data["cand_true_std_syn1"].append(tmp[0])
        self.data["cand_true_std_syn2"].append(tmp[1])
        ############### candidate true additive genetic variance ###############
        tmp = gpmod.var_A(genome["cand"])
        self.data["cand_true_var_A_syn1"].append(tmp[0])
        self.data["cand_true_var_A_syn2"].append(tmp[1])
        ################ candidate true additive genic variance ################
        tmp = gpmod.var_a(genome["cand"])
        self.data["cand_true_var_a_syn1"].append(tmp[0])
        self.data["cand_true_var_a_syn2"].append(tmp[1])
        ##################### candidate true bulmer ratio ######################
        tmp = gpmod.bulmer(genome["cand"])
        self.data["cand_true_bulmer_syn1"].append(tmp[0])
        self.data["cand_true_bulmer_syn2"].append(tmp[1])
        ################# candidate true upper selection limit #################
        tmp = gpmod.usl(genome["cand"], unscale = True)
        self.data["cand_true_usl_syn1"].append(tmp[0])
        self.data["cand_true_usl_syn2"].append(tmp[1])
        ################# candidate true lower selection limit #################
        tmp = gpmod.lsl(genome["cand"], unscale = True)
        self.data["cand_true_lsl_syn1"].append(tmp[0])
        self.data["cand_true_lsl_syn2"].append(tmp[1])
        ########################################################################
        ################## main mean expected heterozygosity ###################
        self.data["main_meh"].append(genome["main"].meh())
        ############################## main means ##############################
        tmp = bval["main"].tmean(unscale = True)
        self.data["main_mean_syn1"].append(tmp[0])
        self.data["main_mean_syn2"].append(tmp[1])
        ########################### main true means ############################
        tmp = main_bval_true.tmean(unscale = True)
        self.data["main_true_mean_syn1"].append(tmp[0])
        self.data["main_true_mean_syn2"].append(tmp[1])
        ####################### main standard deviations #######################
        tmp = bval["main"].tstd(unscale = True)
        self.data["main_std_syn1"].append(tmp[0])
        self.data["main_std_syn2"].append(tmp[1])
        #################### main true standard deviations #####################
        tmp = main_bval_true.tstd(unscale = True)
        self.data["main_true_std_syn1"].append(tmp[0])
        self.data["main_true_std_syn2"].append(tmp[1])
        ##################### main true genetic variances ######################
        tmp = gpmod.var_A(genome["main"])
        self.data["main_true_var_A_syn1"].append(tmp[0])
        self.data["main_true_var_A_syn2"].append(tmp[1])
        ###################### main true genic variances #######################
        tmp = gpmod.var_a(genome["main"])
        self.data["main_true_var_a_syn1"].append(tmp[0])
        self.data["main_true_var_a_syn2"].append(tmp[1])
        ####################### main true bulmer ratios ########################
        tmp = gpmod.bulmer(genome["main"])
        self.data["main_true_bulmer_syn1"].append(tmp[0])
        self.data["main_true_bulmer_syn2"].append(tmp[1])
        ################### main true lower selection limits ###################
        tmp = gpmod.usl(genome["main"], unscale = True)
        self.data["main_true_usl_syn1"].append(tmp[0])
        self.data["main_true_usl_syn2"].append(tmp[1])
        ################### main true lower selection limits ###################
        tmp = gpmod.lsl(genome["main"], unscale = True)
        self.data["main_true_lsl_syn1"].append(tmp[0])
        self.data["main_true_lsl_syn2"].append(tmp[1])
    def reset(self):
        self.data = {
            "rep": [],
            "t_cur": [],
            "cand_meh": [],
            "cand_mean_syn1": [],
            "cand_mean_syn2": [],
            "cand_true_mean_syn1": [],
            "cand_true_mean_syn2": [],
            "cand_std_syn1": [],
            "cand_std_syn2": [],
            "cand_true_std_syn1": [],
            "cand_true_std_syn2": [],
            "cand_true_var_A_syn1": [],
            "cand_true_var_A_syn2": [],
            "cand_true_var_a_syn1": [],
            "cand_true_var_a_syn2": [],
            "cand_true_bulmer_syn1": [],
            "cand_true_bulmer_syn2": [],
            "cand_true_usl_syn1": [],
            "cand_true_usl_syn2": [],
            "cand_true_lsl_syn1": [],
            "cand_true_lsl_syn2": [],
            "main_meh": [],
            "main_mean_syn1": [],
            "main_mean_syn2": [],
            "main_true_mean_syn1": [],
            "main_true_mean_syn2": [],
            "main_std_syn1": [],
            "main_std_syn2": [],
            "main_true_std_syn1": [],
            "main_true_std_syn2": [],
            "main_true_var_A_syn1": [],
            "main_true_var_A_syn2": [],
            "main_true_var_a_syn1": [],
            "main_true_var_a_syn2": [],
            "main_true_bulmer_syn1": [],
            "main_true_bulmer_syn2": [],
            "main_true_usl_syn1": [],
            "main_true_usl_syn2": [],
            "main_true_lsl_syn1": [],
            "main_true_lsl_syn2": [],
        }
        self.data_frontier = {
            "rep": [],
            "t_cur": [],
            "frontier": [],
        }
        self.rep = 0
    def write(self, filename):
        pandas_df = pandas.DataFrame(self.data)
        pandas_df.to_csv(filename, index = False)
    def write_frontier(self, filename):
        tmp_df_ls = []
        for i in range(len(self.data_frontier["frontier"])):
            tmp_df = pandas.DataFrame(
                data = self.data_frontier["frontier"][i].soln_obj,
                columns = ["syn1", "syn2"],
            )
            tmp_df["t_cur"] = self.data["t_cur"][i]
            tmp_df_ls.append(tmp_df)
        pandas_df = pandas.concat(tmp_df_ls)
        pandas_df.to_csv(filename, index = False)

### create init operators ###
init_pselop = MyInitParentSelectionOperator()
init_mateop = MyInitMatingOperator(mate2waydh.progeny_counter, mate2waydh.family_counter)
init_evalop = MyInitEvaluationOperator(algmod_true, ptprot.var_err)
init_sselop = MyInitSurvivorSelectionOperator(nindiv_fam)

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
    burnin = nburnin
)
pselop = MyParentSelectionOperator()
mateop = MyMatingOperator(mate2waydh.progeny_counter, mate2waydh.family_counter)
evalop = MyEvaluationOperator(algmod_true, ptprot.var_err)
sselop = MySurvivorSelectionOperator(nindiv_fam)
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
rsprog.evolve(nrep = 1, ngen = ngen, lbook = lbook, verbose = True)

lbook.write("multiobjective_genomic_selection_program.csv")
lbook.write_frontier("multiobjective_genomic_selection_program_frontier.csv")
