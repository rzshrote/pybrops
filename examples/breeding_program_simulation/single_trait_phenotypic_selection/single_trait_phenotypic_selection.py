#!/usr/bin/env python3

###
### Single-Trait Phenotypic Selection
### #################################

##
## Loading Required Modules and Seeding the global PRNG
## ====================================================

# import libraries
import numpy
import pandas
import pybrops
from matplotlib import pyplot
from pybrops.breed.prot.bv.MeanPhenotypicBreedingValue import MeanPhenotypicBreedingValue
from pybrops.breed.prot.mate.TwoWayCross import TwoWayCross
from pybrops.breed.prot.mate.TwoWayDHCross import TwoWayDHCross
from pybrops.breed.prot.pt.G_E_Phenotyping import G_E_Phenotyping
from pybrops.breed.prot.sel.EstimatedBreedingValueSelection import EstimatedBreedingValueSubsetSelection
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.opt.algo.SortingSubsetOptimizationAlgorithm import SortingSubsetOptimizationAlgorithm
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.gmap.StandardGeneticMap import StandardGeneticMap
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

# seed python random and numpy random
pybrops.core.random.prng.seed(648201539)

##
## Loading Genetic Map Data from a Text File
## =========================================

# read genetic map
gmap = StandardGeneticMap.from_csv(
    "McMullen_2009_US_NAM.gmap",
    vrnt_chrgrp_col     = "chr",
    vrnt_phypos_col     = "pos",
    vrnt_genpos_col     = "cM",
    vrnt_genpos_units   = "cM",
    auto_group          = True,
    auto_build_spline   = True,
    sep                 = "\t",
    header              = 0,
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
panel_pgmat = DensePhasedGenotypeMatrix.from_vcf(
    "widiv_2000SNPs.vcf.gz", # file name to load
    auto_group_vrnt = True,  # automatically sort and group variants
)

# interpolate genetic map positions
panel_pgmat.interp_xoprob(gmap, gmapfn)

##
## Constructing a Single-Trait Genomic Model
## =========================================

# model intercepts: (1,ntrait)
beta = numpy.array([[0.0]], dtype = float)

# marker effects: (nvrnt,1)
mkreffect = numpy.random.normal(
    loc = 0.0,
    scale = 0.05,
    size = (panel_pgmat.nvrnt,1)
)

# trait names: (ntrait,)
trait = numpy.array(["Syn1"], dtype = object)

# create an additive linear genomic model to model traits
algmod = DenseAdditiveLinearGenomicModel(
    beta        = beta,                 # model intercepts
    u_misc      = None,                 # miscellaneous random effects
    u_a         = mkreffect,            # random marker effects
    trait       = trait,                # trait names
    model_name  = "synthetic_model",    # name of the model
    hyperparams = None                  # model parameters
)

##
## Select founders and randomly intermate for 20 generations
## =========================================================

# founder population parameters
nfndr = 40          # number of random founders to select (must be even)
fndr_nmating = 1    # number of times to perform cross configuration (only needed for 3+ way crosses)
fndr_nprogeny = 80  # number of progenies to derive from cross configuration
fndr_nrandmate = 20 # number of random mating generations

# create 2-way cross object
mate2way = TwoWayCross()

# randomly select and pair ``nfndr`` founders
xconfig = numpy.random.choice(panel_pgmat.ntaxa,nfndr).reshape(nfndr//2,2)

# randomly intermate ``nfndr`` founders to create initial hybrids
fndr_pgmat = mate2way.mate(
    pgmat = panel_pgmat,
    xconfig = xconfig,
    nmating = fndr_nmating,
    nprogeny = fndr_nprogeny,
)

# randomly intermate for ``fndr_nrandmate`` generations
# each individual in the population is randomly mated with another individual
# and creates a single progeny so that the population size is held constant
for _ in range(fndr_nrandmate):
    # get the number of taxa
    ntaxa = fndr_pgmat.ntaxa
    # randomly select and pair ``ntaxa`` parents
    xconfig = numpy.empty((ntaxa,2), dtype = int)
    xconfig[:,0] = numpy.random.choice(ntaxa, ntaxa, replace = False)
    xconfig[:,1] = numpy.random.choice(ntaxa, ntaxa, replace = False)
    # randomly intermate ``ntaxa`` parents
    fndr_pgmat = mate2way.mate(
        pgmat = fndr_pgmat,
        xconfig = xconfig,
        nmating = 1,
        nprogeny = 1,
    )

# create a 2-way DH cross object, use the counters from the 2-way cross object
mate2waydh = TwoWayDHCross(
    progeny_counter = mate2way.progeny_counter,
    family_counter  = mate2way.family_counter,
)

# get the number of taxa
ntaxa = fndr_pgmat.ntaxa

# randomly select and pair ``ntaxa`` parents
xconfig = numpy.empty((ntaxa,2), dtype = int)
xconfig[:,0] = numpy.random.choice(ntaxa, ntaxa, replace = False)
xconfig[:,1] = numpy.random.choice(ntaxa, ntaxa, replace = False)

# DH all individuals in the founder population to create our initial breeding population
pgmat = mate2waydh.mate(
    pgmat = fndr_pgmat,
    xconfig = xconfig,
    nmating = 1,
    nprogeny = 1,
)

##
## Create a Phenotyping Protocol Object
## ====================================

# create a phenotyping protocol object to simulate 4 environments with 1 rep each
ptprot = G_E_Phenotyping(
    gpmod = algmod,
    nenv = 4,
    nrep = 1,
)

# set the trait heritability using the initial population
# initial population fits heritability assumptions of being randomly mated
ptprot.set_h2(0.4, pgmat)

##
## Create a Breeding Value Estimation Protocol Object
## ==================================================

# estimate breeding value using mean across environments for simplicity
bvprot = MeanPhenotypicBreedingValue(
    taxa_col = "taxa",
    trait_cols = "Syn1",
)

##
## Create a Selection Protocol Object
## ==================================

# use a hillclimber for the single-objective optimization algorithm
# this is a very general algorithm and may not be the most efficient for
# all single-objective optimizations
soalgo = SortingSubsetOptimizationAlgorithm()

# create a selection protocol that selects based on EBVs
selprot = EstimatedBreedingValueSubsetSelection(
    ntrait      = 1,        # number of expected traits
    ncross      = 20,       # number of cross configurations
    nparent     = 2,        # number of parents per cross configuration
    nmating     = 1,        # number of matings per cross configuration
    nprogeny    = 80,       # number of progeny per mating event
    nobj        = 1,        # number of objectives == ntrait
    soalgo      = soalgo,   # use hillclimber to solve single-objective problem
)

##
## Simulate disjoint recurrent phenotypic selection for 60 generations
## ===================================================================

# number of generations for which to simulate selection
ngen = 60

# make a dictionary logbook
lbook = {
    "gen"           : [],
    "meh"           : [],
    "lsl"           : [],
    "usl"           : [],
    "tbv_min_Syn1"  : [],
    "tbv_mean_Syn1" : [],
    "tbv_max_Syn1"  : [],
    "tbv_std_Syn1"  : [],
    "ebv_min_Syn1"  : [],
    "ebv_mean_Syn1" : [],
    "ebv_max_Syn1"  : [],
    "ebv_std_Syn1"  : [],
}

#
# Simulation Initialization
# -------------------------

# initial phenotyping
pheno_df = ptprot.phenotype(pgmat)

# initial breeding value estimation
bvmat = bvprot.estimate(ptobj=pheno_df)

# log metrics
lbook["gen"].append(0)
lbook["meh"].append(pgmat.meh())
lbook["lsl"].append(algmod.lsl(pgmat)[0])
lbook["usl"].append(algmod.usl(pgmat)[0])
tbv = algmod.gebv(pgmat).unscale()
lbook["tbv_min_Syn1"].append(tbv.min(0)[0])
lbook["tbv_mean_Syn1"].append(tbv.mean(0)[0])
lbook["tbv_max_Syn1"].append(tbv.max(0)[0])
lbook["tbv_std_Syn1"].append(tbv.std(0)[0])
ebv = bvmat.unscale()
lbook["ebv_min_Syn1"].append(ebv.min(0)[0])
lbook["ebv_mean_Syn1"].append(ebv.mean(0)[0])
lbook["ebv_max_Syn1"].append(ebv.max(0)[0])
lbook["ebv_std_Syn1"].append(ebv.std(0)[0])
print("Gen: {0}".format(0))

#
# Main Simulation Loop
# --------------------

# simulate for ``ngen`` generations
for gen in range(1,ngen+1):
    # select individuals
    selcfg = selprot.select(
        pgmat   = pgmat,    # genomes from which to build SelectionConfiguration
        gmat    = None,     # not required by this selection protocol
        ptdf    = None,     # not required by this selection protocol
        bvmat   = bvmat,    # breeding values (required)
        gpmod   = None,     # not required by this selection protocol
        t_cur   = 0,        # not required by this selection protocol
        t_max   = 0,        # not required by this selection protocol
    )
    # mate individuals
    pgmat = mate2waydh.mate(
        pgmat = selcfg.pgmat,
        xconfig = selcfg.xconfig,
        nmating = selcfg.nmating,
        nprogeny = selcfg.nprogeny,
    )
    # phenotype progenies
    pheno_df = ptprot.phenotype(pgmat)
    # estimate breeding values for progenies
    bvmat = bvprot.estimate(ptobj=pheno_df)
    # log metrics
    lbook["gen"].append(gen)
    lbook["meh"].append(pgmat.meh())
    lbook["lsl"].append(algmod.lsl(pgmat)[0])
    lbook["usl"].append(algmod.usl(pgmat)[0])
    tbv = algmod.gebv(pgmat).unscale()
    lbook["tbv_min_Syn1"].append(tbv.min(0)[0])
    lbook["tbv_mean_Syn1"].append(tbv.mean(0)[0])
    lbook["tbv_max_Syn1"].append(tbv.max(0)[0])
    lbook["tbv_std_Syn1"].append(tbv.std(0)[0])
    ebv = bvmat.unscale()
    lbook["ebv_min_Syn1"].append(ebv.min(0)[0])
    lbook["ebv_mean_Syn1"].append(ebv.mean(0)[0])
    lbook["ebv_max_Syn1"].append(ebv.max(0)[0])
    lbook["ebv_std_Syn1"].append(ebv.std(0)[0])
    print("Gen: {0}".format(gen))

#
# Saving Results to a File
# ------------------------

# create output dataframe and save
lbook_df = pandas.DataFrame(lbook)
lbook_df.to_csv("lbook.csv", sep = ",", index = False)

##
## Visualizing Breeding Program Simulation Results with ``matplotlib``
## ===================================================================

#
# Visualizing True Breeding Values (TBVs)
# ---------------------------------------

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.plot(lbook_df["gen"], lbook_df["tbv_min_Syn1"], label = "Min Population TBV")
ax.plot(lbook_df["gen"], lbook_df["tbv_mean_Syn1"], label = "Mean Population TBV")
ax.plot(lbook_df["gen"], lbook_df["tbv_max_Syn1"], label = "Max Population TBV")
ax.plot(lbook_df["gen"], lbook_df["lsl"], label = "Lower Selection Limit")
ax.plot(lbook_df["gen"], lbook_df["usl"], label = "Upper Selection Limit")
ax.set_title("Single-Trait Recurrent Phenotypic Selection")
ax.set_xlabel("Generation")
ax.set_ylabel("Synthetic Trait Breeding Value")
ax.legend()
pyplot.savefig("true_breeding_values.png", dpi = 300)
pyplot.close(fig)

#
# Visualizing True Breeding Values (TBVs)
# ---------------------------------------

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.plot(lbook_df["gen"], lbook_df["ebv_min_Syn1"], label = "Min Population EBV")
ax.plot(lbook_df["gen"], lbook_df["ebv_mean_Syn1"], label = "Mean Population EBV")
ax.plot(lbook_df["gen"], lbook_df["ebv_max_Syn1"], label = "Max Population EBV")
ax.plot(lbook_df["gen"], lbook_df["lsl"], label = "Lower Selection Limit")
ax.plot(lbook_df["gen"], lbook_df["usl"], label = "Upper Selection Limit")
ax.set_title("Single-Trait Recurrent Phenotypic Selection")
ax.set_xlabel("Generation")
ax.set_ylabel("Synthetic Trait Breeding Value")
ax.legend()
pyplot.savefig("estimated_breeding_values.png", dpi = 300)
pyplot.close(fig)
