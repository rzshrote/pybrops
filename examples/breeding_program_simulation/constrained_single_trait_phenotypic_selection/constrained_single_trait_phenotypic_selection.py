#!/usr/bin/env python3

###
### Constrained Single-Trait Phenotypic Selection
### #############################################

##
## Loading Required Modules and Seeding the global PRNG
## ====================================================

# import libraries
from numbers import Real
import numpy
import pandas
import pybrops
from matplotlib import pyplot
from pybrops.breed.prot.bv.MeanPhenotypicBreedingValue import MeanPhenotypicBreedingValue
from pybrops.breed.prot.mate.TwoWayCross import TwoWayCross
from pybrops.breed.prot.mate.TwoWayDHCross import TwoWayDHCross
from pybrops.breed.prot.pt.G_E_Phenotyping import G_E_Phenotyping
from pybrops.breed.prot.sel.EstimatedBreedingValueSelection import EstimatedBreedingValueSubsetSelection
from pybrops.breed.prot.sel.OptimalContributionSelection import OptimalContributionSubsetSelection
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.opt.algo.SortingSubsetOptimizationAlgorithm import SortingSubsetOptimizationAlgorithm
from pybrops.opt.algo.SteepestDescentSubsetHillClimber import SteepestDescentSubsetHillClimber
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.popgen.cmat.fcty.DenseMolecularCoancestryMatrixFactory import DenseMolecularCoancestryMatrixFactory
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.gmap.StandardGeneticMap import StandardGeneticMap
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

# seed python random and numpy random
pybrops.core.random.prng.seed(52347529)

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

# randomly select and pair 20 parents
xconfig = numpy.random.choice(ntaxa, 40, replace = False).reshape(20,2)

# DH all individuals in the founder population to create our initial breeding population
fndr_pgmat = mate2waydh.mate(
    pgmat = fndr_pgmat,
    xconfig = xconfig,
    nmating = 1,
    nprogeny = 80,
)

##
## Create a Phenotyping Protocol Object
## ====================================

# create a phenotyping protocol object to simulate 4 environments with 1 rep each
ptprot = G_E_Phenotyping(gpmod = algmod, nenv = 4, nrep = 1)

# set the trait heritability using the initial population
# initial population fits heritability assumptions of being randomly mated
ptprot.set_h2(0.4, fndr_pgmat)

##
## Create a Breeding Value Estimation Protocol Object
## ==================================================

# estimate breeding value using mean across environments for simplicity
bvprot = MeanPhenotypicBreedingValue(
    taxa_col = "taxa",
    taxa_grp_col = "taxa_grp",
    trait_cols = "Syn1",
)

##
## Create a Within-Family Selection Function
## =========================================

# define function to do within family selection
def within_family_selection(bvmat: DenseBreedingValueMatrix, nindiv: int):
    order = bvmat.mat.argsort(0)[:,0]
    mask = numpy.full(len(order), False, bool)
    groups = numpy.unique(bvmat.taxa_grp)
    for group in groups:
        tmp = order[bvmat.taxa_grp == group]
        tmp.sort()
        ix = tmp[:nindiv]
        for i in ix:
            mask[order == i] = True
    indices = numpy.flatnonzero(mask)
    return indices

##
## Create a Constrained Selection Protocol Object
## ==============================================

# create a dense molecular coancestry matrix factory
cmatfcty = DenseMolecularCoancestryMatrixFactory()

# define an objective transformation function
def obj_trans(
        decnvec: numpy.ndarray,
        latentvec: numpy.ndarray, 
        **kwargs: dict
    ) -> numpy.ndarray:
    """
    Receive an incoming vector of [MGR,BV1,...,BVn] and transform it to
    [BV1,...,BVn].

    Where::
    
        - MGR is the mean genomic relationship (kinship; in range [0,1]).
        - BVn is the nth mean breeding value for the subset.

    Parameters
    ----------
    decnvec : numpy.ndarray
        A decision space vector of shape (ndecn,)
    latentvec : numpy.ndarray
        A latent space function vector of shape (1+ntrait,)
    
    Returns
    -------
    out : numpy.ndarray
        A vector of shape (ntrait,).
    """
    # extract trait(s) as objective(s)
    return latentvec[1:]

# define an inequality constraint violation function
def ineqcv_trans(
        decnvec: numpy.ndarray,
        latentvec: numpy.ndarray, 
        maxinb: Real,
        **kwargs: dict
    ) -> numpy.ndarray:
    """
    A custom inequality constraint violation function.

    Parameters
    ----------
    decnvec : numpy.ndarray
        A decision space vector of shape (ndecn,)
    latentvec : numpy.ndarray
        A latent space function vector of shape (1+ntrait,)
    minvec : numpy.ndarray
        Vector of minimum values for which the latent vector can take.
    
    Returns
    -------
    out : numpy.ndarray
        An inequality constraint violation vector of shape (1,).
    """
    # calculate constraint violation for inbreeding
    out = numpy.array([max(latentvec[0] - maxinb, 0.0)], dtype = float)
    # return inequality constraint violation array
    return out

# use a hillclimber for the single-objective optimization algorithm
soalgo = SteepestDescentSubsetHillClimber()

# create a selection protocol that selects based on EBVs with an inbreeding constraint
constrained_selprot = OptimalContributionSubsetSelection(
    ntrait       = 1,            # number of expected traits
    cmatfcty     = cmatfcty,     # coancestry/kinship matrix factory
    unscale      = True,         # unscale breeding values to human-readable format
    ncross       = 20,           # number of cross configurations
    nparent      = 2,            # number of parents per cross configuration
    nmating      = 1,            # number of matings per cross configuration
    nprogeny     = 80,           # number of progeny per mating event
    nobj         = 1,            # number of objectives == ntrait
    obj_trans    = obj_trans,    # latent vector transformation to create objective function
    nineqcv      = 1,            # number of inequality constraint violations
    ineqcv_trans = ineqcv_trans, # latent vector transformation to create inequality constraints
    ineqcv_trans_kwargs = {      # keyword arguments
        "maxinb": 1.0
    },
    soalgo       = soalgo,       # use hillclimber to solve single-objective problem
)

##
## Create an Unconstrained Selection Protocol Object
## =================================================

# use a sorting algorithm for the single-objective optimization algorithm
soalgo = SortingSubsetOptimizationAlgorithm()

# create a selection protocol that selects based on EBVs with an inbreeding constraint
unconstrained_selprot = EstimatedBreedingValueSubsetSelection(
    ntrait       = 1,            # number of expected traits
    ncross       = 20,           # number of cross configurations
    nparent      = 2,            # number of parents per cross configuration
    nmating      = 1,            # number of matings per cross configuration
    nprogeny     = 80,           # number of progeny per mating event
    nobj         = 1,            # number of objectives == ntrait
    soalgo       = soalgo,       # use sorting algorithm to solve single-objective problem
)

##
## Simulate Constrained and Unconstrained Phenotypic Selection for 60 Generations
## ==============================================================================

#
# Rudimentary Logbooks
# --------------------

# make a dictionary logbook
constrained_lbook = {
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
# Constrained Simulation Initialization
# -------------------------------------

# copy founder population
pgmat = fndr_pgmat.deepcopy()

# initial phenotyping
pheno_df = ptprot.phenotype(pgmat)

# initial breeding value estimation
bvmat = bvprot.estimate(pheno_df, pgmat)

# get candidate indices using within family selection
indices = within_family_selection(bvmat, 8) # select top 10%

# get parental candidates
cand_pgmat = pgmat.select_taxa(indices)
cand_bvmat = bvmat.select_taxa(indices)

# log metrics
constrained_lbook["gen"].append(0)
constrained_lbook["meh"].append(pgmat.meh())
constrained_lbook["lsl"].append(algmod.lsl(pgmat)[0])
constrained_lbook["usl"].append(algmod.usl(pgmat)[0])
tbv = algmod.gebv(pgmat).unscale()
constrained_lbook["tbv_min_Syn1"].append(tbv.min(0)[0])
constrained_lbook["tbv_mean_Syn1"].append(tbv.mean(0)[0])
constrained_lbook["tbv_max_Syn1"].append(tbv.max(0)[0])
constrained_lbook["tbv_std_Syn1"].append(tbv.std(0)[0])
ebv = bvmat.unscale()
constrained_lbook["ebv_min_Syn1"].append(ebv.min(0)[0])
constrained_lbook["ebv_mean_Syn1"].append(ebv.mean(0)[0])
constrained_lbook["ebv_max_Syn1"].append(ebv.max(0)[0])
constrained_lbook["ebv_std_Syn1"].append(ebv.std(0)[0])
print("Gen: {0}".format(0))

#
# Constrained Simulation Main Loop
# --------------------------------

# number of generations for which to simulate selection
ngen = 60

# create evenly spaced maximum inbreeding allowed across ``ngen`` generations
maxinb = numpy.linspace(0.77, 1.0, ngen+1)

# simulate for ``ngen`` generations
for gen in range(1,ngen+1):
    # get candidate mask using within family selection
    indices = within_family_selection(bvmat, 8) # select top 10%
    # get parental candidates
    cand_pgmat = pgmat.select_taxa(indices)
    cand_bvmat = bvmat.select_taxa(indices)
    # set the inbreeding constraint
    constrained_selprot.ineqcv_trans_kwargs["maxinb"] = maxinb[gen]
    # select individuals
    selcfg = constrained_selprot.select(
        pgmat   = cand_pgmat,   # genomes from which to build SelectionConfiguration
        gmat    = cand_pgmat,   # genotypes (required)
        ptdf    = None,         # not required by this selection protocol
        bvmat   = cand_bvmat,   # breeding values (required)
        gpmod   = None,         # not required by this selection protocol
        t_cur   = 0,            # not required by this selection protocol
        t_max   = 0,            # not required by this selection protocol
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
    # estimate breeding values for progenies and align to pgmat
    bvmat = bvprot.estimate(pheno_df, pgmat)
    # log metrics
    constrained_lbook["gen"].append(gen)
    constrained_lbook["meh"].append(pgmat.meh())
    constrained_lbook["lsl"].append(algmod.lsl(pgmat)[0])
    constrained_lbook["usl"].append(algmod.usl(pgmat)[0])
    tbv = algmod.gebv(pgmat).unscale()
    constrained_lbook["tbv_min_Syn1"].append(tbv.min(0)[0])
    constrained_lbook["tbv_mean_Syn1"].append(tbv.mean(0)[0])
    constrained_lbook["tbv_max_Syn1"].append(tbv.max(0)[0])
    constrained_lbook["tbv_std_Syn1"].append(tbv.std(0)[0])
    ebv = bvmat.unscale()
    constrained_lbook["ebv_min_Syn1"].append(ebv.min(0)[0])
    constrained_lbook["ebv_mean_Syn1"].append(ebv.mean(0)[0])
    constrained_lbook["ebv_max_Syn1"].append(ebv.max(0)[0])
    constrained_lbook["ebv_std_Syn1"].append(ebv.std(0)[0])
    print("Gen: {0}".format(gen))

#
# Saving Results to a File
# ------------------------

# create output dataframe and save
constrained_lbook_df = pandas.DataFrame(constrained_lbook)
constrained_lbook_df.to_csv("constrained_lbook.csv", sep = ",", index = False)

#
# Rudimentary Logbooks
# --------------------

# make a dictionary logbook
unconstrained_lbook = {
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
# Constrained Simulation Initialization
# -------------------------------------

# copy founder population
pgmat = fndr_pgmat.deepcopy()

# initial phenotyping
pheno_df = ptprot.phenotype(pgmat)

# initial breeding value estimation
bvmat = bvprot.estimate(pheno_df, pgmat)

# get candidate indices using within family selection
indices = within_family_selection(bvmat, 8) # select top 10%

# get parental candidates
cand_pgmat = pgmat.select_taxa(indices)
cand_bvmat = bvmat.select_taxa(indices)

# log metrics
unconstrained_lbook["gen"].append(0)
unconstrained_lbook["meh"].append(pgmat.meh())
unconstrained_lbook["lsl"].append(algmod.lsl(pgmat)[0])
unconstrained_lbook["usl"].append(algmod.usl(pgmat)[0])
tbv = algmod.gebv(pgmat).unscale()
unconstrained_lbook["tbv_min_Syn1"].append(tbv.min(0)[0])
unconstrained_lbook["tbv_mean_Syn1"].append(tbv.mean(0)[0])
unconstrained_lbook["tbv_max_Syn1"].append(tbv.max(0)[0])
unconstrained_lbook["tbv_std_Syn1"].append(tbv.std(0)[0])
ebv = bvmat.unscale()
unconstrained_lbook["ebv_min_Syn1"].append(ebv.min(0)[0])
unconstrained_lbook["ebv_mean_Syn1"].append(ebv.mean(0)[0])
unconstrained_lbook["ebv_max_Syn1"].append(ebv.max(0)[0])
unconstrained_lbook["ebv_std_Syn1"].append(ebv.std(0)[0])
print("Gen: {0}".format(0))

#
# Constrained Simulation Main Loop
# --------------------------------

# number of generations for which to simulate selection
ngen = 60

# simulate for ``ngen`` generations
for gen in range(1,ngen+1):
    # get candidate mask using within family selection
    indices = within_family_selection(bvmat, 8) # select top 10%
    # get parental candidates
    cand_pgmat = pgmat.select_taxa(indices)
    cand_bvmat = bvmat.select_taxa(indices)
    # select individuals
    selcfg = unconstrained_selprot.select(
        pgmat   = cand_pgmat,   # genomes from which to build SelectionConfiguration
        gmat    = cand_pgmat,   # genotypes (required)
        ptdf    = None,         # not required by this selection protocol
        bvmat   = cand_bvmat,   # breeding values (required)
        gpmod   = None,         # not required by this selection protocol
        t_cur   = 0,            # not required by this selection protocol
        t_max   = 0,            # not required by this selection protocol
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
    # estimate breeding values for progenies and align to pgmat
    bvmat = bvprot.estimate(pheno_df, pgmat)
    # log metrics
    unconstrained_lbook["gen"].append(gen)
    unconstrained_lbook["meh"].append(pgmat.meh())
    unconstrained_lbook["lsl"].append(algmod.lsl(pgmat)[0])
    unconstrained_lbook["usl"].append(algmod.usl(pgmat)[0])
    tbv = algmod.gebv(pgmat).unscale()
    unconstrained_lbook["tbv_min_Syn1"].append(tbv.min(0)[0])
    unconstrained_lbook["tbv_mean_Syn1"].append(tbv.mean(0)[0])
    unconstrained_lbook["tbv_max_Syn1"].append(tbv.max(0)[0])
    unconstrained_lbook["tbv_std_Syn1"].append(tbv.std(0)[0])
    ebv = bvmat.unscale()
    unconstrained_lbook["ebv_min_Syn1"].append(ebv.min(0)[0])
    unconstrained_lbook["ebv_mean_Syn1"].append(ebv.mean(0)[0])
    unconstrained_lbook["ebv_max_Syn1"].append(ebv.max(0)[0])
    unconstrained_lbook["ebv_std_Syn1"].append(ebv.std(0)[0])
    print("Gen: {0}".format(gen))

#
# Saving Results to a File
# ------------------------

# create output dataframe and save
unconstrained_lbook_df = pandas.DataFrame(unconstrained_lbook)
unconstrained_lbook_df.to_csv("unconstrained_lbook.csv", sep = ",", index = False)

##
## Visualizing Breeding Program Simulation Results with ``matplotlib``
## ===================================================================

#
# Visualizing True Breeding Values (TBVs)
# ---------------------------------------

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.plot(constrained_lbook_df["gen"],   constrained_lbook_df["tbv_mean_Syn1"],   '-b',  label = "Const. Sel.: Mean Pop. TBV")
ax.plot(constrained_lbook_df["gen"],   constrained_lbook_df["lsl"],             ':b',  label = "Const. Sel.: LSL")
ax.plot(constrained_lbook_df["gen"],   constrained_lbook_df["usl"],             '--b', label = "Const. Sel.: USL")
ax.plot(unconstrained_lbook_df["gen"], unconstrained_lbook_df["tbv_mean_Syn1"], '-r',  label = "Unconst. Sel.: Mean Pop. TBV")
ax.plot(unconstrained_lbook_df["gen"], unconstrained_lbook_df["lsl"],           ':r',  label = "Unconst. Sel.: LSL")
ax.plot(unconstrained_lbook_df["gen"], unconstrained_lbook_df["usl"],           '--r', label = "Unconst. Sel.: USL")
ax.set_title("Single-Trait Recurrent Phenotypic Selection")
ax.set_xlabel("Generation")
ax.set_ylabel("Synthetic Trait Breeding Value")
ax.legend()
pyplot.savefig("constrained_single_trait_phenotypic_selection_true_breeding_values.png", dpi = 300)
pyplot.close(fig)

#
# Visualizing Estimated Breeding Values (TBVs)
# --------------------------------------------

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.plot(constrained_lbook_df["gen"],   constrained_lbook_df["ebv_mean_Syn1"],   '-b',  label = "Const. Sel.: Mean Pop. EBV")
ax.plot(constrained_lbook_df["gen"],   constrained_lbook_df["lsl"],             ':b',  label = "Const. Sel.: LSL")
ax.plot(constrained_lbook_df["gen"],   constrained_lbook_df["usl"],             '--b', label = "Const. Sel.: USL")
ax.plot(unconstrained_lbook_df["gen"], unconstrained_lbook_df["ebv_mean_Syn1"], '-r',  label = "Unconst. Sel.: Mean Pop. EBV")
ax.plot(unconstrained_lbook_df["gen"], unconstrained_lbook_df["lsl"],           ':r',  label = "Unconst. Sel.: LSL")
ax.plot(unconstrained_lbook_df["gen"], unconstrained_lbook_df["usl"],           '--r', label = "Unconst. Sel.: USL")
ax.set_title("Single-Trait Recurrent Phenotypic Selection")
ax.set_xlabel("Generation")
ax.set_ylabel("Synthetic Trait Breeding Value")
ax.legend()
pyplot.savefig("constrained_single_trait_phenotypic_selection_estimated_breeding_values.png", dpi = 300)
pyplot.close(fig)

#
# Visualizing Mean Expected Heterozygosity (MEH)
# ==============================================

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.plot(constrained_lbook_df["gen"],   constrained_lbook_df["meh"],   '-b',  label = "Const. Sel.: Pop. MEH")
ax.plot(unconstrained_lbook_df["gen"], unconstrained_lbook_df["meh"], '-r',  label = "Unconst. Sel.: Pop. MEH")
ax.set_title("Single-Trait Recurrent Phenotypic Selection")
ax.set_xlabel("Generation")
ax.set_ylabel("Mean Expected Heterozygosity")
ax.legend()
pyplot.savefig("constrained_single_trait_phenotypic_selection_mean_expected_heterozygosity.png", dpi = 300)
pyplot.close(fig)

#
# Visualizing True Breeding Value Standard Deviations
# ===================================================

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.plot(constrained_lbook_df["gen"],   constrained_lbook_df["tbv_std_Syn1"],   '-b',  label = "Const. Sel.: Pop. TBV SD")
ax.plot(unconstrained_lbook_df["gen"], unconstrained_lbook_df["tbv_std_Syn1"], '-r',  label = "Unconst. Sel.: Pop. TBV SD")
ax.set_title("Single-Trait Recurrent Phenotypic Selection")
ax.set_xlabel("Generation")
ax.set_ylabel("Synthetic Trait Breeding Value Standard Deviation")
ax.legend()
pyplot.savefig("constrained_single_trait_phenotypic_selection_true_breeding_value_standard_deviation.png", dpi = 300)
pyplot.close(fig)

#
# Visualizing Estimated Breeding Value Standard Deviations
# ========================================================

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.plot(constrained_lbook_df["gen"],   constrained_lbook_df["ebv_std_Syn1"],   '-b',  label = "Const. Sel.: Pop. EBV SD")
ax.plot(unconstrained_lbook_df["gen"], unconstrained_lbook_df["ebv_std_Syn1"], '-r',  label = "Unconst. Sel.: Pop. EBV SD")
ax.set_title("Single-Trait Recurrent Phenotypic Selection")
ax.set_xlabel("Generation")
ax.set_ylabel("Synthetic Trait Breeding Value Standard Deviation")
ax.legend()
pyplot.savefig("constrained_single_trait_phenotypic_selection_estimated_breeding_value_standard_deviation.png", dpi = 300)
pyplot.close(fig)
