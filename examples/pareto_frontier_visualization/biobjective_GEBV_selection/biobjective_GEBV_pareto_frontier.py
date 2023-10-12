#!/usr/bin/env python3

###
### Bi-Objective GEBV Selection Pareto Frontier Visualization
### #########################################################

##
## Loading Required Modules and Seeding the global PRNG
## ====================================================

# import libraries
import numpy
from matplotlib import pyplot
import pybrops
from pybrops.breed.prot.sel.GenomicEstimatedBreedingValueSelection import GenomicEstimatedBreedingValueSubsetSelection
from pybrops.opt.algo.NSGA2SubsetGeneticAlgorithm import NSGA2SubsetGeneticAlgorithm
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

# seed python random and numpy random
pybrops.core.random.prng.seed(31621463)

##
## Loading Genotypic Data from a VCF File
## ======================================

# read unphased genetic markers from a vcf file
gmat = DenseGenotypeMatrix.from_vcf(
    "widiv_2000SNPs.vcf.gz", # file name to load
    auto_group_vrnt = True,  # automatically sort and group variants
)

##
## Constructing a 2-Trait Genomic Model
## ====================================

# make marker effects for two traits which are competing in nature
# marker effects array is of shape (nvrnt, 2)
mkreffect = numpy.random.multivariate_normal(
    mean    = numpy.array([0.0, 0.0]), 
    cov     = numpy.array([
                [ 1.0, -0.4],
                [-0.4,  1.0]
            ]), 
    size    = gmat.nvrnt
)

# create an additive linear genomic model to model traits
algmod = DenseAdditiveLinearGenomicModel(
    beta    = numpy.float64([[10.0, 25.0]]), # model intercepts
    u_misc  = None,                          # miscellaneous random effects
    u_a     = mkreffect,                     # random marker effects
    trait   = numpy.array(                   # trait names
                ["syn1","syn2"],
                dtype=object
            ),
    model_name = "synthetic_model",          # name of the model
    hyperparams = None                       # model parameters
)

##
## Constructing a GEBV Subset Selection object
## ===========================================

# create custom multi-objective algorithm for optimization
# use NSGA-II and evolve for 1000 generations
moalgo = NSGA2SubsetGeneticAlgorithm(
    ngen = 1000,    # number of generations to evolve
    pop_size = 100  # number of parents in population
)

# construct a subset selection object for GEBV selection
selprot = GenomicEstimatedBreedingValueSubsetSelection(
    ntrait = 2,         # number of traits to expect from GEBV matrix
    unscale = True,     # unscale GEBVs so that units are in unscaled trait values
    ncross = 10,        # number of breeding crosses to select
    nparent = 2,        # number of parents per breeding cross to select
    nmating = 1,        # number of times parents are mated per cross
    nprogeny = 40,      # number of progenies to derive from a mating event
    nobj = 2,           # number of objectives (ntrait)
    moalgo = moalgo,    # custom multi-objective algorithm
)

##
## Estimating the Pareto Frontier
## ==============================

# estimate pareto frontier using optimization algorithm
selsoln = selprot.mosolve(
    pgmat = None,       # argument not utilized
    gmat = gmat,        # ``gmat`` argument required
    ptdf = None,        # argument not utilized
    bvmat = None,       # argument not utilized
    gpmod = algmod,     # ``gpmod`` argument required
    t_cur = 0,          # argument not utilized
    t_max = 0,          # argument not utilized
)

##
## Visualizing the Pareto Frontier with ``matplotlib``
## ===================================================

# get the pareto frontier
# negate the objectives to get the mean GEBV since optimization problems are always minimizing
xdata = -selsoln.soln_obj[:,0]
ydata = -selsoln.soln_obj[:,1]

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.scatter(xdata, ydata)
ax.set_title("Bi-Objective GEBV Selection Pareto Frontier")
ax.set_xlabel("Synthetic Trait 1 Mean GEBV")
ax.set_ylabel("Synthetic Trait 2 Mean GEBV")
pyplot.savefig("biobjective_GEBV_pareto_frontier.png", dpi = 250)
pyplot.close(fig)
