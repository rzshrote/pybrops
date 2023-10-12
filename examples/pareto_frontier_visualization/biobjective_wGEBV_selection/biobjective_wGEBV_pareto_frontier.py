#!/usr/bin/env python3

import numpy
from matplotlib import pyplot

import pybrops
from pybrops.breed.prot.sel.WeightedGenomicSelection import WeightedGenomicSubsetSelection
from pybrops.opt.algo.NSGA2SubsetGeneticAlgorithm import NSGA2SubsetGeneticAlgorithm
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

# seed python random and numpy random
pybrops.core.random.prng.seed(31621463)

################### Load genotypic data ####################

# read unphased genetic markers from a vcf file
gmat = DenseGenotypeMatrix.from_vcf(
    "widiv_2000SNPs.vcf.gz", # file name to load
    auto_group_vrnt = True,  # automatically sort and group variants
)

################# Construct genomic model ##################

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

############# Construct GEBV Selection object ##############

# create custom multi-objective algorithm for optimization
# use NSGA-II and evolve for 1000 generations
moalgo = NSGA2SubsetGeneticAlgorithm(
    ngen = 1000,    # number of generations to evolve
    pop_size = 100  # number of parents in population
)

# construct a subset selection object for weighted genomic selection
selprot = WeightedGenomicSubsetSelection(
    ntrait = 2,         # number of traits to expect from GEBV matrix
    ncross = 10,        # number of breeding crosses to select
    nparent = 2,        # number of parents per breeding cross to select
    nmating = 1,        # number of times parents are mated per cross
    nprogeny = 40,      # number of progenies to derive from a mating event
    nobj = 2,           # number of objectives (ntrait)
    moalgo = moalgo,    # custom multi-objective algorithm
)

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

# get the pareto frontier
# negate 
xdata = -selsoln.soln_obj[:,0]
ydata = -selsoln.soln_obj[:,1]

# establish axis labels
xlabel = "Syn 1 Trait wGEBV Sum"
ylabel = "Syn 2 Trait wGEBV Sum"

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.scatter(xdata, ydata)
ax.set_title("Bi-Objective Weighted Genomic Selection Pareto Frontier")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
pyplot.savefig("wGEBV_2d_frontier.png", dpi = 250)
pyplot.close(fig)
