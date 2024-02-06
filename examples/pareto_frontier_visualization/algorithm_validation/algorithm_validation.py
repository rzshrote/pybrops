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
from pybrops.opt.algo.SteepestDescentSubsetHillClimber import SteepestDescentSubsetHillClimber
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

# seed python random and numpy random
pybrops.core.random.prng.seed(31621463)

from pymoo.config import Config
Config.warnings['not_compiled'] = False

##
## Loading Genotypic Data from a VCF File
## ======================================

# read unphased genetic markers from a vcf file
gmat = DenseGenotypeMatrix.from_vcf(
    "widiv_2000SNPs.vcf.gz", # file name to load
    auto_group_vrnt = True,  # automatically sort and group variants
)
# select first 100 individuals to keep problem size smaller
gmat = gmat.select_taxa(numpy.arange(100))

##
## Constructing a Bi-Trait Genomic Model
## =====================================

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
    trait   = numpy.array(["syn1","syn2"], dtype=object), # trait names
)

##
## Setup for constrained optimization
## ----------------------------------

# suppose we desire to minimize the negated GEBV of our first trait, 
# using our second trait as a constraint

# define an objective transformation function
def obj_trans(
        decnvec: numpy.ndarray,
        latentvec: numpy.ndarray, 
        **kwargs: dict
    ) -> numpy.ndarray:
    """
    A custom objective transformation function.

    Parameters
    ----------
    decnvec : numpy.ndarray
        A decision space vector of shape (ndecn,)
    latentvec : numpy.ndarray
        A latent space function vector of shape (ntrait,)
    
    Returns
    -------
    out : numpy.ndarray
        A vector of shape (1,).
    """
    # extract the first trait as an objective
    return latentvec[0:1]

def obj_trans_max0(
        decnvec: numpy.ndarray,
        latentvec: numpy.ndarray, 
        **kwargs: dict
    ) -> numpy.ndarray:
    return latentvec[0:1]

def obj_trans_max1(
        decnvec: numpy.ndarray,
        latentvec: numpy.ndarray, 
        **kwargs: dict
    ) -> numpy.ndarray:
    return latentvec[1:2]

# define an inequality constraint violation function
def ineqcv_trans(
        decnvec: numpy.ndarray,
        latentvec: numpy.ndarray, 
        val: float,
        **kwargs: dict
    ) -> numpy.ndarray:
    """
    A custom inequality constraint violation function.

    Parameters
    ----------
    decnvec : numpy.ndarray
        A decision space vector of shape (ndecn,)
    latentvec : numpy.ndarray
        A latent space function vector of shape (ntrait,)
    minvec : numpy.ndarray
        Vector of minimum values for which the latent vector can take.
    
    Returns
    -------
    out : numpy.ndarray
        An inequality constraint violation vector of shape (ntrait,).
    """
    # negate second trait
    negval = -latentvec[1]
    outval = 0.0 if negval > val else val - negval
    out = numpy.array([outval], float)
    return out

##
## Construct optimization algorithms
## =================================

# create custom multi-objective algorithm for optimization
# use NSGA-II and evolve for 1000 generations
moalgo = NSGA2SubsetGeneticAlgorithm(
    ngen = 500,    # number of generations to evolve
    pop_size = 100  # number of parents in population
)

# create a hillclimber algorithm
hillclimber = SteepestDescentSubsetHillClimber()

##
## Constructing a GEBV Subset Selection Object
## ===========================================

# construct a subset selection object for GEBV selection
selprot_unconstrained = GenomicEstimatedBreedingValueSubsetSelection(
    ntrait = 2,         # number of traits to expect from GEBV matrix
    unscale = True,     # unscale GEBVs so that units are in unscaled trait values
    ncross = 10,        # number of breeding crosses to select
    nparent = 2,        # number of parents per breeding cross to select
    nmating = 1,        # number of times parents are mated per cross
    nprogeny = 40,      # number of progenies to derive from a mating event
    nobj = 2,           # number of objectives (ntrait)
    moalgo = moalgo,    # custom multi-objective algorithm
)

selprot_max0 = GenomicEstimatedBreedingValueSubsetSelection(
    ntrait = 2,         # number of traits to expect from GEBV matrix
    unscale = True,     # unscale GEBVs so that units are in unscaled trait values
    ncross = 10,        # number of breeding crosses to select
    nparent = 2,        # number of parents per breeding cross to select
    nmating = 1,        # number of times parents are mated per cross
    nprogeny = 40,      # number of progenies to derive from a mating event
    nobj = 1, # one since due to transformation function
    obj_trans = obj_trans_max0,
    soalgo = hillclimber,
)

selprot_max1 = GenomicEstimatedBreedingValueSubsetSelection(
    ntrait = 2,         # number of traits to expect from GEBV matrix
    unscale = True,     # unscale GEBVs so that units are in unscaled trait values
    ncross = 10,        # number of breeding crosses to select
    nparent = 2,        # number of parents per breeding cross to select
    nmating = 1,        # number of times parents are mated per cross
    nprogeny = 40,      # number of progenies to derive from a mating event
    nobj = 1, # one since due to transformation function
    obj_trans = obj_trans_max1,
    soalgo = hillclimber,
)

# create the constrained selection protocol for 10 two-way crosses
selprot_constrained1 = GenomicEstimatedBreedingValueSubsetSelection(
    ntrait = 2,
    unscale = True,     # unscale GEBVs so that units are in unscaled trait values
    ncross = 10, # ten crosses total
    nparent = 2, # two-way
    nmating = 1,
    nprogeny = 40,
    nobj = 1, # one since due to transformation function
    obj_trans = obj_trans,
    nineqcv = 1,
    ineqcv_trans = ineqcv_trans,
    ineqcv_trans_kwargs = {"val": 50},
    soalgo = hillclimber,
)
selprot_constrained2 = GenomicEstimatedBreedingValueSubsetSelection(
    ntrait = 2,
    unscale = True,     # unscale GEBVs so that units are in unscaled trait values
    ncross = 10, # ten crosses total
    nparent = 2, # two-way
    nmating = 1,
    nprogeny = 40,
    nobj = 1, # one since due to transformation function
    obj_trans = obj_trans,
    nineqcv = 1,
    ineqcv_trans = ineqcv_trans,
    ineqcv_trans_kwargs = {"val": 40},
    soalgo = hillclimber,
)
selprot_constrained3 = GenomicEstimatedBreedingValueSubsetSelection(
    ntrait = 2,
    unscale = True,     # unscale GEBVs so that units are in unscaled trait values
    ncross = 10, # ten crosses total
    nparent = 2, # two-way
    nmating = 1,
    nprogeny = 40,
    nobj = 1, # one since due to transformation function
    obj_trans = obj_trans,
    nineqcv = 1,
    ineqcv_trans = ineqcv_trans,
    ineqcv_trans_kwargs = {"val": 30},
    soalgo = hillclimber,
)
selprot_constrained4 = GenomicEstimatedBreedingValueSubsetSelection(
    ntrait = 2,
    unscale = True,     # unscale GEBVs so that units are in unscaled trait values
    ncross = 10, # ten crosses total
    nparent = 2, # two-way
    nmating = 1,
    nprogeny = 40,
    nobj = 1, # one since due to transformation function
    obj_trans = obj_trans,
    nineqcv = 1,
    ineqcv_trans = ineqcv_trans,
    ineqcv_trans_kwargs = {"val": 20},
    soalgo = hillclimber,
)

##
## Estimating the Pareto Frontier
## ==============================

selprob_unconstrained = selprot_unconstrained.problem(
    pgmat = None,       # argument not utilized
    gmat = gmat,        # ``gmat`` argument required
    ptdf = None,        # argument not utilized
    bvmat = None,       # argument not utilized
    gpmod = algmod,     # ``gpmod`` argument required
    t_cur = 0,          # argument not utilized
    t_max = 0,          # argument not utilized
)

# estimate pareto frontier using optimization algorithm
print("Estimating Pareto Frontier ... ", end = "")
selsoln_unconstrained = selprot_unconstrained.mosolve(
    pgmat = None,       # argument not utilized
    gmat = gmat,        # ``gmat`` argument required
    ptdf = None,        # argument not utilized
    bvmat = None,       # argument not utilized
    gpmod = algmod,     # ``gpmod`` argument required
    t_cur = 0,          # argument not utilized
    t_max = 0,          # argument not utilized
)
print("Done")

# use epsilon constraint method to estimate single solution point
print("Validating Pareto Frontier using Hillclimber ... ", end = "")
selsoln_max0 = selprot_max0.sosolve(
    pgmat = None,       # argument not utilized
    gmat = gmat,        # ``gmat`` argument required
    ptdf = None,        # argument not utilized
    bvmat = None,       # argument not utilized
    gpmod = algmod,     # ``gpmod`` argument required
    t_cur = 0,          # argument not utilized
    t_max = 0,          # argument not utilized
)
selsoln_max1 = selprot_max1.sosolve(
    pgmat = None,       # argument not utilized
    gmat = gmat,        # ``gmat`` argument required
    ptdf = None,        # argument not utilized
    bvmat = None,       # argument not utilized
    gpmod = algmod,     # ``gpmod`` argument required
    t_cur = 0,          # argument not utilized
    t_max = 0,          # argument not utilized
)
selsoln_constrained1 = selprot_constrained1.sosolve(
    pgmat = None,       # argument not utilized
    gmat = gmat,        # ``gmat`` argument required
    ptdf = None,        # argument not utilized
    bvmat = None,       # argument not utilized
    gpmod = algmod,     # ``gpmod`` argument required
    t_cur = 0,          # argument not utilized
    t_max = 0,          # argument not utilized
)
selsoln_constrained2 = selprot_constrained2.sosolve(
    pgmat = None,       # argument not utilized
    gmat = gmat,        # ``gmat`` argument required
    ptdf = None,        # argument not utilized
    bvmat = None,       # argument not utilized
    gpmod = algmod,     # ``gpmod`` argument required
    t_cur = 0,          # argument not utilized
    t_max = 0,          # argument not utilized
)
selsoln_constrained3 = selprot_constrained3.sosolve(
    pgmat = None,       # argument not utilized
    gmat = gmat,        # ``gmat`` argument required
    ptdf = None,        # argument not utilized
    bvmat = None,       # argument not utilized
    gpmod = algmod,     # ``gpmod`` argument required
    t_cur = 0,          # argument not utilized
    t_max = 0,          # argument not utilized
)
selsoln_constrained4 = selprot_constrained4.sosolve(
    pgmat = None,       # argument not utilized
    gmat = gmat,        # ``gmat`` argument required
    ptdf = None,        # argument not utilized
    bvmat = None,       # argument not utilized
    gpmod = algmod,     # ``gpmod`` argument required
    t_cur = 0,          # argument not utilized
    t_max = 0,          # argument not utilized
)
print("Done")

# evaluate constrained selection
obj1, ineq1, eq1 = selprob_unconstrained.evalfn(selsoln_constrained1.soln_decn[0])
obj2, ineq2, eq2 = selprob_unconstrained.evalfn(selsoln_constrained2.soln_decn[0])
obj3, ineq3, eq3 = selprob_unconstrained.evalfn(selsoln_constrained3.soln_decn[0])
obj4, ineq4, eq4 = selprob_unconstrained.evalfn(selsoln_constrained4.soln_decn[0])
obj5, ineq5, eq5 = selprob_unconstrained.evalfn(selsoln_max0.soln_decn[0])
obj6, ineq6, eq6 = selprob_unconstrained.evalfn(selsoln_max1.soln_decn[0])

##
## Visualizing the Pareto Frontier with ``matplotlib``
## ===================================================

# get the pareto frontier
# negate the objectives to get the mean GEBV since optimization problems are always minimizing
xdata_unconstrained = -selsoln_unconstrained.soln_obj[:,0]
ydata_unconstrained = -selsoln_unconstrained.soln_obj[:,1]
xdata_constrained = -numpy.array([obj1[0],obj2[0],obj3[0],obj4[0],obj5[0],obj6[0]])
ydata_constrained = -numpy.array([obj1[1],obj2[1],obj3[1],obj4[1],obj5[1],obj6[1]])

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.scatter(xdata_unconstrained, ydata_unconstrained, label = "NSGA-II")
ax.scatter(xdata_constrained, ydata_constrained, label = "Hillclimber")
ax.set_title("Bi-Objective GEBV Selection Pareto Frontier")
ax.set_xlabel("Synthetic Trait 1 Mean GEBV")
ax.set_ylabel("Synthetic Trait 2 Mean GEBV")
pyplot.savefig("biobjective_GEBV_pareto_frontier.png", dpi = 250)
pyplot.close(fig)
