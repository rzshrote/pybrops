#!/usr/bin/env python3

import os
import numpy
import pandas
from matplotlib import pyplot
from matplotlib import animation

import pybrops
from pybrops.breed.prot.sel.OptimalContributionSelection import OptimalContributionSubsetSelection
from pybrops.opt.algo.NSGA2SubsetGeneticAlgorithm import NSGA2SubsetGeneticAlgorithm
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.gmap.StandardGeneticMap import StandardGeneticMap
from pybrops.popgen.cmat.fcty.DenseMolecularCoancestryMatrixFactory import DenseMolecularCoancestryMatrixFactory

# seed python random and numpy random
pybrops.core.random.prng.seed(31621463)

##################### Read genetic map #####################

# create a standard genetic map object by reading from a csv-type file
gmap = StandardGeneticMap.from_csv(
    "McMullen_2009_US_NAM.gmap", # genetic map file name
    vrnt_chrgrp_col   = "chr",   # chromosome/linkage group column name
    vrnt_phypos_col   = "pos",   # chromosome physical position column name
    vrnt_genpos_col   = "cM",    # chromosome genetic position column
    vrnt_genpos_units = "cM",    # chromosome genetic position column units
    auto_group        = True,    # automatically sort and group markers within chromosome groups
    auto_build_spline = True,    # automatically build interpolation spline
    sep               = '\t',    # field delimiter for csv file format
)

############### Create genetic map function ################

# create a Haldane mapping function for a basic scenario
gmapfn = HaldaneMapFunction()

################### Load genotypic data ####################

# read unphased genetic markers from a vcf file
gmat = DenseGenotypeMatrix.from_vcf(
    "widiv_2000SNPs.vcf.gz", # file name to load
    auto_group_vrnt = True,  # automatically sort and group variants
)

# interpolate crossover probabilies for each locus using the genetic map and 
# genetic map function created above
gmat.interp_xoprob(gmap = gmap, gmapfn = gmapfn)

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

#### Calculate breeding value matrix from genomic model ####

# calculate the GEBVs from the genotype matrix
bvmat = algmod.gebv(gtobj = gmat)

################### Construct OCS object ###################

# create coancestry matrix factory object for creating
# identity by state coancestry matrices
ibscmatfcty = DenseMolecularCoancestryMatrixFactory()

# create custom multi-objective algorithm for optimization
algo = NSGA2SubsetGeneticAlgorithm(
    ngen = 1000,    # number of generations to evolve
    pop_size = 100  # number of parents in population
)
# construct a subset selection problem for OCS
selprot = OptimalContributionSubsetSelection(
    ntrait = 2,             # number of expected traits
    cmatfcty = ibscmatfcty, # identity by state
    unscale = True,         # whether to unscale breeding values
    ncross = 10,            # number of breeding crosses to select
    nparent = 2,            # number of parents per breeding cross to select
    nmating = 1,            # number of times parents are mated per cross
    nprogeny = 40,          # number of progenies to derive from a mating event
    nobj = 3,               # number of objectives (1+ntrait)
    moalgo = algo,          # custom multi-objective algorithm
    # leave all other arguments as their default values
)

# estimate pareto frontier
selsoln = selprot.mosolve(
    pgmat = None,
    gmat = gmat,
    ptdf = None,
    bvmat = bvmat,
    gpmod = algmod,
    t_cur = 0,
    t_max = 0
)

# get the pareto frontier
frontier = selsoln.soln_obj

# get axis data
xdata = frontier[:,0]
ydata = frontier[:,1]
zdata = frontier[:,2]

# establish axis labels
xlabel = "Inbreeding"
ylabel = "-Syn 1 Trait"
zlabel = "-Syn 2 Trait"

# create static figure
fig = pyplot.figure()
ax = pyplot.axes(projection = '3d')
ax.scatter3D(xdata, ydata, zdata)
ax.set_title("Multi-Objective Optimal Contribution Selection Pareto Frontier")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_zlabel(zlabel)
pyplot.savefig("ocs_3d_frontier.png", dpi = 250)
pyplot.close(fig)

# create animation frames output directory
outdir = "frames"
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# create animation frames
for i in range(100):
    fig = pyplot.figure()
    ax = pyplot.axes(projection = '3d')
    ax.scatter3D(xdata, ydata, zdata)
    ax.set_title("Multi-Objective Optimal Contribution Selection Pareto Frontier")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.view_init(elev = 30., azim = 3.6 * i)
    s = outdir + "/" + "ocs_3d_frontier_" + str(i).zfill(3) + ".png"
    pyplot.savefig(s, dpi = 250)
    pyplot.close(fig)
