#!/usr/bin/env python3

import os
import numpy
import pandas
from matplotlib import pyplot
from matplotlib import animation

import pybrops
from pybrops.algo.opt.NSGA3UnityConstraintGeneticAlgorithm import NSGA3UnityConstraintGeneticAlgorithm
from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.breed.prot.sel.transfn import trans_ndpt_to_vec_dist
from pybrops.breed.prot.sel.OptimalContributionSelection import OptimalContributionSelection
from pybrops.core.random import global_prng
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import DenseMolecularCoancestryMatrix

# seed python random and numpy random
pybrops.core.random.seed(31621463)

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

################### Calculate genotypes ####################
gtprot = DenseUnphasedGenotyping()                              # genotyping protocols
dgmat = gtprot.genotype(pgmat = dpgmat)                         # get genotypes (add across phases)

################# Construct genomic model ##################
# trait means
mkrmean = numpy.array([0.0, 0.0])

# make two traits that are competing in nature
mkrcov = numpy.array([
    [ 1.0, -0.4],
    [-0.4,  1.0]
])

# make marker effects of shape (nvrnt, 2)
mkreffect = global_prng.multivariate_normal(mkrmean, mkrcov, dpgmat.nvrnt)

dalgmod = DenseAdditiveLinearGenomicModel(  # create model
    beta = numpy.float64([[10.0, 25.0]]),   # model intercepts
    u_misc = None,                          # miscellaneous random effects
    u_a = mkreffect,                        # random marker effects
    trait = numpy.object_(["syn1","syn2"]), # trait names
    model_name = "synthetic_model",         # name of the model
    params = None                           # model parameters
)

################### Construct OCS object ###################
def inbfn(t_cur, t_max):
    return 0.75

# make selection object
sel = OptimalContributionSelection(
    nparent = 20,
    ncross = 1,
    nprogeny = 40,
    inbfn = inbfn,
    cmatcls = DenseMolecularCoancestryMatrix,
    bvtype = "gebv",
    method = "pareto",
    objfn_trans = None, # sum of two traits
    objfn_trans_kwargs = None, # no kwargs
    objfn_wt = numpy.array([1.0, 1.0, 1.0]), # maximizing function
    ndset_trans = trans_ndpt_to_vec_dist,
    ndset_trans_kwargs = {
        "objfn_wt": numpy.array([1.0, 1.0, 1.0]),   # all objectives maximizing
        "wt": numpy.array([1./3., 1./3., 1./3.])    # 1/3 equal weight to all
    },
    moalgo = NSGA3UnityConstraintGeneticAlgorithm(
        ngen = 1500,            # number of generations to evolve
        mu = 100,               # number of parents in population
        lamb = 100,             # number of progeny to produce
        cxeta = 30.0,           # crossover variance parameter
        muteta = 20.0,          # mutation crossover parameter
        refpnts = None,         # hyperplane reference points
        save_logbook = False,   # whether to save logs or not
        rng = global_prng       # PRNG source
    ),
    rng = global_prng
)

# get pareto frontier
frontier, sel_config = sel.pareto(
    pgmat = None,
    gmat = dgmat,
    ptdf = None,
    bvmat = None,
    gpmod = dalgmod,
    t_cur = 0,
    t_max = 20
)

xdata = frontier[:,0]
ydata = frontier[:,1]
zdata = frontier[:,2]

xlabel = "inbreeding"
ylabel = dalgmod.trait[0]
zlabel = dalgmod.trait[1]

# create static figure
fig = pyplot.figure()
ax = pyplot.axes(projection = '3d')
ax.scatter3D(xdata, ydata, zdata)
ax.set_title("Multi-Objective Optimal Contribution Selection Pareto Frontier")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_zlabel(zlabel)
pyplot.savefig("ocs_3d_frontier.png", dpi = 250)

# create animation
fig = pyplot.figure()
ax = pyplot.axes(projection = '3d')

def init():
    ax.scatter3D(xdata, ydata, zdata)
    ax.set_title("Multi-Objective Optimal Contribution Selection Pareto Frontier")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    return fig,

def animate(i):
    ax.view_init(elev = 30., azim = 3.6 * i)
    # ax.view_init(elev = 30., azim = i)
    return fig,

# create and same animation
outdir = "frames"
if not os.path.isdir(outdir):
    os.mkdir(outdir)
init()
for i in range(100):
# for i in range(360):
    animate(i)
    s = outdir + "/" + "ocs_3d_frontier_" + str(i).zfill(3) + ".png"
    pyplot.savefig(s, dpi = 250)
