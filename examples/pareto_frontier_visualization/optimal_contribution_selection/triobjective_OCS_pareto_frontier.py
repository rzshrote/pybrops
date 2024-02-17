#!/usr/bin/env python3

###
### Tri-Objective Optimal Contribution Selection Pareto Frontier Visualization
### ##########################################################################

##
## Loading Required Modules and Seeding the global PRNG
## ====================================================

# import libraries
import os
import numpy
from matplotlib import pyplot
from matplotlib import rcParams
rcParams['font.family'] = 'Liberation Serif' # set default font
from PIL import Image
import pybrops
from pybrops.breed.prot.sel.OptimalContributionSelection import OptimalContributionSubsetSelection
from pybrops.opt.algo.NSGA3SubsetGeneticAlgorithm import NSGA3SubsetGeneticAlgorithm
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.cmat.fcty.DenseMolecularCoancestryMatrixFactory import DenseMolecularCoancestryMatrixFactory

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
    trait   = numpy.array(                   # trait names
                ["syn1","syn2"],
                dtype=object
            ),
    model_name = "synthetic_model",          # name of the model
    hyperparams = None                       # model parameters
)

##
## Constructing a Breeding Value Matrix
## ====================================

# calculate the GEBVs from the genotype matrix
bvmat = algmod.gebv(gtobj = gmat)

##
## Constructing an Optimal Contribution Subset Selection Object
## ============================================================

# create coancestry matrix factory object for creating
# identity by state coancestry matrices
ibscmatfcty = DenseMolecularCoancestryMatrixFactory()

# create custom multi-objective algorithm for optimization
algo = NSGA3SubsetGeneticAlgorithm(
    ngen = 1500,    # number of generations to evolve
    pop_size = 100, # number of parents in population
    nrefpts = 91,   # number of reference points for optimization
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

##
## Estimating the Pareto Frontier
## ==============================

# estimate pareto frontier
selsoln = selprot.mosolve(
    pgmat = None,   # argument not utilized
    gmat = gmat,    # ``gmat`` argument required
    ptdf = None,    # argument not utilized
    bvmat = bvmat,  # ``bvmat`` argument required
    gpmod = None,   # argument not utilized
    t_cur = 0,      # argument not utilized
    t_max = 0,      # argument not utilized
)

##
## Visualizing the Pareto Frontier with ``matplotlib``
## ===================================================

#
# Creating a static image
# -----------------------

# image base name
basename = "triobjective_OCS_pareto_frontier"

# get axis data
x =  selsoln.soln_obj[:,0] # 2 * mean kinship (additive relationship/inbreeding)
y = -selsoln.soln_obj[:,1] # negate to get EBV
z = -selsoln.soln_obj[:,2] # negate to get EBV
z2 = numpy.ones(shape = x.shape) * min(z)

# create static figure
fig = pyplot.figure(figsize=(5,5))
ax = pyplot.axes(projection = "3d")
ax.scatter(x, y, z, color = "blue")
for i,j,k,h in zip(x,y,z,z2):
    a = 0.5 * (i - min(x)) / (max(x) - min(x)) + 0.5
    ax.plot([i,i],[j,j],[k,h], color="cornflowerblue", alpha = a)

ax.set_title("Multi-Objective Optimal Contribution Selection Pareto Frontier")
ax.set_xlabel("Inbreeding")
ax.set_ylabel("Synthetic Trait 1 Mean EBV")
ax.set_zlabel("Synthetic Trait 2 Mean EBV")
ax.view_init(elev = 30., azim = 32)
pyplot.savefig(basename + ".png", dpi = 250)
pyplot.close(fig)

#
# Creating an animation
# ---------------------

# image base name
basename = "triobjective_OCS_pareto_frontier"

# create animation frames output directory
outdir = "frames"
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# create animation frames
for i in range(360):
    fig = pyplot.figure()
    ax = pyplot.axes(projection = '3d')
    ax.scatter3D(x, y, z)
    ax.set_title("Multi-Objective Optimal Contribution Selection Pareto Frontier")
    ax.set_xlabel("Inbreeding")
    ax.set_ylabel("Synthetic Trait 1 Mean EBV")
    ax.set_zlabel("Synthetic Trait 2 Mean EBV")
    ax.view_init(elev = 30., azim = i)
    pyplot.savefig(outdir + "/" + basename + "_" + str(i).zfill(3) + ".png", dpi = 250)
    pyplot.close(fig)

# construct filenames from which to read
filenames = [outdir + "/" + basename + "_" + str(i).zfill(3) + ".png" for i in range(360)]

# read image files from which to create animation using PIL
images = [Image.open(filename) for filename in filenames]

# resize images to 50% size using PIL
images_resize = [img.resize(tuple(px // 2 for px in img.size)) for img in images]

# get first image
img = images_resize[0]

# create gif by appending remaining images to end of first image
img.save(
    basename + ".gif", 
    save_all = True, 
    append_images = images_resize[1:], 
    optimize = True, 
    duration = 55,      # inverse of speed
    loop = 0,           # loop indefinitely
)
