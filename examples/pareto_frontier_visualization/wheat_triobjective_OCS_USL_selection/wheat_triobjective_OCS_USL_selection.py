#!/usr/bin/env python3

###
### Custom Tri-Objective OCS/USL Selection Wheat
### ############################################

##
## Loading Required Modules and Seeding the global PRNG
## ====================================================
import os
from typing import Tuple
import numpy
import pandas
from matplotlib import pyplot
from PIL import Image
from matplotlib import rcParams

from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import DenseMolecularCoancestryMatrix
rcParams['font.family'] = 'Liberation Serif' # set default font

import pybrops
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.opt.algo.NSGA3SubsetGeneticAlgorithm import NSGA3SubsetGeneticAlgorithm
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.opt.prob.SubsetProblem import SubsetProblem

# seed python random and numpy random
pybrops.core.random.prng.seed(78582734)

##
## Reading External Data
## =====================

#
# Loading Genotypic Data from a CSV File
# --------------------------------------

# read marker matrix from file
wheat_markers_df = pandas.read_csv("wheat_markers.csv.gz", index_col=0)

# construct a genotype matrix from read inputs
gmat = DenseGenotypeMatrix(
    mat = wheat_markers_df.to_numpy('int8'),
    taxa = wheat_markers_df.index.to_numpy(object),
    vrnt_name = wheat_markers_df.columns.to_numpy(object),
    ploidy = 2,
)

#
# Loading Genomic Model Data from CSV Files
# -----------------------------------------

# read intercepts from file
wheat_intercepts_df = pandas.read_csv("wheat_intercepts.csv.gz", index_col=0)

# read marker effects from file
wheat_marker_effects_df = pandas.read_csv("wheat_marker_effects.csv.gz", index_col=0)

# construct an additive linear genomic model from read inputs
# use only the first two traits to create the genomic model
algmod = DenseAdditiveLinearGenomicModel(
    beta = wheat_intercepts_df.to_numpy(float)[:,0:1],
    u_misc = None,
    u_a = wheat_marker_effects_df.to_numpy(float)[:,0:1],
    model_name = "uni-trait",
)

##
## Defining a Custom Problem Class for OCS/USL Selection
## =====================================================

# define class
class OCSUSLSubsetProblem(SubsetProblem):
    ### class constructor ###
    def __init__(
            self,
            nselindiv: int,
            ploidy: int,
            gebvvec: numpy.ndarray,
            C: numpy.ndarray,
            haplomat: numpy.ndarray,
            **kwargs: dict
        ) -> None:
        """
        Constructor for OCSUSLSubsetProblem.
        
        Parameters
        ----------
        nselindiv : int
            Number of individuals to select. Must be less than the number of taxa.
        ploidy : int
            Ploidy level of the organism.
        gebvvec : numpy.ndarray
            Breeding value vector of shape (n,) for a single trait.
        C : numpy.ndarray
            Cholesky decomposition (upper triangle matrix) of a kinship matrix of shape (n,n).
        haplomat : numpy.ndarray
            A haplotype effect matrix of shape ``(m,n,h)``.

            Where:

            - ``m`` is the number of chromosome phases (2 for diploid, etc.).
            - ``n`` is the number of individuals.
            - ``h`` is the number of haplotype blocks.
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        self.nselindiv = nselindiv
        self.ploidy = ploidy
        self.gebvvec = gebvvec
        self.C = C
        self.haplomat = haplomat
        ntaxa = len(gebvvec)
        decn_space = numpy.arange(ntaxa)
        decn_space_lower = numpy.repeat(0, self.nselindiv)
        decn_space_upper = numpy.repeat(ntaxa-1, self.nselindiv)
        super().__init__(
            ndecn = nselindiv, 
            decn_space = decn_space, 
            decn_space_lower = decn_space_lower, 
            decn_space_upper = decn_space_upper, 
            nobj = 3, 
            **kwargs
        )
    ### method required by PyBrOpS interface ###
    def evalfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]:
        """
        Evaluate a candidate solution for the given Problem.

        This calculates three vectors which are to be minimized:

        .. math::

            \\mathbf{v_{obj}} = \\mathbf{w_{obj} \\odot F_{obj}(x)} \\
            \\mathbf{v_{ineqcv}} = \\mathbf{w_{ineqcv} \\odot G_{ineqcv}(x)} \\
            \\mathbf{v_{eqcv}} = \\mathbf{w_{eqcv} \\odot H_{eqcv}(x)}
        
        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,)``.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : tuple
            A tuple ``(obj, ineqcv, eqcv)``.
            
            Where:
            
            - ``obj`` is a numpy.ndarray of shape ``(nobj,)`` that contains 
                objective function evaluations.
            - ``ineqcv`` is a numpy.ndarray of shape ``(nineqcv,)`` that contains 
                inequality constraint violation values.
            - ``eqcv`` is a numpy.ndarray of shape ``(neqcv,)`` that contains 
                equality constraint violation values.
        """
        # calculate individual contribution
        # scalar
        indcontrib = 1.0 / len(x)
        # calculate negated genetic gain (minimizing objective)
        # (n,)[(k,)] -> (k,)
        # (k,).sum(0, keepdims=True) -> (1,)
        # scalar * (1,) -> (1,)
        gain = -indcontrib * self.gebvvec[x].sum(0, keepdims = True)
        # calculate negated upper selection limit from haplotype values (minimizing objective)
        # (m,n,h)[:,(k,),:] -> (m,k,h)
        # (m,k,h).max((0,1)) -> (h,)
        # (h,).sum(0, keepdims=True) -> (1,)
        # scalar * (1,) -> (1,)
        usl = -self.ploidy * self.haplomat[:,x,:].max((0,1)).sum(0, keepdims = True)
        # calculate mean genomic relationship (minimizing objective)
        # (n,n)[:,(k,)] -> (n,k)
        # scalar * (n,k).sum(1) -> (n,)
        # norm2( (n,), keepdims=True ) -> (1,)
        mgr = numpy.linalg.norm(indcontrib * self.C[:,x].sum(1), ord = 2, keepdims = True)
        # concatenate vectors to get objective function vector
        # (1,) cat (1,) cat (1,) -> (3,)
        obj = self.obj_wt * numpy.concatenate([gain, usl, mgr])
        # calculate inequality constraint violations
        # (0,)
        ineqcv = self.ineqcv_wt * numpy.zeros(self.nineqcv)
        # calculate equality constraint violations
        eqcv = self.eqcv_wt * numpy.zeros(self.neqcv)
        # return (3,), (0,), (0,)
        return obj, ineqcv, eqcv
    ### method required by PyMOO interface ###
    # use default ``_evaluate`` method which uses the ``evalfn`` method

##
## Constructing a Breeding Value Matrix
## ====================================

# calculate the GEBVs from the genotype matrix
bvmat = algmod.gebv(gtobj = gmat)

# extract first trait in unscaled terms
gebvvec = bvmat.unscale()[:,0]

##
## Calculating a Haplotype Value Matrix
##

# (n,p)
M = gmat.mat
# (p,1)
B = algmod.u_a

# take half because we're dealing with homozygotes and we're diploid
# (n,p) * (1,p) -> (n,p)
# scalar * (n,p) -> (n,p)
V = 0.5 * M * B.T

# stack
# (n,p) -> (m,n,p)
haplomat = numpy.stack([V,V])

##
## Calculating an Identity-By-State Kinship Matrix and its Decomposition
## =====================================================================

# calculate identity by state coancestry
cmat = DenseMolecularCoancestryMatrix.from_gmat(gmat)

# extract kinship as numpy.ndarray
K = cmat.mat_asformat("kinship")

# calculate Cholesky decomposition of kinship matrix
# convert to upper triangle matrix
# (n,n).T -> (n,n)
C = numpy.linalg.cholesky(K).T

##
## Constructing an OCS/USL Problem Object
## ======================================

# construct optimization problem
prob = OCSUSLSubsetProblem(
    nselindiv = 20,
    ploidy = 2,
    gebvvec = gebvvec,
    C = C,
    haplomat = haplomat,
)

##
## Constructing a Custom Genetic Algorithm Object
## ==============================================

# create custom multi-objective algorithm for optimization
moea = NSGA3SubsetGeneticAlgorithm(
    ngen = 2000,    # number of generations to evolve
    pop_size = 100, # number of parents in population
    nrefpts = 91,   # number of reference points for optimization in 3d space
)

##
## Estimating the Pareto Frontier
## ==============================

print("optimizing...")
# minimize the optimization problem
soln = moea.minimize(prob)
print("...finished optimization")

##
## Visualizing the Pareto Frontier with ``matplotlib``
## ===================================================

#
# Creating an animation
# ---------------------

print("building animation...")
# set default font size
rcParams['font.size'] = 10

# image base name
basename = "wheat_triobjective_OCS_USL_pareto_frontier"

# create animation frames output directory
outdir = "frames"
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# get axis data
x = -soln.soln_obj[:,0] # negate to get mean genetic gain
y = -soln.soln_obj[:,1] # negate to get upper selection limit
z =  soln.soln_obj[:,2] # get mean genomic relationship

# create animation frames
for i in range(360):
    fig = pyplot.figure()
    ax = pyplot.axes(projection = '3d')
    ax.scatter3D(x, y, z)
    ax.set_title("Optimal Contribution + Upper Selection Limit\nPareto Frontier in Wheat")
    ax.set_xlabel("Yield Mean GEBV")
    ax.set_ylabel("Upper Selection Limit")
    ax.set_zlabel("Mean IBS Kinship")
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
print("...animation built")
