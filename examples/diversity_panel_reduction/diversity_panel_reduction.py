#!/usr/bin/env python3

import numpy
from matplotlib import pyplot
from matplotlib import animation

import pybrops
from pybrops.opt.algo.UnconstrainedNSGA2BinaryGeneticAlgorithm import UnconstrainedNSGA2BinaryGeneticAlgorithm
from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import DenseMolecularCoancestryMatrix

# seed python random and numpy random
pybrops.core.random.seed(194711822)

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
X = dpgmat.mat_asformat("{0,1,2}")                              # get genotype matrix

########## Calculate molecular coancestry matrix ###########
dcmat = DenseMolecularCoancestryMatrix.from_gmat(dpgmat)        # calculate coancestry
K = dcmat.mat_asformat("kinship")                               # extract kinship as numpy.ndarray

# define an objective function
def objfn1(x):
    x = numpy.float64(x)
    f1 = x.sum()                # minimizing function
    if f1 == 0.0:
        return numpy.array([numpy.inf, numpy.inf])
    y = x / x.sum()
    f2 = y.dot(K).dot(y) # minimizing function
    out = numpy.array([f1, f2])
    return out

# construct NSGA-II object
moea = UnconstrainedNSGA2BinaryGeneticAlgorithm(
    ngen = 400,
    mu = 400,
    lamb = 400
)

# optimize objectives using NSGA-II
frontier, sel_config, misc = moea.optimize(
    objfn = objfn1,
    k = dcmat.ntaxa,
    sspace = None,
    objfn_wt = [-1.0, -1.0]
)

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.scatter(frontier[:,0], frontier[:,1])
ax.set_title("Pareto frontier for number of selected individuals and molecular diversity")
ax.set_xlabel("Number of individuals")
ax.set_ylabel("Mean molecular coancestry")
pyplot.savefig("pareto_frontier1.png", dpi = 250)

# define objective function
X = dpgmat.mat_asformat("{0,1,2}")                              # get genotype matrix
def objfn2(x):
    mask = numpy.array(x, dtype = "bool")
    x = numpy.array(x, dtype = "float64")
    f1 = x.sum()            # minimizing function
    if f1 == 0.0:
        return numpy.array([numpy.inf, numpy.inf])
    p = X[mask,:].sum(0)
    f2 = (p == 0).sum() + (p == (2*mask.sum())).sum()
    out = numpy.array([f1, f2])
    return out

# construct NSGA-II object
moea = UnconstrainedNSGA2BinaryGeneticAlgorithm(
    ngen = 550,
    mu = 200,
    lamb = 200
)

# optimize objectives using NSGA-II
frontier, sel_config, misc = moea.optimize(
    objfn = objfn2,
    k = dcmat.ntaxa,
    sspace = None,
    objfn_wt = [-1.0, -1.0]
)

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.scatter(frontier[:,0], frontier[:,1])
ax.set_title("Pareto frontier for number of selected individuals and number of fixed alleles")
ax.set_xlabel("Number of individuals")
ax.set_ylabel("Number of alleles fixed")
pyplot.savefig("pareto_frontier2.png", dpi = 250)
