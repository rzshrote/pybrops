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
from pybrops.breed.prot.sel.OptimalContributionSelection import OptimalContributionSubsetSelection
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.opt.algo.NSGA2SubsetGeneticAlgorithm import NSGA2SubsetGeneticAlgorithm
from pybrops.opt.algo.SortingSteepestDescentSubsetHillClimber import SortingSteepestDescentSubsetHillClimber
from pybrops.opt.algo.SteepestDescentSubsetHillClimber import SteepestDescentSubsetHillClimber
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
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
pgmat = mate2waydh.mate(
    pgmat = fndr_pgmat,
    xconfig = xconfig,
    nmating = 1,
    nprogeny = 80,
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
    taxa_grp_col = "taxa_grp",
    trait_cols = "Syn1",
)

##
## Create a Selection Protocol Object
## ==================================

# create a dense molecular coancestry matrix factory
cmatfcty = DenseMolecularCoancestryMatrixFactory()

# create custom multi-objective algorithm for optimization
# use NSGA-II and evolve for 1000 generations
moalgo = NSGA2SubsetGeneticAlgorithm(
    ngen = 2000,    # number of generations to evolve
    pop_size = 100  # number of parents in population
)

# create a selection protocol that selects based on EBVs
selprot = OptimalContributionSubsetSelection(
    ntrait       = 1,            # number of expected traits
    cmatfcty     = cmatfcty,     # coancestry/kinship matrix factory
    unscale      = True,         # unscale breeding values to human-readable format
    ncross       = 20,           # number of cross configurations
    nparent      = 2,            # number of parents per cross configuration
    nmating      = 1,            # number of matings per cross configuration
    nprogeny     = 80,           # number of progeny per mating event
    nobj         = 2,            # number of objectives == 1 + ntrait
    moalgo       = moalgo,
)

#
# Perform Within-Family Selection
# ===============================

def within_family_selection(bvmat: DenseBreedingValueMatrix, nindiv: int) -> numpy.ndarray:
    order = numpy.arange(bvmat.ntaxa)
    value = bvmat.mat[:,0] # get trait breeding values
    indices = []
    groups = numpy.unique(bvmat.taxa_grp)
    for group in groups:
        mask = bvmat.taxa_grp == group
        tmp_order = order[mask]
        tmp_value = value[mask]
        value_argsort = tmp_value.argsort()
        ix = value_argsort[::-1][:nindiv]
        indices.append(tmp_order[ix])
    indices = numpy.concatenate(indices)
    return indices

# initial phenotyping
pheno_df = ptprot.phenotype(pgmat)

# initial breeding value estimation
bvmat = bvprot.estimate(pheno_df, pgmat)

# get candidate mask using within family selection
indices = within_family_selection(bvmat, 8)

# get parental candidates
cand_pgmat = pgmat.select_taxa(indices)
cand_bvmat = bvmat.select_taxa(indices)

# construct an initial coancestry matrix from whole population
cmat = cmatfcty.from_gmat(cand_pgmat)

# extract the mean inbreeding in population
meaninb = cmat.mean("kinship")

mosoln = selprot.mosolve(
    pgmat   = cand_pgmat,
    gmat    = cand_pgmat,
    ptdf    = None,
    bvmat   = cand_bvmat,
    gpmod   = None,
    t_cur   = 0,
    t_max   = 0,
)

##
## Visualizing the Pareto Frontier with ``matplotlib``
## ===================================================

# get the pareto frontier
# negate the objectives to get the mean GEBV since optimization problems are always minimizing
xdata =  mosoln.soln_obj[:,0]
ydata = -mosoln.soln_obj[:,1]

# create static figure
fig = pyplot.figure()
ax = pyplot.axes()
ax.scatter(xdata, ydata)
ax.set_title("Single-Trait OCS Pareto Frontier")
ax.set_xlabel("Inbreeding")
ax.set_ylabel("Synthetic Trait Mean EBV")
pyplot.savefig("single_trait_OCS_pareto_frontier.png", dpi = 300)
pyplot.close(fig)
