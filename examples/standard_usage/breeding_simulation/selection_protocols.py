#!/usr/bin/env python3

import numpy
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix


###
### Loading Class Modules
### =====================

###
### abstract interface classes
### --------------------------

# interface for all selection protocols
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol

# interfaces for all individual selection protocols
from pybrops.breed.prot.sel.BinarySelectionProtocol import BinarySelectionProtocol
from pybrops.breed.prot.sel.IntegerSelectionProtocol import IntegerSelectionProtocol
from pybrops.breed.prot.sel.RealSelectionProtocol import RealSelectionProtocol
from pybrops.breed.prot.sel.SubsetSelectionProtocol import SubsetSelectionProtocol

# interfaces for all mate selection protocols
from pybrops.breed.prot.sel.BinaryMateSelectionProtocol import BinaryMateSelectionProtocol
from pybrops.breed.prot.sel.IntegerMateSelectionProtocol import IntegerMateSelectionProtocol
from pybrops.breed.prot.sel.RealMateSelectionProtocol import RealMateSelectionProtocol
from pybrops.breed.prot.sel.SubsetMateSelectionProtocol import SubsetMateSelectionProtocol

###
### concrete implementation classes
### -------------------------------

# EBV selection
from pybrops.breed.prot.sel.EstimatedBreedingValueSelection import EstimatedBreedingValueBinarySelection
from pybrops.breed.prot.sel.EstimatedBreedingValueSelection import EstimatedBreedingValueIntegerSelection
from pybrops.breed.prot.sel.EstimatedBreedingValueSelection import EstimatedBreedingValueRealSelection
from pybrops.breed.prot.sel.EstimatedBreedingValueSelection import EstimatedBreedingValueSubsetSelection

# EMBV selection
from pybrops.breed.prot.sel.ExpectedMaximumBreedingValueSelection import ExpectedMaximumBreedingValueBinarySelection
from pybrops.breed.prot.sel.ExpectedMaximumBreedingValueSelection import ExpectedMaximumBreedingValueIntegerSelection
from pybrops.breed.prot.sel.ExpectedMaximumBreedingValueSelection import ExpectedMaximumBreedingValueRealSelection
from pybrops.breed.prot.sel.ExpectedMaximumBreedingValueSelection import ExpectedMaximumBreedingValueSubsetSelection

# within family EBV selection
from pybrops.breed.prot.sel.FamilyEstimatedBreedingValueSelection import FamilyEstimatedBreedingValueBinarySelection
from pybrops.breed.prot.sel.FamilyEstimatedBreedingValueSelection import FamilyEstimatedBreedingValueIntegerSelection
from pybrops.breed.prot.sel.FamilyEstimatedBreedingValueSelection import FamilyEstimatedBreedingValueRealSelection
from pybrops.breed.prot.sel.FamilyEstimatedBreedingValueSelection import FamilyEstimatedBreedingValueSubsetSelection

# GEBV selection
from pybrops.breed.prot.sel.GenomicEstimatedBreedingValueSelection import GenomicEstimatedBreedingValueBinarySelection
from pybrops.breed.prot.sel.GenomicEstimatedBreedingValueSelection import GenomicEstimatedBreedingValueIntegerSelection
from pybrops.breed.prot.sel.GenomicEstimatedBreedingValueSelection import GenomicEstimatedBreedingValueRealSelection
from pybrops.breed.prot.sel.GenomicEstimatedBreedingValueSelection import GenomicEstimatedBreedingValueSubsetSelection

# GB selection
from pybrops.breed.prot.sel.GenotypeBuilderSelection import GenotypeBuilderSubsetSelection

# optimal contribution selection
from pybrops.breed.prot.sel.OptimalContributionSelection import OptimalContributionBinarySelection
from pybrops.breed.prot.sel.OptimalContributionSelection import OptimalContributionIntegerSelection
from pybrops.breed.prot.sel.OptimalContributionSelection import OptimalContributionRealSelection
from pybrops.breed.prot.sel.OptimalContributionSelection import OptimalContributionSubsetSelection

# OHV selection
from pybrops.breed.prot.sel.OptimalHaploidValueSelection import OptimalHaploidValueBinarySelection
from pybrops.breed.prot.sel.OptimalHaploidValueSelection import OptimalHaploidValueIntegerSelection
from pybrops.breed.prot.sel.OptimalHaploidValueSelection import OptimalHaploidValueRealSelection
from pybrops.breed.prot.sel.OptimalHaploidValueSelection import OptimalHaploidValueSubsetSelection

# OPV selection
from pybrops.breed.prot.sel.OptimalPopulationValueSelection import OptimalPopulationValueSubsetSelection

# random selection
from pybrops.breed.prot.sel.RandomSelection import RandomBinarySelection
from pybrops.breed.prot.sel.RandomSelection import RandomIntegerSelection
from pybrops.breed.prot.sel.RandomSelection import RandomRealSelection
from pybrops.breed.prot.sel.RandomSelection import RandomSubsetSelection

# UC selection
from pybrops.breed.prot.sel.UsefulnessCriterionSelection import UsefulnessCriterionBinarySelection
from pybrops.breed.prot.sel.UsefulnessCriterionSelection import UsefulnessCriterionIntegerSelection
from pybrops.breed.prot.sel.UsefulnessCriterionSelection import UsefulnessCriterionRealSelection
from pybrops.breed.prot.sel.UsefulnessCriterionSelection import UsefulnessCriterionSubsetSelection

# weighted GEBV selection
from pybrops.breed.prot.sel.WeightedGenomicSelection import WeightedGenomicBinarySelection
from pybrops.breed.prot.sel.WeightedGenomicSelection import WeightedGenomicIntegerSelection
from pybrops.breed.prot.sel.WeightedGenomicSelection import WeightedGenomicRealSelection
from pybrops.breed.prot.sel.WeightedGenomicSelection import WeightedGenomicSubsetSelection

###
### Creating Selection Protocol Classes
### ===================================

#
# Construction from NumPy arrays
# ------------------------------

# create standard GEBV selection protocol in subset decision space
selprot = GenomicEstimatedBreedingValueSubsetSelection(
    ntrait = 2,
    ncross = 10,
    nparent = 2,
    nmating = 1,
    nprogeny = 40,
    nobj = 2,
)

#
# Creating a genomic model
#

# model parameters
nfixed = 1      # number of fixed effects
ntrait = 2      # number of traits
nmisc = 0       # number of miscellaneous random effects
nadditive = 50  # number of additive marker effects

# create dummy values
beta = numpy.random.random((nfixed,ntrait))
u_misc = numpy.random.random((nmisc,ntrait))
u_a = numpy.random.random((nadditive,ntrait))
trait = numpy.array(["Trait"+str(i+1).zfill(2) for i in range(ntrait)], dtype = object)

# create additive linear genomic model
algmod = DenseAdditiveLinearGenomicModel(
    beta = beta,
    u_misc = u_misc,
    u_a = u_a,
    trait = trait,
    model_name = "example",
    params = None
)

#
# Construct random genomes
#

# shape parameters for random genomes
ntaxa = 100
nvrnt = nadditive
ngroup = 20
nchrom = 10
nphase = 2

# create random genotypes
mat = numpy.random.randint(0, 2, size = (nphase,ntaxa,nvrnt)).astype("int8")

# create taxa names
taxa = numpy.array(["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], dtype = object)

# create taxa groups
taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
taxa_grp.sort()

# create marker variant chromsome assignments
vrnt_chrgrp = numpy.random.randint(1, nchrom+1, nvrnt)
vrnt_chrgrp.sort()

# create marker physical positions
vrnt_phypos = numpy.random.choice(1000000, size = nvrnt, replace = False)
vrnt_phypos.sort()

# create marker variant names
vrnt_name = numpy.array(["SNP"+str(i+1).zfill(4) for i in range(nvrnt)], dtype = object)

# create a phased genotype matrix from scratch using NumPy arrays
pgmat = DensePhasedGenotypeMatrix(
    mat = mat,
    taxa = taxa,
    taxa_grp = taxa_grp, 
    vrnt_chrgrp = vrnt_chrgrp,
    vrnt_phypos = vrnt_phypos, 
    vrnt_name = vrnt_name, 
    ploidy = nphase
)

# create a genotype matrix from scratch using NumPy arrays
gmat = DenseGenotypeMatrix(
    mat = mat.sum(0, dtype="int8"),
    taxa = taxa,
    taxa_grp = taxa_grp, 
    vrnt_chrgrp = vrnt_chrgrp,
    vrnt_phypos = vrnt_phypos, 
    vrnt_name = vrnt_name, 
    ploidy = nphase
)

###
### Generating Selection Problems for Optimization
### ==============================================

##
## Setup for unconstrained optimization
## ------------------------------------

# nothing special to do!!!

# create the selection protocol for 10 two-way crosses
selprot_unconstrained = GenomicEstimatedBreedingValueSubsetSelection(
    ntrait = 2,
    ncross = 10, # ten crosses total
    nparent = 2, # two-way
    nmating = 1,
    nprogeny = 40,
    nobj = 2,
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
        maskvec: numpy.ndarray,
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
    maskvec : numpy.ndarray
        A mask vector of shape (ntrait,)
    
    Returns
    -------
    out : numpy.ndarray
        A vector of shape (sum(maskvec),).
    """
    # extract trait(s) as objective(s)
    return latentvec[maskvec]

# define an inequality constraint violation function
def ineqcv_trans(
        decnvec: numpy.ndarray,
        latentvec: numpy.ndarray, 
        minvec: numpy.ndarray,
        **kwargs: dict
    ) -> numpy.ndarray:
    """
    A custom inequality constraint violation function.

    Parameters
    ----------
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
    # calculate constraint violations
    out = minvec - latentvec
    # where constraint violation is negative (no constraint violation), set to zero
    out[out < 0] = 0
    # return inequality constraint violation array
    return out

# for constrained selection, make keyword arguments for our custom 
# transformation functions
obj_trans_kwargs = {
    # we want to select the first trait
    "maskvec": numpy.array([True, False], dtype=bool)
}
ineqcv_trans_kwargs = {
    # we don't care about the first trait's minimum value (negated maximum), so set to -Inf
    # we do care about the second trait's minimum value (negated maximum), so set to -1.0
    "minvec": numpy.array([-numpy.inf, -1.0], dtype=float)
}

# create the constrained selection protocol for 10 two-way crosses
selprot_constrained = GenomicEstimatedBreedingValueSubsetSelection(
    ntrait = 2,
    ncross = 10, # ten crosses total
    nparent = 2, # two-way
    nmating = 1,
    nprogeny = 40,
    nobj = 1, # one since sum(maskvec) == 1
    obj_trans = obj_trans,
    obj_trans_kwargs = obj_trans_kwargs,
    nineqcv = 2,
    ineqcv_trans = ineqcv_trans,
    ineqcv_trans_kwargs = ineqcv_trans_kwargs
)

##
## Generating the selection problem
## --------------------------------

# generate an unconstrained GEBV subset selection problem
# for this selection protocol type, we only need the genotype matrix and a genomic prediction model
prob_unconstrained = selprot_unconstrained.problem(
    pgmat = None,
    gmat = gmat,
    ptdf = None,
    bvmat = None,
    gpmod = algmod,
    t_cur = None,
    t_max = None,
)

# generate a constrained GEBV subset selection problem
# for this selection protocol type, we only need the genotype matrix and a genomic prediction model
prob_constrained = selprot_constrained.problem(
    pgmat = None,
    gmat = gmat,
    ptdf = None,
    bvmat = None,
    gpmod = algmod,
    t_cur = None,
    t_max = None,
)

# generate a random solution to test
soln = numpy.random.choice(prob_constrained.decn_space, prob_constrained.ndecn)

# evaluate the solution in the unconstrained problem
eval_unconstrained = prob_unconstrained.evalfn(soln)

# evaluate the solution in the constrained problem
eval_constrained = prob_constrained.evalfn(soln)

###
### Single-Objective Optimization
### =============================

# perform single-objective optimization using 
soln_constrained = selprot_constrained.sosolve(
    pgmat = None,
    gmat = gmat,
    ptdf = None,
    bvmat = None,
    gpmod = algmod,
    t_cur = None,
    t_max = None,
)

# examine the solution decision vector(s)
soln_constrained.soln_decn

# examine the solution objective function vector(s)
soln_constrained.soln_obj

# examine the solution inequality constraint violation vector(s)
soln_constrained.soln_ineqcv

# examine the soltuion equality constraint violation vector(s)
soln_constrained.soln_eqcv

###
### Multi-Objective Optimization
### ============================

# perform single-objective optimization using 
soln_unconstrained = selprot_unconstrained.mosolve(
    pgmat = None,
    gmat = gmat,
    ptdf = None,
    bvmat = None,
    gpmod = algmod,
    t_cur = None,
    t_max = None,
)

# examine the solution decision vector(s)
soln_unconstrained.soln_decn

# examine the solution objective function vector(s)
soln_unconstrained.soln_obj

# examine the solution inequality constraint violation vector(s)
soln_unconstrained.soln_ineqcv

# examine the soltuion equality constraint violation vector(s)
soln_unconstrained.soln_eqcv

###
### Selection
### =========

selcfg_constrained = selprot_constrained.select(
    pgmat = pgmat,
    gmat = gmat,
    ptdf = None,
    bvmat = None,
    gpmod = algmod,
    t_cur = None,
    t_max = None,
)
