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
    ntrait=None,
    ncross=None,
    nparent=None,
    nmating=None,
    nprogeny=None,
    nobj=None
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
mat = numpy.random.randint(0, nphase+1, size = (ntaxa,nvrnt)).astype("int8")

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
pgmat = DenseGenotypeMatrix(
    mat = mat,
    taxa = taxa,
    taxa_grp = taxa_grp, 
    vrnt_chrgrp = vrnt_chrgrp,
    vrnt_phypos = vrnt_phypos, 
    vrnt_name = vrnt_name, 
    ploidy = nphase
)


selprot.problem()
selprot.sosolve()
selprot.mosolve()
selprot.select()
