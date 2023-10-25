#!/usr/bin/env python3

import numpy
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.breed.prot.pt.G_E_Phenotyping import G_E_Phenotyping

###
### Loading Class Modules
### =====================

# import the BreedingValueProtocol class (an abstract interface class)
from pybrops.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol

# import the TrueBreedingValue class (a concrete implemented class)
from pybrops.breed.prot.bv.TrueBreedingValue import TrueBreedingValue

# import the MeanPhenotypicBreedingValue class (a concrete implemented class)
from pybrops.breed.prot.bv.MeanPhenotypicBreedingValue import MeanPhenotypicBreedingValue

###
### Creating Breeding Value Estimation Protocols
### ============================================

bvprot = MeanPhenotypicBreedingValue(
    taxa_col = "taxa",
    taxa_grp_col = "taxa_grp",
    trait_cols = ["Trait01","Trait02"]
)

###
### Estimating Breeding Values
### ==========================

#
# Creating a true genomic model
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
    hyperparams = None
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
taxa = numpy.array(["Taxon"+str(i+1).zfill(3) for i in range(ntaxa)], dtype = object)

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

#
# Creating a phenotyping object
#

# phenotyping parameters
nenv = 3    # number of environments
nrep = 2    # number of replicates within each environment

# construct phenotyping object
ptprot = G_E_Phenotyping(
    gpmod = algmod,
    nenv = nenv,
    nrep = nrep
)

# set the narrow sense heritability
ptprot.set_h2(
    h2 = numpy.array([0.4, 0.7]),
    pgmat = pgmat
)

#
# Creating phenotypes for mean estimation
#

# phenotype individuals
pheno_df = ptprot.phenotype(pgmat)

#
# Calculating the mean values
#

# without a reference genotype matrix for alignment
bvmat1 = bvprot.estimate(pheno_df)

# with a reference genotype matrix for alignment
bvmat2 = bvprot.estimate(pheno_df, pgmat)
