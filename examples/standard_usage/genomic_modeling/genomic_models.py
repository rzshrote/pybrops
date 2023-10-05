#!/usr/bin/env python3

import copy
import numpy
import pandas

# import GenomicModel classes (abstract interface classes)
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.gmod.NonlinearGenomicModel import NonlinearGenomicModel
from pybrops.model.gmod.LinearGenomicModel import LinearGenomicModel
from pybrops.model.gmod.CoancestryLinearGenomicModel import CoancestryLinearGenomicModel
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.AdditiveDominanceLinearGenomicModel import AdditiveDominanceLinearGenomicModel
from pybrops.model.gmod.AdditiveDominanceEpistaticLinearGenomicModel import AdditiveDominanceEpistaticLinearGenomicModel

# import dense genomic models (concrete implementation classes)
from pybrops.model.gmod.DenseLinearGenomicModel import DenseLinearGenomicModel
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

###
### Pre-demonstration setup: Construct dummy genotypes
### ==================================================

# shape parameters for random genotypes
ntaxa = 100
nvrnt = 50
ngroup = 20
nchrom = 10
ploidy = 2

# create random genotypes
mat = numpy.random.randint(0, ploidy+1, size = (ntaxa,nvrnt)).astype("int8")

# create taxa names
taxa = numpy.array(
    ["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], 
    dtype = object
)

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
vrnt_name = numpy.array(
    ["SNP"+str(i+1).zfill(4) for i in range(nvrnt)],
    dtype = object
)

# create a genotype matrix from scratch using NumPy arrays
gmat = DenseGenotypeMatrix(
    mat = mat,
    taxa = taxa,
    taxa_grp = taxa_grp, 
    vrnt_chrgrp = vrnt_chrgrp,
    vrnt_phypos = vrnt_phypos, 
    vrnt_name = vrnt_name, 
    vrnt_genpos = None,
    vrnt_xoprob = None, 
    vrnt_hapgrp = None, 
    vrnt_hapalt = None,
    vrnt_hapref = None, 
    vrnt_mask = None,
    ploidy = ploidy
)

###
### Creating Genomic Models
### =======================

## 
## Creating genomic models from NumPy arrays
## -----------------------------------------

# model parameters
nfixed = 1      # number of fixed effects
ntrait = 2      # number of traits
nmisc = 0       # number of miscellaneous random effects
nadditive = 50  # number of additive marker effects

# create dummy values
beta = numpy.random.random((nfixed,ntrait))
u_misc = numpy.random.random((nmisc,ntrait))
u_a = numpy.random.random((nadditive,ntrait))
trait = numpy.array(
    ["Trait"+str(i+1).zfill(2) for i in range(ntrait)],
    dtype = object
)

# create additive linear genomic model
algmod = DenseAdditiveLinearGenomicModel(
    beta = beta,
    u_misc = u_misc,
    u_a = u_a,
    trait = trait,
    model_name = "example",
    hyperparams = None
)

# export to pandas, csv, hdf5 to demonstrate code below
algmod.to_csv_dict({"beta":"beta.csv","u_misc":"u_misc.csv","u_a":"u_a.csv"})
algmod.to_hdf5("gmod.h5")

##
## Loading genomic models from dictionaries of Pandas DataFrames
## -------------------------------------------------------------

# for LinearGenomicModel classes only!

nfixed = 1      # number of fixed effects
ntrait = 2      # number of traits
nmisc = 0       # number of miscellaneous random effects
nadditive = 50  # number of additive marker effects

# create dummy values
beta = numpy.random.random((nfixed,ntrait))
u_misc = numpy.random.random((nmisc,ntrait))
u_a = numpy.random.random((nadditive,ntrait))
trait = numpy.array(["Trait"+str(i+1) for i in range(ntrait)], dtype=object)

# create dictionary with dataframes
# need required fields of "beta", "u_misc", "u_a" for this class
df_dict = {
    "beta": pandas.DataFrame(beta, columns = trait),
    "u_misc": pandas.DataFrame(u_misc, columns = trait),
    "u_a": pandas.DataFrame(u_a, columns = trait)
}

# construct model from dictionary of pandas dataframes
mod = DenseAdditiveLinearGenomicModel.from_pandas_dict(df_dict)

##
## Loading genomic models from dictionaries of CSV file names
## ----------------------------------------------------------

# for LinearGenomicModel classes only!

# create dictionary with filenames
# need required fields of "beta", "u_misc", "u_a" for this class
dic = {
    "beta": "beta.csv",
    "u_misc": "u_misc.csv",
    "u_a": "u_a.csv",
}

# construct model from dictionary of pandas dataframes
mod = DenseAdditiveLinearGenomicModel.from_csv_dict(dic)

##
## Loading genomic models HDF5 files
## ---------------------------------

# load from HDF5 file
mod = DenseAdditiveLinearGenomicModel.from_hdf5("gmod.h5")

### 
### Genomic Model Properties
### ========================

# general properties
tmp = algmod.model_name
tmp = algmod.hyperparams

# model properties
tmp = algmod.nexplan
tmp = algmod.nparam

# trait properties
tmp = algmod.ntrait
tmp = algmod.trait

###
### Linear Genomic Model Properties
### ===============================

# linear model properties
tmp = algmod.nexplan_beta
tmp = algmod.nexplan_u
tmp = algmod.nparam_beta
tmp = algmod.nparam_u
tmp = algmod.beta
tmp = algmod.u

###
### Additive Linear Genomic Model Properties
### ========================================

# additive linear model properties
tmp = algmod.nexplan_u_a
tmp = algmod.nexplan_u_misc
tmp = algmod.nparam_u_a
tmp = algmod.nparam_u_misc
tmp = algmod.u_a
tmp = algmod.u_misc

###
### Additive Dominance Linear Genomic Model Properties
### ==================================================

# additive dominance linear model properties
# tmp = algmod.nexplan_u_d
# tmp = algmod.nparam_u_d
# tmp = algmod.u_d

###
### Additive Dominance Epistatic Linear Genomic Model Properties
### ============================================================

# additive dominance epistatic linear model properties
# tmp = algmod.nexplan_u_i
# tmp = algmod.nparam_u_i
# tmp = algmod.u_i

############################################################
################## Method Demonstrations ###################
############################################################

###
### Copying models
###

#
# Shallow copying
#

# copy a genomic model
tmp = copy.copy(algmod)
tmp = algmod.copy()

#
# Deep copying
#

# deep copy a genomic model
tmp = copy.deepcopy(algmod)
tmp = algmod.deepcopy()

###
### Model fitting methods
###

algmod.fit_numpy()
algmod.fit()

###
### Model prediction methods
### ========================

# create random genotypes to test
X = numpy.ones((ntaxa,1))
Z = numpy.random.randint(0, ploidy+1, size = (ntaxa,nvrnt)).astype("int8")

# predict genotypic values using numpy arrays
out = algmod.predict_numpy(X, Z)

# predict genotypic values using objects
out = algmod.predict(cvobj = X, gtobj = gmat)

###
### Score model prediction accuracy
### ===============================

# create some dummy matrices for input
X = numpy.ones((ntaxa,1))
B = algmod.beta
Z = gmat.mat
U = algmod.u_a
Y = X@B + Z@U
e = numpy.random.normal(size = Y.shape)
Y += e

# score predictions using numpy arrays
out = algmod.score_numpy(Y, X, Z)

# score predictions using objects
out = algmod.score(ptobj = Y, cvobj = X, gtobj = gmat)

###
### Predicting genomic estimated breeding values
###

# predict GEBVs using numpy arrays
out = algmod.gebv_numpy(Z)

# predict GEBVs using objects
out = algmod.gebv(gmat)

###
### Calculating population genetic variance terms
###

#
# Predicting genetic variance
#

# predict Var(G) using numpy arrays
out = algmod.var_G_numpy(Z)

# predict Var(G) using objects
out = algmod.var_G(gmat)

#
# Predicting additive genetic variance
#

# predict Var(A) using numpy arrays
out = algmod.var_A_numpy(Z)

# predict Var(A) using objects
out = algmod.var_A(gmat)

#
# Predicting additive genic variance
#

# predict Var(a) using numpy arrays
out = algmod.var_a_numpy(
    p = gmat.afreq(),
    ploidy = gmat.ploidy
)

# predict Var(a) using objects
out = algmod.var_a(gmat)

#
# Predicting the Bulmer effect
#

# predict Bulmer effect using numpy arrays
out = algmod.bulmer_numpy(
    Z,
    p = gmat.afreq(),
    ploidy = gmat.ploidy
)

# predict Bulmer effect using objects
out = algmod.bulmer(gmat)

###
### Calculating population selection limits
###

#
# Upper selection limit
#

# upper selection limit using numpy arrays
out = algmod.usl_numpy(
    p = gmat.afreq(),
    ploidy = gmat.ploidy
)

# upper selection limit using objects
out = algmod.usl(gtobj = gmat)

#
# Lower selection limit
#

# lower selection limit using numpy arrays
out = algmod.lsl_numpy(
    p = gmat.afreq(),
    ploidy = gmat.ploidy
)

# lower selection limit using objects
out = algmod.lsl(gtobj = gmat)

###
### Calculating favorable allele metrics
###

# calculate favorable allele counts
out = algmod.facount(gmat)

# calculate favorable allele frequencies
out = algmod.fafreq(gmat)

# calculate favorable allele availability at loci in a population
out = algmod.faavail(gmat)

# calculate favorable allele fixation at loci in a population
out = algmod.fafixed(gmat)

###
### Calculating deleterious allele metrics
###

# calculate deleterious allele counts
out = algmod.dacount(gmat)

# calculate deleterious allele frequencies
out = algmod.dafreq(gmat)

# calculate deleterious allele availability at loci in a population
out = algmod.daavail(gmat)

# calculate deleterious allele fixation at loci in a population
out = algmod.dafixed(gmat)

###
### Exporting Genomic Models
### ========================

## 
## Exporting to dictionaries of Pandas DataFrames
## ----------------------------------------------

# for LinearGenomicModel classes only!

# export all trait columns
out = algmod.to_pandas_dict()

## 
## Exporting to CSV files
## ----------------------

# for LinearGenomicModel classes only!

# create dictionary with filenames
# need required fields of "beta", "u_misc", "u_a" for this class
dic = {
    "beta": "beta.csv",
    "u_misc": "u_misc.csv",
    "u_a": "u_a.csv",
}

# export all trait columns to csv files
algmod.to_csv_dict(dic)

## 
## Exporting to HDF5
## -----------------

# writing a genomic model to a file
algmod.to_hdf5("saved_algmod.h5")


algmod.bulmer()
algmod.bulmer_numpy()
algmod.copy()
algmod.daavail()
algmod.dacount()
algmod.dafixed()
algmod.dafreq()
algmod.deepcopy()
algmod.faavail()
algmod.facount()
algmod.fafixed()
algmod.fafreq()
algmod.fit()
algmod.fit_numpy()
algmod.from_csv_dict()
algmod.from_hdf5()
algmod.from_pandas_dict()
algmod.gebv()
algmod.gebv_numpy()
algmod.gegv()
algmod.gegv_numpy()
algmod.lsl()
algmod.lsl_numpy()
algmod.predict()
algmod.predict_numpy()
algmod.score()
algmod.score_numpy()
algmod.to_csv_dict()
algmod.to_hdf5()
algmod.to_pandas_dict()
algmod.usl()
algmod.usl_numpy()
algmod.var_a()
algmod.var_A()
algmod.var_a_numpy()
algmod.var_A_numpy()
algmod.var_G()
algmod.var_G_numpy()
