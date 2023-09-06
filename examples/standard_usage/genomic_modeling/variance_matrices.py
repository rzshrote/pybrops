#!/usr/bin/env python3

import copy
import os
import numpy
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

###
### Loading Variance Matrix Modules
### ===============================

#
# Loading genetic variance matrix modules
# ---------------------------------------

# import abstract interface classes
from pybrops.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix
from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix

# import semi-abstract classes
from pybrops.model.vmat.DenseGeneticVarianceMatrix import DenseGeneticVarianceMatrix
from pybrops.model.vmat.DenseAdditiveGeneticVarianceMatrix import DenseAdditiveGeneticVarianceMatrix

# import concrete implemented classes
from pybrops.model.vmat.DenseTwoWayDHAdditiveGeneticVarianceMatrix import DenseTwoWayDHAdditiveGeneticVarianceMatrix
from pybrops.model.vmat.DenseThreeWayDHAdditiveGeneticVarianceMatrix import DenseThreeWayDHAdditiveGeneticVarianceMatrix
from pybrops.model.vmat.DenseFourWayDHAdditiveGeneticVarianceMatrix import DenseFourWayDHAdditiveGeneticVarianceMatrix
from pybrops.model.vmat.DenseDihybridDHAdditiveGeneticVarianceMatrix import DenseDihybridDHAdditiveGeneticVarianceMatrix

#
# Loading genic variance matrix modules
# -------------------------------------

# import abstract interface classes
from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix
from pybrops.model.vmat.AdditiveGenicVarianceMatrix import AdditiveGenicVarianceMatrix

# import semi-abstract classes
from pybrops.model.vmat.DenseGenicVarianceMatrix import DenseGenicVarianceMatrix
from pybrops.model.vmat.DenseAdditiveGenicVarianceMatrix import DenseAdditiveGenicVarianceMatrix

# import concrete implemented classes
from pybrops.model.vmat.DenseTwoWayDHAdditiveGenicVarianceMatrix import DenseTwoWayDHAdditiveGenicVarianceMatrix
from pybrops.model.vmat.DenseThreeWayDHAdditiveGenicVarianceMatrix import DenseThreeWayDHAdditiveGenicVarianceMatrix
from pybrops.model.vmat.DenseFourWayDHAdditiveGenicVarianceMatrix import DenseFourWayDHAdditiveGenicVarianceMatrix
from pybrops.model.vmat.DenseDihybridDHAdditiveGenicVarianceMatrix import DenseDihybridDHAdditiveGenicVarianceMatrix

###
### Creating Variance Matrices
### ==========================

#
# Creating variance matrices from NumPy arrays
# --------------------------------------------

# shape parameters for random genotypes
ntaxa = 100
ntrait = 2
ngroup = 20

# create random variance values
mat = numpy.random.uniform(0, 1, size = (ntaxa,ntaxa,ntrait))

# create taxa names
taxa = numpy.array(["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], dtype = object)

# create taxa groups
taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
taxa_grp.sort()

# create trait names
trait = numpy.array(["trait"+str(i+1).zfill(2) for i in range(ntrait)], dtype = object)

# create genetic variance matrix
vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix(
    mat = mat,
    taxa = taxa,
    taxa_grp = taxa_grp,
    trait = trait
)

# create genic variance matrix
gvmat = DenseTwoWayDHAdditiveGenicVarianceMatrix(
    mat = mat,
    taxa = taxa,
    taxa_grp = taxa_grp,
    trait = trait
)

#
# Creating variance matrices from genomic models
# ----------------------------------------------

# create a dummy genomic model
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

# shape parameters for random genotypes
ntaxa = 100
nvrnt = nadditive
ngroup = 20
nchrom = 10
ploidy = 2

# create random genotypes
mat = numpy.random.randint(0, 2, size = (ploidy,ntaxa,nvrnt)).astype("int8")

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

# create marker genetic positions
vrnt_genpos = numpy.random.random(nvrnt)
vrnt_genpos.sort()

# create marker variant names
vrnt_name = numpy.array(["SNP"+str(i+1).zfill(4) for i in range(nvrnt)], dtype = object)

# create a genotype matrix from scratch using NumPy arrays
pgmat = DensePhasedGenotypeMatrix(
    mat = mat,
    taxa = taxa,
    taxa_grp = taxa_grp, 
    vrnt_chrgrp = vrnt_chrgrp,
    vrnt_phypos = vrnt_phypos, 
    vrnt_genpos = vrnt_genpos,
    vrnt_name = vrnt_name, 
    ploidy = ploidy
)
pgmat.group_vrnt()

# calculate genetic variance matrix from GenomicModel
vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_gmod(
    gmod = algmod,
    pgmat = pgmat,
    ncross = 1,
    nprogeny = 10,
    nself = 0,
    gmapfn = HaldaneMapFunction()
)

# calculate genetic variance matrix from AdditiveLinearGenomicModel
vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_algmod(
    algmod = algmod,
    pgmat = pgmat,
    ncross = 1,
    nprogeny = 10,
    nself = 0,
    gmapfn = HaldaneMapFunction()
)

# calculate genic variance matrix from GenomicModel
gvmat = DenseTwoWayDHAdditiveGenicVarianceMatrix.from_gmod(
    gmod = algmod,
    pgmat = pgmat,
    nprogeny = 10
)

# calculate genic variance matrix from AdditiveLinearGenomicModel
gvmat = DenseTwoWayDHAdditiveGenicVarianceMatrix.from_algmod(
    algmod = algmod,
    pgmat = pgmat,
    nprogeny = 10
)

# save example variance matrices to HDF5 files
# vmat.to_hdf5("saved_vmat.h5")
# gvmat.to_hdf5("saved_gvmat.h5")

###
### Loading variance matrices from HDF5 files
### -----------------------------------------

# read genetic variance matrix from HDF5 file
vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_hdf5("saved_vmat.h5")
gvmat = DenseTwoWayDHAdditiveGenicVarianceMatrix.from_hdf5("saved_gvmat.h5")

###
### Variance matrix properties 
### ==========================

###
### General properties 
### ------------------

# get the raw variance matrix pointer
tmp = vmat.mat
tmp = gvmat.mat

# get the number of dimensions for the variance matrix
tmp = vmat.mat_ndim
tmp = gvmat.mat_ndim

# get the shape of the variance matrix
tmp = vmat.mat_shape
tmp = gvmat.mat_shape

# get the expected parental genomic contribution
tmp = vmat.epgc
tmp = gvmat.epgc

###
### Square properties
### -----------------

# get the number of square axes for the variance matrix
tmp = vmat.nsquare
tmp = gvmat.nsquare

# get the axes indices for the square axes for the variance matrix
tmp = vmat.square_axes
tmp = gvmat.square_axes

# get the lengths of the square axes for the variance matrix
tmp = vmat.square_axes_len
tmp = gvmat.square_axes_len

###
### Taxa properties
### ---------------

# get the number of taxa for the variance matrix
tmp = vmat.ntaxa
tmp = gvmat.ntaxa

# get the names of the taxa for the variance matrix
tmp = vmat.taxa
tmp = gvmat.taxa

# get the taxa axis for the variance matrix
tmp = vmat.taxa_axis
tmp = gvmat.taxa_axis

# get the taxa groups for the variance matrix
tmp = vmat.taxa_grp
tmp = gvmat.taxa_grp

# get the taxa group lengths for the variance matrix
tmp = vmat.taxa_grp_len
tmp = gvmat.taxa_grp_len

# get the taxa group names for the variance matrix
tmp = vmat.taxa_grp_name
tmp = gvmat.taxa_grp_name

# get the taxa group start indices for the variance matrix
tmp = vmat.taxa_grp_stix
tmp = gvmat.taxa_grp_stix

# get the taxa group stop indices for the variance matrix
tmp = vmat.taxa_grp_spix
tmp = gvmat.taxa_grp_spix

###
### Trait properties
### ----------------

# get the number of traits represented by the variance matrix
tmp = vmat.ntrait
tmp = gvmat.ntrait

# get the names of the traits represented by the variance matrix
tmp = vmat.trait
tmp = gvmat.trait

# get the trait axis for the variance matrix
tmp = vmat.trait_axis
tmp = gvmat.trait_axis

###
### Copying
### =======

# copy a genetic variance matrix
tmp = copy.copy(vmat)
tmp = vmat.copy()

# copy a genic variance matrix
tmp = copy.copy(gvmat)
tmp = gvmat.copy()

# deep copy a genetic variance matrix
tmp = copy.deepcopy(vmat)
tmp = vmat.deepcopy()

# deep copy a genic variance matrix
tmp = copy.deepcopy(gvmat)
tmp = gvmat.deepcopy()

###
### Variance Matrix Element Copy-On-Manipulation
### ============================================

##
## ``adjoin`` examples
## -------------------

# create a new variance matrices to demonstrate
newvmat = vmat.deepcopy()
newgvmat = gvmat.deepcopy()

# adjoin variance matrices along the taxa axis
tmp = vmat.adjoin(newvmat, axis = vmat.taxa_axis)
tmp = vmat.adjoin_taxa(newvmat)
tmp = gvmat.adjoin(newgvmat, axis = gvmat.taxa_axis)
tmp = gvmat.adjoin_taxa(newgvmat)

# adjoin variance matrices along the trait axis
tmp = vmat.adjoin(newvmat, axis = vmat.trait_axis)
tmp = vmat.adjoin_trait(newvmat)
tmp = gvmat.adjoin(newgvmat, axis = gvmat.trait_axis)
tmp = gvmat.adjoin_trait(newgvmat)

##
## ``delete`` examples
## -------------------

#
# ``delete`` taxa examples
# ++++++++++++++++++++++++

# delete first taxon using an integer
tmp = vmat.delete(0, axis = vmat.taxa_axis)
tmp = vmat.delete_taxa(0)
tmp = gvmat.delete(0, axis = gvmat.taxa_axis)
tmp = gvmat.delete_taxa(0)

# delete first five taxa using a slice
tmp = vmat.delete(slice(0,5), axis = vmat.taxa_axis)
tmp = vmat.delete_taxa(slice(0,5))
tmp = gvmat.delete(slice(0,5), axis = gvmat.taxa_axis)
tmp = gvmat.delete_taxa(slice(0,5))

# delete first five taxa using a Sequence
tmp = vmat.delete([0,1,2,3,4], axis = vmat.taxa_axis)
tmp = vmat.delete_taxa([0,1,2,3,4])
tmp = gvmat.delete([0,1,2,3,4], axis = gvmat.taxa_axis)
tmp = gvmat.delete_taxa([0,1,2,3,4])

#
# ``delete`` traits examples
# ++++++++++++++++++++++++++

# delete first trait using an integer
tmp = vmat.delete(0, axis = vmat.trait_axis)
tmp = vmat.delete_trait(0)
tmp = gvmat.delete(0, axis = gvmat.trait_axis)
tmp = gvmat.delete_trait(0)

# delete first two traits using a slice
tmp = vmat.delete(slice(0,2), axis = vmat.trait_axis)
tmp = vmat.delete_trait(slice(0,2))
tmp = gvmat.delete(slice(0,2), axis = gvmat.trait_axis)
tmp = gvmat.delete_trait(slice(0,2))

# delete first two traits using a Sequence
tmp = vmat.delete([0,1], axis = vmat.trait_axis)
tmp = vmat.delete_trait([0,1])
tmp = gvmat.delete([0,1], axis = gvmat.trait_axis)
tmp = gvmat.delete_trait([0,1])

##
## ``insert`` examples
## -------------------

# # create a new variance matrix to demonstrate
# newvmat = vmat.deepcopy()
# newgvmat = gvmat.deepcopy()

# # insert variance matrix along the taxa axis before index 0
# tmp = vmat.insert(0, newvmat, axis = vmat.taxa_axis)
# tmp = vmat.insert_taxa(0, newvmat)
# tmp = gvmat.insert(0, newgvmat, axis = gvmat.taxa_axis)
# tmp = gvmat.insert_taxa(0, newgvmat)

# # insert variance matrix along the trait axis before index 0
# tmp = vmat.insert(0, newvmat, axis = vmat.trait_axis)
# tmp = vmat.insert_trait(0, newvmat)
# tmp = gvmat.insert(0, newgvmat, axis = gvmat.trait_axis)
# tmp = gvmat.insert_trait(0, newgvmat)

##
## ``select`` examples
## -------------------

# select first five taxa using a Sequence
tmp = vmat.select([0,1,2,3,4], axis = vmat.taxa_axis)
tmp = vmat.select_taxa([0,1,2,3,4])
tmp = gvmat.select([0,1,2,3,4], axis = gvmat.taxa_axis)
tmp = gvmat.select_taxa([0,1,2,3,4])

# select first two traits using a Sequence
tmp = vmat.select([0,1], axis = vmat.trait_axis)
tmp = vmat.select_trait([0,1])
tmp = gvmat.select([0,1], axis = gvmat.trait_axis)
tmp = gvmat.select_trait([0,1])

###
### Variance Matrix Element In-Place-Manipulation
### =============================================

##
## ``append`` examples
## -------------------

# append variance matrices along the taxa axis
tmp = vmat.deepcopy()                   # copy original
tmp.append(vmat, axis = tmp.taxa_axis)  # append original to copy
tmp = gvmat.deepcopy()                  # copy original
tmp.append(gvmat, axis = tmp.taxa_axis) # append original to copy

tmp = vmat.deepcopy()                   # copy original
tmp.append_taxa(vmat)                   # append original to copy
tmp = gvmat.deepcopy()                  # copy original
tmp.append_taxa(gvmat)                  # append original to copy

# append variance matrices along the trait axis
tmp = vmat.deepcopy()                   # copy original
tmp.append(vmat, axis = tmp.trait_axis) # append original to copy
tmp = gvmat.deepcopy()                  # copy original
tmp.append(gvmat, axis = tmp.trait_axis)# append original to copy

tmp = vmat.deepcopy()                   # copy original
tmp.append_trait(vmat)                  # append original to copy
tmp = gvmat.deepcopy()                  # copy original
tmp.append_trait(gvmat)                 # append original to copy

##
## ``remove`` examples
## -------------------

#
# ``remove`` taxa examples
#

# remove first taxon using an integer
tmp = vmat.deepcopy()                           # copy original
tmp.remove(0, axis = vmat.taxa_axis)            # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove(0, axis = gvmat.taxa_axis)           # remove from copy

tmp = vmat.deepcopy()                           # copy original
tmp.remove_taxa(0)                              # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove_taxa(0)                              # remove from copy

# remove first five taxa using a slice
tmp = vmat.deepcopy()                           # copy original
tmp.remove(slice(0,5), axis = vmat.taxa_axis)   # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove(slice(0,5), axis = gvmat.taxa_axis)  # remove from copy

tmp = vmat.deepcopy()                           # copy original
tmp.remove_taxa(slice(0,5))                     # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove_taxa(slice(0,5))                     # remove from copy

# remove first five taxa using a Sequence
tmp = vmat.deepcopy()                           # copy original
tmp.remove([0,1,2,3,4], axis = vmat.taxa_axis)  # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove([0,1,2,3,4], axis = gvmat.taxa_axis) # remove from copy

tmp = vmat.deepcopy()                           # copy original
tmp.remove_taxa([0,1,2,3,4])                    # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove_taxa([0,1,2,3,4])                    # remove from copy

#
# ``remove`` traits examples
#

# remove first trait using an integer
tmp = vmat.deepcopy()                           # copy original
tmp.remove(0, axis = vmat.trait_axis)           # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove(0, axis = gvmat.trait_axis)          # remove from copy

tmp = vmat.deepcopy()                           # copy original
tmp.remove_trait(0)                             # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove_trait(0)                             # remove from copy

# remove first trait using a slice
tmp = vmat.deepcopy()                           # copy original
tmp.remove(slice(0,1), axis = vmat.trait_axis)  # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove(slice(0,1), axis = gvmat.trait_axis) # remove from copy

tmp = vmat.deepcopy()                           # copy original
tmp.remove_trait(slice(0,1))                    # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove_trait(slice(0,1))                    # remove from copy

# remove first trait using a Sequence
tmp = vmat.deepcopy()                           # copy original
tmp.remove([0], axis = vmat.trait_axis)         # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove([0], axis = gvmat.trait_axis)        # remove from copy

tmp = vmat.deepcopy()                           # copy original
tmp.remove_trait([0])                           # remove from copy
tmp = gvmat.deepcopy()                          # copy original
tmp.remove_trait([0])                           # remove from copy

##
## ``incorp`` examples
## -------------------

# incorp variance matrix along the taxa axis before index 0
tmp = vmat.deepcopy()                           # copy original
tmp.incorp(0, vmat, axis = vmat.taxa_axis)      # incorporate into copy
tmp = gvmat.deepcopy()                          # copy original
tmp.incorp(0, gvmat, axis = gvmat.taxa_axis)    # incorporate into copy

tmp = vmat.deepcopy()                           # copy original
tmp.incorp_taxa(0, vmat)                        # incorporate into copy
tmp = gvmat.deepcopy()                          # copy original
tmp.incorp_taxa(0, gvmat)                       # incorporate into copy

# incorp variance matrix along the trait axis before index 0
tmp = vmat.deepcopy()                           # copy original
tmp.incorp(0, vmat, axis = vmat.trait_axis)     # incorporate into copy
tmp = gvmat.deepcopy()                          # copy original
tmp.incorp(0, gvmat, axis = gvmat.trait_axis)   # incorporate into copy

tmp = vmat.deepcopy()                           # copy original
tmp.incorp_trait(0, vmat)                       # incorporate into copy
tmp = gvmat.deepcopy()                          # copy original
tmp.incorp_trait(0, gvmat)                      # incorporate into copy

##
## ``concat`` examples
## -------------------

# concatenate along the taxa axis
tmp = vmat.concat([vmat, vmat], axis = vmat.taxa_axis)
tmp = vmat.concat_taxa([vmat, vmat])
tmp = gvmat.concat([gvmat, gvmat], axis = gvmat.taxa_axis)
tmp = gvmat.concat_taxa([gvmat, gvmat])

# concatenate along the trait axis
tmp = vmat.concat([vmat, vmat], axis = vmat.trait_axis)
tmp = vmat.concat_trait([vmat, vmat])
tmp = gvmat.concat([gvmat, gvmat], axis = gvmat.trait_axis)
tmp = gvmat.concat_trait([gvmat, gvmat])

###
### Grouping and sorting
### ====================

##
## Reordering
## ----------

#
# taxa ``reorder``ing example
# +++++++++++++++++++++++++++

# create reordering indices
indices = numpy.arange(vmat.ntaxa)
numpy.random.shuffle(indices)

# reorder values along the taxa axis
tmp = vmat.deepcopy()
tmp.reorder(indices, axis = tmp.taxa_axis)
tmp.reorder_taxa(indices)
tmp = gvmat.deepcopy()
tmp.reorder(indices, axis = tmp.taxa_axis)
tmp.reorder_taxa(indices)

#
# trait ``reorder``ing example
# ++++++++++++++++++++++++++++

# create reordering indices
indices = numpy.arange(vmat.ntrait)
numpy.random.shuffle(indices)

# reorder values along the trait axis
tmp = vmat.deepcopy()
tmp.reorder(indices, axis = tmp.trait_axis)
tmp.reorder_trait(indices)
tmp = gvmat.deepcopy()
tmp.reorder(indices, axis = tmp.trait_axis)
tmp.reorder_trait(indices)

##
## Lexsorting
## ----------

#
# taxa ``lexsort`` example
# ++++++++++++++++++++++++

# create lexsort keys for taxa
key1 = numpy.random.randint(0, 10, vmat.ntaxa)
key2 = numpy.arange(vmat.ntaxa)
numpy.random.shuffle(key2)

# lexsort along the taxa axis
vmat.lexsort((key2,key1), axis = vmat.taxa_axis)
vmat.lexsort_taxa((key2,key1))
gvmat.lexsort((key2,key1), axis = gvmat.taxa_axis)
gvmat.lexsort_taxa((key2,key1))

#
# trait ``lexsort`` example
# +++++++++++++++++++++++++

# create lexsort keys for trait
key1 = numpy.random.randint(0, 10, vmat.ntaxa)
key2 = numpy.arange(vmat.ntaxa)
numpy.random.shuffle(key2)

# lexsort along the trait axis
vmat.lexsort((key2,key1), axis = vmat.taxa_axis)
vmat.lexsort_taxa((key2,key1))
gvmat.lexsort((key2,key1), axis = gvmat.taxa_axis)
gvmat.lexsort_taxa((key2,key1))

##
## Sorting
## -------

#
# taxa ``sort``ing example
# ++++++++++++++++++++++++

# sort along taxa axis
tmp = vmat.deepcopy()
tmp.sort(axis = tmp.taxa_axis)
tmp.sort_taxa()
tmp = gvmat.deepcopy()
tmp.sort(axis = tmp.taxa_axis)
tmp.sort_taxa()

#
# trait ``sort``ing example
# +++++++++++++++++++++++++

# sort along trait axis
tmp = vmat.deepcopy()
tmp.sort(axis = tmp.trait_axis)
tmp.sort_trait()
tmp = gvmat.deepcopy()
tmp.sort(axis = tmp.trait_axis)
tmp.sort_trait()

##
## Grouping
## --------

#
# taxa ``group``ing example
# +++++++++++++++++++++++++

# sort genetic along taxa axis
tmp = vmat.deepcopy()
tmp.group(axis = tmp.taxa_axis)
tmp.group_taxa()
# determine whether grouping has occurred along the taxa axis
tmp.is_grouped(axis = tmp.taxa_axis)
tmp.is_grouped_taxa()

# sort genic variance matrix along taxa axis
tmp = gvmat.deepcopy()
tmp.group(axis = tmp.taxa_axis)
tmp.group_taxa()
# determine whether grouping has occurred along the taxa axis
tmp.is_grouped(axis = tmp.taxa_axis)
tmp.is_grouped_taxa()

###
### Saving Breeding Value Matrices
### ==============================

#
# write to HDF5
# -------------

# remove exported file if it exists
if os.path.exists("saved_vmat.h5"):
    os.remove("saved_vmat.h5")

# write a breeding value matrix to an HDF5 file
vmat.to_hdf5("saved_vmat.h5")



vmat.is_square()
vmat.to_csv()
vmat.to_hdf5()
