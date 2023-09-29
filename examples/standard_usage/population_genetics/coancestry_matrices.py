#!/usr/bin/env python3

import os
import copy
import numpy

### 
### Loading Coancestry Matrix Modules
### =================================

# import the CoancestryMatrix class (an abstract interface class)
from pybrops.popgen.cmat.CoancestryMatrix import CoancestryMatrix

# import the DenseCoancestryMatrix class (a semi-abstract class)
from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix

# import the DenseMolecularCoancestryMatrix class (a concrete implemented class)
from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import DenseMolecularCoancestryMatrix

# import the DenseVanRadenCoancestryMatrix class (a concrete implemented class)
from pybrops.popgen.cmat.DenseVanRadenCoancestryMatrix import DenseVanRadenCoancestryMatrix

# import the DenseYangCoancestryMatrix class (a concrete implemented class)
from pybrops.popgen.cmat.DenseYangCoancestryMatrix import DenseYangCoancestryMatrix
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

### 
### Creating Coancestry Matrices
### ============================

## 
## Creating coancestry matrices from NumPy arrays
## ----------------------------------------------

# shape parameters
ntaxa = 100
ngroup = 20

# create random coancestries
mat = numpy.random.uniform(0.0, 1.0, size = (ntaxa,ntaxa))

# create taxa names
taxa = numpy.array(
    ["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], 
    dtype = object
)

# create taxa groups
taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
taxa_grp.sort()

# create a coancestry matrix from NumPy arrays
cmat = DenseMolecularCoancestryMatrix(
    mat = mat,
    taxa = taxa,
    taxa_grp = taxa_grp
)

### 
### Creating coancestry matrices from GenotypeMatrix objects
### --------------------------------------------------------

# shape parameters for random genotypes
ntaxa = 100
nvrnt = 1000
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

# group taxa and variants
gmat.group_taxa()
gmat.group_vrnt()

# construct Coancestry Matrix from a Genotype Matrix
cmat = DenseMolecularCoancestryMatrix.from_gmat(gmat = gmat)

# save for use in examples below
df = cmat.to_pandas()
cmat.to_csv("saved_coancestry_matrix.csv")

## 
## Creating coancestry matrices from Pandas DataFrames
## ---------------------------------------------------

# load from pandas.DataFrame
tmp = DenseMolecularCoancestryMatrix.from_pandas(
    df = df,
    taxa_col = "taxa",          # column from which to load taxa
    taxa_grp_col = "taxa_grp",  # column from which to load taxa groups
    taxa = "all",               # load all taxa
)

### 
### Loading coancestry matrices from CSV files
### ------------------------------------------

# load from pandas.DataFrame
tmp = DenseMolecularCoancestryMatrix.from_csv(
    filename = "saved_coancestry_matrix.csv",
    taxa_col = "taxa",          # column from which to load taxa
    taxa_grp_col = "taxa_grp",  # column from which to load taxa groups
    taxa = "all",               # load all taxa
)

### 
### Loading coancestry matrices from HDF5 files
### -------------------------------------------

# read from file
cmat = DenseMolecularCoancestryMatrix.from_hdf5("sample_coancestry_matrix.h5")

###
### Coancestry matrix general properties
###

tmp = cmat.mat              # get the raw coancestry matrix pointer
tmp = cmat.mat_ndim         # get the number of dimensions for the coancestry matrix
tmp = cmat.mat_shape        # get the shape of the coancestry matrix

###
### Coancestry matrix taxa properties
###

tmp = cmat.ntaxa           # get the number of taxa represented by the coancestry matrix
tmp = cmat.taxa            # get the names of the taxa
tmp = cmat.taxa_axis       # get the matrix axis along which taxa are stored
tmp = cmat.taxa_grp        # get an optional taxa group label
tmp = cmat.taxa_grp_name   # if taxa are sorted by group: get the names of the groups
tmp = cmat.taxa_grp_stix   # if taxa are sorted by group: get the start indices (inclusive) for each group
tmp = cmat.taxa_grp_spix   # if taxa are sorted by group: get the stop indices (exclusive) for each group
tmp = cmat.taxa_grp_len    # if taxa are sorted by group: get the length of each group

###
### Coancestry matrix square properties
###

tmp = cmat.nsquare          # get the number of square axes for the coancestry matrix
tmp = cmat.square_axes      # get the axes indices for the square axes for the coancestry matrix
tmp = cmat.square_axes_len  # get the lengths of the square axes for the coancestry matrix

###
### Copying
###

# copy a coancestry matrix
tmp = copy.copy(cmat)
tmp = cmat.copy()

# deep copy a coancestry matrix
tmp = copy.deepcopy(cmat)
tmp = cmat.deepcopy()

###
### Genotype Matrix Element Copy-On-Manipulation
###

##
## adjoin examples
##

# create a new coancestry matrix to demonstrate
new = cmat.deepcopy()

# adjoin coancestry matrices along the taxa axis
tmp = cmat.adjoin(new, axis = cmat.taxa_axis)
tmp = cmat.adjoin_taxa(new)

##
## delete examples
##

#
# delete taxa examples
#

# delete first taxon using an integer
tmp = cmat.delete(0, axis = cmat.taxa_axis)
tmp = cmat.delete_taxa(0)

# delete first five taxa using a slice
tmp = cmat.delete(slice(0,5), axis = cmat.taxa_axis)
tmp = cmat.delete_taxa(slice(0,5))

# delete first five taxa using a Sequence
tmp = cmat.delete([0,1,2,3,4], axis = cmat.taxa_axis)
tmp = cmat.delete_taxa([0,1,2,3,4])

##
## insert examples
##

# # create a new coancestry matrix to demonstrate
# new = cmat.deepcopy()

# # insert coancestry matrix along the taxa axis before index 0
# tmp = cmat.insert(0, new, axis = cmat.taxa_axis)
# tmp = cmat.insert_taxa(0, new)

# insert coancestry matrix along the trait axis before index 0
# tmp = cmat.insert(0, new, axis = cmat.trait_axis)
# tmp = cmat.insert_trait(0, new)

##
## select examples
##

# select first five taxa using a Sequence
tmp = cmat.select([0,1,2,3,4], axis = cmat.taxa_axis)
tmp = cmat.select_taxa([0,1,2,3,4])

###
### Genotype Matrix Element In-Place-Manipulation
###

##
## append examples
##

# append coancestry matrices along the taxa axis
tmp = cmat.deepcopy()                   # copy original
tmp.append(cmat, axis = tmp.taxa_axis)  # append original to copy

tmp = cmat.deepcopy()                   # copy original
tmp.append_taxa(cmat)                   # append original to copy

##
## remove examples
##

#
# remove taxa examples
#

# remove first taxon using an integer
tmp = cmat.deepcopy()                           # copy original
tmp.remove(0, axis = cmat.taxa_axis)            # remove from copy

tmp = cmat.deepcopy()                           # copy original
tmp.remove_taxa(0)                               # remove from copy

# remove first five taxa using a slice
tmp = cmat.deepcopy()                           # copy original
tmp.remove(slice(0,5), axis = cmat.taxa_axis)   # remove from copy

tmp = cmat.deepcopy()                           # copy original
tmp.remove_taxa(slice(0,5))                      # remove from copy

# remove first five taxa using a Sequence
tmp = cmat.deepcopy()                           # copy original
tmp.remove([0,1,2,3,4], axis = cmat.taxa_axis)  # remove from copy

tmp = cmat.deepcopy()                           # copy original
tmp.remove_taxa([0,1,2,3,4])                     # remove from copy

##
## incorp examples
##

# incorp coancestry matrix along the taxa axis before index 0
tmp = cmat.deepcopy()                           # copy original
tmp.incorp(0, cmat, axis = cmat.taxa_axis)     # incorporate into copy

tmp = cmat.deepcopy()                           # copy original
tmp.incorp_taxa(0, cmat)                        # incorporate into copy

##
## concat examples
##

# # concatenate along the taxa axis
# tmp = cmat.concat([cmat, cmat], axis = cmat.taxa_axis)
# tmp = cmat.concat_taxa([cmat, cmat])

###
### Grouping and sorting
###

##
## Reordering
##

#
# taxa reordering example
#

# create reordering indices
indices = numpy.arange(cmat.ntaxa)
numpy.random.shuffle(indices)
tmp = cmat.deepcopy()

# reorder values along the taxa axis
tmp.reorder(indices, axis = tmp.taxa_axis)
tmp.reorder_taxa(indices)

##
## Lexsorting
##

#
# taxa lexsort example
#

# create lexsort keys for taxa
key1 = numpy.random.randint(0, 10, cmat.ntaxa)
key2 = numpy.arange(cmat.ntaxa)
numpy.random.shuffle(key2)

# lexsort along the taxa axis
cmat.lexsort((key2,key1), axis = cmat.taxa_axis)
cmat.lexsort_taxa((key2,key1))

##
## Sorting
##

# make copy
tmp = cmat.deepcopy()

#
# taxa sorting example
#

# sort along taxa axis
tmp.sort(axis = tmp.taxa_axis)
tmp.sort_taxa()

##
## Grouping
##

# make copy
tmp = cmat.deepcopy()

#
# taxa grouping example
#

# sort along taxa axis
tmp.group(axis = tmp.taxa_axis)
tmp.group_taxa()

# determine whether grouping has occurred along the taxa axis
out = tmp.is_grouped(axis = tmp.taxa_axis)
out = tmp.is_grouped_taxa()

###
### Coancestry/kinship Methods
###

#
# Get the coancestry at a specific matrix coordinate
#

out = cmat.coancestry(0,0)
out = cmat[0,0] # NOT guaranteed to be in correct format

#
# Get the kinship at a specific matrix coordinate
#

out = cmat.kinship(0,0)
out = 0.5 * cmat[0,0] # NOT guaranteed to be in correct format

#
# Get the coancestry matrix as a specific format
#

out = cmat.mat_asformat(format = "kinship")
out = cmat.mat_asformat(format = "coancestry")

#
# Determine if the coancestry matrix is positive semidefinite (convex)
#

out = cmat.is_positive_semidefinite()

#
# Apply a jitter along the diagonal to try to make the matrix positive semidefinite
#

out = cmat.apply_jitter()

#
# Calculate the maximum attainable inbreeding after 1 generation
#

out = cmat.max_inbreeding()
out = cmat.max_inbreeding(format = "kinship")
out = cmat.max_inbreeding(format = "coancestry")

#
# Calculate the minimum attainable inbreeding after 1 generation
#

out = cmat.min_inbreeding()
out = cmat.min_inbreeding(format = "kinship")
out = cmat.min_inbreeding(format = "coancestry")

#
# Calculate the inverse of the coancestry matrix
#

out = cmat.inverse()
out = cmat.inverse(format = "kinship")
out = cmat.inverse(format = "coancestry")

###
### Summary Statistics
###

# get the max for the whole coancestry matrix
out = cmat.max()

# get the mean for the whole coancestry matrix
out = cmat.mean()

# get the min for the whole coancestry matrix
out = cmat.min()

###
### Saving Breeding Value Matrices
###

## 
## Exporting to Pandas DataFrame
## -----------------------------

# export to a pandas.DataFrame
# use default column names to export
df = cmat.to_pandas()

## 
## Exporting to CSV
## ----------------

# export to a CSV
# use default column names to export
cmat.to_csv("saved_coancestry_matrix.csv")

## 
## Exporting to HDF5
## -----------------

# remove exported file if it exists
if os.path.exists("saved_coancestry_matrix.h5"):
    os.remove("saved_coancestry_matrix.h5")

# write a coancestry matrix to an HDF5 file
cmat.to_hdf5("saved_coancestry_matrix.h5")



cmat.is_square