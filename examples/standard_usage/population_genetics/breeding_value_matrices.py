#!/usr/bin/env python3

import copy
import numpy

# import the BreedingValueMatrix class (an abstract interface class)
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix

# import the DenseBreedingValueMatrix class (a concrete implemented class)
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix

###
### Breeding Value Matrix Object Creation
###

#
# construct from NumPy
#

# shape parameters
ntaxa = 100
ntrait = 3
ngroup = 20

# create random breeding values
mat = numpy.random.normal(size = (ntaxa,ntrait))

# create taxa names
taxa = numpy.array(
    ["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], 
    dtype = object
)

# create taxa groups
taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
taxa_grp.sort()

# create trait names
trait = numpy.array(
    ["trait"+str(i+1).zfill(2) for i in range(ntrait)],
    dtype = object
)

# create a breeding value matrix from NumPy arrays
bvmat = DenseBreedingValueMatrix(
    mat = mat,
    location = 0.0,
    scale = 1.0,
    taxa = taxa,
    taxa_grp = taxa_grp,
    trait = trait
)

#
# read from HDF5
#

# read a breeding value matrix from an HDF5 file
bvmat = DenseBreedingValueMatrix.from_hdf5("sample_breeding_values.h5")

###
### Breeding value matrix general properties
###

tmp = bvmat.mat         # get the raw breeding value matrix pointer
tmp = bvmat.mat_ndim    # get the number of dimensions for the breeding value matrix
tmp = bvmat.mat_shape   # get the breeding value matrix shape
tmp = bvmat.ntaxa       # get the number of taxa represented by the breeding value matrix
tmp = bvmat.ntrait      # get the number of traits represented by the breeding value matrix
tmp = bvmat.location    # get the location of the breeding value matrix if it has been transformed
tmp = bvmat.scale       # get the scale of the breeding value matrix if it has been transformed

###
### Breeding value matrix taxa properties
###

tmp = bvmat.taxa            # get the names of the taxa
tmp = bvmat.taxa_axis       # get the matrix axis along which taxa are stored
tmp = bvmat.taxa_grp        # get an optional taxa group label
tmp = bvmat.taxa_grp_name   # if taxa are sorted by group: get the names of the groups
tmp = bvmat.taxa_grp_stix   # if taxa are sorted by group: get the start indices (inclusive) for each group
tmp = bvmat.taxa_grp_spix   # if taxa are sorted by group: get the stop indices (exclusive) for each group
tmp = bvmat.taxa_grp_len    # if taxa are sorted by group: get the length of each group

###
### Breeding value matrix trait properties
###

tmp = bvmat.trait       # get the names of the traits
tmp = bvmat.trait_axis  # get the matrix axis along which traits are stored

###
### Copying
###

# copy a breeding value matrix
tmp = copy.copy(bvmat)
tmp = bvmat.copy()

# deep copy a breeding value matrix
tmp = copy.deepcopy(bvmat)
tmp = bvmat.deepcopy()

###
### Genotype Matrix Element Copy-On-Manipulation
###

##
## adjoin examples
##

# create a new genotype matrix to demonstrate
new = bvmat.deepcopy()

# adjoin genotype matrices along the taxa axis
tmp = bvmat.adjoin(new, axis = bvmat.taxa_axis)
tmp = bvmat.adjoin_taxa(new)

# adjoin genotype matrices along the trait axis
tmp = bvmat.adjoin(new, axis = bvmat.trait_axis)
tmp = bvmat.adjoin_trait(new)

##
## delete examples
##

#
# delete taxa examples
#

# delete first taxon using an integer
tmp = bvmat.delete(0, axis = bvmat.taxa_axis)
tmp = bvmat.delete_taxa(0)

# delete first five taxa using a slice
tmp = bvmat.delete(slice(0,5), axis = bvmat.taxa_axis)
tmp = bvmat.delete_taxa(slice(0,5))

# delete first five taxa using a Sequence
tmp = bvmat.delete([0,1,2,3,4], axis = bvmat.taxa_axis)
tmp = bvmat.delete_taxa([0,1,2,3,4])

#
# delete traits examples
#

# delete first trait using an integer
tmp = bvmat.delete(0, axis = bvmat.trait_axis)
tmp = bvmat.delete_trait(0)

# delete first two traits using a slice
tmp = bvmat.delete(slice(0,2), axis = bvmat.trait_axis)
tmp = bvmat.delete_trait(slice(0,2))

# delete first two traits using a Sequence
tmp = bvmat.delete([0,1], axis = bvmat.trait_axis)
tmp = bvmat.delete_trait([0,1])

##
## insert examples
##

# create a new genotype matrix to demonstrate
new = bvmat.deepcopy()

# insert genotype matrix along the taxa axis before index 0
tmp = bvmat.insert(0, new, axis = bvmat.taxa_axis)
tmp = bvmat.insert_taxa(0, new)

# insert genotype matrix along the trait axis before index 0
# tmp = bvmat.insert(0, new, axis = bvmat.trait_axis)
# tmp = bvmat.insert_trait(0, new)

##
## select examples
##

# select first five taxa using a Sequence
tmp = bvmat.select([0,1,2,3,4], axis = bvmat.taxa_axis)
tmp = bvmat.select_taxa([0,1,2,3,4])

# select first two traits using a Sequence
tmp = bvmat.select([0,1], axis = bvmat.trait_axis)
tmp = bvmat.select_trait([0,1])

###
### Genotype Matrix Element In-Place-Manipulation
###

##
## append examples
##

# append genotype matrices along the taxa axis
tmp = bvmat.deepcopy()                   # copy original
tmp.append(bvmat, axis = tmp.taxa_axis)  # append original to copy

tmp = bvmat.deepcopy()                   # copy original
tmp.append_taxa(bvmat)                   # append original to copy

# append genotype matrices along the trait axis
tmp = bvmat.deepcopy()                   # copy original
tmp.append(bvmat, axis = tmp.trait_axis) # append original to copy

tmp = bvmat.deepcopy()                   # copy original
tmp.append_trait(bvmat)                  # append original to copy

##
## remove examples
##

#
# remove taxa examples
#

# remove first taxon using an integer
tmp = bvmat.deepcopy()                           # copy original
tmp.remove(0, axis = bvmat.taxa_axis)            # remove from copy

tmp = bvmat.deepcopy()                           # copy original
tmp.remove_taxa(0)                               # remove from copy

# remove first five taxa using a slice
tmp = bvmat.deepcopy()                           # copy original
tmp.remove(slice(0,5), axis = bvmat.taxa_axis)   # remove from copy

tmp = bvmat.deepcopy()                           # copy original
tmp.remove_taxa(slice(0,5))                      # remove from copy

# remove first five taxa using a Sequence
tmp = bvmat.deepcopy()                           # copy original
tmp.remove([0,1,2,3,4], axis = bvmat.taxa_axis)  # remove from copy

tmp = bvmat.deepcopy()                           # copy original
tmp.remove_taxa([0,1,2,3,4])                     # remove from copy

#
# remove traits examples
#

# remove first trait using an integer
tmp = bvmat.deepcopy()                           # copy original
tmp.remove(0, axis = bvmat.trait_axis)           # remove from copy

tmp = bvmat.deepcopy()                           # copy original
tmp.remove_trait(0)                              # remove from copy

# remove first two traits using a slice
tmp = bvmat.deepcopy()                           # copy original
tmp.remove(slice(0,2), axis = bvmat.trait_axis)  # remove from copy

tmp = bvmat.deepcopy()                           # copy original
tmp.remove_trait(slice(0,2))                     # remove from copy

# remove first two traits using a Sequence
tmp = bvmat.deepcopy()                           # copy original
tmp.remove([0,1], axis = bvmat.trait_axis)       # remove from copy

tmp = bvmat.deepcopy()                           # copy original
tmp.remove_trait([0,1])                          # remove from copy

##
## incorp examples
##

# incorp genotype matrix along the taxa axis before index 0
tmp = bvmat.deepcopy()                           # copy original
tmp.incorp(0, bvmat, axis = bvmat.taxa_axis)     # incorporate into copy

tmp = bvmat.deepcopy()                           # copy original
tmp.incorp_taxa(0, bvmat)                        # incorporate into copy

# incorp genotype matrix along the trait axis before index 0
# tmp = bvmat.deepcopy()                           # copy original
# tmp.incorp(0, bvmat, axis = bvmat.trait_axis)    # incorporate into copy

# tmp = bvmat.deepcopy()                           # copy original
# tmp.incorp_trait(0, bvmat)                       # incorporate into copy

##
## concat examples
##

# concatenate along the taxa axis
tmp = bvmat.concat([bvmat, bvmat], axis = bvmat.taxa_axis)
tmp = bvmat.concat_taxa([bvmat, bvmat])

# concatenate along the trait axis
tmp = bvmat.concat([bvmat, bvmat], axis = bvmat.trait_axis)
tmp = bvmat.concat_trait([bvmat, bvmat])

###
### Summary Statistics
###

# get the indices of the taxa having the maximum values for each trait
out = bvmat.targmax()

# get the indices of the taxa having the minimum values for each trait
out = bvmat.targmin()

# get the maximum breeding values for each trait
out = bvmat.tmax()

# get the mean breeding values for each trait
out = bvmat.tmean()

# get the minimum breeding values for each trait
out = bvmat.tmin()

# get the breeding value ranges for each trait
out = bvmat.trange()

# get the breeding value standard deviations for each trait
out = bvmat.tstd()

# get the breeding value variances for each trait
out = bvmat.tvar()

# de-transform a breeding value matrix 
out = bvmat.descale()



bvmat.from_hdf5
bvmat.from_numpy
bvmat.group
bvmat.group_taxa
bvmat.is_grouped
bvmat.is_grouped_taxa
bvmat.lexsort
bvmat.lexsort_taxa
bvmat.lexsort_trait
bvmat.reorder
bvmat.reorder_taxa
bvmat.reorder_trait
bvmat.sort
bvmat.sort_taxa
bvmat.sort_trait
bvmat.to_hdf5
