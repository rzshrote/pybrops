#!/usr/bin/env python3

import copy
import numpy

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

###
### Coancestry Matrix Object Creation
###

#
# construct from NumPy
#

# shape parameters
ntaxa = 100
ngroup = 20

# create random breeding values
mat = numpy.random.uniform(0.0, 1.0, size = (ntaxa,ntaxa))

# create taxa names
taxa = numpy.array(
    ["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], 
    dtype = object
)

# create taxa groups
taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
taxa_grp.sort()

# create a breeding value matrix from NumPy arrays
cmat = DenseMolecularCoancestryMatrix(
    mat = mat,
    taxa = taxa,
    taxa_grp = taxa_grp
)

###
### Coancestry matrix general properties
###

tmp = cmat.mat              # get the raw coancestry matrix pointer
tmp = cmat.mat_ndim         # get the number of dimensions for the coancestry matrix
tmp = cmat.mat_shape        # get the shape of the coancestry matrix

###
### Coancestry matrix taxa properties
###

tmp = cmat.ntaxa           # get the number of taxa represented by the breeding value matrix
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

# copy a breeding value matrix
tmp = copy.copy(cmat)
tmp = cmat.copy()

# deep copy a breeding value matrix
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

# create a new coancestry matrix to demonstrate
new = cmat.deepcopy()

# insert coancestry matrix along the taxa axis before index 0
tmp = cmat.insert(0, new, axis = cmat.taxa_axis)
tmp = cmat.insert_taxa(0, new)

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

# concatenate along the taxa axis
tmp = cmat.concat([cmat, cmat], axis = cmat.taxa_axis)
tmp = cmat.concat_taxa([cmat, cmat])

###
### Summary Statistics
###


cmat.apply_jitter
cmat.coancestry
cmat.copy
cmat.deepcopy
cmat.from_gmat
cmat.from_hdf5
cmat.group
cmat.group_taxa
cmat.inverse
cmat.is_grouped
cmat.is_grouped_taxa
cmat.is_positive_semidefinite
cmat.is_square
cmat.kinship
cmat.lexsort
cmat.lexsort_taxa
cmat.mat_asformat
cmat.max
cmat.max_inbreeding
cmat.mean
cmat.min
cmat.min_inbreeding
cmat.reorder
cmat.reorder_taxa
cmat.select
cmat.select_taxa
cmat.sort
cmat.sort_taxa
cmat.to_hdf5