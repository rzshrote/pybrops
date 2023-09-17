#!/usr/bin/env python3

import os
import copy
import numpy

###
### Loading Class Modules
### =====================

# import the PhasedGenotypeMatrix class (an abstract interface class)
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

# import the DensePhasedGenotypeMatrix class (a concrete implemented class)
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

###
### Creating Phased Genotype Matrices
### =================================

#
# Creating phased genotype matrices from NumPy arrays
# ---------------------------------------------------

# shape parameters for random genotypes
ntaxa = 100
nvrnt = 1000
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

#
# Loading phased genotype matrices from VCF files
# -----------------------------------------------

# read a phased genotype matrix from VCF file
pgmat = DensePhasedGenotypeMatrix.from_vcf("widiv_2000SNPs.vcf.gz")

#
# Loading phased genotype matrices from HDF5 files
# ------------------------------------------------

# read a genotype matrix from HDF5 file
pgmat = DensePhasedGenotypeMatrix.from_hdf5("widiv_2000SNPs.h5")

###
### Phased Genotype Matrix Properties
### =================================

#
# General properties
# ------------------

tmp = pgmat.mat              # gain access to raw genotype matrix pointer
tmp = pgmat.mat_ndim         # get the number of dimensions for the genotype matrix
tmp = pgmat.mat_shape        # get genotype matrix shape
tmp = pgmat.mat_format       # get genotype matrix format
tmp = pgmat.ploidy           # get the ploidy of the taxa represented by the genotype matrix

#
# Phase properties
# ----------------

tmp = pgmat.nphase           # get the number of chromosome phases represented by the genotype matrix
tmp = pgmat.phase_axis       # get the axis along which phases are stored

#
# Taxa properties
# ---------------

tmp = pgmat.ntaxa            # get the number of taxa represented by the genotype matrix
tmp = pgmat.taxa             # get the names of the taxa
tmp = pgmat.taxa_axis        # get the matrix axis along which taxa are stored
tmp = pgmat.taxa_grp         # get an optional taxa group label
tmp = pgmat.taxa_grp_name    # if taxa are sorted by group: get the names of the groups
tmp = pgmat.taxa_grp_stix    # if taxa are sorted by group: get the start indices (inclusive) for each group
tmp = pgmat.taxa_grp_spix    # if taxa are sorted by group: get the stop indices (exclusive) for each group
tmp = pgmat.taxa_grp_len     # if taxa are sorted by group: get the length of each group

#
# Marker variant properties
# -------------------------

tmp = pgmat.nvrnt            # get the number of genotype variants represented by the genotype matrix
tmp = pgmat.vrnt_name        # get the names of the marker variants
tmp = pgmat.vrnt_axis        # get the axis along which marker variants are stored
tmp = pgmat.vrnt_chrgrp      # get the chromosome to which a marker variant belongs
tmp = pgmat.vrnt_phypos      # get the physical position of a marker variant
tmp = pgmat.vrnt_genpos      # get the genetic position of a marker variant
tmp = pgmat.vrnt_xoprob      # get the crossover probability between the previous marker
tmp = pgmat.vrnt_hapref      # get the reference haplotype for the marker variant
tmp = pgmat.vrnt_hapalt      # get the alternative haplotype for the marker variant
tmp = pgmat.vrnt_hapgrp      # get the haplotype grouping for the marker variant
tmp = pgmat.vrnt_mask        # get a mask associated with the marker variants
tmp = pgmat.vrnt_chrgrp_name # get the names of the chromosomes if sorted
tmp = pgmat.vrnt_chrgrp_stix # get the start indices (inclusive) for each chromosome if sorted
tmp = pgmat.vrnt_chrgrp_spix # get the stop indices (exclusive) for each chromosome if sorted
tmp = pgmat.vrnt_chrgrp_len  # get the length of each chromosome if sorted

###
### Copying Phased Genotype Matrices
### ================================

# copy a phased genotype matrix
tmp = copy.copy(pgmat)
tmp = pgmat.copy()

# deep copy a phased genotype matrix
tmp = copy.deepcopy(pgmat)
tmp = pgmat.deepcopy()

###
### Phased Genotype Matrix Element Copy-On-Manipulation
### ===================================================

##
## Adjoining elements
## ------------------

# create a new genotype matrix to demonstrate
new = pgmat.deepcopy()

# adjoin genotype matrices along the taxa axis
tmp = pgmat.adjoin(new, axis = pgmat.taxa_axis)
tmp = pgmat.adjoin_taxa(new)

# adjoin genotype matrices along the variant axis
tmp = pgmat.adjoin(new, axis = pgmat.vrnt_axis)
tmp = pgmat.adjoin_vrnt(new)

##
## Deleting elements
## -----------------

#
# ``delete`` taxa
# +++++++++++++++

# delete first taxon using an integer
tmp = pgmat.delete(0, axis = pgmat.taxa_axis)
tmp = pgmat.delete_taxa(0)

# delete first five taxa using a slice
tmp = pgmat.delete(slice(0,5), axis = pgmat.taxa_axis)
tmp = pgmat.delete_taxa(slice(0,5))

# delete first five taxa using a Sequence
tmp = pgmat.delete([0,1,2,3,4], axis = pgmat.taxa_axis)
tmp = pgmat.delete_taxa([0,1,2,3,4])

#
# ``delete`` marker variants
# ++++++++++++++++++++++++++

# delete first marker variant using an integer
tmp = pgmat.delete(0, axis = pgmat.vrnt_axis)
tmp = pgmat.delete_vrnt(0)

# delete first five marker variants using a slice
tmp = pgmat.delete(slice(0,5), axis = pgmat.vrnt_axis)
tmp = pgmat.delete_vrnt(slice(0,5))

# delete first five marker variants using a Sequence
tmp = pgmat.delete([0,1,2,3,4], axis = pgmat.vrnt_axis)
tmp = pgmat.delete_vrnt([0,1,2,3,4])

##
## Inserting elements
## ------------------

# create a new genotype matrix to demonstrate
new = pgmat.deepcopy()

# insert genotype matrix along the taxa axis before index 0
tmp = pgmat.insert(0, new, axis = pgmat.taxa_axis)
tmp = pgmat.insert_taxa(0, new)

# insert genotype matrix along the variant axis before index 0
tmp = pgmat.insert(0, new, axis = pgmat.vrnt_axis)
tmp = pgmat.insert_vrnt(0, new)

##
## Selecting elements
## ------------------

# select first five taxa using a Sequence
tmp = pgmat.select([0,1,2,3,4], axis = pgmat.taxa_axis)
tmp = pgmat.select_taxa([0,1,2,3,4])

# select first five marker variants using a Sequence
tmp = pgmat.select([0,1,2,3,4], axis = pgmat.vrnt_axis)
tmp = pgmat.select_vrnt([0,1,2,3,4])

###
### Phased Genotype Matrix Element In-Place-Manipulation
### ====================================================

##
## Appending elements
## ------------------

# append genotype matrices along the taxa axis
tmp = pgmat.deepcopy()                   # copy original
tmp.append(pgmat, axis = tmp.taxa_axis)  # append original to copy

tmp = pgmat.deepcopy()                   # copy original
tmp.append_taxa(pgmat)                   # append original to copy

# append genotype matrices along the variant axis
tmp = pgmat.deepcopy()                   # copy original
tmp.append(pgmat, axis = tmp.vrnt_axis)  # append original to copy

tmp = pgmat.deepcopy()                   # copy original
tmp.append_vrnt(pgmat)                   # append original to copy

##
## Removing elements
## -----------------

#
# ``remove`` taxa
# +++++++++++++++

# remove first taxon using an integer
tmp = pgmat.deepcopy()                           # copy original
tmp.remove(0, axis = pgmat.taxa_axis)            # remove from copy

tmp = pgmat.deepcopy()                           # copy original
tmp.remove_taxa(0)                               # remove from copy

# remove first five taxa using a slice
tmp = pgmat.deepcopy()                           # copy original
tmp.remove(slice(0,5), axis = pgmat.taxa_axis)   # remove from copy

tmp = pgmat.deepcopy()                           # copy original
tmp.remove_taxa(slice(0,5))                      # remove from copy

# remove first five taxa using a Sequence
tmp = pgmat.deepcopy()                           # copy original
tmp.remove([0,1,2,3,4], axis = pgmat.taxa_axis)  # remove from copy

tmp = pgmat.deepcopy()                           # copy original
tmp.remove_taxa([0,1,2,3,4])                     # remove from copy

#
# ``remove`` marker variants
# ++++++++++++++++++++++++++

# remove first marker variant using an integer
tmp = pgmat.deepcopy()                           # copy original
tmp.remove(0, axis = pgmat.vrnt_axis)            # remove from copy

tmp = pgmat.deepcopy()                           # copy original
tmp.remove_vrnt(0)                               # remove from copy

# remove first five marker variants using a slice
tmp = pgmat.deepcopy()                           # copy original
tmp.remove(slice(0,5), axis = pgmat.vrnt_axis)   # remove from copy

tmp = pgmat.deepcopy()                           # copy original
tmp.remove_vrnt(slice(0,5))                      # remove from copy

# remove first five marker variants using a Sequence
tmp = pgmat.deepcopy()                           # copy original
tmp.remove([0,1,2,3,4], axis = pgmat.vrnt_axis)  # remove from copy

tmp = pgmat.deepcopy()                           # copy original
tmp.remove_vrnt([0,1,2,3,4])                     # remove from copy

##
## Incorporating elements
## ----------------------

# incorp genotype matrix along the taxa axis before index 0
tmp = pgmat.deepcopy()                           # copy original
tmp.incorp(0, pgmat, axis = pgmat.taxa_axis)     # incorporate into copy

tmp = pgmat.deepcopy()                           # copy original
tmp.incorp_taxa(0, pgmat)                        # incorporate into copy

# incorp genotype matrix along the variant axis before index 0
tmp = pgmat.deepcopy()                           # copy original
tmp.incorp(0, pgmat, axis = pgmat.vrnt_axis)     # incorporate into copy

tmp = pgmat.deepcopy()                           # copy original
tmp.incorp_vrnt(0, pgmat)                        # incorporate into copy

##
## Concatenating matrices
## ----------------------

# concatenate along the taxa axis
tmp = pgmat.concat([pgmat, pgmat], axis = pgmat.taxa_axis)
tmp = pgmat.concat_taxa([pgmat, pgmat])

# concatenate along the variant axis
tmp = pgmat.concat([pgmat, pgmat], axis = pgmat.vrnt_axis)
tmp = pgmat.concat_vrnt([pgmat, pgmat])

###
### Grouping and Sorting
### ====================

##
## Reordering
## ----------

#
# ``reorder`` taxa
# ++++++++++++++++

# create reordering indices
indices = numpy.arange(pgmat.ntaxa)
numpy.random.shuffle(indices)
tmp = pgmat.deepcopy()

# reorder values along the taxa axis
tmp.reorder(indices, axis = tmp.taxa_axis)
tmp.reorder_taxa(indices)

#
# ``reorder`` marker variants
# +++++++++++++++++++++++++++

# create reordering indices
indices = numpy.arange(pgmat.nvrnt)
numpy.random.shuffle(indices)
tmp = pgmat.deepcopy()

# reorder values along the marker variant axis
tmp = pgmat.deepcopy()
tmp.reorder(indices, axis = tmp.vrnt_axis)
tmp.reorder_vrnt(indices)

##
## Lexsorting
## ----------

#
# ``lexsort`` taxa
# ++++++++++++++++

# create lexsort keys for taxa
key1 = numpy.random.randint(0, 10, pgmat.ntaxa)
key2 = numpy.arange(pgmat.ntaxa)
numpy.random.shuffle(key2)

# lexsort along the taxa axis
pgmat.lexsort((key2,key1), axis = pgmat.taxa_axis)
pgmat.lexsort_taxa((key2,key1))

#
# ``lexsort`` marker variants
# +++++++++++++++++++++++++++

# create lexsort keys for marker variants
key1 = numpy.random.randint(0, 10, pgmat.nvrnt)
key2 = numpy.arange(pgmat.nvrnt)
numpy.random.shuffle(key2)

# lexsort along the marker variant axis
pgmat.lexsort((key2,key1), axis = pgmat.vrnt_axis)
pgmat.lexsort_vrnt((key2,key1))

##
## Sorting
## -------

#
# ``sort`` taxa
# +++++++++++++

# sort along taxa axis
tmp = pgmat.deepcopy()
tmp.sort(axis = tmp.taxa_axis)
tmp.sort_taxa()

#
# ``sort`` marker variants
# ++++++++++++++++++++++++

# sort along marker variant axis
tmp = pgmat.deepcopy()
tmp.sort(axis = tmp.vrnt_axis)
tmp.sort_vrnt()

##
## Grouping
## --------

#
# ``group`` taxa
# ++++++++++++++

# sort along taxa axis
tmp = pgmat.deepcopy()
tmp.group(axis = tmp.taxa_axis)
tmp.group_taxa()
# determine whether grouping has occurred along the taxa axis
tmp.is_grouped(axis = tmp.taxa_axis)
tmp.is_grouped_taxa()

#
# ``group`` marker variants
# +++++++++++++++++++++++++

# sort along vrnt axis
tmp = pgmat.deepcopy()
tmp.group(axis = tmp.vrnt_axis)
tmp.group_vrnt()
# determine whether grouping has occurred along the vrnt axis
tmp.is_grouped(axis = tmp.vrnt_axis)
tmp.is_grouped_vrnt()

###
### Summary Statistics
### ==================

# count the number of major alleles across all taxa
out = pgmat.acount()
out = pgmat.acount(dtype = "int32")

# calculate the allele frequency across all taxa
out = pgmat.afreq()
out = pgmat.afreq(dtype = "float32")

# calculate whether a locus is polymorphic across all taxa 
out = pgmat.apoly()
out = pgmat.apoly(dtype = int)

# count the number of genotypes across all taxa
out = pgmat.gtcount()
out = pgmat.gtcount(dtype = "int32")

# calculate the genotype frequency across all taxa
out = pgmat.gtfreq()
out = pgmat.gtfreq(dtype = "float32")

# calculate the minor allele frequency across all taxa
out = pgmat.maf()
out = pgmat.maf(dtype = "float32")

# calculate the mean expected heterozygosity for the population
out = pgmat.meh()
out = pgmat.meh(dtype = "float32")

# count the number of major alleles individually within taxa
out = pgmat.tacount()
out = pgmat.tacount(dtype = "int32")

# calculate the allele frequency individually within taxa
out = pgmat.tafreq()
out = pgmat.tafreq(dtype = "float32")

###
### Saving Genotype Matrices
### ========================

#
# Write to HDF5
# -------------

# remove exported file if it exists
if os.path.exists("saved_genotypes.h5"):
    os.remove("saved_genotypes.h5")

# write a breeding value matrix to an HDF5 file
pgmat.to_hdf5("saved_genotypes.h5")


pgmat.interp_genpos
pgmat.interp_xoprob
pgmat.mat_asformat
