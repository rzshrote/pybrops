#!/usr/bin/env python3

# import the GenotypeMatrix class (an abstract interface class)
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

# import the DenseGenotypeMatrix class (a concrete implemented class)
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

###
### Genotype Matrix Object Creation
###

# read a genotype matrix from file
gmat = DenseGenotypeMatrix.from_vcf("widiv_2000SNPs.vcf.gz")

###
### Genotype matrix general properties
###

# gain access to raw genotype matrix pointer
gmat_mat_ptr = gmat.mat

# get the number of dimensions for the genotype matrix
gmat_ndim = gmat.mat_ndim

# get genotype matrix shape
gmat_shape = gmat.mat_shape

# get genotype matrix format
gmat_format = gmat.mat_format

# get the ploidy of the taxa represented by the genotype matrix
gmat_ploidy = gmat.ploidy

# get the number of chromosome phases represented by the genotype matrix
gmat_nphase = gmat.nphase

# get the number of taxa represented by the genotype matrix
gmat_ntaxa = gmat.ntaxa

# get the number of genotype variants represented by the genotype matrix
gmat_nvrnt = gmat.nvrnt

###
### Genotype matrix taxa properties
###

# get the names of the taxa
gmat_taxa = gmat.taxa

# get the matrix axis along which taxa are stored
gmat_taxa_axis = gmat.taxa_axis

# get an optional taxa group label
gmat_taxa_grp = gmat.taxa_grp

# if taxa are sorted by group: get the names of the groups
gmat_taxa_grp_name = gmat.taxa_grp_name

# if taxa are sorted by group: get the start indices (inclusive) for each group
gmat_taxa_grp_stix = gmat.taxa_grp_stix

# if taxa are sorted by group: get the stop indices (exclusive) for each group
gmat_taxa_grp_spix = gmat.taxa_grp_spix

# if taxa are sorted by group: get the length of each group
gmat_taxa_grp_len = gmat.taxa_grp_len

###
### Genotype matrix marker variant properties
###

# get the names of the marker variants
gmat_vrnt_name = gmat.vrnt_name

# get the axis along which marker variants are stored
gmat_vrnt_axis = gmat.vrnt_axis

# get the chromosome to which a marker variant belongs
gmat_vrnt_chrgrp = gmat.vrnt_chrgrp

# get the physical position of a marker variant
gmat_vrnt_phypos = gmat.vrnt_phypos

# get the genetic position of a marker variant
gmat_vrnt_genpos = gmat.vrnt_genpos

# get the crossover probability between the previous marker
gmat_vrnt_xoprob = gmat.vrnt_xoprob

# get the reference haplotype for the marker variant
gmat_vrnt_hapref = gmat.vrnt_hapref

# get the alternative haplotype for the marker variant
gmat_vrnt_hapalt = gmat.vrnt_hapalt

# get the haplotype grouping for the marker variant
gmat_vrnt_hapgrp = gmat.vrnt_hapgrp

# get a mask associated with the marker variants
gmat_vrnt_mask = gmat.vrnt_mask

# if marker variants are sorted by chromosome: 
# get the names of the chromosomes
gmat_vrnt_chrgrp_name = gmat.vrnt_chrgrp_name

# if marker variants are sorted by chromosome: 
# get the start indices (inclusive) for each chromosome
gmat_vrnt_chrgrp_stix = gmat.vrnt_chrgrp_stix

# if marker variants are sorted by chromosome: 
# get the stop indices (exclusive) for each chromosome
gmat_vrnt_chrgrp_spix = gmat.vrnt_chrgrp_spix

# if marker variants are sorted by chromosome: 
# get the length of each chromosome
gmat_vrnt_chrgrp_len = gmat.vrnt_chrgrp_len

###
### Copying
###

# copy a genotype matrix
gmat_copy = gmat.copy()

# deep copy a genotype matrix
gmat_deepcopy = gmat.deepcopy()

###
### Genotype Matrix Element Copy-On-Manipulation
###

##
## adjoin examples
##

# create a new genotype matrix to demonstrate
new = gmat.deepcopy()

# adjoin genotype matrices along the taxa axis
tmp = gmat.adjoin(new, axis = gmat.taxa_axis)
tmp = gmat.adjoin_taxa(new)

# adjoin genotype matrices along the variant axis
tmp = gmat.adjoin(new, axis = gmat.vrnt_axis)
tmp = gmat.adjoin_vrnt(new)

##
## delete examples
##

#
# delete taxa examples
#

# delete first taxon using an integer
tmp = gmat.delete(0, axis = gmat.taxa_axis)
tmp = gmat.delete_taxa(0)

# delete first five taxa using a slice
tmp = gmat.delete(slice(0,5), axis = gmat.taxa_axis)
tmp = gmat.delete_taxa(slice(0,5))

# delete first five taxa using a Sequence
tmp = gmat.delete([0,1,2,3,4], axis = gmat.taxa_axis)
tmp = gmat.delete_taxa([0,1,2,3,4])

#
# delete marker variants examples
#

# delete first marker variant using an integer
tmp = gmat.delete(0, axis = gmat.vrnt_axis)
tmp = gmat.delete_vrnt(0)

# delete first five marker variants using a slice
tmp = gmat.delete(slice(0,5), axis = gmat.vrnt_axis)
tmp = gmat.delete_vrnt(slice(0,5))

# delete first five marker variants using a Sequence
tmp = gmat.delete([0,1,2,3,4], axis = gmat.vrnt_axis)
tmp = gmat.delete_vrnt([0,1,2,3,4])

##
## insert examples
##

# create a new genotype matrix to demonstrate
new = gmat.deepcopy()

# insert genotype matrix along the taxa axis before index 0
tmp = gmat.insert(0, new, axis = gmat.taxa_axis)
tmp = gmat.insert_taxa(0, new)

# insert genotype matrix along the variant axis before index 0
# tmp = gmat.insert(0, new, axis = gmat.vrnt_axis)
# tmp = gmat.insert_vrnt(0, new)

##
## select examples
##

# select first five taxa using a Sequence
tmp = gmat.select([0,1,2,3,4], axis = gmat.taxa_axis)
tmp = gmat.select_taxa([0,1,2,3,4])

# select first five marker variants using a Sequence
tmp = gmat.select([0,1,2,3,4], axis = gmat.vrnt_axis)
tmp = gmat.select_vrnt([0,1,2,3,4])

###
### Genotype Matrix Element In-Place-Manipulation
###

##
## append examples
##

# append genotype matrices along the taxa axis
tmp = gmat.deepcopy()                   # copy original
tmp.append(gmat, axis = tmp.taxa_axis)  # append original to copy

tmp = gmat.deepcopy()                   # copy original
tmp.append_taxa(gmat)                   # append original to copy

# append genotype matrices along the variant axis
tmp = gmat.deepcopy()                   # copy original
tmp.append(gmat, axis = tmp.vrnt_axis)  # append original to copy

tmp = gmat.deepcopy()                   # copy original
tmp.append_vrnt(gmat)                   # append original to copy

##
## remove examples
##

#
# remove taxa examples
#

# remove first taxon using an integer
tmp = gmat.deepcopy()                           # copy original
tmp.remove(0, axis = gmat.taxa_axis)            # remove from copy

tmp = gmat.deepcopy()                           # copy original
tmp.remove_taxa(0)                              # remove from copy

# remove first five taxa using a slice
tmp = gmat.deepcopy()                           # copy original
tmp.remove(slice(0,5), axis = gmat.taxa_axis)   # remove from copy

tmp = gmat.deepcopy()                           # copy original
tmp.remove_taxa(slice(0,5))                     # remove from copy

# remove first five taxa using a Sequence
tmp = gmat.deepcopy()                           # copy original
tmp.remove([0,1,2,3,4], axis = gmat.taxa_axis)  # remove from copy

tmp = gmat.deepcopy()                           # copy original
tmp.remove_taxa([0,1,2,3,4])                    # remove from copy

#
# remove marker variants examples
#

# remove first marker variant using an integer
tmp = gmat.deepcopy()                           # copy original
tmp.remove(0, axis = gmat.vrnt_axis)            # remove from copy

tmp = gmat.deepcopy()                           # copy original
tmp.remove_vrnt(0)                              # remove from copy

# remove first five marker variants using a slice
tmp = gmat.deepcopy()                           # copy original
tmp.remove(slice(0,5), axis = gmat.vrnt_axis)   # remove from copy

tmp = gmat.deepcopy()                           # copy original
tmp.remove_vrnt(slice(0,5))                     # remove from copy

# remove first five marker variants using a Sequence
tmp = gmat.deepcopy()                           # copy original
tmp.remove([0,1,2,3,4], axis = gmat.vrnt_axis)  # remove from copy

tmp = gmat.deepcopy()                           # copy original
tmp.remove_vrnt([0,1,2,3,4])                    # remove from copy

##
## incorp examples
##

# incorp genotype matrix along the taxa axis before index 0
tmp = gmat.deepcopy()                           # copy original
tmp.incorp(0, gmat, axis = gmat.taxa_axis)      # incorporate into copy

tmp = gmat.deepcopy()                           # copy original
tmp.incorp_taxa(0, gmat)                        # incorporate into copy

# incorp genotype matrix along the variant axis before index 0
# tmp = gmat.deepcopy()                           # copy original
# tmp.incorp(0, gmat, axis = gmat.vrnt_axis)      # incorporate into copy

# tmp = gmat.deepcopy()                           # copy original
# tmp.incorp_vrnt(0, gmat)                        # incorporate into copy

##
## concat examples
##

# concatenate along the taxa axis
tmp = gmat.concat([gmat, gmat], axis = gmat.taxa_axis)
tmp = gmat.concat_taxa([gmat, gmat])

# concatenate along the variant axis
tmp = gmat.concat([gmat, gmat], axis = gmat.vrnt_axis)
tmp = gmat.concat_vrnt([gmat, gmat])

###
### Summary Statistics
###

#
#
#

# count the number of major alleles across all taxa
out = gmat.acount()
out = gmat.acount(out = "int32")

# calculate the allele frequency across all taxa
out = gmat.afreq()
out = gmat.afreq(dtype = "float32")

# calculate whether a locus is polymorphic across all taxa 
out = gmat.apoly()
out = gmat.apoly(dtype = int)

# count the number of genotypes across all taxa
out = gmat.gtcount()
out = gmat.gtcount(dtype = "int32")

# calculate the genotype frequency across all taxa
out = gmat.gtfreq()
out = gmat.gtfreq(dtype = "float32")

# calculate the minor allele frequency across all taxa
out = gmat.maf()
out = gmat.maf(dtype = "float32")

# calculate the mean expected heterozygosity for the population
out = gmat.meh()
out = gmat.meh(dtype = "float32")

# count the number of major alleles individually within taxa
out = gmat.tacount()
out = gmat.tacount(dtype = "int32")

# calculate the allele frequency individually within taxa
out = gmat.tafreq()
out = gmat.tafreq(dtype = "float32")



gmat.assign_hapgrp
gmat.from_hdf5
gmat.from_vcf
gmat.group
gmat.group_taxa
gmat.group_vrnt
gmat.interp_genpos
gmat.interp_xoprob
gmat.is_grouped
gmat.is_grouped_taxa
gmat.is_grouped_vrnt
gmat.lexsort
gmat.lexsort_taxa
gmat.lexsort_vrnt
gmat.mat_asformat
gmat.reorder
gmat.reorder_taxa
gmat.reorder_vrnt
gmat.sort
gmat.sort_taxa
gmat.sort_vrnt
gmat.to_hdf5
