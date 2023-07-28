#!/usr/bin/env python3

# import the GenotypeMatrix class (an abstract interface class)
import copy
import numpy
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

# import the DenseGenotypeMatrix class (a concrete implemented class)
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

###
### Genotype Matrix Object Creation
###

#
# construct from NumPy
#

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

#
# read from VCF
#

# read a genotype matrix from VCF file
gmat = DenseGenotypeMatrix.from_vcf("widiv_2000SNPs.vcf.gz")

#
# read from HDF5
#

# read a genotype matrix from HDF5 file
gmat = DenseGenotypeMatrix.from_hdf5("widiv_2000SNPs.h5")

###
### Genotype matrix general properties
###

tmp = gmat.mat          # gain access to raw genotype matrix pointer
tmp = gmat.mat_ndim     # get the number of dimensions for the genotype matrix
tmp = gmat.mat_shape    # get genotype matrix shape
tmp = gmat.mat_format   # get genotype matrix format
tmp = gmat.ploidy       # get the ploidy of the taxa represented by the genotype matrix
tmp = gmat.nphase       # get the number of chromosome phases represented by the genotype matrix
tmp = gmat.ntaxa        # get the number of taxa represented by the genotype matrix
tmp = gmat.nvrnt        # get the number of genotype variants represented by the genotype matrix

###
### Genotype matrix taxa properties
###

tmp = gmat.taxa             # get the names of the taxa
tmp = gmat.taxa_axis        # get the matrix axis along which taxa are stored
tmp = gmat.taxa_grp         # get an optional taxa group label
tmp = gmat.taxa_grp_name    # if taxa are sorted by group: get the names of the groups
tmp = gmat.taxa_grp_stix    # if taxa are sorted by group: get the start indices (inclusive) for each group
tmp = gmat.taxa_grp_spix    # if taxa are sorted by group: get the stop indices (exclusive) for each group
tmp = gmat.taxa_grp_len     # if taxa are sorted by group: get the length of each group

###
### Genotype matrix marker variant properties
###

tmp = gmat.vrnt_name        # get the names of the marker variants
tmp = gmat.vrnt_axis        # get the axis along which marker variants are stored
tmp = gmat.vrnt_chrgrp      # get the chromosome to which a marker variant belongs
tmp = gmat.vrnt_phypos      # get the physical position of a marker variant
tmp = gmat.vrnt_genpos      # get the genetic position of a marker variant
tmp = gmat.vrnt_xoprob      # get the crossover probability between the previous marker
tmp = gmat.vrnt_hapref      # get the reference haplotype for the marker variant
tmp = gmat.vrnt_hapalt      # get the alternative haplotype for the marker variant
tmp = gmat.vrnt_hapgrp      # get the haplotype grouping for the marker variant
tmp = gmat.vrnt_mask        # get a mask associated with the marker variants
tmp = gmat.vrnt_chrgrp_name # get the names of the chromosomes if sorted
tmp = gmat.vrnt_chrgrp_stix # get the start indices (inclusive) for each chromosome if sorted
tmp = gmat.vrnt_chrgrp_spix # get the stop indices (exclusive) for each chromosome if sorted
tmp = gmat.vrnt_chrgrp_len  # get the length of each chromosome if sorted

###
### Copying
###

# copy a genotype matrix
tmp = copy.copy(gmat)
tmp = gmat.copy()

# deep copy a genotype matrix
tmp = copy.deepcopy(gmat)
tmp = gmat.deepcopy()

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

# count the number of major alleles across all taxa
out = gmat.acount()
out = gmat.acount(dtype = "int32")

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
