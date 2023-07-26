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

# get the crossover probability between the 
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

###
### Summary Statistics
###

# count the number of major alleles across all taxa
gmat.acount()




gmat.acount
gmat.adjoin
gmat.adjoin_taxa
gmat.adjoin_vrnt
gmat.afreq
gmat.apoly
gmat.append
gmat.append_taxa
gmat.append_vrnt
gmat.assign_hapgrp
gmat.concat
gmat.concat_taxa
gmat.concat_vrnt
gmat.delete
gmat.delete_taxa
gmat.delete_vrnt
gmat.from_hdf5
gmat.from_vcf
gmat.group
gmat.group_taxa
gmat.group_vrnt
gmat.gtcount
gmat.gtfreq
gmat.incorp
gmat.incorp_taxa
gmat.incorp_vrnt
gmat.insert
gmat.insert_taxa
gmat.insert_vrnt
gmat.interp_genpos
gmat.interp_xoprob
gmat.is_grouped
gmat.is_grouped_taxa
gmat.is_grouped_vrnt
gmat.lexsort
gmat.lexsort_taxa
gmat.lexsort_vrnt
gmat.maf
gmat.mat_asformat
gmat.meh
gmat.remove
gmat.remove_taxa
gmat.remove_vrnt
gmat.reorder
gmat.reorder_taxa
gmat.reorder_vrnt
gmat.select
gmat.select_taxa
gmat.select_vrnt
gmat.sort
gmat.sort_taxa
gmat.sort_vrnt
gmat.tacount
gmat.tafreq
gmat.to_hdf5
