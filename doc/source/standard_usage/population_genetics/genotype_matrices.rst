Genotype Matrices
#################

Loading Genotype Matrix Modules
===============================

Genotype matrix support in PyBrOpS is found in the ``pybrops.popgen.gmat`` submodule. The ``GenotypeMatrix`` abstract class is the basal interface for all genotype matrices. Classes can be imported as follows:

.. code-block:: python

    # import the GenotypeMatrix class (an abstract interface class)
    from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

    # import the DenseGenotypeMatrix class (a concrete implemented class)
    from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

Loading Genotype Matrices from VCF Files
========================================

VCF files can be loaded into the Python environment as genotype matrices using the ``from_vcf`` method.

.. code-block:: python

    # read a genotype matrix from file
    gmat = DenseGenotypeMatrix.from_vcf("widiv_2000SNPs.vcf.gz")

Genotype Matrix General Properties
==================================

Access to raw matrix pointer
----------------------------

.. code-block:: python

    # gain access to raw genotype matrix pointer
    gmat_mat_ptr = gmat.mat

Number of genotype matrix dimensions
------------------------------------

.. code-block:: python

    # get the number of dimensions for the genotype matrix
    gmat_ndim = gmat.mat_ndim

Shape of the genotype matrix
----------------------------

.. code-block:: python

    # get genotype matrix shape
    gmat_shape = gmat.mat_shape

Format of the genotype matrix
-----------------------------

.. code-block:: python

    # get genotype matrix format
    gmat_format = gmat.mat_format

Ploidy of taxa represented in a genotype matrix
-----------------------------------------------

.. code-block:: python

    # get the ploidy of the taxa represented by the genotype matrix
    gmat_ploidy = gmat.ploidy


Number of chromosome phases in a genotype matrix
------------------------------------------------

.. code-block:: python

    # get the number of chromosome phases represented by the genotype matrix
    gmat_nphase = gmat.nphase

Number of taxa in a genotype matrix
-----------------------------------

.. code-block:: python

    # get the number of taxa represented by the genotype matrix
    gmat_ntaxa = gmat.ntaxa

Number of genetic markers in a genotype matrix
----------------------------------------------

.. code-block:: python

    # get the number of genotype variants represented by the genotype matrix
    gmat_nvrnt = gmat.nvrnt

Genotype Matrix Taxa Properties
===============================

Taxa names
----------

.. code-block:: python

    # get the names of the taxa
    gmat_taxa = gmat.taxa

Axis along which taxa are stored
--------------------------------

.. code-block:: python

    # get the matrix axis along which taxa are stored
    gmat_taxa_axis = gmat.taxa_axis

Taxa group labels
-----------------

.. code-block:: python

    # get an optional taxa group label
    gmat_taxa_grp = gmat.taxa_grp

Unique taxa group label names
-----------------------------

.. code-block:: python

    # if taxa are sorted by group: get the names of the groups
    gmat_taxa_grp_name = gmat.taxa_grp_name

Taxa group start indices
------------------------

.. code-block:: python

    # if taxa are sorted by group: get the start indices (inclusive) for each group
    gmat_taxa_grp_stix = gmat.taxa_grp_stix

Taxa group stop indices
-----------------------

.. code-block:: python

    # if taxa are sorted by group: get the stop indices (exclusive) for each group
    gmat_taxa_grp_spix = gmat.taxa_grp_spix

Taxa group length
-----------------

.. code-block:: python

    # if taxa are sorted by group: get the length of each group
    gmat_taxa_grp_len = gmat.taxa_grp_len

Genotype Matrix Marker Variant Properties
=========================================

Marker variant names
--------------------

.. code-block:: python

    # get the names of the marker variants
    gmat_vrnt_name = gmat.vrnt_name

Axis along which marker variants are stored
-------------------------------------------

.. code-block:: python

    # get the axis along which marker variants are stored
    gmat_vrnt_axis = gmat.vrnt_axis

Marker variant chromosome group designations
--------------------------------------------

.. code-block:: python

    # get the chromosome to which a marker variant belongs
    gmat_vrnt_chrgrp = gmat.vrnt_chrgrp

Marker variant chromosome physical positions
--------------------------------------------

.. code-block:: python

    # get the physical position of a marker variant
    gmat_vrnt_phypos = gmat.vrnt_phypos

Marker variant chromosome genetic map positions
-----------------------------------------------

.. code-block:: python

    # get the genetic position of a marker variant
    gmat_vrnt_genpos = gmat.vrnt_genpos

Marker variant sequential crossover probabilities
-------------------------------------------------

.. code-block:: python

    # get the crossover probability between the previous marker
    gmat_vrnt_xoprob = gmat.vrnt_xoprob

Marker variant reference haplotype alleles
------------------------------------------

.. code-block:: python

    # get the reference haplotype for the marker variant
    gmat_vrnt_hapref = gmat.vrnt_hapref

Marker variant alternative haplotype alleles
--------------------------------------------

.. code-block:: python

    # get the alternative haplotype for the marker variant
    gmat_vrnt_hapalt = gmat.vrnt_hapalt

Marker variant haplotype group designations
-------------------------------------------

.. code-block:: python

    # get the haplotype grouping for the marker variant
    gmat_vrnt_hapgrp = gmat.vrnt_hapgrp

Marker variant mask
-------------------

.. code-block:: python

    # get a mask associated with the marker variants
    gmat_vrnt_mask = gmat.vrnt_mask

Unique marker variant group label names
---------------------------------------

.. code-block:: python

    # if marker variants are sorted by chromosome: 
    # get the names of the chromosomes
    gmat_vrnt_chrgrp_name = gmat.vrnt_chrgrp_name

Marker variant group start indices
----------------------------------

.. code-block:: python

    # if marker variants are sorted by chromosome: 
    # get the start indices (inclusive) for each chromosome
    gmat_vrnt_chrgrp_stix = gmat.vrnt_chrgrp_stix

Marker variant group stop indices
---------------------------------

.. code-block:: python

    # if marker variants are sorted by chromosome: 
    # get the stop indices (exclusive) for each chromosome
    gmat_vrnt_chrgrp_spix = gmat.vrnt_chrgrp_spix

Marker variant group length
---------------------------

.. code-block:: python

    # if marker variants are sorted by chromosome: 
    # get the length of each chromosome
    gmat_vrnt_chrgrp_len = gmat.vrnt_chrgrp_len


Copying Genotype Matrices
=========================

Shallow copying
---------------

.. code-block:: python

    # copy a genotype matrix
    gmat_copy = gmat.copy()

Deep copying
------------

.. code-block:: python

    # deep copy a genotype matrix
    gmat_deepcopy = gmat.deepcopy()

Genotype Matrix Element Copy-On-Manipulation
============================================

Adjoin elements
---------------

.. code-block:: python

    # create a new genotype matrix to demonstrate
    new = gmat.deepcopy()

    # adjoin genotype matrices along the taxa axis
    tmp = gmat.adjoin(new, axis = gmat.taxa_axis)
    tmp = gmat.adjoin_taxa(new)

    # adjoin genotype matrices along the variant axis
    tmp = gmat.adjoin(new, axis = gmat.vrnt_axis)
    tmp = gmat.adjoin_vrnt(new)

Delete elements
---------------

.. code-block:: python

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

Insert elements
---------------

.. code-block:: python

    # create a new genotype matrix to demonstrate
    new = gmat.deepcopy()

    # insert genotype matrix along the taxa axis before index 0
    tmp = gmat.insert(0, new, axis = gmat.taxa_axis)
    tmp = gmat.insert_taxa(0, new)

    # insert genotype matrix along the variant axis before index 0
    tmp = gmat.insert(0, new, axis = gmat.vrnt_axis)
    tmp = gmat.insert_vrnt(0, new)

Select elements
---------------

.. code-block:: python

    # select first five taxa using a Sequence
    tmp = gmat.select([0,1,2,3,4], axis = gmat.taxa_axis)
    tmp = gmat.select_taxa([0,1,2,3,4])

    # select first five marker variants using a Sequence
    tmp = gmat.select([0,1,2,3,4], axis = gmat.vrnt_axis)
    tmp = gmat.select_vrnt([0,1,2,3,4])

Genotype Matrix Element In-Place-Manipulation
=============================================

Append elements
---------------

.. code-block:: python

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

Remove elements
---------------

.. code-block:: python

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

Incorporate elements
--------------------

.. code-block:: python

    # incorp genotype matrix along the taxa axis before index 0
    tmp = gmat.deepcopy()                           # copy original
    tmp.incorp(0, gmat, axis = gmat.taxa_axis)      # incorporate into copy

    tmp = gmat.deepcopy()                           # copy original
    tmp.incorp_taxa(0, gmat)                        # incorporate into copy

    # incorp genotype matrix along the variant axis before index 0
    tmp = gmat.deepcopy()                           # copy original
    tmp.incorp(0, gmat, axis = gmat.vrnt_axis)      # incorporate into copy

    tmp = gmat.deepcopy()                           # copy original
    tmp.incorp_vrnt(0, gmat)                        # incorporate into copy

Concatenate elements
--------------------

.. code-block:: python

    # concatenate along the taxa axis
    tmp = gmat.concat([gmat, gmat], axis = gmat.taxa_axis)
    tmp = gmat.concat_taxa([gmat, gmat])

    # concatenate along the variant axis
    tmp = gmat.concat([gmat, gmat], axis = gmat.vrnt_axis)
    tmp = gmat.concat_vrnt([gmat, gmat])


Summary Statistics
==================

Population allele counts
------------------------

.. code-block:: python

    # count the number of major alleles across all taxa
    out = gmat.acount()
    out = gmat.acount(out = "int32")

Population allele frequencies
-----------------------------

.. code-block:: python

    # calculate the allele frequency across all taxa
    out = gmat.afreq()
    out = gmat.afreq(dtype = "float32")

Population allele polymorphism presence
---------------------------------------

.. code-block:: python

    # calculate whether a locus is polymorphic across all taxa 
    out = gmat.apoly()
    out = gmat.apoly(dtype = int)

Population genotype counts
--------------------------

.. code-block:: python

    # count the number of genotypes across all taxa
    out = gmat.gtcount()
    out = gmat.gtcount(dtype = "int32")

Population genotype frequencies
-------------------------------

.. code-block:: python

    # calculate the genotype frequency across all taxa
    out = gmat.gtfreq()
    out = gmat.gtfreq(dtype = "float32")

Population minor allele frequencies
-----------------------------------

.. code-block:: python

    # calculate the minor allele frequency across all taxa
    out = gmat.maf()
    out = gmat.maf(dtype = "float32")

Population mean expected heterozygosity
---------------------------------------

.. code-block:: python

    # calculate the mean expected heterozygosity for the population
    out = gmat.meh()
    out = gmat.meh(dtype = "float32")

Taxa allele counts
------------------

.. code-block:: python

    # count the number of major alleles individually within taxa
    out = gmat.tacount()
    out = gmat.tacount(dtype = "int32")

Taxa allele frequencies
-----------------------

.. code-block:: python

    # calculate the allele frequency individually within taxa
    out = gmat.tafreq()
    out = gmat.tafreq(dtype = "float32")
