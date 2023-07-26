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

Names of the taxa in a genotype matrix
--------------------------------------

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

.. code-block:: python

    # get the names of the marker variants
    gmat_vrnt_name = gmat.vrnt_name

.. code-block:: python

    # get the axis along which marker variants are stored
    gmat_vrnt_axis = gmat.vrnt_axis

.. code-block:: python

    # get the chromosome to which a marker variant belongs
    gmat_vrnt_chrgrp = gmat.vrnt_chrgrp

.. code-block:: python

    # get the physical position of a marker variant
    gmat_vrnt_phypos = gmat.vrnt_phypos

.. code-block:: python

    # get the genetic position of a marker variant
    gmat_vrnt_genpos = gmat.vrnt_genpos

.. code-block:: python

    # get the crossover probability between the 
    gmat_vrnt_xoprob = gmat.vrnt_xoprob

.. code-block:: python

    # get the reference haplotype for the marker variant
    gmat_vrnt_hapref = gmat.vrnt_hapref

.. code-block:: python

    # get the alternative haplotype for the marker variant
    gmat_vrnt_hapalt = gmat.vrnt_hapalt

.. code-block:: python

    # get the haplotype grouping for the marker variant
    gmat_vrnt_hapgrp = gmat.vrnt_hapgrp

.. code-block:: python

    # get a mask associated with the marker variants
    gmat_vrnt_mask = gmat.vrnt_mask

.. code-block:: python

    # if marker variants are sorted by chromosome: 
    # get the names of the chromosomes
    gmat_vrnt_chrgrp_name = gmat.vrnt_chrgrp_name

.. code-block:: python

    # if marker variants are sorted by chromosome: 
    # get the start indices (inclusive) for each chromosome
    gmat_vrnt_chrgrp_stix = gmat.vrnt_chrgrp_stix

.. code-block:: python

    # if marker variants are sorted by chromosome: 
    # get the stop indices (exclusive) for each chromosome
    gmat_vrnt_chrgrp_spix = gmat.vrnt_chrgrp_spix

.. code-block:: python

    # if marker variants are sorted by chromosome: 
    # get the length of each chromosome
    gmat_vrnt_chrgrp_len = gmat.vrnt_chrgrp_len


Copying Genotype Matrices
=========================

.. code-block:: python

    # copy a genotype matrix
    gmat_copy = gmat.copy()

    # deep copy a genotype matrix
    gmat_deepcopy = gmat.deepcopy()

Genotype Matrix Element Copy-On-Manipulation
============================================

Genotype Matrix Element In-Place-Manipulation
=============================================

Summary Statistics
==================

Allele Counts
-------------




