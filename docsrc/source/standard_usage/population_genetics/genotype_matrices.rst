Genotype Matrices
#################

Class Family Overview
=====================

Perhaps the most important family of classes in PyBrOpS is the ``GenotypeMatrix`` object family. The purpose of this family of classes is to store and represent genotypic data. ``GenotypeMatrix`` classes can be used in the estimation of genomic prediction models, the estimation of breeding values, the calculation of genomic relationship matrices, and to make selection decisions. Within PyBrOpS it is possible to read genotypic data from VCF files to create ``GenotypeMatrix`` objects, allowing for real-world data to be used in breeding program simulations.

Summary of Genotype Matrix Classes
==================================

Genotype matrix classes in PyBrOpS are found in the ``pybrops.popgen.gmat`` module. Contained in this moduel are several ``GenotypeMatrix`` class type definitions which are summarized in the table below.

.. list-table:: Summary of genotype matrix classes in the ``pybrops.popgen.gmat`` module
    :widths: 25 15 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``GenotypeMatrix``
      - Abstract
      - Interface for all genotype matrix child classes.
    * - ``DenseGenotypeMatrix``
      - Concrete
      - Class representing dense, genotype matrices.

Loading Class Modules
=====================

Genotype matrix classes can be imported as follows:

.. code-block:: python

    # import the GenotypeMatrix class (an abstract interface class)
    from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

    # import the DenseGenotypeMatrix class (a concrete implemented class)
    from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

Creating Genotype Matrices
==========================

Genotype matrices can be created using several methods. Genotype matrices may be constructed from raw NumPy arrays, by reading data from VCF files, or by reading data from an HDF5 file containing a saved genotype matrix.

Creating genotype matrices from NumPy arrays
--------------------------------------------

Genotype matrices can be created from raw NumPy arrays. The example below illustrates creating a ``DenseGenotypeMatrix``. The only requirement to construct a ``DenseGenotypeMatrix`` is a ``(n,p)`` numerical matrix containing genotypic codings. Here, ``n`` is the number of taxa and ``p`` is the number of marker variants. The matrix must be coded in the ``{0,1,2}`` (number of alternative alleles) format and must have an ``int8`` data type.

Additional optional metadata may be stored along with a ``DenseGenotypeMatrix`` including taxa names (``taxa``), taxa groups (``taxa_grp``), marker variant chromosome assignments (``vrnt_chrgrp``), marker variant chromosome physical positions (``vrnt_phypos``), marker variant names (``vrnt_name``), marker variant genetic map positions (``vrnt_genpos``), sequential recombination probabilities between markers (``vrnt_xoprob``), marker variant haplotype group assignments (``vrnt_hapgrp``), reference haplotype (``vrnt_hapref``), alternative haplotype (``vrnt_hapalt``), a variant mask (``vrnt_mask``), and finally the ploidy of the genotypes (``ploid``).

.. code-block:: python

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

Loading genotype matrices from VCF files
----------------------------------------

VCF files can be loaded into Python scopes as genotype matrices using the ``from_vcf`` method.

.. code-block:: python

    # read a genotype matrix from file
    gmat = DenseGenotypeMatrix.from_vcf("widiv_2000SNPs.vcf.gz")

Loading genotype matrices from HDF5 files
-----------------------------------------

Genotype matrices in PyBrOpS can be exported to HDF5 files via the ``to_hdf5`` method. These files can later be read into PyBrOpS via the ``from_hdf5`` method. The example below illustrates loading a ``DenseGenotypeMatrix`` into memory from an HDF5 file:

.. code-block:: python

    # read a genotype matrix from HDF5 file
    gmat = DenseGenotypeMatrix.from_hdf5("widiv_2000SNPs.h5")

Genotype Matrix Properties
==========================

``GenotypeMatrix`` objects share a set of common object properties which can be grouped into three groups: general properties, taxa properties, and marker variant properties.

General properties
------------------

General properties of a genotype matrix primarily relate to the shape and format of the genotype matrix itself. These properties are summarized in the table below.

.. list-table:: Summary of ``GenotypeMatrix`` general properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``mat``
      - Pointer to the raw genotype matrix pointer
    * - ``mat_ndim``
      - The number of dimensions for the genotype matrix
    * - ``mat_shape``
      - Genotype matrix shape
    * - ``mat_format``
      - Genotype matrix format
    * - ``ploidy``
      - The ploidy of the taxa represented by the genotype matrix
    * - ``nphase``
      - The number of chromosome phases represented by the genotype matrix

Taxa properties
---------------

Genotype matrices can store basic information on the taxa represented by the matrix. These data include taxa names and any grouping labels. Grouping labels may be used to indicate a family of individuals. Taxa and taxa grouping data can be accessed using the object properties summarized in the table below.

.. list-table:: Summary of ``GenotypeMatrix`` taxa properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ntaxa``
      - The number of taxa represented by the genotype matrix
    * - ``taxa``
      - The names of the taxa
    * - ``taxa_axis``
      - The matrix axis along which taxa are stored
    * - ``taxa_grp``
      - Taxa group label
    * - ``taxa_grp_name``
      - The names of the taxa groups
    * - ``taxa_grp_stix``
      - The start indices (inclusive) for each taxa group, post sorting and grouping
    * - ``taxa_grp_spix``
      - The stop indices (exclusive) for each taxa group, post sorting and grouping
    * - ``taxa_grp_len``
      - The length of each taxa group, post sorting and grouping

Marker variant properties
-------------------------

Genotype matrices can also store basic information about the genetic markers which are represented by the matrix. Key information about a marker variant's name, chromosome assignment, physical and genetic positions, and sequential crossover probabilities can be stored. The table below summarizes marker variant properties in genotype matrices.

.. list-table:: Summary of ``GenotypeMatrix`` marker variant properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nvrnt``
      - The number of genotype variants represented by the genotype matrix
    * - ``vrnt_name``
      - The names of the marker variants
    * - ``vrnt_axis``
      - The axis along which marker variants are stored
    * - ``vrnt_chrgrp``
      - The chromosome to which a marker variant belongs
    * - ``vrnt_phypos``
      - The physical position of a marker variant
    * - ``vrnt_genpos``
      - The genetic position of a marker variant
    * - ``vrnt_xoprob``
      - The crossover probability between the previous marker
    * - ``vrnt_hapref``
      - The reference haplotype for the marker variant
    * - ``vrnt_hapalt``
      - The alternative haplotype for the marker variant
    * - ``vrnt_hapgrp``
      - The haplotype grouping for the marker variant
    * - ``vrnt_mask``
      - A mask associated with the marker variants
    * - ``vrnt_chrgrp_name``
      - The names of the chromosomes
    * - ``vrnt_chrgrp_stix``
      - The start indices (inclusive) for each chromosome, post sorting and grouping
    * - ``vrnt_chrgrp_spix``
      - The stop indices (exclusive) for each chromosome, post sorting and grouping
    * - ``vrnt_chrgrp_len``
      - The length of each chromosome, post sorting and grouping

Copying Genotype Matrices
=========================

Genotype matrices can be copied. There are two methods that can be used to copy a genotype matrix: shallow copying and deep copying.

Shallow copying
---------------

In shallow copying, only references to a ``GenotypeMatrix``'s data are copied to a new genotype matrix. Copying is only one level deep so any changes to the data in the original object will be reflected in the copied object.

.. code-block:: python

    # copy a genotype matrix
    tmp = copy.copy(gmat)
    tmp = gmat.copy()

Deep copying
------------

In deep copying, a ``GenotypeMatrix``'s data is recursively copied to the deepest level. Changes to data values in the original genotype matrix will not affect the data values in the copy. 

.. code-block:: python

    # deep copy a genotype matrix
    tmp = copy.deepcopy(gmat)
    tmp = gmat.deepcopy()

Copy-On Element Manipulation
============================

Genotype matrices have several methods which create a modified copy of the original matrix, leaving the original genotype matrix unaltered.

Adjoin elements
---------------

The ``adjoin`` family of methods allows for a genotype matrix to be adjoined to another genotype matrix, creating a new matrix in the process. Use of the ``adjoin`` method family is demonstrated in the code below.

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

The ``delete`` family of methods allows for rows (taxa) and columns (marker variants) of the genotype matrix to be removed in a copy of the original. Use of the ``delete`` method family is demonstrated in the code below.

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

The ``insert`` family of methods allows for rows (taxa) and columns (marker variants) of the genotype matrix to be inserted into a copy of the original matrix. Use of the ``insert`` method family is demonstrated in the code below.

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

The ``select`` family of methods allows for rows (taxa) and columns (marker variants) of the genotype matrix to be extracted to a copy of the original matrix. The ``select`` family is the inverse of the ``delete`` family. Use of the ``select`` method family is demonstrated in the code below.

.. code-block:: python

    # select first five taxa using a Sequence
    tmp = gmat.select([0,1,2,3,4], axis = gmat.taxa_axis)
    tmp = gmat.select_taxa([0,1,2,3,4])

    # select first five marker variants using a Sequence
    tmp = gmat.select([0,1,2,3,4], axis = gmat.vrnt_axis)
    tmp = gmat.select_vrnt([0,1,2,3,4])

In-Place Element Manipulation
=============================

Genotype matrices have several methods which create a modify a matrix in-place.

Append elements
---------------

The ``append`` family of methods allows for new rows (taxa) and columns (marker variants) to be appended to the genotype matrix. The code segment below demonstrates their use. 

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

The ``remove`` family of methods allows for rows (taxa) and columns (marker variants) to be removed from a genotype matrix. A demonstration of their use can be seen below. 

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

The ``incorp`` family of methods allows for new rows (taxa) and columns (marker variants) to be inserted at specific locations a genotype matrix. Use of the ``incorp`` family is demonstrated in the code segment below below. 

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

The ``concat`` family of methods allows for multiple genotype matrices to be concatenated to each other. The code segment below demonstrates their use. 

.. code-block:: python

    # concatenate along the taxa axis
    tmp = gmat.concat([gmat, gmat], axis = gmat.taxa_axis)
    tmp = gmat.concat_taxa([gmat, gmat])

    # concatenate along the variant axis
    tmp = gmat.concat([gmat, gmat], axis = gmat.vrnt_axis)
    tmp = gmat.concat_vrnt([gmat, gmat])

Grouping and Sorting
====================

Genotype matrices in PyBrOpS have several sorting and grouping focused methods. Sorting methods can be used to sort taxa alphanumerically and marker variants according to their chromosome positions. Grouping methods sort taxa and marker variants, and calculate metadata to allow for indexing of taxa groups and chromosomes groups.

Reordering elements
-------------------

Taxa and marker variants in a genotype matrix can be reordered using the ``reorder`` family of methods.

.. code-block:: python

    #
    # taxa reordering example
    #

    # create reordering indices
    indices = numpy.arange(gmat.ntaxa)
    numpy.random.shuffle(indices)
    tmp = gmat.deepcopy()

    # reorder values along the taxa axis
    tmp.reorder(indices, axis = tmp.taxa_axis)
    tmp.reorder_taxa(indices)

    #
    # marker variant reordering example
    #

    # create reordering indices
    indices = numpy.arange(gmat.nvrnt)
    numpy.random.shuffle(indices)
    tmp = gmat.deepcopy()

    # reorder values along the marker variant axis
    tmp = gmat.deepcopy()
    tmp.reorder(indices, axis = tmp.vrnt_axis)
    tmp.reorder_vrnt(indices)

Lexsorting elements
-------------------

An indirect sort for the taxa and marker variants axes can be performed using the ``lexsort`` family of methods.

.. code-block:: python

    #
    # taxa lexsort example
    #

    # create lexsort keys for taxa
    key1 = numpy.random.randint(0, 10, gmat.ntaxa)
    key2 = numpy.arange(gmat.ntaxa)
    numpy.random.shuffle(key2)

    # lexsort along the taxa axis
    gmat.lexsort((key2,key1), axis = gmat.taxa_axis)
    gmat.lexsort_taxa((key2,key1))

    #
    # marker variant lexsort example
    #

    # create lexsort keys for marker variants
    key1 = numpy.random.randint(0, 10, gmat.nvrnt)
    key2 = numpy.arange(gmat.nvrnt)
    numpy.random.shuffle(key2)

    # lexsort along the marker variant axis
    gmat.lexsort((key2,key1), axis = gmat.vrnt_axis)
    gmat.lexsort_vrnt((key2,key1))

Sorting elements
----------------

Sorting along taxa and marker variant axes can be done using the ``sort`` family of methods.

.. code-block:: python

    # make copy
    tmp = gmat.deepcopy()

    #
    # taxa sorting example
    #

    # sort along taxa axis
    tmp.sort(axis = tmp.taxa_axis)
    tmp.sort_taxa()

    #
    # marker variant sorting example
    #

    # sort along marker variant axis
    tmp.sort(axis = tmp.vrnt_axis)
    tmp.sort_vrnt()

Grouping elements
-----------------

Grouping along taxa and marker variant axes can be done using the ``sort`` family of methods.

.. code-block:: python

    # make copy
    tmp = gmat.deepcopy()

    #
    # taxa grouping example
    #

    # sort along taxa axis
    tmp.group(axis = tmp.taxa_axis)
    tmp.group_taxa()

    # determine whether grouping has occurred along the taxa axis
    tmp.is_grouped(axis = tmp.taxa_axis)
    tmp.is_grouped_taxa()

    #
    # marker variant grouping example
    #

    # sort along vrnt axis
    tmp.group(axis = tmp.vrnt_axis)
    tmp.group_vrnt()

    # determine whether grouping has occurred along the vrnt axis
    tmp.is_grouped(axis = tmp.vrnt_axis)
    tmp.is_grouped_vrnt()

Summary Statistics
==================

Various summary statistics can be calculated from genotype matrices. PyBrOpS offers several common routines which are summarized in the sections below.

Population allele counts
------------------------

Allele counts of the dominant allele (the allele coded as 1) at each locus may be calculated using the ``acount`` method. The code below demonstrates this method's use.

.. code-block:: python

    # count the number of major alleles across all taxa
    out = gmat.acount()
    out = gmat.acount(dtype = "int32")

Population allele frequencies
-----------------------------

Allele frequencies of the dominant allele (the allele coded as 1) at each locus may be calculated using the ``afreq`` method. The code below demonstrates this method's use.

.. code-block:: python

    # calculate the allele frequency across all taxa
    out = gmat.afreq()
    out = gmat.afreq(dtype = "float32")

Population allele polymorphism presence
---------------------------------------

The presense of allele polymorphisms at each locus can be determined using the ``apoly`` method. The code below demonstrates this method's use.

.. code-block:: python

    # calculate whether a locus is polymorphic across all taxa 
    out = gmat.apoly()
    out = gmat.apoly(dtype = int)

Population genotype counts
--------------------------

Genotype counts at each locus may be calculated using the ``gtcount`` method. This method counts the number of homozygous recessive, heterozygous, and homozygous dominant individuals. The code below demonstrates this method's use.

.. code-block:: python

    # count the number of genotypes across all taxa
    out = gmat.gtcount()
    out = gmat.gtcount(dtype = "int32")

Population genotype frequencies
-------------------------------

Genotype frequencies at each locus may be calculated using the ``gtfreq`` method. The code below demonstrates this method's use.

.. code-block:: python

    # calculate the genotype frequency across all taxa
    out = gmat.gtfreq()
    out = gmat.gtfreq(dtype = "float32")

Population minor allele frequencies
-----------------------------------

The minor allele frequency at each locus may be calculated using the ``maf`` method. The code below demonstrates the use of this method.

.. code-block:: python

    # calculate the minor allele frequency across all taxa
    out = gmat.maf()
    out = gmat.maf(dtype = "float32")

Population mean expected heterozygosity
---------------------------------------

The mean expected heterozygosity of the genotype matrix as a whole may be calculated using the ``meh`` method. Use of the ``meh`` method can be seen in the code below.

.. code-block:: python

    # calculate the mean expected heterozygosity for the population
    out = gmat.meh()
    out = gmat.meh(dtype = "float32")

Taxa allele counts
------------------

Allele counts of the dominant allele within each individual taxon may be calculated using the ``tacount`` method. A demonstration of this method's use is below.

.. code-block:: python

    # count the number of major alleles individually within taxa
    out = gmat.tacount()
    out = gmat.tacount(dtype = "int32")

Taxa allele frequencies
-----------------------

Allele frequencies of the dominant allele within each individual taxon may be calculated using the ``tafreq`` method. The ``tafreq`` method can be used like so.

.. code-block:: python

    # calculate the allele frequency individually within taxa
    out = gmat.tafreq()
    out = gmat.tafreq(dtype = "float32")

Saving Genotype Matrices
========================

Genotype matrices can be exported and saved to disk. Currently, only one export format is available: HDF5.

Exporting to HDF5
-----------------

Genotype matrices can be exported to the `HDF5 format <https://www.hdfgroup.org/>`_. The code below demonstrates how to export a genotype matrix to an HDF5 file.

.. code-block:: python

    # write a breeding value matrix to an HDF5 file
    gmat.to_hdf5("saved_genotypes.h5")

