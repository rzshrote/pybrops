Genotype Matrices
#################

Class Family Overview
=====================

Perhaps the most important family of classes in PyBrOpS is the ``GenotypeMatrix`` object family. The purpose of this family of classes is to store and represent genotypic data. ``GenotypeMatrix`` classes can be used in the estimation of genomic prediction models, the estimation of breeding values, the calculation of genomic relationship matrices, and to make selection decisions. Within PyBrOpS it is possible to read genotypic data from VCF files to create ``GenotypeMatrix`` objects, allowing for real-world data to be used in breeding program simulations.

Loading Genotype Matrix Modules
===============================

Genotype matrix support in PyBrOpS is found in the ``pybrops.popgen.gmat`` submodule. The ``GenotypeMatrix`` abstract class is the basal interface for all genotype matrices. Classes can be imported as follows:

.. code-block:: python

    # import the GenotypeMatrix class (an abstract interface class)
    from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

    # import the DenseGenotypeMatrix class (a concrete implemented class)
    from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

Creating Genotype Matrices
==========================

Creating genotype matrices from NumPy arrays
--------------------------------------------

Genotype matrices can be created from raw ``numpy.ndarray``'s. The example below illustrates creating a ``DenseGenotypeMatrix``. The only requirement to construct a ``DenseGenotypeMatrix`` is a ``(n,p)`` numerical matrix containing genotypic codings. Here, ``n`` is the number of taxa and ``p`` is the number of marker variants. The matrix must be coded in the ``{0,1,2}`` (number of alternative alleles) format and must have an ``int8`` data type.

Additional optional metadata may be stored along with a ``DenseGenotypeMatrix`` including taxa names (``taxa``), taxa groups (``taxa_grp``), marker variant chromosome assignments (``vrnt_chrgrp``), marker variant chromosme physical positions (``vrnt_phypos``), marker variant names (``vrnt_name``), marker variant genetic map positions (``vrnt_genpos``), sequential recombination probabilities between markers (``vrnt_xoprob``), marker variant haplotype group assignments (``vrnt_hapgrp``), reference haplotype (``vrnt_hapref``), alternative haplotype (``vrnt_hapalt``), a variant mask (``vrnt_mask``), and finally the ploidy of the genotypes (``ploid``).

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

VCF files can be loaded into the Python environment as genotype matrices using the ``from_vcf`` method.

.. code-block:: python

    # read a genotype matrix from file
    gmat = DenseGenotypeMatrix.from_vcf("widiv_2000SNPs.vcf.gz")

Loading genotype matrices from HDF5 files
-----------------------------------------

Genotype matrices in PyBrOpS can be exported to HDF5 files via the ``to_hdf5`` method. These files can later be read into PyBrOpS via the ``from_hdf5`` method. The example below illustrates loading a ``GenotypeMatrix`` into memory from an HDF5 file:

.. code-block:: python

    # read a genotype matrix from HDF5 file
    gmat = DenseGenotypeMatrix.from_hdf5("widiv_2000SNPs.h5")

Genotype Matrix General Properties
==================================

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
    * - ``ntaxa``
      - The number of taxa represented by the genotype matrix
    * - ``nvrnt``
      - The number of genotype variants represented by the genotype matrix

Genotype Matrix Taxa Properties
===============================

.. list-table:: Summary of ``GenotypeMatrix`` taxa properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
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

Genotype Matrix Marker Variant Properties
=========================================

.. list-table:: Summary of ``GenotypeMatrix`` marker variant properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
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

Shallow copying
---------------

.. code-block:: python

    # copy a genotype matrix
    tmp = copy.copy(gmat)
    tmp = gmat.copy()

Deep copying
------------

.. code-block:: python

    # deep copy a genotype matrix
    tmp = copy.deepcopy(gmat)
    tmp = gmat.deepcopy()

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
    out = gmat.acount(dtype = "int32")

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
