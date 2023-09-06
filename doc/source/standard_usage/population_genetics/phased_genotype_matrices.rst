Phased Genotype Matrices
########################

Class Family Overview
=====================

Loading Class Modules
=====================

.. code-block:: python

    # import the PhasedGenotypeMatrix class (an abstract interface class)
    from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

    # import the DensePhasedGenotypeMatrix class (a concrete implemented class)
    from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

Creating Phased Genotype Matrices
=================================

Creating phased genotype matrices from NumPy arrays
---------------------------------------------------

.. code-block:: python

    # shape parameters for random genotypes
    ntaxa = 100
    nvrnt = 1000
    ngroup = 20
    nchrom = 10
    ploidy = 2

    # create random genotypes
    mat = numpy.random.randint(0, ploidy+1, size = (ntaxa,nvrnt)).astype("int8")

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
        vrnt_genpos = None,
        vrnt_xoprob = None, 
        vrnt_hapgrp = None, 
        vrnt_hapalt = None,
        vrnt_hapref = None, 
        vrnt_mask = None,
        ploidy = ploidy
    )

Loading phased genotype matrices from VCF files
-----------------------------------------------

.. code-block:: python

    # read a phased genotype matrix from VCF file
    pgmat = DensePhasedGenotypeMatrix.from_vcf("widiv_2000SNPs.vcf.gz")

Loading phased genotype matrices from HDF5 files
------------------------------------------------

.. code-block:: python

    # read a genotype matrix from HDF5 file
    pgmat = DensePhasedGenotypeMatrix.from_hdf5("widiv_2000SNPs.h5")

Phased Genotype Matrix Properties
=================================

General properties
------------------

.. list-table:: Summary of ``PhasedGenotypeMatrix`` general properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``mat``
      - Pointer to the raw phased genotype matrix pointer
    * - ``mat_ndim``
      - The number of dimensions for the phased genotype matrix
    * - ``mat_shape``
      - Genotype matrix shape
    * - ``mat_format``
      - Genotype matrix format
    * - ``ploidy``
      - The ploidy of the taxa represented by the phased genotype matrix

Phase properties
----------------

.. list-table:: Summary of ``PhasedGenotypeMatrix`` phase properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nphase``
      - The number of chromosome phases represented by the phased genotype matrix
    * - ``phase_axis``
      - The matrix axis along which phases are stored

Taxa properties
---------------

.. list-table:: Summary of ``PhasedGenotypeMatrix`` taxa properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ntaxa``
      - The number of taxa represented by the phased genotype matrix
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

.. list-table:: Summary of ``PhasedGenotypeMatrix`` marker variant properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nvrnt``
      - The number of genotype variants represented by the phased genotype matrix
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

Copying Phased Genotype Matrices
================================

.. code-block:: python

    # copy a phased genotype matrix
    tmp = copy.copy(pgmat)
    tmp = pgmat.copy()

    # deep copy a phased genotype matrix
    tmp = copy.deepcopy(pgmat)
    tmp = pgmat.deepcopy()

Phased Genotype Matrix Element Copy-On-Manipulation
===================================================

Adjoining elements
------------------

.. code-block:: python

    # create a new genotype matrix to demonstrate
    new = pgmat.deepcopy()

    # adjoin genotype matrices along the taxa axis
    tmp = pgmat.adjoin(new, axis = pgmat.taxa_axis)
    tmp = pgmat.adjoin_taxa(new)

    # adjoin genotype matrices along the variant axis
    tmp = pgmat.adjoin(new, axis = pgmat.vrnt_axis)
    tmp = pgmat.adjoin_vrnt(new)

Deleting elements
-----------------

``delete`` taxa
+++++++++++++++

.. code-block:: python

    # delete first taxon using an integer
    tmp = pgmat.delete(0, axis = pgmat.taxa_axis)
    tmp = pgmat.delete_taxa(0)

    # delete first five taxa using a slice
    tmp = pgmat.delete(slice(0,5), axis = pgmat.taxa_axis)
    tmp = pgmat.delete_taxa(slice(0,5))

    # delete first five taxa using a Sequence
    tmp = pgmat.delete([0,1,2,3,4], axis = pgmat.taxa_axis)
    tmp = pgmat.delete_taxa([0,1,2,3,4])

``delete`` marker variants
++++++++++++++++++++++++++

.. code-block:: python

    # delete first marker variant using an integer
    tmp = pgmat.delete(0, axis = pgmat.vrnt_axis)
    tmp = pgmat.delete_vrnt(0)

    # delete first five marker variants using a slice
    tmp = pgmat.delete(slice(0,5), axis = pgmat.vrnt_axis)
    tmp = pgmat.delete_vrnt(slice(0,5))

    # delete first five marker variants using a Sequence
    tmp = pgmat.delete([0,1,2,3,4], axis = pgmat.vrnt_axis)
    tmp = pgmat.delete_vrnt([0,1,2,3,4])

Inserting elements
------------------

.. code-block:: python

    # create a new genotype matrix to demonstrate
    new = pgmat.deepcopy()

    # insert genotype matrix along the taxa axis before index 0
    tmp = pgmat.insert(0, new, axis = pgmat.taxa_axis)
    tmp = pgmat.insert_taxa(0, new)

    # insert genotype matrix along the variant axis before index 0
    tmp = pgmat.insert(0, new, axis = pgmat.vrnt_axis)
    tmp = pgmat.insert_vrnt(0, new)

Selecting elements
------------------

.. code-block:: python

    # select first five taxa using a Sequence
    tmp = pgmat.select([0,1,2,3,4], axis = pgmat.taxa_axis)
    tmp = pgmat.select_taxa([0,1,2,3,4])

    # select first five marker variants using a Sequence
    tmp = pgmat.select([0,1,2,3,4], axis = pgmat.vrnt_axis)
    tmp = pgmat.select_vrnt([0,1,2,3,4])

Phased Genotype Matrix Element In-Place-Manipulation
====================================================

Appending elements
------------------

.. code-block:: python

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

Removing elements
-----------------

``remove`` taxa
+++++++++++++++

.. code-block:: python

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

``remove`` marker variants
++++++++++++++++++++++++++

.. code-block:: python

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

Incorporating elements
----------------------

.. code-block:: python

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

Concatenating matrices
----------------------

.. code-block:: python

    # concatenate along the taxa axis
    tmp = pgmat.concat([pgmat, pgmat], axis = pgmat.taxa_axis)
    tmp = pgmat.concat_taxa([pgmat, pgmat])

    # concatenate along the variant axis
    tmp = pgmat.concat([pgmat, pgmat], axis = pgmat.vrnt_axis)
    tmp = pgmat.concat_vrnt([pgmat, pgmat])

Grouping and Sorting
====================

Reordering
----------

``reorder`` taxa
++++++++++++++++

.. code-block:: python

    # create reordering indices
    indices = numpy.arange(pgmat.ntaxa)
    numpy.random.shuffle(indices)
    tmp = pgmat.deepcopy()

    # reorder values along the taxa axis
    tmp.reorder(indices, axis = tmp.taxa_axis)
    tmp.reorder_taxa(indices)

``reorder`` marker variants
+++++++++++++++++++++++++++

.. code-block:: python

    # create reordering indices
    indices = numpy.arange(pgmat.nvrnt)
    numpy.random.shuffle(indices)
    tmp = pgmat.deepcopy()
    
    # reorder values along the marker variant axis
    tmp = pgmat.deepcopy()
    tmp.reorder(indices, axis = tmp.vrnt_axis)
    tmp.reorder_vrnt(indices)

Lexsorting
----------

``lexsort`` taxa
++++++++++++++++

.. code-block:: python

    # create lexsort keys for taxa
    key1 = numpy.random.randint(0, 10, pgmat.ntaxa)
    key2 = numpy.arange(pgmat.ntaxa)
    numpy.random.shuffle(key2)

    # lexsort along the taxa axis
    pgmat.lexsort((key2,key1), axis = pgmat.taxa_axis)
    pgmat.lexsort_taxa((key2,key1))

``lexsort`` marker variants
+++++++++++++++++++++++++++

.. code-block:: python

    # create lexsort keys for marker variants
    key1 = numpy.random.randint(0, 10, pgmat.nvrnt)
    key2 = numpy.arange(pgmat.nvrnt)
    numpy.random.shuffle(key2)

    # lexsort along the marker variant axis
    pgmat.lexsort((key2,key1), axis = pgmat.vrnt_axis)
    pgmat.lexsort_vrnt((key2,key1))

Sorting
-------

``sort`` taxa
+++++++++++++

.. code-block:: python

    # sort along taxa axis
    tmp = pgmat.deepcopy()
    tmp.sort(axis = tmp.taxa_axis)
    tmp.sort_taxa()

``sort`` marker variants
++++++++++++++++++++++++

.. code-block:: python

    # sort along marker variant axis
    tmp = pgmat.deepcopy()
    tmp.sort(axis = tmp.vrnt_axis)
    tmp.sort_vrnt()

Grouping
--------

``group`` taxa
++++++++++++++

.. code-block:: python

    # sort along taxa axis
    tmp = pgmat.deepcopy()
    tmp.group(axis = tmp.taxa_axis)
    tmp.group_taxa()
    # determine whether grouping has occurred along the taxa axis
    tmp.is_grouped(axis = tmp.taxa_axis)
    tmp.is_grouped_taxa()

``group`` marker variants
+++++++++++++++++++++++++

.. code-block:: python

    # sort along vrnt axis
    tmp = pgmat.deepcopy()
    tmp.group(axis = tmp.vrnt_axis)
    tmp.group_vrnt()
    # determine whether grouping has occurred along the vrnt axis
    tmp.is_grouped(axis = tmp.vrnt_axis)
    tmp.is_grouped_vrnt()

Summary Statistics
==================

.. code-block:: python

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

Saving Genotype Matrices
========================

Write to HDF5
-------------

.. code-block:: python

    # write a breeding value matrix to an HDF5 file
    pgmat.to_hdf5("saved_genotypes.h5")
