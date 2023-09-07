Genotyping Protocols
####################

Class Family Overview
=====================

Summary of Genotyping Protocol Classes
======================================

.. list-table:: Summary of classes in the ``pybrops.breed.prot.gt`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``GenotypingProtocol``
      - Abstract
      - Interface for all coancestry matrix child classes.
    * - ``DenseUnphasedGenotyping``
      - Concrete
      - Class representing whole-genome genotyping which produces a dense unphased genotype matrix.

Loading Class Modules
=====================

.. code-block:: python

    # import the GenotypingProtocol class (an abstract interface class)
    from pybrops.breed.prot.gt.GenotypingProtocol import GenotypingProtocol

    # import the DenseUnphasedGenotyping class (a concrete implemented class)
    from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping

Creating Genotyping Protocols
=============================

.. code-block:: python

    # object creation takes no arguments
    gtprot = DenseUnphasedGenotyping()

Genotyping Individuals from their Genomes
=========================================

.. code-block:: python

    #
    # Construct random genomes
    #

    # shape parameters for random genomes
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
    # Genotype individuals
    #

    # get unphased genotypes from phased genotypes
    gmat = gtprot.genotype(pgmat=pgmat)

    # check if output is unphased
    out = isinstance(gmat, DenseGenotypeMatrix) and not isinstance(gmat, DensePhasedGenotypeMatrix)
