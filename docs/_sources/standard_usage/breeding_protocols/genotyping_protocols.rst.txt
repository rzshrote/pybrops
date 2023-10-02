Genotyping Protocols
####################

Class Family Overview
=====================

The ``GenotypingProtocol`` family of classes is used to simulate the genotyping of individuals. ``GenotypeProtocol`` classes take a set of phased genomes (represented by a ``PhasedGenotypeMatrix``) and return a (typically unphased) genotype matrix of the individuals in the input.

Summary of Genotyping Protocol Classes
======================================

Genotyping protocols can be found in the ``pybrops.breed.prot.gt`` module. PyBrOpS currently only has a single implemented class in this family, the ``DenseUnphasedGenotyping`` class, which represents flawless, whole-genome unphased genotyping. The genotyping protocols module is deliberately left open for users to implement their own classes if so desired.

.. list-table:: Summary of classes in the ``pybrops.breed.prot.gt`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``GenotypingProtocol``
      - Abstract
      - Interface for all genotyping protocol classes.
    * - ``DenseUnphasedGenotyping``
      - Concrete
      - Class representing whole-genome genotyping which produces a dense unphased genotype matrix.

Genotyping Protocol Properties
==============================

Genotyping protocols do not have any required properties defined in their interface.

Loading Class Modules
=====================

Genotyping protocols may be imported as illustrated below.

.. code-block:: python

    # import the GenotypingProtocol class (an abstract interface class)
    from pybrops.breed.prot.gt.GenotypingProtocol import GenotypingProtocol

    # import the DenseUnphasedGenotyping class (a concrete implemented class)
    from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping

Creating Genotyping Protocols
=============================

Genotyping protocol class construction is entirely implementation dependent. Below is an example of how to construct a ``DenseUnphasedGenotyping`` object, which takes no arguments for its construction.

.. code-block:: python

    # object creation takes no arguments
    gtprot = DenseUnphasedGenotyping()

Genotyping Individuals from their Genomes
=========================================

Individuals may be genotyped using the ``genotype`` method. This method takes a set of individual genomes (represented using a ``PhasedGenotypeMatrix``) and returns a ``GenotypeMatrix`` representing genotyped individuals. A demonstration of how to use the ``genotype`` method is shown below.

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
