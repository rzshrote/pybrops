Variance Matrices
#################

Class Family Overview
=====================

The ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` families of classes allow for representation of the expected progeny trait variances from a cross between individuals. The purpose of these families of classes is to calculate expected progeny trait variances assuming linkage for the ``GeneticVarianceMatrix`` family of classe and no linkage for the ``GenicVarianceMatrix`` family of classes. Both families of variance matrices utilize genomic models to calculate variances and are designed to be agnostic of ``GenomicModel`` type. For additive linear genomic models, there exist deterministic equations to calculate progeny variance for two-, three-, and four-way crosses. ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` families which assume an additive linear genomic model have interfaces of ``AdditiveGeneticVarianceMatrix`` and ``AdditiveGenicVarianceMatrix``, respectfully.

Deriving from the ``AdditiveGeneticVarianceMatrix`` interface are several implemented classes which are useful. They are summarized below:

Summary of Variance Matrix Classes
==================================

Genetic variance matrix classes
-------------------------------

.. list-table:: Summary of genetic variance classes in the ``pybrops.model.vmat`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``GeneticVarianceMatrix``
      - Abstract
      - Interface for all genetic variance matrix child classes.
    * - ``AdditiveGeneticVarianceMatrix``
      - Abstract
      - Interface for all additive genetic variance matrix child classes.
    * - ``DenseGeneticVarianceMatrix``
      - Semi-Abstract
      - Semi-implemented class for deriving new dense genetic variance matrix child classes.
    * - ``DenseAdditiveGeneticVarianceMatrix``
      - Semi-Abstract
      - Semi-implemented class for deriving new dense additive genetic variance matrix child classes.
    * - ``DenseTwoWayDHAdditiveGeneticVarianceMatrix``
      - Concrete
      - Class representing genetic variance matrices calculated from two-way crosses.
    * - ``DenseThreeWayDHAdditiveGeneticVarianceMatrix``
      - Concrete
      - Class representing genetic variance matrices calculated from three-way crosses.
    * - ``DenseFourWayDHAdditiveGeneticVarianceMatrix``
      - Concrete
      - Class representing genetic variance matrices calculated from four-way crosses.
    * - ``DenseDihybridDHAdditiveGeneticVarianceMatrix``
      - Concrete
      - Class representing genetic variance matrices calculated from dihybrid crosses.

Genic variance matrix classes
-----------------------------

.. list-table:: Summary of genic variance classes in the ``pybrops.model.vmat`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``GenicVarianceMatrix``
      - Abstract
      - Interface for all genic variance matrix child classes.
    * - ``AdditiveGenicVarianceMatrix``
      - Abstract
      - Interface for all additive genic variance matrix child classes.
    * - ``DenseGenicVarianceMatrix``
      - Semi-Abstract
      - Semi-implemented class for deriving new dense genic variance matrix child classes.
    * - ``DenseAdditiveGenicVarianceMatrix``
      - Semi-Abstract
      - Semi-implemented class for deriving new dense additive genic variance matrix child classes.
    * - ``DenseTwoWayDHAdditiveGenicVarianceMatrix``
      - Concrete
      - Class representing genic variance matrices calculated from two-way crosses.
    * - ``DenseThreeWayDHAdditiveGenicVarianceMatrix``
      - Concrete
      - Class representing genic variance matrices calculated from three-way crosses.
    * - ``DenseFourWayDHAdditiveGenicVarianceMatrix``
      - Concrete
      - Class representing genic variance matrices calculated from four-way crosses.
    * - ``DenseDihybridDHAdditiveGenicVarianceMatrix``
      - Concrete
      - Class representing genic variance matrices calculated from dihybrid crosses.

Loading Variance Matrix Modules
===============================

Loading genetic variance matrix modules
---------------------------------------

.. code-block:: python

    # import abstract interface classes
    from pybrops.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix
    from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix

    # import semi-abstract classes
    from pybrops.model.vmat.DenseGeneticVarianceMatrix import DenseGeneticVarianceMatrix
    from pybrops.model.vmat.DenseAdditiveGeneticVarianceMatrix import DenseAdditiveGeneticVarianceMatrix

    # import concrete implemented classes
    from pybrops.model.vmat.DenseTwoWayDHAdditiveGeneticVarianceMatrix import DenseTwoWayDHAdditiveGeneticVarianceMatrix
    from pybrops.model.vmat.DenseThreeWayDHAdditiveGeneticVarianceMatrix import DenseThreeWayDHAdditiveGeneticVarianceMatrix
    from pybrops.model.vmat.DenseFourWayDHAdditiveGeneticVarianceMatrix import DenseFourWayDHAdditiveGeneticVarianceMatrix
    from pybrops.model.vmat.DenseDihybridDHAdditiveGeneticVarianceMatrix import DenseDihybridDHAdditiveGeneticVarianceMatrix

Loading genic variance matrix modules
-------------------------------------

.. code-block:: python

    # import abstract interface classes
    from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix
    from pybrops.model.vmat.AdditiveGenicVarianceMatrix import AdditiveGenicVarianceMatrix

    # import semi-abstract classes
    from pybrops.model.vmat.DenseGenicVarianceMatrix import DenseGenicVarianceMatrix
    from pybrops.model.vmat.DenseAdditiveGenicVarianceMatrix import DenseAdditiveGenicVarianceMatrix

    # import concrete implemented classes
    from pybrops.model.vmat.DenseTwoWayDHAdditiveGenicVarianceMatrix import DenseTwoWayDHAdditiveGenicVarianceMatrix
    from pybrops.model.vmat.DenseThreeWayDHAdditiveGenicVarianceMatrix import DenseThreeWayDHAdditiveGenicVarianceMatrix
    from pybrops.model.vmat.DenseFourWayDHAdditiveGenicVarianceMatrix import DenseFourWayDHAdditiveGenicVarianceMatrix
    from pybrops.model.vmat.DenseDihybridDHAdditiveGenicVarianceMatrix import DenseDihybridDHAdditiveGenicVarianceMatrix

Creating Variance Matrices
==========================

Creating variance matrices from NumPy arrays
--------------------------------------------

.. code-block:: python

    # shape parameters for random genotypes
    ntaxa = 100
    ntrait = 2
    ngroup = 20

    # create random variance values
    mat = numpy.random.uniform(0, 1, size = (ntaxa,ntaxa,ntrait))

    # create taxa names
    taxa = numpy.array(["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], dtype = object)

    # create taxa groups
    taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
    taxa_grp.sort()

    # create trait names
    trait = numpy.array(["trait"+str(i+1).zfill(2) for i in range(ntrait)], dtype = object)

    # create genetic variance matrix
    vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix(
        mat = mat,
        taxa = taxa,
        taxa_grp = taxa_grp,
        trait = trait
    )

    # create genic variance matrix
    gvmat = DenseTwoWayDHAdditiveGenicVarianceMatrix(
        mat = mat,
        taxa = taxa,
        taxa_grp = taxa_grp,
        trait = trait
    )

Creating variance matrices from genomic models
----------------------------------------------

.. code-block:: python

    # create a dummy genomic model
    # model parameters
    nfixed = 1      # number of fixed effects
    ntrait = 2      # number of traits
    nmisc = 0       # number of miscellaneous random effects
    nadditive = 50  # number of additive marker effects

    # create dummy values
    beta = numpy.random.random((nfixed,ntrait))
    u_misc = numpy.random.random((nmisc,ntrait))
    u_a = numpy.random.random((nadditive,ntrait))
    trait = numpy.array(["Trait"+str(i+1).zfill(2) for i in range(ntrait)], dtype = object)

    # create additive linear genomic model
    algmod = DenseAdditiveLinearGenomicModel(
        beta = beta,
        u_misc = u_misc,
        u_a = u_a,
        trait = trait,
        model_name = "example",
        params = None
    )

    # shape parameters for random genotypes
    ntaxa = 100
    nvrnt = nadditive
    ngroup = 20
    nchrom = 10
    ploidy = 2

    # create random genotypes
    mat = numpy.random.randint(0, 2, size = (ploidy,ntaxa,nvrnt)).astype("int8")

    # create taxa names
    taxa = numpy.array(["Taxon"+str(i+1).zfill(3) for i in range(ntaxa)], dtype = object)

    # create taxa groups
    taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
    taxa_grp.sort()

    # create marker variant chromsome assignments
    vrnt_chrgrp = numpy.random.randint(1, nchrom+1, nvrnt)
    vrnt_chrgrp.sort()

    # create marker physical positions
    vrnt_phypos = numpy.random.choice(1000000, size = nvrnt, replace = False)
    vrnt_phypos.sort()

    # create marker genetic positions
    vrnt_genpos = numpy.random.random(nvrnt)
    vrnt_genpos.sort()

    # create marker variant names
    vrnt_name = numpy.array(["SNP"+str(i+1).zfill(4) for i in range(nvrnt)], dtype = object)

    # create a genotype matrix from scratch using NumPy arrays
    pgmat = DensePhasedGenotypeMatrix(
        mat = mat,
        taxa = taxa,
        taxa_grp = taxa_grp, 
        vrnt_chrgrp = vrnt_chrgrp,
        vrnt_phypos = vrnt_phypos, 
        vrnt_genpos = vrnt_genpos,
        vrnt_name = vrnt_name, 
        ploidy = ploidy
    )
    pgmat.group_vrnt()

    # calculate genetic variance matrix from GenomicModel
    vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_gmod(
        gmod = algmod,
        pgmat = pgmat,
        ncross = 1,
        nprogeny = 10,
        nself = 0,
        gmapfn = HaldaneMapFunction()
    )

    # calculate genetic variance matrix from AdditiveLinearGenomicModel
    vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_algmod(
        algmod = algmod,
        pgmat = pgmat,
        ncross = 1,
        nprogeny = 10,
        nself = 0,
        gmapfn = HaldaneMapFunction()
    )

    # calculate genic variance matrix from GenomicModel
    gvmat = DenseTwoWayDHAdditiveGenicVarianceMatrix.from_gmod(
        gmod = algmod,
        pgmat = pgmat,
        nprogeny = 10
    )

    # calculate genic variance matrix from AdditiveLinearGenomicModel
    gvmat = DenseTwoWayDHAdditiveGenicVarianceMatrix.from_algmod(
        algmod = algmod,
        pgmat = pgmat,
        nprogeny = 10
    )

Loading variance matrices from HDF5 files
-----------------------------------------

.. code-block:: python

    # read genetic variance matrix from HDF5 file
    vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_hdf5("saved_vmat.h5")
    gvmat = DenseTwoWayDHAdditiveGenicVarianceMatrix.from_hdf5("saved_gvmat.h5")

Variance matrix properties 
==========================

General properties 
------------------

.. list-table:: Summary of ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` general properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``mat``
      - The raw variance matrix pointer
    * - ``mat_ndim``
      - The number of dimensions for the variance matrix
    * - ``mat_shape``
      - The variance matrix shape
    * - ``epgc``
      - The expected parental genomic contribution for each parental axis

Square properties
-----------------

.. list-table:: Summary of ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` square matrix properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nsquare``
      - The number of square axes for the variance matrix
    * - ``square_axes``
      - The axes indices for the square axes for the variance matrix
    * - ``square_axes_len``
      - The lengths of the square axes for the variance matrix

Taxa properties
---------------

.. list-table:: Summary of ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` taxa properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ntaxa``
      - The number of taxa represented by the variance matrix
    * - ``taxa``
      - The names of the taxa
    * - ``taxa_axis``
      - The matrix axis along which taxa are stored
    * - ``taxa_grp``
      - An optional taxa group label
    * - ``taxa_grp_name``
      - If taxa are sorted by group: get the names of the groups
    * - ``taxa_grp_stix``
      - If taxa are sorted by group: get the start indices (inclusive) for each group
    * - ``taxa_grp_spix``
      - If taxa are sorted by group: get the stop indices (exclusive) for each group
    * - ``taxa_grp_len``
      - If taxa are sorted by group: get the length of each group

Trait properties
----------------

.. list-table:: Summary of ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` square matrix properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ntrait``
      - The number of square axes for the variance matrix
    * - ``trait``
      - The names of the traits represented by the variance matrix
    * - ``trait_axis``
      - The trait axis for the variance matrix

Copying
=======

.. code-block:: python

    # copy a genetic variance matrix
    tmp = copy.copy(vmat)
    tmp = vmat.copy()

    # copy a genic variance matrix
    tmp = copy.copy(gvmat)
    tmp = gvmat.copy()

    # deep copy a genetic variance matrix
    tmp = copy.deepcopy(vmat)
    tmp = vmat.deepcopy()

    # deep copy a genic variance matrix
    tmp = copy.deepcopy(gvmat)
    tmp = gvmat.deepcopy()
