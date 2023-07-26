Genomic Modeling Data Types
###########################

The ``pybrops.gmod`` module houses submodules, classes, and functions related to genomic modeling. In this context, genomic modeling refers to anything which takes genomic information as an input to create an output.

Genomic Models
**************

The ``GenomicModel`` family of classes allow for the representation of any type of genomic model. The purpose of this family of classes is primarily to provide functionality for fitting a genomic prediction model and estimating breeding values. ``GenomicModel`` classes also provide functionality for estimating population genetic, and genic variances, allele value metrics, and upper and lower selection limits.

There are many interfaces deriving from the main ``GenomicModel`` interface. Briefly, the ``GenomicModel`` family can be divided into ``LinearGenomicModel`` classes and ``NonlinearGenomicModel`` classes. As their names suggest, the former defines genomic models which are linear in nature, while the latter defines models which are non-linear in nature. The ``NonlinearGenomicModel`` interface may be useful for defining machine learning models.

Deriving from the ``LinearGenomicModel`` interface are several subtypes of linear genomic models. Derivative interfaces are summarized below:

.. list-table::
    :widths: 25 50
    :header-rows: 1

    * - Interface
      - Description
    * - ``AdditiveLinearGenomicModel``
      - Defines genomic models which assume strictly additive allelic effects.
    * - ``CoancestryLinearGenomicModel``
      - Defines genomic models which calculate breeding values from coancestry relationships.
    * - ``AdditiveDominanceLinearGenomicModel``
      - Defines genomic models which assume additive and dominance allelic effects.
    * - ``AdditiveDominanceEpistaticLinearGenomicModel``
      - Defines genomic models which assume additive, dominance, and epistatic allelic effects.

Variance Matrices
*****************
The ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` families of classes allow for representation of the expected progeny trait variances from a cross between individuals. The purpose of these families of classes is to calculate expected progeny trait variances assuming linkage for the ``GeneticVarianceMatrix`` family of classe and no linkage for the ``GenicVarianceMatrix`` family of classes. Both families of variance matrices utilize genomic models to calculate variances and are designed to be agnostic of ``GenomicModel`` type. For additive linear genomic models, there exist deterministic equations to calculate progeny variance for two-, three-, and four-way crosses. ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` families which assume an additive linear genomic model have interfaces of ``AdditiveGeneticVarianceMatrix`` and ``AdditiveGenicVarianceMatrix``, respectfully.

Deriving from the ``AdditiveGeneticVarianceMatrix`` interface are several implemented classes which are useful. They are summarized below:

.. list-table::
    :widths: 25 50
    :header-rows: 1

    * - Implemented Class
      - Description
    * - ``DenseTwoWayDHAdditiveGeneticVarianceMatrix``
      - Implements genetic variance matrices calculated from two-way crosses.
    * - ``DenseThreeWayDHAdditiveGeneticVarianceMatrix``
      - Implements genetic variance matrices calculated from three-way crosses.
    * - ``DenseFourWayDHAdditiveGeneticVarianceMatrix``
      - Implements genetic variance matrices calculated from four-way crosses.
    * - ``DenseDihybridDHAdditiveGeneticVarianceMatrix``
      - Implements genetic variance matrices calculated from dihybrid crosses.
