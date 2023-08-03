Genomic Models
##############

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
