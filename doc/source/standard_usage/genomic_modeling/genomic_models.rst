Genomic Models
##############

Class Family Overview
=====================

The ``GenomicModel`` family of classes allow for the representation of any type of genomic model. The purpose of this family of classes is primarily to provide functionality for fitting a genomic prediction model and estimating breeding values. ``GenomicModel`` classes also provide functionality for estimating population genetic, and genic variances, allele value metrics, and upper and lower selection limits.

Summary of Genomic Model Classes
================================

There are many interfaces deriving from the main ``GenomicModel`` interface. Briefly, the ``GenomicModel`` family can be divided into ``LinearGenomicModel`` classes and ``NonlinearGenomicModel`` classes. As their names suggest, the former defines genomic models which are linear in nature, while the latter defines models which are non-linear in nature. The ``NonlinearGenomicModel`` interface may be useful for defining machine learning models.

Deriving from the ``LinearGenomicModel`` interface are several subtypes of linear genomic models. Derivative interfaces are summarized below:

.. list-table:: Summary of classes in ``pybrops.model.gmod`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``GenomicModel``
      - Abstract
      - Interface for all genomic models.
    * - ``NonlinearGenomicModel``
      - Abstract
      - Interface for non-linear genomic models.
    * - ``LinearGenomicModel``
      - Abstract
      - Interface for linear genomic models.
    * - ``CoancestryLinearGenomicModel``
      - Abstract
      - Interface for genomic models which calculate breeding values from coancestry relationships.
    * - ``AdditiveLinearGenomicModel``
      - Abstract
      - Interface for genomic models which assume strictly additive allelic effects.
    * - ``AdditiveDominanceLinearGenomicModel``
      - Abstract
      - Interface for genomic models which assume additive and dominance allelic effects.
    * - ``AdditiveDominanceEpistaticLinearGenomicModel``
      - Abstract
      - Interface for genomic models which assume additive, dominance, and epistatic allelic effects.
    * - ``DenseLinearGenomicModel``
      - Concrete
      - Class representing a generic linear genomic model.
    * - ``DenseAdditiveLinearGenomicModel``
      - Concrete
      - Class representing a generic additive linear genomic model.

Loading Genomic Model Modules
=============================

.. code-block:: python

    # import GenomicModel classes (abstract interface classes)
    from pybrops.model.gmod.GenomicModel import GenomicModel
    from pybrops.model.gmod.NonlinearGenomicModel import NonlinearGenomicModel
    from pybrops.model.gmod.LinearGenomicModel import LinearGenomicModel
    from pybrops.model.gmod.CoancestryLinearGenomicModel import CoancestryLinearGenomicModel
    from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
    from pybrops.model.gmod.AdditiveDominanceLinearGenomicModel import AdditiveDominanceLinearGenomicModel
    from pybrops.model.gmod.AdditiveDominanceEpistaticLinearGenomicModel import AdditiveDominanceEpistaticLinearGenomicModel

    # import dense genomic models (concrete implementation classes)
    from pybrops.model.gmod.DenseLinearGenomicModel import DenseLinearGenomicModel
    from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel

Creating Genomic Models
=======================

Creating genomic models from NumPy arrays
-----------------------------------------

