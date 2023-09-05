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

.. code-block:: python

    # model parameters
    nfixed = 1      # number of fixed effects
    ntrait = 2      # number of traits
    nmisc = 0       # number of miscellaneous random effects
    nadditive = 50  # number of additive marker effects

    # create dummy values
    beta = numpy.random.random((nfixed,ntrait))
    u_misc = numpy.random.random((nmisc,ntrait))
    u_a = numpy.random.random((nadditive,ntrait))
    trait = numpy.array(
        ["Trait"+str(i+1).zfill(2) for i in range(ntrait)],
        dtype = object
    )

    # create additive linear genomic model
    algmod = DenseAdditiveLinearGenomicModel(
        beta = beta,
        u_misc = u_misc,
        u_a = u_a,
        trait = trait,
        model_name = "example",
        params = None
    )

Genomic Model Properties
========================

General properties
------------------

.. list-table:: Summary of ``GenotypeMatrix`` general properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``model_name``
      - Name assigned to the genomic model
    * - ``params``
      - Model hyperparameters (e.g. for fitting)

Model coefficient properties
----------------------------

.. list-table:: Summary of ``GenotypeMatrix`` taxa properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``beta``
      - Coefficients for fixed effects
    * - ``u``
      - Coefficients for all random effects
    * - ``u_a``
      - Coefficients for additive marker random effects
    * - ``u_misc``
      - Coefficients for miscellaneous random effects

Trait properties
----------------

.. list-table:: Summary of ``GenotypeMatrix`` marker variant properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ntrait``
      - Number of traits for which the model predicts
    * - ``trait``
      - Names of traits for which the model predicts

Copying Genomic Models
======================

Shallow copying
---------------

.. code-block:: python

    # copy a genomic model
    tmp = copy.copy(algmod)
    tmp = algmod.copy()

Deep copying
------------

.. code-block:: python

    # deep copy a genomic model
    tmp = copy.deepcopy(algmod)
    tmp = algmod.deepcopy()


Model prediction methods
========================

.. code-block:: python

    # create random genotypes to test
    X = numpy.ones((ntaxa,1))
    Z = numpy.random.randint(0, ploidy+1, size = (ntaxa,nvrnt)).astype("int8")

    # predict genotypic values using numpy arrays
    out = algmod.predict_numpy(X, Z)

    # predict genotypic values using objects
    out = algmod.predict(cvobj = X, gtobj = gmat)

Score model prediction accuracy
===============================

.. code-block:: python

    # create some dummy matrices for input
    X = numpy.ones((ntaxa,1))
    B = algmod.beta
    Z = gmat.mat
    U = algmod.u_a
    Y = X@B + Z@U
    e = numpy.random.normal(size = Y.shape)
    Y += e

    # score predictions using numpy arrays
    out = algmod.score_numpy(Y, X, Z)

    # score predictions using objects
    out = algmod.score(ptobj = Y, cvobj = X, gtobj = gmat)

Predicting genomic estimated breeding values
============================================

.. code-block:: python

    # predict GEBVs using numpy arrays
    out = algmod.gebv_numpy(Z)

    # predict GEBVs using objects
    out = algmod.gebv(gmat)

Calculating population genetic variance terms
=============================================

Predicting genetic variance
---------------------------

.. code-block:: python

    # predict Var(G) using numpy arrays
    out = algmod.var_G_numpy(Z)

    # predict Var(G) using objects
    out = algmod.var_G(gmat)

Predicting additive genetic variance
------------------------------------

.. code-block:: python

    # predict Var(A) using numpy arrays
    out = algmod.var_A_numpy(Z)

    # predict Var(A) using objects
    out = algmod.var_A(gmat)

Predicting additive genic variance
----------------------------------

.. code-block:: python

    # predict Var(a) using numpy arrays
    out = algmod.var_a_numpy(
        p = gmat.afreq(),
        ploidy = gmat.ploidy
    )

    # predict Var(a) using objects
    out = algmod.var_a(gmat)

Predicting the Bulmer effect
----------------------------

.. code-block:: python

    # predict Bulmer effect using numpy arrays
    out = algmod.bulmer_numpy(
        Z,
        p = gmat.afreq(),
        ploidy = gmat.ploidy
    )

    # predict Bulmer effect using objects
    out = algmod.bulmer(gmat)

Calculating population selection limits
=======================================

Upper selection limit
---------------------

.. code-block:: python

    # upper selection limit using numpy arrays
    out = algmod.usl_numpy(
        p = gmat.afreq(),
        ploidy = gmat.ploidy
    )

    # upper selection limit using objects
    out = algmod.usl(gtobj = gmat)

Lower selection limit
---------------------

.. code-block:: python

    # lower selection limit using numpy arrays
    out = algmod.lsl_numpy(
        p = gmat.afreq(),
        ploidy = gmat.ploidy
    )

    # lower selection limit using objects
    out = algmod.lsl(gtobj = gmat)

Calculating favorable allele metrics
====================================

.. code-block:: python

    # calculate favorable allele counts
    out = algmod.facount(gmat)

    # calculate favorable allele frequencies
    out = algmod.fafreq(gmat)

    # calculate favorable allele availability at loci in a population
    out = algmod.faavail(gmat)

    # calculate favorable allele fixation at loci in a population
    out = algmod.fafixed(gmat)

Calculating deleterious allele metrics
======================================

.. code-block:: python

    # calculate deleterious allele counts
    out = algmod.dacount(gmat)

    # calculate deleterious allele frequencies
    out = algmod.dafreq(gmat)

    # calculate deleterious allele availability at loci in a population
    out = algmod.daavail(gmat)

    # calculate deleterious allele fixation at loci in a population
    out = algmod.dafixed(gmat)

Reading and writing a genomic model
===================================

.. code-block:: python

    # writing a genomic model to a file
    algmod.to_hdf5("saved_algmod.h5")

    # reading a genomic model from a file
    out = DenseAdditiveLinearGenomicModel.from_hdf5("saved_algmod.h5")
