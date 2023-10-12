Genomic Models
##############

Class Family Overview
=====================

The ``GenomicModel`` family of classes allow for the representation of any type of genomic model. The purpose of this family of classes to provide functionality for fitting genomic prediction models and estimating breeding values. ``GenomicModel`` classes also provide functionality for estimating population genetic, and genic variances, allele value metrics, and upper and lower selection limits.

Summary of Genomic Model Classes
================================

The central ``GenomicModel`` interface defines all basic genomic models operations. The ``GenomicModel`` interface is intended to be extremely broad as there are many modeling methods by which genetic values may be estimated and modeled. Deriving from this interface are two broad categories: linear genomic models and non-linear genomic models. Linear and non-linear genomic model behaviors are defined by the ``LinearGenomicModel`` and ``NonlinearGenomicModel`` interfaces, respectively. As their names suggest, the former defines genomic models which are linear in nature, while the latter defines models which are non-linear in nature. The ``NonlinearGenomicModel`` interface may be useful for defining machine learning models.

Deriving from the ``LinearGenomicModel`` interface are several subtypes of linear genomic models. Derivatives include additive genomic models, additive-dominance genomic models, and additive-dominance-epistasis genomic models. Additive genomic models are models where genomic marker effects are assumed to be strictly additive in nature. The behavior or additive genomic models is defined by the ``AddititiveLinearGenomicModel`` interface. Additive-dominance genomic models are models where genomic markers exhibit additive and dominance effects. The behavior of these types of models is defined by the ``AdditiveDominanceLinearGenomicModel``. Finally, additive-dominance-epistasis genomic models are models where genomic markers exhibit additive, dominance, and epistatic effects. The behavior of these types of models is defined by the ``AdditiveDominanceEpistaticLinearGenomicModel``.

Below is a summary of the genomic model abstract interfaces present in PyBrOpS.

.. list-table:: Summary of abstract classes in the ``pybrops.model.gmod`` module
    :widths: 25 15 50
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
    * - | ``AdditiveDominance``
        | ``LinearGenomicModel``
      - Abstract
      - Interface for genomic models which assume additive and dominance allelic effects.
    * - | ``AdditiveDominanceEpistatic``
        | ``LinearGenomicModel``
      - Abstract
      - Interface for genomic models which assume additive, dominance, and epistatic allelic effects.

PyBrOpS implements several generic linear genomic model classes since these model types are common. Below is a summary of implemented genomic model classes in PyBrOpS.

.. list-table:: Summary of concrete classes in the ``pybrops.model.gmod`` module
    :widths: 25 15 50
    :header-rows: 1

    * - ``DenseLinearGenomicModel``
      - Concrete
      - Class representing a generic linear genomic model.
    * - ``DenseAdditiveLinearGenomicModel``
      - Concrete
      - Class representing a generic additive linear genomic model.

Genomic Model Properties
========================

Genomic models share several properties in common. These properties can be categorized into three groups: general properties, model properties, and trait properties. These property groupings are summarized in the following three subsections. 

In addition to base ``GenomicModel`` properties, derivatives of the ``GenomicModel`` interface may define extra properties as in the case of the ``LinearGenomicModel`` family of classes. The subsections following the ``GenomicModel`` subsections detail properties which have been added to classes in the ``LinearGenomicModel`` class family.

Genomic Models: General properties
----------------------------------

All genomic models share a set of general, basic properties. These properties are the name of the model and any hyperparameters that the model may have. General genomic model properties are summarized below.

.. list-table:: Summary of ``GenomicModel`` general properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``model_name``
      - Name assigned to the genomic model
    * - ``hyperparams``
      - Model hyperparameters (e.g. for fitting)

Genomic Models: Model properties
--------------------------------

All genomic models also share a set of model properties pertaining to the dimensionality of the model. These include the number of explanatory variables that the model uses as predictors and the number of parameters of the model. These model properties are intensionally very basic because of the broad diversity of modeling techniques which may be used to model genotypic values.

.. list-table:: Summary of ``GenomicModel`` model properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nexplan``
      - Number of explanatory variables required by the model
    * - ``nparam``
      - Number of model parameters

Genomic Models: Trait properties
--------------------------------

Genomic models also contain the names of traits for which the model predicts. Trait related properties are summarized in the table below.

.. list-table:: Summary of ``GenomicModel`` trait properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ntrait``
      - Number of traits for which the model predicts
    * - ``trait``
      - Names of traits for which the model predicts

Linear Genomic Models: Additional model properties
--------------------------------------------------

For ``LinearGenomicModel`` s, additional fixed and random effect model properties are added on top of the properties in the base ``GenomicModel`` class. These additional properties pertain to the number of explanatory variables, number of model parameters, and regression coefficients for the linear mixed model. Below summarizes the additional model properties.

.. list-table:: Summary of ``LinearGenomicModel`` additional model properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nexplan_beta``
      - Number of fixed effect explanatory variables required by the model
    * - ``nparam_beta``
      - Number of fixed effect parameters
    * - ``beta``
      - Coefficients for fixed effects
    * - ``nexplan_u``
      - Number of random effect explanatory variables required by the model
    * - ``nparam_u``
      - Number of random effect parameters
    * - ``u``
      - Coefficients for random effects

Additive Linear Genomic Models: Additional model properties
-----------------------------------------------------------

For ``AdditiveLinearGenomicModel`` s, the random effects properties of the ``LinearGenomicModel`` class are subdivided into miscellaneous random effects and additive genomic marker effects. The ``AdditiveLinearGenomicModel`` interface adds additional properties to access these random effect subdivisions. Below summarizes model properties which are added.

.. list-table:: Summary of ``AdditiveLinearGenomicModel`` additional model properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nexplan_u_misc``
      - Number of miscellaneous random effect explanatory variables required by the model
    * - ``nparam_u_misc``
      - Number of miscellaneous random effect parameters
    * - ``u_misc``
      - Coefficients for miscellaneous random effects
    * - ``nexplan_u_a``
      - Number of additive genomic marker explanatory variables required by the model
    * - ``nparam_u_a``
      - Number of additive genomic marker parameters
    * - ``u_a``
      - Coefficients for additive genomic marker effects

Additive Dominance Linear Genomic Models: Additional model properties
---------------------------------------------------------------------

For the ``AdditiveDominanceLinearGenomicModel`` class family, an additional set of dominance random effect properties is added on top of the base ``AdditiveLinearGenomicModel`` class. Below summarizes the dominance model properties which are added.

.. list-table:: Summary of ``AdditiveDominanceLinearGenomicModel`` additional model properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nexplan_u_d``
      - Number of dominance genomic marker explanatory variables required by the model
    * - ``nparam_u_d``
      - Number of dominance genomic marker parameters
    * - ``u_d``
      - Coefficients for dominance genomic marker effects

Additive Dominance Epistatic Linear Genomic Models: Model coefficient properties
--------------------------------------------------------------------------------

For the ``AdditiveDominanceEpistaticLinearGenomicModel`` class family, an additional set of epistatic random effect properties is added on top of the base ``AdditiveDominanceLinearGenomicModel`` class. Below summarizes the epistatic model properties which are added.

.. list-table:: Summary of ``AdditiveDominanceEpistaticLinearGenomicModel`` model coefficient properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nexplan_u_i``
      - Number of epistatic genomic marker explanatory variables required by the model
    * - ``nparam_u_i``
      - Number of epistatic genomic marker parameters
    * - ``u_i``
      - Coefficients for epistatic genomic marker effects

Loading Genomic Model Modules
=============================

Genomic model classes can be loaded from the ``pybrops.model.gmod`` module using the import statements in the code below.

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

Genomic models can be created using multiple methods. All genomic models can be created from raw NumPy arrays using the class constructor or loaded from HDF5 files. In addition, ``LinearGenomicModel`` classes can be loaded from Pandas DataFrames or from CSV files. The sections below demonstrate how to load a ``DenseAdditiveLinearGenomicModel`` using all four of the aforementioned creation mechanisms.

Creating genomic models from raw NumPy arrays
---------------------------------------------

The code example below demonstrates how to create a dense additive linear genomic model from raw NumPy arrays using the class constructor. Most genomic model class constructors accept NumPy arrays for their arguments since this is the lowest common denominator. This said, the genomic model interface defined by PyBrOpS does not require that a class adhere to this norm. A genomic model constructor is implementation dependent, allowing for freedom in cases where construction from NumPy arrays is impractical or undesirable.

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

Loading linear genomic models from dictionaries of Pandas DataFrames
--------------------------------------------------------------------

``LinearGenomicModel`` s can be loaded from dictionaries of Pandas DataFrames. Below demonstrates how the ``from_pandas_dict`` can be used to load a model from multiple Pandas DataFrames.

.. code-block:: python

    nfixed = 1      # number of fixed effects
    ntrait = 2      # number of traits
    nmisc = 0       # number of miscellaneous random effects
    nadditive = 50  # number of additive marker effects

    # create dummy values
    beta = numpy.random.random((nfixed,ntrait))
    u_misc = numpy.random.random((nmisc,ntrait))
    u_a = numpy.random.random((nadditive,ntrait))
    trait = numpy.array(["Trait"+str(i+1) for i in range(ntrait)], dtype=object)

    # create dictionary with dataframes
    # need required fields of "beta", "u_misc", "u_a" for this class
    df_dict = {
        "beta": pandas.DataFrame(beta, columns = trait),
        "u_misc": pandas.DataFrame(u_misc, columns = trait),
        "u_a": pandas.DataFrame(u_a, columns = trait)
    }

    # construct model from dictionary of pandas dataframes
    mod = DenseAdditiveLinearGenomicModel.from_pandas_dict(df_dict)

Loading linear genomic models from dictionaries of CSV file names
-----------------------------------------------------------------

CSV files may be used to load linear genomic models in a similar manner to loading from Pandas DataFrames. The code below demonstrates how to load a linear genomic model from multiple CSV files containing different sets of marker coefficients.

.. code-block:: python

    # create dictionary with filenames
    # need required fields of "beta", "u_misc", "u_a" for this class
    dic = {
        "beta": "beta.csv",
        "u_misc": "u_misc.csv",
        "u_a": "u_a.csv",
    }

    # construct model from dictionary of pandas dataframes
    mod = DenseAdditiveLinearGenomicModel.from_csv_dict(dic)

Loading genomic models HDF5 files
---------------------------------

All genomic models may be saved and reloaded from HDF5 files. The code below demonstrates how to read a ``GenomicModel`` from an HDF5 file.

.. code-block:: python

    # load from HDF5 file
    mod = DenseAdditiveLinearGenomicModel.from_hdf5("gmod.h5")

Copying Genomic Models
======================

Genomic models may be copied using either shallow copying or deep copying. The subsections below demonstrate both copying methods.

Shallow copying
---------------

.. |link_copy_copy| replace:: ``copy.copy``
.. _link_copy_copy: https://docs.python.org/3/library/copy.html#copy.copy

In shallow copying, references for a ``GenomicModel`` s internal data are copied to a new object. Copying may be accomplished using the ``copy`` method or the base Python function |link_copy_copy|_.

.. code-block:: python

    # copy a genomic model
    tmp = copy.copy(algmod)
    tmp = algmod.copy()

Deep copying
------------

.. |link_copy_deepcopy| replace:: ``copy.deepcopy``
.. _link_copy_deepcopy: https://docs.python.org/3/library/copy.html#copy.deepcopy

In deep copying, a ``GenomicModel``'s internal data is recursively copied to a new object. Deep copying may be accomplished using the ``deepcopy`` method or the base Python function |link_copy_deepcopy|_.

.. code-block:: python

    # deep copy a genomic model
    tmp = copy.deepcopy(algmod)
    tmp = algmod.deepcopy()

Model prediction methods
========================

Phenotypic predictions can be made with ``GenomicModel`` objects using the ``predict_numpy`` and ``predict`` methods. The former method accepts raw NumPy arrays, while the latter method accepts raw NumPy arrays and other object classes such as ``GenotypeMatrix`` classes. The code below demonstrates the use of these methods.

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

The prediction accuracy of a ``GenomicModel`` can be assessed using the ``predict_numpy`` and ``predict`` methods. The former method accepts raw NumPy arrays, while the latter method accepts raw NumPy arrays and other object classes such as ``GenotypeMatrix`` classes. The code below demonstrates the use of these methods.

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

Genomic estimated breeding values (GEBVs) can be estimated by a ``GenomicModel`` object using the ``gebv_numpy`` and ``gebv`` methods. The former method accepts raw NumPy arrays, while the latter method accepts raw NumPy arrays and ``GenotypeMatrix`` classes. The code below demonstrates the use of these methods.

.. code-block:: python

    # predict GEBVs using numpy arrays
    out = algmod.gebv_numpy(Z)

    # predict GEBVs using objects
    out = algmod.gebv(gmat)

Predicting genomic estimated genotypic values
=============================================

Genomic estimated genotypic values (GEGVs) can be estimated by a ``GenomicModel`` object using the ``gegv_numpy`` and ``gegv`` methods. The former method accepts raw NumPy arrays, while the latter method accepts raw NumPy arrays and ``GenotypeMatrix`` classes. The code below demonstrates the use of these methods.

.. code-block:: python

    # predict GEGVs using numpy arrays
    out = algmod.gegv_numpy(Z)

    # predict GEGVs using objects
    out = algmod.gegv(gmat)


Calculating population genetic variance terms
=============================================

Genomic models may be used to calculate various population genetic variance statistics. These statistics include genetic variance, additive genetic variance, additive genic variance, and the Bulmer effect for a population. The subsections detail the calculation of these metrics.

Predicting genetic variance
---------------------------

Population genetic variance can be calculated using the ``var_G_numpy`` and ``var_G`` methods. The code below demonstrates their use.

.. code-block:: python

    # predict Var(G) using numpy arrays
    out = algmod.var_G_numpy(Z)

    # predict Var(G) using objects
    out = algmod.var_G(gmat)

Predicting additive genetic variance
------------------------------------

Population additive genetic variance can be calculated using the ``var_A_numpy`` and ``var_A`` methods. The code below demonstrates their use.

.. code-block:: python

    # predict Var(A) using numpy arrays
    out = algmod.var_A_numpy(Z)

    # predict Var(A) using objects
    out = algmod.var_A(gmat)

Predicting additive genic variance
----------------------------------

Population additive genic variance (the genetic variance assuming complete linkage equilibrium) can be calculated using the ``var_a_numpy`` and ``var_a`` methods. The code below demonstrates their use.

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

The Bulmer effect for a population can be calculated using the ``bulmer_numpy`` and ``bulmer`` methods. The code below demonstrates their use.

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

Population selection limits are an important measure of how much valuable genetic variation can be extracted from a given population. ``GenomicModel`` classes provide both upper and lower selection limit calculations which are detailed below.

Upper selection limit
---------------------

Upper selection limits (USLs) can be calculated using the ``usl_numpy`` and ``usl`` methods. The upper selection limit represents the maximum breeding value achieved by the fixation of all favorable alleles in a given population. The code below demonstrates how to calculate the USL.

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

Lower selection limits (USLs) can be calculated using the ``lsl_numpy`` and ``lsl`` methods. The lower selection limit represents the minimum breeding value achieved by the fixation of all deleterious alleles in a given population. The code below demonstrates how to calculate the LSL.

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

Metrics on the favorable alleles within a population may be of importance when conducting breeding simulations. The favorable allele count across each locus for a population can be calculated using the ``facount`` method. The favorable allele frequency across each locus for a population can be calculated using the ``fafreq`` method. A mask of the favorable allele availability across each locus for a population can be calculated using the ``faavail`` method. Finally, a mask of the favorable allele fixation across each locus for a population can be calculated using the ``fafixed`` method. The code below demonstrates the use of these four methods.

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

Metrics on the deleterious alleles within a population may be of importance when conducting breeding simulations. The deleterious allele count across each locus for a population can be calculated using the ``facount`` method. The deleterious allele frequency across each locus for a population can be calculated using the ``fafreq`` method. A mask of the deleterious allele availability across each locus for a population can be calculated using the ``faavail`` method. Finally, a mask of the deleterious allele fixation across each locus for a population can be calculated using the ``fafixed`` method. The code below demonstrates the use of these four methods.

.. code-block:: python

    # calculate deleterious allele counts
    out = algmod.dacount(gmat)

    # calculate deleterious allele frequencies
    out = algmod.dafreq(gmat)

    # calculate deleterious allele availability at loci in a population
    out = algmod.daavail(gmat)

    # calculate deleterious allele fixation at loci in a population
    out = algmod.dafixed(gmat)

Exporting Genomic Models
========================

``GenomicModel`` objects can be exported for later use or import into another program. All ``GenomicModel`` objects can be saved as an HDF5 file. ``LinearGenomicModel`` objects may also be exported to Pandas DataFrames and CSV files.

Exporting to dictionaries of Pandas DataFrames
----------------------------------------------

``LinearGenomicModel`` objects can be exported to a dictionary of Pandas DataFrames using the ``to_pandas_dict`` method. The code below demonstrates the use of this function.

.. code-block:: python

    # export all trait columns
    out = algmod.to_pandas_dict()

Exporting to CSV files
----------------------

``LinearGenomicModel`` objects can be also be exported to a set of CSV files using the ``to_csv_dict`` method. The code below demonstrates the use of this function.

.. code-block:: python

    # create dictionary with filenames
    # need required fields of "beta", "u_misc", "u_a" for this class
    dic = {
        "beta": "beta.csv",
        "u_misc": "u_misc.csv",
        "u_a": "u_a.csv",
    }

    # export all trait columns to csv files
    algmod.to_csv_dict(dic)

Exporting to HDF5
-----------------

Finally, all ``GenomicModel`` objects can be saved to an HDF5 file for later reloading into PyBrOpS. The code below demonstrates how to export a ``GenomicModel`` to an HDF5 file.

.. code-block:: python

    # writing a genomic model to a file
    algmod.to_hdf5("saved_algmod.h5")
