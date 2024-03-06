Variance Matrices
#################

Class Family Overview
=====================

The ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` families of classes allow for representation of the expected progeny trait variances from a cross between individuals. The purpose of these families of classes is to calculate expected progeny trait variances assuming linkage for the ``GeneticVarianceMatrix`` family of classe and no linkage for the ``GenicVarianceMatrix`` family of classes. Both families of variance matrices utilize genomic models to calculate variances and are designed to be agnostic of ``GenomicModel`` type. For additive linear genomic models, there exist deterministic equations to calculate progeny variance for two-, three-, and four-way crosses. ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` families which assume an additive linear genomic model have interfaces of ``AdditiveGeneticVarianceMatrix`` and ``AdditiveGenicVarianceMatrix``, respectfully.

Deriving from the ``AdditiveGeneticVarianceMatrix`` interface are several implemented classes which are useful. They are summarized below:

Summary of Variance Matrix Classes
==================================

Genetic and genic variance matrices are found in the ``pybrops.model.vmat`` module. The two subsections below summarize the genetic and genic variance matrices

Genetic variance matrix classes
-------------------------------

Genetic variance matrices, which represent the expected progeny trait variances assuming linkage disequilibrium, are summarized in the table below.

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

Genic variance matrices, which represent the expected progeny trait variances assuming linkage equilibrium, are summarized in the table below.

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

Variance Matrix Properties
==========================

Both genetic and genic variance matrix classes share a common set of properties. These properties can be grouped into four categories: general properties, taxa properties, trait properties, and square properties. The following four subsections summarize these property categories.

General properties
------------------

Genetic and genic variance matrices share several shape properties that are common to all ``Matrix`` classes. In addition, genetic and genic variance matrices have an ``epgc`` property which is used to indicate the expected parental genomic contribution of each parent corresponding to each parental axis. The table below summarizes variance matrix general properties.

.. list-table:: Summary of ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` general properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``mat``
      - The raw genetic or genic variance matrix pointer
    * - ``mat_ndim``
      - The number of dimensions for the genetic or genic variance matrix
    * - ``mat_shape``
      - The genetic or genic variance matrix shape
    * - ``epgc``
      - The expected parental genomic contribution corresponding to each parental axis.

Taxa properties
---------------

Genetic and genic variance matrices have several taxa related properties including taxa names, taxa group identities, and sorting metadata, which can be used for quick group access and sorting. The table below summarizes variance matrix taxa properties.

.. list-table:: Summary of ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` taxa properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ntaxa``
      - The number of taxa represented by the genetic or genic variance matrix
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

Genetic and genic variance matrices have several trait related properties of which the most important is the trait names.

.. list-table:: Summary of ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` trait properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ntrait``
      - The number of traits represented by the genetic or genic variance matrix
    * - ``trait``
      - The names of the traits
    * - ``trait_axis``
      - The matrix axis along which traits are stored

Square matrix properties
------------------------

Since genetic and genic matrices are square by nature, they also have several properties which extract data regarding their squareness. These properties are summarized below.

.. list-table:: Summary of ``GeneticVarianceMatrix`` and ``GenicVarianceMatrix`` square matrix properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nsquare``
      - The number of square axes for the genetic or genic variance matrix
    * - ``square_axes``
      - The axes indices for the square axes for the genetic or genic variance matrix
    * - ``square_axes_len``
      - The lengths of the square axes for the genetic or genic variance matrix

Loading Variance Matrix Modules
===============================

Loading genetic variance matrix modules
---------------------------------------

Importing genetic variance matrix classes can be accomplished using the following import statements:

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

Importing genic variance matrix classes can be accomplished using the following import statements:

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

Genetic and genic variance matrices can be created using several method including from raw NumPy arrays, from genotype matrices and genomic models, from Pandas DataFrames, from CSV files, and from HDF5 files. The following subsections detail the creation or loading of genetic or genic variance matrices from their corresponding sources.

Creating variance matrices from NumPy arrays
--------------------------------------------

Using the constructor of a ``GeneticVarianceMatrix`` or ``GenicVarianceMatrix`` class, one can create genetic and genic variance matrices, respectively, from NumPy arrays. The example below demonstrates the creation of a ``DenseTwoWayDHAdditiveGeneticVarianceMatrix`` and ``DenseTwoWayDHAdditiveGenicVarianceMatrix``` objects from raw NumPy arrays.

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

Calculating variance matrices from genomic models
-------------------------------------------------

Genetic and genic variance matrices may also be calculated from a combination of ``PhasedGenotypeMatrix`` and ``GenomicModel`` objects. This can be accomplished using the ``from_gmod`` or ``from_algmod`` class methods. The code below demonstrates how to use this method to accomplish this task.

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
        nmating = 1,
        nprogeny = 10,
        nself = 0,
        gmapfn = HaldaneMapFunction()
    )

    # calculate genetic variance matrix from AdditiveLinearGenomicModel
    vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_algmod(
        algmod = algmod,
        pgmat = pgmat,
        nmating = 1,
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

Creating variance matrices from Pandas DataFrames
-------------------------------------------------

Genetic and genic variance matrices may be read from Pandas DataFrames. The ``from_pandas`` class method may be used to read a ``GeneticVarianceMatrix`` or ``GenicVarianceMatrix`` from a Pandas DataFrame. The code example below demonstrates this method's usage.

.. code-block:: python

    # create dummy pandas dataframe
    df = pandas.DataFrame({
        "Female": ["A","A","B","B",],
        "Female Group": [0,0,1,1,],
        "Male": ["A","B","A","B",],
        "Male Group": [0,0,1,1,],
        "Trait": ["Yield","Yield","Yield","Yield",],
        "Variance": [0.2,0.3,0.5,0.7],
    })

    # construct genetic variance matrix from pandas dataframe
    tmp_vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_pandas(
        df = df,
        female_col = "Female",
        female_grp_col = "Female Group",
        male_col = "Male",
        male_grp_col = "Male Group",
        trait_col = "Trait",
        variance_col = "Variance",
    )

    # construct genic variance matrix from pandas dataframe
    tmp_gvmat = DenseTwoWayDHAdditiveGenicVarianceMatrix.from_pandas(
        df = df,
        female_col = "Female",
        female_grp_col = "Female Group",
        male_col = "Male",
        male_grp_col = "Male Group",
        trait_col = "Trait",
        variance_col = "Variance",
    )

Loading variance matrices from CSV files
----------------------------------------

Genetic and genic variance matrices may also be read from CSV files in a manner similar to Pandas DataFrames. The ``from_csv`` class method may be used to load genetic or genic variance matrices from CSV files. The following code block demonstrates the usage of this method.

.. code-block:: python

    # create dummy pandas dataframe and export as CSV
    df = pandas.DataFrame({
        "Female": ["A","A","B","B",],
        "Female Group": [0,0,1,1,],
        "Male": ["A","B","A","B",],
        "Male Group": [0,0,1,1,],
        "Trait": ["Yield","Yield","Yield","Yield",],
        "Variance": [0.2,0.3,0.5,0.7],
    })
    df.to_csv("saved_df.csv", index = False)

    # construct genetic variance matrix from csv
    tmp_vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_csv(
        filename = "saved_df.csv",
        female_col = "Female",
        female_grp_col = "Female Group",
        male_col = "Male",
        male_grp_col = "Male Group",
        trait_col = "Trait",
        variance_col = "Variance",
    )

    # construct genic variance matrix from csv
    tmp_gvmat = DenseTwoWayDHAdditiveGenicVarianceMatrix.from_csv(
        filename = "saved_df.csv",
        female_col = "Female",
        female_grp_col = "Female Group",
        male_col = "Male",
        male_grp_col = "Male Group",
        trait_col = "Trait",
        variance_col = "Variance",
    )

Loading variance matrices from HDF5 files
-----------------------------------------

As with all classes in the ``Matrix`` family, ``GeneticVarianceMatrix`` or ``GenicVarianceMatrix`` objects may be imported and exported to an HDF5 format. To read saved genetic or genic variance matrices from an HDF5 file, use the ``from_hdf5`` class method. The code below demonstrates the use of this method.

.. code-block:: python

    # read genetic variance matrix from HDF5 file
    vmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_hdf5("saved_vmat.h5")
    gvmat = DenseTwoWayDHAdditiveGenicVarianceMatrix.from_hdf5("saved_gvmat.h5")

Copying Variance Matrices
=========================

Genetic and genic variance matrices may be copied using two methods: shallow copying and deep copying.

Shallow copying
---------------

.. |link_copy_copy| replace:: ``copy.copy``
.. _link_copy_copy: https://docs.python.org/3/library/copy.html#copy.copy

In shallow copying, references to a ``GeneticVarianceMatrix`` or ``GenicVarianceMatrix``'s data are copied to a new genetic or genic variance matrix object. Copying is only one level deep which means that changes to the original object may affect data values in the copied object. The code below illustrates the use of the ``copy`` method bound to ``GeneticVarianceMatrix`` or ``GenicVarianceMatrix`` objects and the base Python function |link_copy_copy|_ which can both be used to shallow copy a genetic or genic variance matrix object.

.. code-block:: python

    # copy a genetic variance matrix
    tmp = copy.copy(vmat)
    tmp = vmat.copy()

    # copy a genic variance matrix
    tmp = copy.copy(gvmat)
    tmp = gvmat.copy()

Deep copying
------------

.. |link_copy_deepcopy| replace:: ``copy.deepcopy``
.. _link_copy_deepcopy: https://docs.python.org/3/library/copy.html#copy.deepcopy

In deep copying, data in a ``GeneticVarianceMatrix`` or ``GenicVarianceMatrix`` is recursively copied to a new genetic or genic variance matrix object. Copying occurs down to the deepest levels so that changes to the original object will not affect data values in the copied object. The code below illustrates the use of the ``deepcopy`` method bound to ``GeneticVarianceMatrix`` or ``GenicVarianceMatrix`` objects and the base Python function |link_copy_deepcopy|_ which can both be used to deep copy a genetic or genic variance matrix object.

.. code-block:: python

    # deep copy a genetic variance matrix
    tmp = copy.deepcopy(vmat)
    tmp = vmat.deepcopy()

    # deep copy a genic variance matrix
    tmp = copy.deepcopy(gvmat)
    tmp = gvmat.deepcopy()

Copy-On Element Manipulation
============================

Genetic or genic variance matrices have several methods by which modifed copies of the original matrix can be made. These are called copy-on element manipulation methods. Matrices may have taxa and trait axes adjoined, deleted, inserted, or selected. The following sections demonstrate the use of these method families.

Adjoining elements
------------------

The ``adjoin`` family of methods allows for taxa and trait axes of a genetic or genic variance matrix to be adjoined together, creating a new matrix in the process. Use of the ``adjoin`` method family is demonstrated in the code below.

.. code-block:: python

    # create a new variance matrices to demonstrate
    newvmat = vmat.deepcopy()
    newgvmat = gvmat.deepcopy()

    # adjoin variance matrices along the taxa axis
    tmp = vmat.adjoin(newvmat, axis = vmat.taxa_axis)
    tmp = vmat.adjoin_taxa(newvmat)
    tmp = gvmat.adjoin(newgvmat, axis = gvmat.taxa_axis)
    tmp = gvmat.adjoin_taxa(newgvmat)

    # adjoin variance matrices along the trait axis
    tmp = vmat.adjoin(newvmat, axis = vmat.trait_axis)
    tmp = vmat.adjoin_trait(newvmat)
    tmp = gvmat.adjoin(newgvmat, axis = gvmat.trait_axis)
    tmp = gvmat.adjoin_trait(newgvmat)

Deleting elements
-----------------

The ``delete`` family of methods allows for taxa and trait axes of a genetic or genic variance matrix to be removed in a copy of the original. Use of the ``delete`` method family is demonstrated in the code below.


``delete`` taxa
+++++++++++++++

.. code-block:: python

    # delete first taxon using an integer
    tmp = vmat.delete(0, axis = vmat.taxa_axis)
    tmp = vmat.delete_taxa(0)
    tmp = gvmat.delete(0, axis = gvmat.taxa_axis)
    tmp = gvmat.delete_taxa(0)

    # delete first five taxa using a slice
    tmp = vmat.delete(slice(0,5), axis = vmat.taxa_axis)
    tmp = vmat.delete_taxa(slice(0,5))
    tmp = gvmat.delete(slice(0,5), axis = gvmat.taxa_axis)
    tmp = gvmat.delete_taxa(slice(0,5))

    # delete first five taxa using a Sequence
    tmp = vmat.delete([0,1,2,3,4], axis = vmat.taxa_axis)
    tmp = vmat.delete_taxa([0,1,2,3,4])
    tmp = gvmat.delete([0,1,2,3,4], axis = gvmat.taxa_axis)
    tmp = gvmat.delete_taxa([0,1,2,3,4])

``delete`` traits
+++++++++++++++++

.. code-block:: python

    # delete first trait using an integer
    tmp = vmat.delete(0, axis = vmat.trait_axis)
    tmp = vmat.delete_trait(0)
    tmp = gvmat.delete(0, axis = gvmat.trait_axis)
    tmp = gvmat.delete_trait(0)

    # delete first two traits using a slice
    tmp = vmat.delete(slice(0,2), axis = vmat.trait_axis)
    tmp = vmat.delete_trait(slice(0,2))
    tmp = gvmat.delete(slice(0,2), axis = gvmat.trait_axis)
    tmp = gvmat.delete_trait(slice(0,2))

    # delete first two traits using a Sequence
    tmp = vmat.delete([0,1], axis = vmat.trait_axis)
    tmp = vmat.delete_trait([0,1])
    tmp = gvmat.delete([0,1], axis = gvmat.trait_axis)
    tmp = gvmat.delete_trait([0,1])

Inserting elements
------------------

The ``insert`` family of methods allows for taxa and trait axes of a genetic or genic variance matrix to be inserted into a copy of the original matrix. Use of the ``insert`` method family is demonstrated in the code below.

.. code-block:: python

    # create a new variance matrix to demonstrate
    newvmat = vmat.deepcopy()
    newgvmat = gvmat.deepcopy()

    # insert variance matrix along the taxa axis before index 0
    tmp = vmat.insert(0, newvmat, axis = vmat.taxa_axis)
    tmp = vmat.insert_taxa(0, newvmat)
    tmp = gvmat.insert(0, newgvmat, axis = gvmat.taxa_axis)
    tmp = gvmat.insert_taxa(0, newgvmat)

    # insert variance matrix along the trait axis before index 0
    tmp = vmat.insert(0, newvmat, axis = vmat.trait_axis)
    tmp = vmat.insert_trait(0, newvmat)
    tmp = gvmat.insert(0, newgvmat, axis = gvmat.trait_axis)
    tmp = gvmat.insert_trait(0, newgvmat)

Selecting elements
------------------

The ``select`` family of methods allows for taxa and trait axes of the genetic or genic variance matrix to be selected and extracted to a copy of the original matrix. Use of the ``select`` method family is demonstrated in the code below.

.. code-block:: python

    # select first five taxa using a Sequence
    tmp = vmat.select([0,1,2,3,4], axis = vmat.taxa_axis)
    tmp = vmat.select_taxa([0,1,2,3,4])
    tmp = gvmat.select([0,1,2,3,4], axis = gvmat.taxa_axis)
    tmp = gvmat.select_taxa([0,1,2,3,4])

    # select first two traits using a Sequence
    tmp = vmat.select([0,1], axis = vmat.trait_axis)
    tmp = vmat.select_trait([0,1])
    tmp = gvmat.select([0,1], axis = gvmat.trait_axis)
    tmp = gvmat.select_trait([0,1])

In-Place Element Manipulation
=============================

Genetic or genic variance matrices have several methods which execute in-place element manipulations. These are called in-place element manipulation methods. Genetic or genic variance matrices may have taxa and trait axes appended, removed, incorporated, or concatenated. The following sections demonstrate the use of these method families.

Appending elements
------------------

The ``append`` family of methods allows for new taxa and trait axes to be appended to the genetic or genic variance matrix. The code segment below demonstrates their use. 

.. code-block:: python

    # append variance matrices along the taxa axis
    tmp = vmat.deepcopy()                   # copy original
    tmp.append(vmat, axis = tmp.taxa_axis)  # append original to copy
    tmp = gvmat.deepcopy()                  # copy original
    tmp.append(gvmat, axis = tmp.taxa_axis) # append original to copy

    tmp = vmat.deepcopy()                   # copy original
    tmp.append_taxa(vmat)                   # append original to copy
    tmp = gvmat.deepcopy()                  # copy original
    tmp.append_taxa(gvmat)                  # append original to copy

    # append variance matrices along the trait axis
    tmp = vmat.deepcopy()                   # copy original
    tmp.append(vmat, axis = tmp.trait_axis) # append original to copy
    tmp = gvmat.deepcopy()                  # copy original
    tmp.append(gvmat, axis = tmp.trait_axis)# append original to copy

    tmp = vmat.deepcopy()                   # copy original
    tmp.append_trait(vmat)                  # append original to copy
    tmp = gvmat.deepcopy()                  # copy original
    tmp.append_trait(gvmat)                 # append original to copy

Removing elements
-----------------

The ``remove`` family of methods allows for taxa and trait axes to be removed from a genetic or genic variance matrix. A demonstration of their use can be seen below. 

``remove`` taxa
+++++++++++++++

.. code-block:: python

    # remove first taxon using an integer
    tmp = vmat.deepcopy()                           # copy original
    tmp.remove(0, axis = vmat.taxa_axis)            # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove(0, axis = gvmat.taxa_axis)           # remove from copy

    tmp = vmat.deepcopy()                           # copy original
    tmp.remove_taxa(0)                              # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove_taxa(0)                              # remove from copy

    # remove first five taxa using a slice
    tmp = vmat.deepcopy()                           # copy original
    tmp.remove(slice(0,5), axis = vmat.taxa_axis)   # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove(slice(0,5), axis = gvmat.taxa_axis)  # remove from copy

    tmp = vmat.deepcopy()                           # copy original
    tmp.remove_taxa(slice(0,5))                     # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove_taxa(slice(0,5))                     # remove from copy

    # remove first five taxa using a Sequence
    tmp = vmat.deepcopy()                           # copy original
    tmp.remove([0,1,2,3,4], axis = vmat.taxa_axis)  # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove([0,1,2,3,4], axis = gvmat.taxa_axis) # remove from copy

    tmp = vmat.deepcopy()                           # copy original
    tmp.remove_taxa([0,1,2,3,4])                    # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove_taxa([0,1,2,3,4])                    # remove from copy

``remove`` traits
+++++++++++++++++

.. code-block:: python

    # remove first trait using an integer
    tmp = vmat.deepcopy()                           # copy original
    tmp.remove(0, axis = vmat.trait_axis)           # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove(0, axis = gvmat.trait_axis)          # remove from copy

    tmp = vmat.deepcopy()                           # copy original
    tmp.remove_trait(0)                             # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove_trait(0)                             # remove from copy

    # remove first trait using a slice
    tmp = vmat.deepcopy()                           # copy original
    tmp.remove(slice(0,1), axis = vmat.trait_axis)  # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove(slice(0,1), axis = gvmat.trait_axis) # remove from copy

    tmp = vmat.deepcopy()                           # copy original
    tmp.remove_trait(slice(0,1))                    # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove_trait(slice(0,1))                    # remove from copy

    # remove first trait using a Sequence
    tmp = vmat.deepcopy()                           # copy original
    tmp.remove([0], axis = vmat.trait_axis)         # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove([0], axis = gvmat.trait_axis)        # remove from copy

    tmp = vmat.deepcopy()                           # copy original
    tmp.remove_trait([0])                           # remove from copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.remove_trait([0])                           # remove from copy

Incorporating elements
----------------------

The ``incorp`` family of methods allows for new taxa and trait axes to be inserted at specific locations in a genetic or genic variance matrix. Use of the ``incorp`` family is demonstrated in the code segment below below. 

.. code-block:: python

    # incorp variance matrix along the taxa axis before index 0
    tmp = vmat.deepcopy()                           # copy original
    tmp.incorp(0, vmat, axis = vmat.taxa_axis)      # incorporate into copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.incorp(0, gvmat, axis = gvmat.taxa_axis)    # incorporate into copy

    tmp = vmat.deepcopy()                           # copy original
    tmp.incorp_taxa(0, vmat)                        # incorporate into copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.incorp_taxa(0, gvmat)                       # incorporate into copy

    # incorp variance matrix along the trait axis before index 0
    tmp = vmat.deepcopy()                           # copy original
    tmp.incorp(0, vmat, axis = vmat.trait_axis)     # incorporate into copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.incorp(0, gvmat, axis = gvmat.trait_axis)   # incorporate into copy

    tmp = vmat.deepcopy()                           # copy original
    tmp.incorp_trait(0, vmat)                       # incorporate into copy
    tmp = gvmat.deepcopy()                          # copy original
    tmp.incorp_trait(0, gvmat)                      # incorporate into copy

Concatenating matrices
----------------------

The ``concat`` family of methods allows for multiple genetic or genic variance matrices to be concatenated to each other along taxa or trait axes. The code segment below demonstrates their use. 

.. code-block:: python

    # concatenate along the taxa axis
    tmp = vmat.concat([vmat, vmat], axis = vmat.taxa_axis)
    tmp = vmat.concat_taxa([vmat, vmat])
    tmp = gvmat.concat([gvmat, gvmat], axis = gvmat.taxa_axis)
    tmp = gvmat.concat_taxa([gvmat, gvmat])

    # concatenate along the trait axis
    tmp = vmat.concat([vmat, vmat], axis = vmat.trait_axis)
    tmp = vmat.concat_trait([vmat, vmat])
    tmp = gvmat.concat([gvmat, gvmat], axis = gvmat.trait_axis)
    tmp = gvmat.concat_trait([gvmat, gvmat])

Grouping and Sorting
====================

Genetic or genic variance matrices in PyBrOpS have several sorting and grouping focused methods. Sorting methods can be used to reorder, sort, and group taxa axes alphanumerically and reorder and sort trait axes.. The following sections demonstrate the use of the ``reorder``, ``lexsort``, ``sort``, and ``group`` method families.

Reordering elements
-------------------

Taxa and trait axes in a genetic or genic variance matrix can be reordered using the ``reorder`` family of methods. Demonstrations of this method family are below.

``reorder`` taxa
++++++++++++++++

.. code-block:: python

    # create reordering indices
    indices = numpy.arange(vmat.ntaxa)
    numpy.random.shuffle(indices)

    # reorder values along the taxa axis
    tmp = vmat.deepcopy()
    tmp.reorder(indices, axis = tmp.taxa_axis)
    tmp.reorder_taxa(indices)
    tmp = gvmat.deepcopy()
    tmp.reorder(indices, axis = tmp.taxa_axis)
    tmp.reorder_taxa(indices)

``reorder`` traits
++++++++++++++++++

.. code-block:: python

    # create reordering indices
    indices = numpy.arange(vmat.ntrait)
    numpy.random.shuffle(indices)

    # reorder values along the trait axis
    tmp = vmat.deepcopy()
    tmp.reorder(indices, axis = tmp.trait_axis)
    tmp.reorder_trait(indices)
    tmp = gvmat.deepcopy()
    tmp.reorder(indices, axis = tmp.trait_axis)
    tmp.reorder_trait(indices)

Lexsorting elements
-------------------

An indirect stable sort - or lexsort - for taxa or trait axes can be performed using the ``lexsort`` family of methods. The code segment below illustrates the use of this family of methods.

``lexsort`` taxa
++++++++++++++++

.. code-block:: python

    # create lexsort keys for taxa
    key1 = numpy.random.randint(0, 10, vmat.ntaxa)
    key2 = numpy.arange(vmat.ntaxa)
    numpy.random.shuffle(key2)

    # lexsort along the taxa axis
    vmat.lexsort((key2,key1), axis = vmat.taxa_axis)
    vmat.lexsort_taxa((key2,key1))
    gvmat.lexsort((key2,key1), axis = gvmat.taxa_axis)
    gvmat.lexsort_taxa((key2,key1))

``lexsort`` traits
++++++++++++++++++

.. code-block:: python

    # create lexsort keys for trait
    key1 = numpy.random.randint(0, 10, vmat.ntaxa)
    key2 = numpy.arange(vmat.ntaxa)
    numpy.random.shuffle(key2)

    # lexsort along the trait axis
    tmp = vmat.lexsort((key2,key1), axis = vmat.taxa_axis)
    tmp = vmat.lexsort_taxa((key2,key1))
    tmp = gvmat.lexsort((key2,key1), axis = gvmat.taxa_axis)
    tmp = gvmat.lexsort_taxa((key2,key1))

Sorting elements
----------------

Alphanumeric sorting along taxa and trait axes can be done using the ``sort`` family of methods. Sorting examples are illustrated below.

``sort`` taxa
+++++++++++++

.. code-block:: python

    # sort along taxa axis
    tmp = vmat.deepcopy()
    tmp.sort(axis = tmp.taxa_axis)
    tmp.sort_taxa()
    tmp = gvmat.deepcopy()
    tmp.sort(axis = tmp.taxa_axis)
    tmp.sort_taxa()

``sort`` traits
+++++++++++++++

.. code-block:: python

    # sort along trait axis
    tmp = vmat.deepcopy()
    tmp.sort(axis = tmp.trait_axis)
    tmp.sort_trait()
    tmp = gvmat.deepcopy()
    tmp.sort(axis = tmp.trait_axis)
    tmp.sort_trait()

Grouping elements
-----------------

Grouping along taxa axes can be done using the ``group`` family of methods. The following code illustrates the use of the ``group`` method family along the taxa axes of a genetic or genic variance matrix.

.. code-block:: python

    #
    # group taxa
    # ++++++++++

    # sort genetic along taxa axis
    tmp = vmat.deepcopy()
    tmp.group(axis = tmp.taxa_axis)
    tmp.group_taxa()
    # determine whether grouping has occurred along the taxa axis
    out = tmp.is_grouped(axis = tmp.taxa_axis)
    out = tmp.is_grouped_taxa()

    # sort genic variance matrix along taxa axis
    tmp = gvmat.deepcopy()
    tmp.group(axis = tmp.taxa_axis)
    tmp.group_taxa()
    # determine whether grouping has occurred along the taxa axis
    out = tmp.is_grouped(axis = tmp.taxa_axis)
    out = tmp.is_grouped_taxa()

Square matrix functions
=======================

Since variance matrices are also inherently square, there are a methods dealing with the squareness of a variance matrix.

Determine whether all square axes are of equal length
-----------------------------------------------------

The ``is_square`` method for genetic and genic variance matrices can be used to determine if a variance matrix is square.

.. code-block:: python

    # boolean value
    out = vmat.is_square()
    out = gvmat.is_square()

Exporting Variance Matrices
===========================

Genetic and genic variance matrices may be exported to multiple formats including Pandas DataFrames, CSV files, and HDF5 files. The following subsections provide export examples.

Exporting to Pandas DataFrame
-----------------------------

The ``to_pandas`` method can be used to export a genetic or genic variance matrix to a Pandas DataFrame. Column names may be optionally provided to override default column names.

.. code-block:: python

    # write variance matrices to a pandas dataframe
    out_df = vmat.to_pandas()
    out_df = gvmat.to_pandas()

Exporting to CSV
----------------

The ``to_csv`` method can be used to export a genetic or genic variance matrix to a CSV file. Like the ``to_pandas`` method, column names may be optionally provided to override default column names.

.. code-block:: python

    # write variance matrices to a CSV file
    vmat.to_csv("saved_vmat.csv")
    gvmat.to_csv("saved_gvmat.csv")

Exporting to HDF5
------------------

Genetic and genic variance matrices can be written to an HDF5 file using the ``to_hdf5`` method. The code below demonstrates the use of this method.

.. code-block:: python

    # write variance matrices to an HDF5 file
    vmat.to_hdf5("saved_vmat.h5")
    gvmat.to_hdf5("saved_gvmat.h5")
