Breeding Value Matrices
#######################

Class Family Overview
=====================

The ``BreedingValueMatrix`` family of classes is used to represent breeding values as its name implies. ``BreedingValueMatrix`` objects can be used in the estimation of genomic prediction models and to make selection decisions. Since breeding values are typically mean-centered and sometimes scaled, breeding value matrices have ``location`` and ``scale`` properties to reconstitute un-scaled values. ``BreedingValueMatrix`` objects store additional taxa and trait metadata which serve as labels for rows and columns, respectively.

Summary of Breeding Value Matrix Classes
========================================

Breeding value matrix classes in PyBrOpS are found in the ``pybrops.popgen.bvmat`` module. Within this module are several ``BreedingValueMatrix`` class definitions which are summarized in the table below.

.. list-table:: Summary of breeding value matrix classes in the ``pybrops.popgen.bvmat`` module
    :widths: 25 15 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``BreedingValueMatrix``
      - Abstract
      - Interface for all breeding value matrix child classes.
    * - ``DenseBreedingValueMatrix``
      - Concrete
      - Class representing dense, breeding value matrices.

Breeding Value Matrix Properties
================================

Breeding value matrices have numerous properties. These properties can be grouped into three main groupings: general properties, taxa properties, and trait properties. Tables summarizing these properties can be read below.

Breeding value matrix general properties
----------------------------------------

Breeding value matrices share several shape properties that are common to all ``Matrix`` classes. In addition, breeding value matrices have ``location`` and ``scale`` properties which specify the trait mean and standard deviation, respectively. These properties are helpful if a breeding value matrix has been centered around zero and scaled to a unit standard deviation.

.. list-table:: Summary of ``BreedingValueMatrix`` general properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``mat``
      - The raw breeding value matrix pointer
    * - ``mat_ndim``
      - The number of dimensions for the breeding value matrix
    * - ``mat_shape``
      - The breeding value matrix shape
    * - ``location``
      - The location of the breeding value matrix if it has been transformed
    * - ``scale``
      - The scale of the breeding value matrix if it has been transformed

Breeding value matrix taxa properties
-------------------------------------

Breeding value matrices have several taxa related properties including taxa names, taxa group identities, and sorting metadata, which can be used for quick group access and sorting.

.. list-table:: Summary of ``BreedingValueMatrix`` taxa properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ntaxa``
      - The number of taxa represented by the breeding value matrix
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

Breeding value matrix trait properties
--------------------------------------

Breeding value matrices have several trait related properties of which the most important is the trait names.

.. list-table:: Summary of ``BreedingValueMatrix`` trait properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ntrait``
      - The number of traits represented by the breeding value matrix
    * - ``trait``
      - The names of the traits
    * - ``trait_axis``
      - The matrix axis along which traits are stored

Loading Breeding Value Matrix Modules
=====================================

Breeding value matrix classes can be imported as demonstrated in the code chunk below:

.. code-block:: python

    # import the BreedingValueMatrix class (an abstract interface class)
    from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix

    # import the DenseBreedingValueMatrix class (a concrete implemented class)
    from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix

Creating Breeding Value Matrices
================================

Breeding value matrices can be created using several method including from raw NumPy arrays, from Pandas DataFrames, from CSV files, and from HDF5 files. The following subsections detail the creation or loading of breeding value matrices from these sources.

Creating breeding value matrices from NumPy arrays
--------------------------------------------------

Using the ``DenseBreedingValueMatrix`` constructor, one can create a breeding value matrix from NumPy arrays.

.. code-block:: python

    # shape parameters
    ntaxa = 100
    ntrait = 3
    ngroup = 20

    # create random breeding values
    mat = numpy.random.normal(size = (ntaxa,ntrait))

    # create taxa names
    taxa = numpy.array(
        ["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], 
        dtype = object
    )

    # create taxa groups
    taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
    taxa_grp.sort()

    # create trait names
    trait = numpy.array(
        ["trait"+str(i+1).zfill(2) for i in range(ntrait)],
        dtype = object
    )

    # create a breeding value matrix from NumPy arrays
    bvmat = DenseBreedingValueMatrix(
        mat = mat,
        location = 0.0,
        scale = 1.0,
        taxa = taxa,
        taxa_grp = taxa_grp,
        trait = trait
    )

Using the ``from_numpy`` class method, one can also create a breeding value matrix from NumPy arrays. The difference between using this method and using the constructor is that this class method will automatically scale the input matrix to have zero mean and unit variance. Location and scale information will be stored in the ``location`` and ``scale`` properties of the created breeding value matrix.

.. code-block:: python

    # shape parameters
    ntaxa = 100
    ntrait = 3
    ngroup = 20

    # create random breeding values
    mat = numpy.random.normal(size = (ntaxa,ntrait))

    # create taxa names
    taxa = numpy.array(
        ["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], 
        dtype = object
    )

    # create taxa groups
    taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
    taxa_grp.sort()

    # create trait names
    trait = numpy.array(
        ["trait"+str(i+1).zfill(2) for i in range(ntrait)],
        dtype = object
    )

    bvmat = DenseBreedingValueMatrix.from_numpy(
        a = mat,
        taxa = taxa,
        taxa_grp = taxa_grp,
        trait = trait
    )

Creating breeding value matrices from Pandas DataFrames
-------------------------------------------------------

Breeding value matrices can be created from Pandas DataFrames. To do this, use the ``from_pandas`` class method. The code block below demonstrates how to use the ``from_pandas`` method to accomplish this.

.. code-block:: python

    # create dummy pandas dataframe
    df = pandas.DataFrame({
        "taxa": ["Taxon"+str(i).zfill(3) for i in range(1,101)],
        "taxa_grp": numpy.repeat([1,2,3,4,5], 20),
        "Trait1": numpy.random.random(100),
        "Trait2": numpy.random.random(100),
        "Trait3": numpy.random.random(100),
    })

    # construct breeding value matrix from pandas dataframe
    # use explicit column name identifiers as method arguments
    bvmat = DenseBreedingValueMatrix.from_pandas(
        df = df,
        location = 0.0,
        scale = 1.0,
        taxa_col = "taxa",
        taxa_grp_col = "taxa_grp",
        trait_cols = ["Trait1","Trait2","Trait3"],
    )

Loading breeding value matrices from CSV files
----------------------------------------------

Breeding value matrices can be read from CSV files. To read a breeding value matrix from a CSV file, use the ``from_csv`` class method. The following code illustrates the use of this method.

.. code-block:: python

    # read from a CSV file
    # use explicit column name identifiers as method arguments
    bvmat = DenseBreedingValueMatrix.from_csv(
        filename = "sample_breeding_values.csv",
        location = 0.0,
        scale = 1.0,
        taxa_col = "taxa",
        taxa_grp_col = "taxa_grp",
        trait_cols = ["Trait1","Trait2","Trait3"],
    )

Loading breeding value matrices from HDF5 files
-----------------------------------------------

Most matrix object types in PyBrOpS allow for both the import and export of matrices into an HDF5 format. To read saved breeding value matrices from an HDF5 file, use the ``from_hdf5`` class method. The code below demonstrates the use of this method to load a breeding value matrix from an HDF5 file.

.. code-block:: python

    # read a breeding value matrix from an HDF5 file
    bvmat = DenseBreedingValueMatrix.from_hdf5("sample_breeding_values.h5")


Copying Breeding Value Matrices
===============================

Copying breeding value matrices can be accomplished using two different methods: by shallow copying or by deep copying.

Shallow copying
---------------

.. |link_copy_copy| replace:: ``copy.copy``
.. _link_copy_copy: https://docs.python.org/3/library/copy.html#copy.copy

In shallow copying, references to a ``BreedingValueMatrix``'s data are copied to a new breeding value matrix object. Copying is only one level deep which means that changes to the original object may affect data values in the copied object. The code below illustrates the use of the ``copy`` method bound to ``BreedingValueMatrix`` objects and the base Python function |link_copy_copy|_ which can both be used to shallow copy a breeding value matrix object.

.. code-block:: python

    # copy a breeding value matrix
    tmp = copy.copy(bvmat)
    tmp = bvmat.copy()

Deep copying
------------

.. |link_copy_deepcopy| replace:: ``copy.deepcopy``
.. _link_copy_deepcopy: https://docs.python.org/3/library/copy.html#copy.deepcopy

In deep copying, data in a ``BreedingValueMatrix`` is recursively copied to a new breeding value matrix object. Copying occurs down to the deepest levels so that changes to the original object will not affect data values in the copied object. The code below illustrates the use of the ``deepcopy`` method bound to ``BreedingValueMatrix`` objects and the base Python function |link_copy_deepcopy|_ which can both be used to deep copy a breeding value matrix object.

.. code-block:: python

    # deep copy a breeding value matrix
    tmp = copy.deepcopy(bvmat)
    tmp = bvmat.deepcopy()


Copy-On Element Manipulation
============================

Breeding value matrices have several methods by which modifed copies of the original matrix can be made. These are called copy-on element manipulation methods. Matrices may have rows and/or columns adjoined, deleted, inserted, or selected. The following sections demonstrate the use of these method families.

Adjoin elements
---------------

The ``adjoin`` family of methods allows for rows (taxa) and columns (traits) of a breeding value matrix to be adjoined together, creating a new matrix in the process. Use of the ``adjoin`` method family is demonstrated in the code below.

.. code-block:: python

    # create a new breeding value matrix to demonstrate
    new = bvmat.deepcopy()

    # adjoin breeding value matrices along the taxa axis
    tmp = bvmat.adjoin(new, axis = bvmat.taxa_axis)
    tmp = bvmat.adjoin_taxa(new)

    # adjoin breeding value matrices along the trait axis
    tmp = bvmat.adjoin(new, axis = bvmat.trait_axis)
    tmp = bvmat.adjoin_trait(new)

Delete elements
---------------

The ``delete`` family of methods allows for rows (taxa) and columns (traits) of a breeding value matrix to be removed in a copy of the original. Use of the ``delete`` method family is demonstrated in the code below.

.. code-block:: python

    #
    # delete taxa examples
    #

    # delete first taxon using an integer
    tmp = bvmat.delete(0, axis = bvmat.taxa_axis)
    tmp = bvmat.delete_taxa(0)

    # delete first five taxa using a slice
    tmp = bvmat.delete(slice(0,5), axis = bvmat.taxa_axis)
    tmp = bvmat.delete_taxa(slice(0,5))

    # delete first five taxa using a Sequence
    tmp = bvmat.delete([0,1,2,3,4], axis = bvmat.taxa_axis)
    tmp = bvmat.delete_taxa([0,1,2,3,4])

    #
    # delete traits examples
    #

    # delete first trait using an integer
    tmp = bvmat.delete(0, axis = bvmat.trait_axis)
    tmp = bvmat.delete_trait(0)

    # delete first two traits using a slice
    tmp = bvmat.delete(slice(0,2), axis = bvmat.trait_axis)
    tmp = bvmat.delete_trait(slice(0,2))

    # delete first two traits using a Sequence
    tmp = bvmat.delete([0,1], axis = bvmat.trait_axis)
    tmp = bvmat.delete_trait([0,1])

Insert elements
---------------

The ``insert`` family of methods allows for rows (taxa) and columns (traits) of a breeding value matrix to be inserted into a copy of the original matrix. Use of the ``insert`` method family is demonstrated in the code below.

.. code-block:: python

    # create a new breeding value matrix to demonstrate
    new = bvmat.deepcopy()

    # insert breeding value matrix along the taxa axis before index 0
    tmp = bvmat.insert(0, new, axis = bvmat.taxa_axis)
    tmp = bvmat.insert_taxa(0, new)

    # insert breeding value matrix along the trait axis before index 0
    tmp = bvmat.insert(0, new, axis = bvmat.trait_axis)
    tmp = bvmat.insert_trait(0, new)

Select elements
---------------

The ``select`` family of methods allows for rows (taxa) and columns (traits) of the breeding value matrix to be selected and extracted to a copy of the original matrix. Use of the ``select`` method family is demonstrated in the code below.

.. code-block:: python

    # select first five taxa using a Sequence
    tmp = bvmat.select([0,1,2,3,4], axis = bvmat.taxa_axis)
    tmp = bvmat.select_taxa([0,1,2,3,4])

    # select first two traits using a Sequence
    tmp = bvmat.select([0,1], axis = bvmat.trait_axis)
    tmp = bvmat.select_trait([0,1])

In-Place Element Manipulation
=============================

Breeding value matrices have several methods which execute in-place element manipulations. These are called in-place element manipulation methods. Breeding value matrices may have taxa rows and/or trait columns appended, removed, incorporated, or concatenated. The following sections demonstrate the use of these method families.

Append elements
---------------

The ``append`` family of methods allows for new rows (taxa) and columns (traits) to be appended to the breeding value matrix. The code segment below demonstrates their use. 

.. code-block:: python

    # append breeding value matrices along the taxa axis
    tmp = bvmat.deepcopy()                   # copy original
    tmp.append(bvmat, axis = tmp.taxa_axis)  # append original to copy

    tmp = bvmat.deepcopy()                   # copy original
    tmp.append_taxa(bvmat)                   # append original to copy

    # append breeding value matrices along the trait axis
    tmp = bvmat.deepcopy()                   # copy original
    tmp.append(bvmat, axis = tmp.trait_axis) # append original to copy

    tmp = bvmat.deepcopy()                   # copy original
    tmp.append_trait(bvmat)                  # append original to copy

Remove elements
---------------

The ``remove`` family of methods allows for rows (taxa) and columns (traits) to be removed from a breeding value matrix. A demonstration of their use can be seen below. 

.. code-block:: python

    #
    # remove taxa examples
    #

    # remove first taxon using an integer
    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove(0, axis = bvmat.taxa_axis)            # remove from copy

    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove_taxa(0)                               # remove from copy

    # remove first five taxa using a slice
    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove(slice(0,5), axis = bvmat.taxa_axis)   # remove from copy

    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove_taxa(slice(0,5))                      # remove from copy

    # remove first five taxa using a Sequence
    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove([0,1,2,3,4], axis = bvmat.taxa_axis)  # remove from copy

    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove_taxa([0,1,2,3,4])                     # remove from copy

    #
    # remove traits examples
    #

    # remove first trait using an integer
    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove(0, axis = bvmat.trait_axis)           # remove from copy

    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove_trait(0)                              # remove from copy

    # remove first two traits using a slice
    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove(slice(0,2), axis = bvmat.trait_axis)  # remove from copy

    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove_trait(slice(0,2))                     # remove from copy

    # remove first two traits using a Sequence
    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove([0,1], axis = bvmat.trait_axis)       # remove from copy

    tmp = bvmat.deepcopy()                           # copy original
    tmp.remove_trait([0,1])                          # remove from copy

Incorporate elements
--------------------

The ``incorp`` family of methods allows for new rows (taxa) and columns (traits) to be inserted at specific locations a breeding value matrix. Use of the ``incorp`` family is demonstrated in the code segment below below. 

.. code-block:: python

    # incorp breeding value matrix along the taxa axis before index 0
    tmp = bvmat.deepcopy()                           # copy original
    tmp.incorp(0, bvmat, axis = bvmat.taxa_axis)     # incorporate into copy

    tmp = bvmat.deepcopy()                           # copy original
    tmp.incorp_taxa(0, bvmat)                        # incorporate into copy

    # incorp breeding value matrix along the trait axis before index 0
    tmp = bvmat.deepcopy()                           # copy original
    tmp.incorp(0, bvmat, axis = bvmat.trait_axis)    # incorporate into copy

    tmp = bvmat.deepcopy()                           # copy original
    tmp.incorp_trait(0, bvmat)                       # incorporate into copy

Concatenate elements
--------------------

The ``concat`` family of methods allows for multiple breeding value matrices to be concatenated to each other. The code segment below demonstrates their use. 

.. code-block:: python

    # concatenate along the taxa axis
    tmp = bvmat.concat([bvmat, bvmat], axis = bvmat.taxa_axis)
    tmp = bvmat.concat_taxa([bvmat, bvmat])

    # concatenate along the trait axis
    tmp = bvmat.concat([bvmat, bvmat], axis = bvmat.trait_axis)
    tmp = bvmat.concat_trait([bvmat, bvmat])

Grouping and Sorting
====================

Breeding value matrices in PyBrOpS have several sorting and grouping focused methods. Sorting methods can be used to reorder, sort, and group taxa alphanumerically, and reorder and sort traits alphanumerically. The following sections demonstrate the use of the ``reorder``, ``lexsort``, ``sort``, and ``group`` method families.

Reordering elements
-------------------

Taxa and traits in a breeding value matrix can be reordered using the ``reorder`` family of methods. Demonstrations of this method family are below.

.. code-block:: python

    #
    # taxa reordering example
    #

    # create reordering indices
    indices = numpy.arange(bvmat.ntaxa)
    numpy.random.shuffle(indices)
    tmp = bvmat.deepcopy()

    # reorder values along the taxa axis
    tmp.reorder(indices, axis = tmp.taxa_axis)
    tmp.reorder_taxa(indices)

    #
    # trait reordering example
    #

    # create reordering indices
    indices = numpy.arange(bvmat.ntrait)
    numpy.random.shuffle(indices)
    tmp = bvmat.deepcopy()

    # reorder values along the trait axis
    tmp = bvmat.deepcopy()
    tmp.reorder(indices, axis = tmp.trait_axis)
    tmp.reorder_trait(indices)

Lexsorting elements
-------------------

An indirect sort - or lexsort - for taxa and trait axes can be performed using the ``lexsort`` family of methods. The code segment below illustrates the use of this family of methods.

.. code-block:: python

    #
    # taxa lexsort example
    #

    # create lexsort keys for taxa
    key1 = numpy.random.randint(0, 10, bvmat.ntaxa)
    key2 = numpy.arange(bvmat.ntaxa)
    numpy.random.shuffle(key2)

    # lexsort along the taxa axis
    bvmat.lexsort((key2,key1), axis = bvmat.taxa_axis)
    bvmat.lexsort_taxa((key2,key1))

    #
    # trait lexsort example
    #

    # create lexsort keys for trait
    key1 = numpy.random.randint(0, 10, bvmat.ntaxa)
    key2 = numpy.arange(bvmat.ntaxa)
    numpy.random.shuffle(key2)

    # lexsort along the trait axis
    bvmat.lexsort((key2,key1), axis = bvmat.taxa_axis)
    bvmat.lexsort_taxa((key2,key1))

Sorting elements
----------------

Alphanumeric sorting along taxa and trait axes can be done using the ``sort`` family of methods. Sorting examples are illustrated below.

.. code-block:: python

    # make copy
    tmp = bvmat.deepcopy()

    #
    # taxa sorting example
    #

    # sort along taxa axis
    tmp.sort(axis = tmp.taxa_axis)
    tmp.sort_taxa()

    #
    # trait sorting example
    #

    # sort along trait axis
    tmp.sort(axis = tmp.trait_axis)
    tmp.sort_trait()

Grouping elements
-----------------

Grouping along only the taxa axis can be done using the ``group`` family of methods. The following code illustrates the use of the ``group`` method family along the taxa axis of a breeding value matrix.

.. code-block:: python

    # make copy
    tmp = bvmat.deepcopy()

    #
    # taxa grouping example
    #

    # sort and group along taxa axis
    tmp.group(axis = tmp.taxa_axis)
    tmp.group_taxa()

    # determine whether grouping has occurred along the taxa axis
    tmp.is_grouped(axis = tmp.taxa_axis)
    tmp.is_grouped_taxa()


Summary Statistics
==================

Various summary statistics can be calculated from breeding value matrices. PyBrOpS offers several common statistical routines which are described in the subsections below.

Maximum breeding values for each trait
--------------------------------------

The maximum breeding value for each trait may be calculated using the ``tmax`` method. The code below illustrates this method's use.

.. code-block:: python

    # get the maximum breeding values for each trait
    out = bvmat.tmax()

Row (taxa) indices of the individuals with the largest breeding values for each trait can be calculated using the ``targmax`` method. The code below illustrates this method's use.

.. code-block:: python

    # get the indices of the taxa having the maximum values for each trait
    out = bvmat.targmax()

Minimum breeding values for each trait
--------------------------------------

The minimum breeding values for each trait may be calculated using the ``tmax`` method. The code below illustrates this method's use.

.. code-block:: python

    # get the minimum breeding values for each trait
    out = bvmat.tmin()

Row (taxa) indices of the individuals with the smallest breeding values for each trait can be calculated using the ``targmax`` method. The code below illustrates this method's use.

.. code-block:: python

    # get the indices of the taxa having the minimum values for each trait
    out = bvmat.targmin()

Mean breeding values for each trait
-----------------------------------

The mean breeding value for each trait may be calculated using the ``tmean`` method. The code below illustrates the use of this method.

.. code-block:: python

    # get the mean breeding values for each trait
    out = bvmat.tmean()

Breeding value ranges for each trait
------------------------------------

The breeding value range for each trait may be calculated using the ``trange`` method. The code below illustrates the use of this method.

.. code-block:: python

    # get the breeding value ranges for each trait
    out = bvmat.trange()

Breeding value standard deviations for each trait
-------------------------------------------------

The breeding value standard deviation for each trait may be calculated using the ``trange`` method. The code below illustrates the use of this method.

.. code-block:: python

    # get the breeding value standard deviations for each trait
    out = bvmat.tstd()

Breeding value variances for each trait
---------------------------------------

The breeding value variance for each trait may be calculated using the ``trange`` method. The code below illustrates the use of this method.

.. code-block:: python

    # get the breeding value variances for each trait
    out = bvmat.tvar()

Unscaling and de-centering breeding values
------------------------------------------

A de-transformed (unscaled and de-centered) breeding value matrix may be calculated using the ``unscale`` method. The code below illustrates the use of this method.

.. code-block:: python

    # de-transform a breeding value matrix 
    out = bvmat.unscale()

Saving Breeding Value Matrices
==============================

Breeding value matrices may be exported to multiple formats including Pandas DataFrames, CSV files, and HDF5 files. The following subsections provide export examples.

Exporting to Pandas DataFrame
-----------------------------

The ``to_pandas`` method can be used to export a breeding value matrix to a Pandas DataFrame. Column names may be optionally provided to override default column names.

.. code-block:: python

    # export to a pandas.DataFrame
    # use default column names to export
    df = bvmat.to_pandas()


Exporting to CSV
----------------

The ``to_csv`` method can be used to export a breeding value matrix to a CSV file. Column names may be optionally provided to override default column names.

.. code-block:: python

    # export to a CSV
    # use default column names to export
    bvmat.to_csv("saved_breeding_values.csv")


Exporting to HDF5
-----------------

Most matrix object types in PyBrOpS allow for the export of matrices into an HDF5 format. To write breeding value matrices to an HDF5 file, use the ``to_hdf5`` method. The code below demonstrates the use of this method.

.. code-block:: python

    # write a breeding value matrix to an HDF5 file
    bvmat.to_hdf5("saved_breeding_values.h5")
