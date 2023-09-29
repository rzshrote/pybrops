Coancestry Matrices
###################

Class Family Overview
=====================

The ``CoancestryMatrix`` family of classes is used to represent coancestry relationships between individuals. This includes additive relationship matrices and genomic relationship matrices. Coancestry matrices can be used in the estimation of genomic prediction models and to make selection decisions. ``CoancestryMatrix`` objects store additional taxa metadata which serve as labels for rows and columns in the matrix.

Summary of Coancestry Matrix Classes
====================================

Coancestry matrix classes can be found in the ``pybrops.popgen.cmat`` module in PyBrOpS. Within this maodule are several interfaces and implemented classes, which are summarized in the table below.

.. list-table:: Summary of classes in the ``pybrops.popgen.cmat`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``CoancestryMatrix``
      - Abstract
      - Interface for all coancestry matrix child classes.
    * - ``DenseCoancestryMatrix``
      - Semi-Abstract
      - Semi-implemented class for deriving new dense coancestry matrix child classes.
    * - ``DenseMolecularCoancestryMatrix``
      - Concrete
      - Class representing dense molecular coancestry matrices.
    * - ``DenseVanRadenCoancestryMatrix``
      - Concrete
      - Class representing a genomic relationship matrix defined by VanRaden (2008).
    * - ``DenseYangCoancestryMatrix``
      - Concrete
      - Class representing a genomic relationship matrix defined by Yang.

Coancestry Matrix Properties
============================

Coancestry matrices have numerous properties which can be grouped into three main groupings: general properties, taxa properties, and square properties. These properties are summarized in the tables below.

General properties
------------------

Coancestry matrices share several shape properties common to the ``Matrix`` family of classes. These common properties are summarized below.

.. list-table:: Summary of ``CoancestryMatrix`` general properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``mat``
      - The raw coancestry matrix pointer
    * - ``mat_ndim``
      - The number of dimensions for the coancestry matrix
    * - ``mat_shape``
      - The coancestry matrix shape

Taxa properties
---------------

Coancestry matrices have several taxa related properties including taxa names, taxa group identities, and sorting metadata, which can be used for quick group access and sorting. These taxa related properties are summarized below.

.. list-table:: Summary of ``CoancestryMatrix`` taxa properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ntaxa``
      - The number of taxa represented by the coancestry matrix
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

Square matrix properties
------------------------

Since coancestry matrices are square by nature, they also have several properties which extract data regarding their squareness. These properties are summarized below.

.. list-table:: Summary of ``CoancestryMatrix`` square matrix properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nsquare``
      - The number of square axes for the coancestry matrix
    * - ``square_axes``
      - The axes indices for the square axes for the coancestry matrix
    * - ``square_axes_len``
      - The lengths of the square axes for the coancestry matrix

Loading Coancestry Matrix Modules
=================================

Importing coancestry matrix classes can be accomplished using the following import statements:

.. code-block:: python

    # import the CoancestryMatrix class (an abstract interface class)
    from pybrops.popgen.cmat.CoancestryMatrix import CoancestryMatrix

    # import the DenseCoancestryMatrix class (a semi-abstract class)
    from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix

    # import the DenseMolecularCoancestryMatrix class (a concrete implemented class)
    from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import DenseMolecularCoancestryMatrix

    # import the DenseVanRadenCoancestryMatrix class (a concrete implemented class)
    from pybrops.popgen.cmat.DenseVanRadenCoancestryMatrix import DenseVanRadenCoancestryMatrix

    # import the DenseYangCoancestryMatrix class (a concrete implemented class)
    from pybrops.popgen.cmat.DenseYangCoancestryMatrix import DenseYangCoancestryMatrix

Creating Coancestry Matrices
============================

Coancestry matrices can be created using several methods including from raw NumPy arrays, from ``GenotypeMatrix`` objects, from Pandas DataFrames, from CSV files, and from HDF5 fiels. The following subsections detail the creation or loading of coancestry matrices from their corresponding sources.

Creating coancestry matrices from NumPy arrays
----------------------------------------------

Using the constructor of a ``CoancestryMatrix`` class, one can create coancestry matrices from NumPy arrays. The example below demonstrates the creation of a ``DenseMolecularCoancestryMatrix`` object from raw NumPy arrays.

.. code-block:: python

    # shape parameters
    ntaxa = 100
    ngroup = 20

    # create random coancestries
    mat = numpy.random.uniform(0.0, 1.0, size = (ntaxa,ntaxa))

    # create taxa names
    taxa = numpy.array(
        ["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], 
        dtype = object
    )

    # create taxa groups
    taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
    taxa_grp.sort()

    # create a coancestry matrix from NumPy arrays
    cmat = DenseMolecularCoancestryMatrix(
        mat = mat,
        taxa = taxa,
        taxa_grp = taxa_grp
    )

Creating coancestry matrices from GenotypeMatrix objects
--------------------------------------------------------

Coancestry matrices may also be constructed from ``GenotypeMatrix`` objects. This can be accomplished using the ``from_gmat`` class emthod. The code below demonstrates how to use this method to accomplish this task.

.. code-block:: python

    # shape parameters for random genotypes
    ntaxa = 100
    nvrnt = 1000
    ngroup = 20
    nchrom = 10
    ploidy = 2

    # create random genotypes
    mat = numpy.random.randint(0, ploidy+1, size = (ntaxa,nvrnt)).astype("int8")

    # create taxa names
    taxa = numpy.array(
        ["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], 
        dtype = object
    )

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
    vrnt_name = numpy.array(
        ["SNP"+str(i+1).zfill(4) for i in range(nvrnt)],
        dtype = object
    )

    # create a genotype matrix from scratch using NumPy arrays
    gmat = DenseGenotypeMatrix(
        mat = mat,
        taxa = taxa,
        taxa_grp = taxa_grp, 
        vrnt_chrgrp = vrnt_chrgrp,
        vrnt_phypos = vrnt_phypos, 
        vrnt_name = vrnt_name, 
        vrnt_genpos = None,
        vrnt_xoprob = None, 
        vrnt_hapgrp = None, 
        vrnt_hapalt = None,
        vrnt_hapref = None, 
        vrnt_mask = None,
        ploidy = ploidy
    )

    # group taxa and variants
    gmat.group_taxa()
    gmat.group_vrnt()

    # construct Coancestry Matrix from a Genotype Matrix
    cmat = DenseMolecularCoancestryMatrix.from_gmat(gmat = gmat)

Creating coancestry matrices from Pandas DataFrames
---------------------------------------------------

Coancestry matrices may be read from Pandas DataFrames. The ``from_pandas`` class method may be used to read a ``CoancestryMatrix`` from a pandas DataFrame. The code example below demonstrates this method's usage.

.. code-block:: python

    # load from pandas.DataFrame
    tmp = DenseMolecularCoancestryMatrix.from_pandas(
        df = df,
        taxa_col = "taxa",          # column from which to load taxa
        taxa_grp_col = "taxa_grp",  # column from which to load taxa groups
        taxa = "all",               # load all taxa
    )

Loading coancestry matrices from CSV files
------------------------------------------

Coancestry matrices may also be read from CSV files in a manner similar to Pandas DataFrames. The ``from_csv`` class method may be used to load coancestry matrices from csv files. The following code block demonstrates the usage of this method.

.. code-block:: python

    # load from pandas.DataFrame
    tmp = DenseMolecularCoancestryMatrix.from_csv(
        filename = "saved_coancestry_matrix.csv",
        taxa_col = "taxa",          # column from which to load taxa
        taxa_grp_col = "taxa_grp",  # column from which to load taxa groups
        taxa = "all",               # load all taxa
    )

Loading coancestry matrices from HDF5 files
-------------------------------------------

As with all classes in the ``Matrix`` family, ``CoancestryMatrix`` objects may be imported and exported to an HDF5 format. To read saved coancestry matrices from an HDF5 file, use the ``from_hdf5`` class method. The code below demonstrates the use of this method.

.. code-block:: python

    # read from file
    cmat = DenseMolecularCoancestryMatrix.from_hdf5("sample_coancestry_matrix.h5")

Copying Coancestry Matrices
===========================

Coancestry matrices may be copied using two methods: shallow copying and deep copying.

Shallow copying
---------------

.. |link_copy_copy| replace:: ``copy.copy``
.. _link_copy_copy: https://docs.python.org/3/library/copy.html#copy.copy

In shallow copying, references to a ``CoancestryMatrix``'s data are copied to a new coancestry matrix object. Copying is only one level deep which means that changes to the original object may affect data values in the copied object. The code below illustrates the use of the ``copy`` method bound to ``CoancestryMatrix`` objects and the base Python function |link_copy_copy|_ which can both be used to shallow copy a coancestry matrix object.

.. code-block:: python

    # copy a coancestry matrix
    tmp = copy.copy(cmat)
    tmp = cmat.copy()

Deep copying
------------

.. |link_copy_deepcopy| replace:: ``copy.deepcopy``
.. _link_copy_deepcopy: https://docs.python.org/3/library/copy.html#copy.deepcopy

In deep copying, data in a ``CoancestryMatrix`` is recursively copied to a new coancestry matrix object. Copying occurs down to the deepest levels so that changes to the original object will not affect data values in the copied object. The code below illustrates the use of the ``deepcopy`` method bound to ``CoancestryMatrix`` objects and the base Python function |link_copy_deepcopy|_ which can both be used to deep copy a coancestry matrix object.

.. code-block:: python

    # deep copy a coancestry matrix
    tmp = copy.deepcopy(cmat)
    tmp = cmat.deepcopy()

Copy-On Element Manipulation
============================

Coancestry matrices have several methods by which modifed copies of the original matrix can be made. These are called copy-on element manipulation methods. Matrices may have rows and/or columns adjoined, deleted, inserted, or selected. The following sections demonstrate the use of these method families.

Adjoin elements
---------------

The ``adjoin`` family of methods allows for taxa rows and columns of a coancestry matrix to be adjoined together, creating a new matrix in the process. Use of the ``adjoin`` method family is demonstrated in the code below.

.. code-block:: python

    # create a new coancestry matrix to demonstrate
    new = cmat.deepcopy()

    # adjoin coancestry matrices along the taxa axis
    tmp = cmat.adjoin(new, axis = cmat.taxa_axis)
    tmp = cmat.adjoin_taxa(new)

Delete elements
---------------

The ``delete`` family of methods allows for taxa rows and columns of a coancestry matrix to be removed in a copy of the original. Use of the ``delete`` method family is demonstrated in the code below.

.. code-block:: python

    # delete first taxon using an integer
    tmp = cmat.delete(0, axis = cmat.taxa_axis)
    tmp = cmat.delete_taxa(0)

    # delete first five taxa using a slice
    tmp = cmat.delete(slice(0,5), axis = cmat.taxa_axis)
    tmp = cmat.delete_taxa(slice(0,5))

    # delete first five taxa using a Sequence
    tmp = cmat.delete([0,1,2,3,4], axis = cmat.taxa_axis)
    tmp = cmat.delete_taxa([0,1,2,3,4])

Insert elements
---------------

The ``insert`` family of methods allows for taxa rows and columns of a coancestry matrix to be inserted into a copy of the original matrix. Use of the ``insert`` method family is demonstrated in the code below.

Select elements
---------------

The ``select`` family of methods allows for taxa rows and columns of the coancestry matrix to be selected and extracted to a copy of the original matrix. Use of the ``select`` method family is demonstrated in the code below.

.. code-block:: python

    # select first five taxa using a Sequence
    tmp = cmat.select([0,1,2,3,4], axis = cmat.taxa_axis)
    tmp = cmat.select_taxa([0,1,2,3,4])

In-Place Element Manipulation
=============================

Coancestry matrices have several methods which execute in-place element manipulations. These are called in-place element manipulation methods. Coancestry matrices may have taxa rows and/or columns appended, removed, incorporated, or concatenated. The following sections demonstrate the use of these method families.

Append elements
---------------

The ``append`` family of methods allows for new taxa rows and columns to be appended to the coancestry matrix. The code segment below demonstrates their use. 

.. code-block:: python

    # append coancestry matrices along the taxa axis
    tmp = cmat.deepcopy()                   # copy original
    tmp.append(cmat, axis = tmp.taxa_axis)  # append original to copy

    tmp = cmat.deepcopy()                   # copy original
    tmp.append_taxa(cmat)                   # append original to copy

Remove elements
---------------

The ``remove`` family of methods allows for taxa rows and columns to be removed from a coancestry matrix. A demonstration of their use can be seen below. 

.. code-block:: python

    # remove first taxon using an integer
    tmp = cmat.deepcopy()                           # copy original
    tmp.remove(0, axis = cmat.taxa_axis)            # remove from copy

    tmp = cmat.deepcopy()                           # copy original
    tmp.remove_taxa(0)                              # remove from copy

    # remove first five taxa using a slice
    tmp = cmat.deepcopy()                           # copy original
    tmp.remove(slice(0,5), axis = cmat.taxa_axis)   # remove from copy

    tmp = cmat.deepcopy()                           # copy original
    tmp.remove_taxa(slice(0,5))                     # remove from copy

    # remove first five taxa using a Sequence
    tmp = cmat.deepcopy()                           # copy original
    tmp.remove([0,1,2,3,4], axis = cmat.taxa_axis)  # remove from copy

    tmp = cmat.deepcopy()                           # copy original
    tmp.remove_taxa([0,1,2,3,4])                    # remove from copy

Incorporate elements
--------------------

The ``incorp`` family of methods allows for new taxa rows and columns to be inserted at specific locations a coancestry matrix. Use of the ``incorp`` family is demonstrated in the code segment below below. 

.. code-block:: python

    # incorp coancestry matrix along the taxa axis before index 0
    tmp = cmat.deepcopy()                           # copy original
    tmp.incorp(0, cmat, axis = cmat.taxa_axis)      # incorporate into copy

    tmp = cmat.deepcopy()                           # copy original
    tmp.incorp_taxa(0, cmat)                        # incorporate into copy

Concatenate elements
--------------------

The ``concat`` family of methods allows for multiple coancestry matrices to be concatenated to each other. The code segment below demonstrates their use. 

Grouping and Sorting
====================

Coancestry matrices in PyBrOpS have several sorting and grouping focused methods. Sorting methods can be used to reorder, sort, and group taxa alphanumerically. The following sections demonstrate the use of the ``reorder``, ``lexsort``, ``sort``, and ``group`` method families.

Reordering elements
-------------------

Taxa in a coancestry matrix can be reordered using the ``reorder`` family of methods. Demonstrations of this method family are below.

.. code-block:: python

    # create reordering indices
    indices = numpy.arange(cmat.ntaxa)
    numpy.random.shuffle(indices)
    tmp = cmat.deepcopy()

    # reorder values along the taxa axis
    tmp.reorder(indices, axis = tmp.taxa_axis)
    tmp.reorder_taxa(indices)

Lexsorting elements
-------------------

An indirect stable sort - or lexsort - for taxa axes can be performed using the ``lexsort`` family of methods. The code segment below illustrates the use of this family of methods.

.. code-block:: python

    # create lexsort keys for taxa
    key1 = numpy.random.randint(0, 10, cmat.ntaxa)
    key2 = numpy.arange(cmat.ntaxa)
    numpy.random.shuffle(key2)

    # lexsort along the taxa axis
    cmat.lexsort((key2,key1), axis = cmat.taxa_axis)
    cmat.lexsort_taxa((key2,key1))

Sorting elements
----------------

Alphanumeric sorting along taxa axes can be done using the ``sort`` family of methods. Sorting examples are illustrated below.

.. code-block:: python

    # make copy
    tmp = cmat.deepcopy()

    # sort along taxa axis
    tmp.sort(axis = tmp.taxa_axis)
    tmp.sort_taxa()

Grouping elements
-----------------

Grouping along taxa axes can be done using the ``group`` family of methods. The following code illustrates the use of the ``group`` method family along the taxa axes of a coancestry matrix.

.. code-block:: python

    # make copy
    tmp = cmat.deepcopy()

    # sort along taxa axis
    tmp.group(axis = tmp.taxa_axis)
    tmp.group_taxa()

    # determine whether grouping has occurred along the taxa axis
    out = tmp.is_grouped(axis = tmp.taxa_axis)
    out = tmp.is_grouped_taxa()

Coancestry and Kinship Methods
==============================

Retrieving coancestry values
----------------------------

Coancestry values may be retrieved by using the ``coancestry`` method. Retrieval of coancestry values may also be done via indexing, but this method is **not** guaranteed to be in the correct format (coancestry or kinship). The return format via indexing is implementation dependent. The code below demonstrates the use of this method.

.. code-block:: python

    # Get the coancestry at a specific matrix coordinate
    out = cmat.coancestry(0,0)
    out = cmat[0,0] # NOT guaranteed to be in correct format

Retrieving kinship values
-------------------------

Kinship values may be retrieved by using the ``kinship`` method. Like coancestry values, kinship values may also be retrieved via indexing, but this method is not guaranteed to be in the correct format since the internal matrix representation is implementation dependent. The code below demonstrates the use of this method.

.. code-block:: python

    # Get the kinship at a specific matrix coordinate
    out = cmat.kinship(0,0)
    out = 0.5 * cmat[0,0] # NOT guaranteed to be in correct format

Retrieving the coancestry matrix as a specific format
-----------------------------------------------------

Coancestry matrices may be extracted as bare-bones NumPy arrays in kinship or coancestry formats using the ``mat_asformat`` method. The code below demonstrates the usage of this method.

.. code-block:: python

    # Get the coancestry matrix as a specific format
    out = cmat.mat_asformat(format = "kinship")
    out = cmat.mat_asformat(format = "coancestry")

Determining if the coancestry matrix is positive semidefinite
-------------------------------------------------------------

For some optimization, it may be necessary for a coancestry matrix to be positive semidefinite. The ``is_positive_semidefinite`` method may be used to determine if a coancestry matrix is positive semidefinite or not. An example of this method's usage is below.

.. code-block:: python

    # Determine if the coancestry matrix is positive semidefinite (convex)
    out = cmat.is_positive_semidefinite()

Applying jitter values along the diagonal
-----------------------------------------

In the event that a coancestry matrix is not positive semidefinite, it may be helpful to apply a small jitter along the diagonal of the matrix. A jitter can be applied using the ``apply_jitter`` method as is demonstrated below.

.. code-block:: python

    # Apply a jitter along the diagonal to try to make the matrix positive semidefinite
    out = cmat.apply_jitter()

Calculating the matrix inverse
------------------------------

The inverse of a coancestry matrix may be calculated using the ``inverse`` method. Varying format arguments may be used to specific if the inverse is for the kinship representation or the coancestry representation. The code below demonstrates this method's usage.

.. code-block:: python

    # Calculate the inverse of the coancestry matrix
    out = cmat.inverse()
    out = cmat.inverse(format = "kinship")
    out = cmat.inverse(format = "coancestry")

Calculating maximum attainable inbreeding
-----------------------------------------

For particular tasks, it may be useful to calculate the maximum attainable level of inbreeding after one generation. This is equivalent to the maximum value along the diagonal of the coancestry matrix. This may be done using the ``max_inbreeding`` method, demonstrated below.

.. code-block:: python

    # Calculate the maximum attainable inbreeding after 1 generation
    out = cmat.max_inbreeding()
    out = cmat.max_inbreeding(format = "kinship")
    out = cmat.max_inbreeding(format = "coancestry")

Calculating minimum attainable inbreeding
-----------------------------------------

For other tasks, it may be useful to calculate the minimum attainable level of inbreeding after one generation. This may be done using the ``min_inbreeding`` method, demonstrated below.

.. code-block:: python

    # Calculate the minimum attainable inbreeding after 1 generation
    out = cmat.min_inbreeding()
    out = cmat.min_inbreeding(format = "kinship")
    out = cmat.min_inbreeding(format = "coancestry")

Summary Statistics
==================

Maximum coancestry
------------------

The maximum coancestry value across the entire coancestry matrix may be calculated using the ``max`` method. Below is a demonstration of this method.

.. code-block:: python

    # get the max for the whole coancestry matrix
    out = cmat.max()

Mean coancestry
---------------

The mean coancestry across the entire coancestry matrix may be calculated using the ``mean`` method. Below is a demonstration of this method.

.. code-block:: python

    # get the mean for the whole coancestry matrix
    out = cmat.mean()

Minimum coancestry
------------------

The minimum coancestry value across the entire coancestry matrix may be calculated using the ``min`` method. Below is a demonstration of this method.

.. code-block:: python

    # get the min for the whole coancestry matrix
    out = cmat.min()

Exporting Coancestry Matrices
=============================

Coancestry matrices may be exported to multiple formats including Pandas DataFrames, CSV files, and HDF5 files. The following subsections demonstrate how to export coancestry matrices.

Exporting to Pandas DataFrame
-----------------------------

The ``to_pandas`` method can be used to export a coancestry matrix to a Pandas DataFrame. Column names may be optionally provided to override default column names.

.. code-block:: python

    # export to a pandas.DataFrame
    # use default column names to export
    df = cmat.to_pandas()

Exporting to CSV
----------------

The ``to_csv`` method can be used to export a coancestry matrix to a CSV file. Column names may be optionally provided to override default column names.

.. code-block:: python

    # export to a CSV
    # use default column names to export
    cmat.to_csv("saved_coancestry_matrix.csv")

Exporting to HDF5
-----------------

To write coancestry matrices to an HDF5 file, use the ``to_hdf5`` method. The code below demonstrates the use of this method.

.. code-block:: python

    # write a coancestry matrix to an HDF5 file
    cmat.to_hdf5("saved_coancestry_matrix.h5")
