Coancestry Matrices
###################

Class Family Overview
=====================

The ``CoancestryMatrix`` family of classes is used to represent coancestry relationships between individuals. ``CoancestryMatrix`` objects can be used in the estimation of genomic prediction models and to make selection decisions.

Summary of Coancestry Matrix Classes
====================================

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

Creating coancestry matrices from NumPy arrays
----------------------------------------------

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

Loading coancestry matrices from HDF5 files
-------------------------------------------

.. code-block:: python

    # read from file
    cmat = DenseMolecularCoancestryMatrix.from_hdf5("sample_coancestry_matrix.h5")

Coancestry Matrix Properties
============================

General properties
------------------

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
    * - ``location``
      - The location of the coancestry matrix if it has been transformed
    * - ``scale``
      - The scale of the coancestry matrix if it has been transformed

Taxa-related properties
-----------------------

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


Square matrix-related properties
--------------------------------

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


Copying Coancestry Matrices
===========================

Shallow copying
---------------

.. code-block:: python

    # copy a coancestry matrix
    tmp = copy.copy(cmat)
    tmp = cmat.copy()

Deep copying
------------

.. code-block:: python

    # deep copy a coancestry matrix
    tmp = copy.deepcopy(cmat)
    tmp = cmat.deepcopy()

Copy-On Element Manipulation
============================

Adjoin elements
---------------

.. code-block:: python

    # create a new coancestry matrix to demonstrate
    new = cmat.deepcopy()

    # adjoin coancestry matrices along the taxa axis
    tmp = cmat.adjoin(new, axis = cmat.taxa_axis)
    tmp = cmat.adjoin_taxa(new)

Delete elements
---------------

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

Select elements
---------------

.. code-block:: python

    # select first five taxa using a Sequence
    tmp = cmat.select([0,1,2,3,4], axis = cmat.taxa_axis)
    tmp = cmat.select_taxa([0,1,2,3,4])

In-Place Element Manipulation
=============================

Append elements
---------------

.. code-block:: python

    # append coancestry matrices along the taxa axis
    tmp = cmat.deepcopy()                   # copy original
    tmp.append(cmat, axis = tmp.taxa_axis)  # append original to copy

    tmp = cmat.deepcopy()                   # copy original
    tmp.append_taxa(cmat)                   # append original to copy

Remove elements
---------------

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

.. code-block:: python

    # incorp coancestry matrix along the taxa axis before index 0
    tmp = cmat.deepcopy()                           # copy original
    tmp.incorp(0, cmat, axis = cmat.taxa_axis)      # incorporate into copy

    tmp = cmat.deepcopy()                           # copy original
    tmp.incorp_taxa(0, cmat)                        # incorporate into copy

Concatenate elements
--------------------

Grouping and Sorting
====================

Reordering elements
-------------------

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

.. code-block:: python

    # make copy
    tmp = cmat.deepcopy()

    # sort along taxa axis
    tmp.sort(axis = tmp.taxa_axis)
    tmp.sort_taxa()

Grouping elements
-----------------

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

.. code-block:: python

    # Get the coancestry at a specific matrix coordinate
    out = cmat.coancestry(0,0)

Retrieving kinship values
-------------------------

.. code-block:: python

    # Get the kinship at a specific matrix coordinate
    out = cmat.kinship(0,0)

Retrieving the coancestry matrix as a specific format
-----------------------------------------------------

.. code-block:: python

    # Get the coancestry matrix as a specific format
    cmat.mat_asformat(format = "kinship")

Determining if the coancestry matrix is positive semidefinite
-------------------------------------------------------------

.. code-block:: python

    # Determine if the coancestry matrix is positive semidefinite (convex)
    out = cmat.is_positive_semidefinite()

Applying jitter values along the diagonal
-----------------------------------------

.. code-block:: python

    # Apply a jitter along the diagonal to try to make the matrix positive semidefinite
    out = cmat.apply_jitter()

Calculating the matrix inverse
------------------------------

.. code-block:: python

    # Calculate the inverse of the coancestry matrix
    out = cmat.inverse()
    out = cmat.inverse(format = "kinship")

Calculating maximum attainable inbreeding
-----------------------------------------

.. code-block:: python

    # Calculate the maximum attainable inbreeding after 1 generation
    out = cmat.max_inbreeding()
    out = cmat.min_inbreeding(format = "kinship")

Calculating minimum attainable inbreeding
-----------------------------------------

.. code-block:: python

    # Calculate the minimum attainable inbreeding after 1 generation
    out = cmat.min_inbreeding()
    out = cmat.min_inbreeding(format = "kinship")

Summary Statistics
==================

Maximum coancestry
------------------

.. code-block:: python

    # get the max for the whole coancestry matrix
    out = cmat.max()

Mean coancestry
---------------

.. code-block:: python

    # get the mean for the whole coancestry matrix
    out = cmat.mean()

Minimum coancestry
------------------

.. code-block:: python

    # get the min for the whole coancestry matrix
    out = cmat.min()

Saving Coancestry Matrices
==========================

Exporting to HDF5
-----------------

.. code-block:: python

    # write a coancestry matrix to an HDF5 file
    cmat.to_hdf5("saved_coancestry_matrix.h5")
