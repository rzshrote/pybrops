Coancestry Matrices
###################

Class Family Overview
=====================

The ``CoancestryMatrix`` family of classes is used to represent coancestry relationships between individuals. ``CoancestryMatrix`` objects can be used in the estimation of genomic prediction models and to make selection decisions.

Loading Coancestry Matrix Modules
=================================

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


Loading coancestry matrices from HDF5 files
-------------------------------------------

Coancestry Matrix Properties
============================

Coancestry matrix general properties
------------------------------------

Breeding value matrix taxa properties
-------------------------------------

Breeding value matrix trait properties
--------------------------------------

Copying Breeding Value Matrices
===============================

Shallow copying
---------------

Deep copying
------------

Coancestry Matrix Element Copy-On-Manipulation
==============================================

Adjoin elements
---------------

Delete elements
---------------

Insert elements
---------------

Select elements
---------------

Coancestry Matrix Element In-Place-Manipulation
===============================================

Append elements
---------------

Remove elements
---------------

Incorporate elements
--------------------

Concatenate elements
--------------------

