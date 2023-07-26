Genetic Maps
############

Loading ``GeneticMap`` modules
==============================

Genetic map support for PyBrOpS is found in the ``pybrops.popgen.gmap`` submodule. Example imports from this submodule are:

.. code-block:: python

    # import the GeneticMap class (an abstract interface class)
    from pybrops.popgen.gmap.GeneticMap import GeneticMap

    # import the ExtendedGeneticMap class (a concrete class)
    from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap


Reading Genetic Maps from a File
================================

``GeneticMap`` objects can be read from a file.

The example below illustrates the reading of an ``ExtendedGeneticMap`` from a file.

.. code-block:: python

    # read genetic map from file
    # for the purpose of this example, do not automatically group markers 
    # and build an interpolation spline after reading genetic map data.
    gmap = ExtendedGeneticMap.from_egmap(
        "McMullen_2009_US_NAM.egmap",
        auto_group = False,
        auto_build_spline = False
    )

Sorting and Grouping Genetic Maps
=================================

On import from a file, genetic map data needs to be sorted and grouped into marker linkage groups. Often, an import method for a ``GeneticMap`` object will automatically sort and group data, but occasionally it may be necessary to manually sort and group marker data. This can be accomplished using the ``group`` method:

.. code-block:: python

    # group markers based on their chromosome/linkage group
    gmap.group()

To test whether a ``GeneticMap``'s data have been sorted and grouped into linkage groups, the ``is_grouped`` method can be used:

.. code-block:: python

    # determine whether a GeneticMap is grouped using the ``is_grouped`` method
    value = gmap.is_grouped()

Checking for the Congruence of a Genetic Map
============================================

Sometimes physical positions and genetic map positions are in disagreement as to their orderings. This may be caused by errors made in genome assemblies and/or genetic map assemblies. Unfortunately, these disagreements cause issues for interpolation spline construction and need to be removed or corrected. An elementwise physical position-genetic position congruence test can be conducted using the ``congruence`` method:

.. code-block:: python

    # elementwise test of marker congruence
    value = gmap.congruence()

If one desires to test whether all of the markers in the ``GeneticMap`` are congruent, one can use the ``is_congruent`` method:

.. code-block:: python

    # whole genetic map congruence test
    value = gmap.is_congruent()

Building Interpolation Splines
==============================

Before using a ``GeneticMap`` to interpolate genetic position data, an interpolation spline must be constructed. Often, an import method for a ``GeneticMap`` object will automatically construct a spline from the provided data. Occasionally, it may be necessary to manually construct an interpolation spline. The ``build_spline`` method can be used to construct an interpolation spline:

.. code-block:: python

    # construct a linear spline to interpolate genetic map positions
    gmap.build_spline()

To test whether a ``GeneticMap`` has an interpolation spline, the ``has_spline`` method can be used:

.. code-block:: python

    # determine whether a GeneticMap has an interpolation spline using the 
    # ``has_spline`` method
    value = gmap.has_spline()

Interpolating Genetic Positions
===============================

Interpolating genetic map positions from physical positions can be done using the ``interp_genpos`` method:

.. code-block:: python

    ### create new positions to interpolate
    # construct linkage group array: everything is on chromosome 1
    new_vrnt_chrgrp = numpy.array(
        [1, 1, 1, 1, 1], 
        dtype = int
    )
    # construct physical position array
    new_vrnt_phypos = numpy.array(
        [18209321, 19296303, 20115034, 20475472, 21396838], 
        dtype = int
    )

    # interpolate new gentic map positions
    new_vrnt_genpos = gmap.interp_genpos(
        vrnt_chrgrp = new_vrnt_chrgrp,
        vrnt_phypos = new_vrnt_phypos
    )
