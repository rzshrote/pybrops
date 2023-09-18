Genetic Maps
############

Class Family Overview
=====================

Genetic maps are used to estimate marker-wise genetic map positions in mating simulations. PyBrOpS provides functionality for the reading and writing of genetic maps and the interpolation of genetic map positions via ``GeneticMap`` classes. The ability to read a genetic map from a file allows for real genetic recombination data to be used in a breeding simulation, adding realism.

Summary of Genetic Map Classes
==============================

Genetic map support for PyBrOpS is found in the ``pybrops.popgen.gmap`` module. Contained in this submodule are several ``GeneticMap`` class type definitions. These classes are summarized in the table below.

.. list-table:: Summary of classes in the ``pybrops.popgen.gmap`` module
    :widths: 25 15 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``GeneticMap``
      - Abstract
      - Interface for all genetic map child classes.
    * - ``ExtendedGeneticMap``
      - Concrete
      - Class representing genetic maps with additional genetic map metadata.

Loading Genetic Map Modules
===========================

The various genetic map classes summarized above can be imported into a Python scope as demonstrated in the code example below.

.. code-block:: python

    # import the GeneticMap class (an abstract interface class)
    from pybrops.popgen.gmap.GeneticMap import GeneticMap

    # import the ExtendedGeneticMap class (a concrete class)
    from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap

Creating Genetic Maps
=====================

Creating genetic maps from NumPy arrays
---------------------------------------

Genetic map objects can be constructed from raw NumPy arrays, which may be useful for creating simulated genetic maps. The code below demonstrates the construction of an ``ExtendedGeneticMap`` from NumPy arrays.

.. code-block:: python

    # define number of variants
    nvrnt = 100

    # create random chromosome groups 1-9
    chroms = list("123456789")
    vrnt_chrgrp = numpy.random.choice(chroms, nvrnt, True).astype(object)
    vrnt_chrgrp.sort()

    # create variant physical positions in range [1, 2**28]
    vrnt_phypos = numpy.random.randint(1, 2**20, nvrnt)
    vrnt_phypos.sort()

    # create variant genetic positions in range [0,1]
    vrnt_genpos = numpy.random.random(nvrnt)
    vrnt_genpos.sort()

    # create variant names
    vrnt_name = numpy.array(["SNP"+str(i+1).zfill(3) for i in range(nvrnt)], dtype=object)

    # construct genetic map
    gmap = ExtendedGeneticMap(
        vrnt_chrgrp=vrnt_chrgrp,
        vrnt_phypos=vrnt_phypos,
        vrnt_stop=vrnt_phypos,
        vrnt_genpos=vrnt_genpos,
        vrnt_name=vrnt_name,
        vrnt_fncode=None, # not needed
    )

Reading genetic maps from a file
--------------------------------

``GeneticMap`` objects can also be read from a file. This may be useful for users who have their own, empirically determined genetic maps from an organism of their choosing. The example below illustrates the reading of an ``ExtendedGeneticMap`` from a file. Specifically, this is the maize genetic map constructed from the US NAM population as published by McMullen et al. (2009).

.. code-block:: python

    # read genetic map from file
    # for the purpose of this example, do not automatically group markers 
    # and build an interpolation spline after reading genetic map data.
    gmap = ExtendedGeneticMap.from_egmap(
        "McMullen_2009_US_NAM.egmap",
        auto_group = False,
        auto_build_spline = False
    )

Genetic Map Properties
======================

``GeneticMap`` objects have a set of properties shared by all genetic maps. These properties can be grouped into two categories: marker variant properties and spline properties. The former set of properties contain information about the marker set constituting the genetic map, while the latter set of properties contain spline model information necessary for the interpolation of genetic map positions. Marker variant and spline properties are summarized in the tables below.

Marker variant properties
-------------------------

.. list-table:: Summary of ``GeneticMap`` marker variant properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nvrnt``           
      - Number of variants in the Genetic Map
    * - ``vrnt_chrgrp``     
      - Marker variant chromosome group labels
    * - ``vrnt_phypos``     
      - Marker variant chromosome physical positions
    * - ``vrnt_genpos``     
      - Marker variant chromosome genetic positions
    * - ``vrnt_name``       
      - Marker variant names
    * - ``vrnt_chrgrp_name``
      - Names of chromosome groups
    * - ``vrnt_chrgrp_stix``
      - Chromosome group start indices
    * - ``vrnt_chrgrp_spix``
      - Chromosome group stop indices
    * - ``vrnt_chrgrp_len`` 
      - Number of marker variants on each chromosome group

Spline properties
-----------------

.. list-table:: Summary of ``GeneticMap`` spline properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``spline``           
      - Interpolation splines
    * - ``spline_kind``      
      - Interpolation spline type
    * - ``spline_fill_value``
      - Interpolation spline default fill value

Copying Genetic Maps
====================

At times, it may be necessary to copy a genetic map. There are two methods of copying genetic maps: shallow copying and deep copying.

Shallow copying
---------------

In shallow copying, references to a ``GeneticMap``'s variant and spline data are copied to a new genetic map object. Copying is only one level deep and changes to the data in the original object may affect data values in the copied object.

.. code-block:: python

    # copy the genetic map
    tmp = copy.copy(gmap)
    tmp = gmap.copy()

Deep copying
------------

In deep copying, data in a ``GeneticMap``'s variant and spline data are recursively copyied. Copying occurs down to the deepest level making it so that changes to the data in the original object will not affect data values in the copied object.

.. code-block:: python

    # deep copy the genetic map
    tmp = copy.deepcopy(gmap)
    tmp = gmap.deepcopy()

Sorting and Grouping Genetic Maps
=================================

Reordering map elements
-----------------------

In some instances, it may be useful to manually reorder genetic map elements. This may be accomplished by providing an array of reordering indices to the ``reorder`` method.

.. code-block:: python

    # create reordering indices
    indices = numpy.arange(gmap.nvrnt)
    numpy.random.shuffle(indices)
    tmp = gmap.deepcopy()

    # reorder values
    tmp.reorder(indices)

Lexsorting map elements
-----------------------

An indirect stable sort may be performed using the ``lexsort`` method. If the ``lexsort`` method is not provided a set of ``keys``, it defaults to utilizing marker variant chromosome group assignments, marker variant chromosome physical positions, and marker variant chromosome genetic positions in that order of priority. 

.. code-block:: python

    # create lexsort keys
    key1 = numpy.random.randint(0, 10, gmap.nvrnt)
    key2 = numpy.random.choice(gmap.nvrnt, gmap.nvrnt, False)

    # lexsort using keys
    out = gmap.lexsort((key2,key1))

Sorting map elements
--------------------

In-place sorting of marker variants in a ``GeneticMap`` object can be accomplished using the ``sort`` method. The ``sort`` method optionally accepts a set of ``keys`` which can be used to sort the marker variants in the genetic map. If a set of ``keys`` is not provided, the keys are the same defaults as those in the ``lexsort`` method.

.. code-block:: python

    # sort the genetic map
    gmap.sort()

Grouping map elements
---------------------

On import from a file, genetic map data needs to be sorted and grouped into marker linkage groups so that an interpolation spline can be built. Often, the constructor or an import method for a ``GeneticMap`` object will automatically sort and group data, but occasionally it may be necessary to manually sort and group marker data. This can be accomplished using the ``group`` method:

.. code-block:: python

    # group markers based on their chromosome/linkage group
    gmap.group()

To test whether a ``GeneticMap``'s data have been sorted and grouped into linkage groups, the ``is_grouped`` method can be used:

.. code-block:: python

    # determine whether a GeneticMap is grouped using the ``is_grouped`` method
    value = gmap.is_grouped()

Genetic Map Congruency
======================

Checking for congruency
-----------------------

Sometimes physical positions and genetic map positions are in disagreement as to their orderings. This may be caused by errors made in genome assemblies and/or genetic map assemblies. Unfortunately, these disagreements cause issues for interpolation spline construction and need to be removed or corrected. An elementwise physical position-genetic position congruence test can be conducted using the ``congruence`` method:

.. code-block:: python

    # elementwise test of marker congruence
    value = gmap.congruence()

If one desires to test whether all of the markers in the ``GeneticMap`` are congruent, one can use the ``is_congruent`` method:

.. code-block:: python

    # whole genetic map congruence test
    value = gmap.is_congruent()

Removing map discrepancies
--------------------------

Loci where the physical positions and genetic positions are not in agreement may be automatically removed using the ``remove_discrepancies``. Be mindful that a manual inspection and correction of a genetic map's discrepancies may be superior to this automatic method.

The code below demonstrates automatic discrepancy removal:

.. code-block:: python

    # automatically remove discrepancies
    gmap.remove_discrepancies()


Building Interpolation Splines
==============================

Before using a ``GeneticMap`` to interpolate genetic position data, an interpolation spline must be constructed. Often, the constructor or an import method for a ``GeneticMap`` object will automatically construct a spline from the provided data. Occasionally, it may be necessary to manually construct an interpolation spline. The ``build_spline`` method can be used to construct an interpolation spline:

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
    chrgrp = numpy.array([1, 1, 1, 1, 1], dtype = int)

    # construct physical position array
    phypos = numpy.array([18203210,19293034,20110347,20474722,21398386], dtype = int)

    # interpolate new gentic map positions
    genpos = gmap.interp_genpos(
        vrnt_chrgrp = chrgrp,
        vrnt_phypos = phypos
    )
