Genetic Map Functions
#####################

Class Family Overview
=====================

The ``GeneticMapFunction`` family of classes allows for the representation of genetic map functions within PyBrOpS. Genetic map functions are used to convert genetic map positions into recombination probabilities, which are useful in the simulation of genetic recombination and the estimation of expected progeny variance. Popular genetic map functions like the Haldane and Kosambi map functions are provided by PyBrOpS out-of-box.

Summary of Genetic Map Function Classes
=======================================

Genetic map function support in PyBrOpS is found in the ``pybrops.popgen.gmap`` module. Contained in this submodule are several ``GeneticMapFunction`` class type definitions. These classes are summarized in the table below.

.. list-table:: Summary of genetic map function classes in the ``pybrops.popgen.gmap`` module
    :widths: 25 15 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``GeneticMapFunction``
      - Abstract
      - Interface for all genetic map function child classes.
    * - ``HaldaneMapFunction``
      - Concrete
      - Class representing the Haldane genetic map function.
    * - ``KosambiMapFunction``
      - Concrete
      - Class representing the Kosambi genetic map function.

Loading Genetic Map Function Modules
====================================

Genetic map function support in PyBrOpS is found in the ``pybrops.popgen.gmap`` submodule. The ``GeneticMapFunction`` abstract class is the basal interface for all genetic map functions. PyBrOpS supports both the Haldane and Kosambi map functions, two of the most popular genetic map functions. Classes can be imported as so:

.. code-block:: python

    # import the GeneticMapFunction class (an abstract interface class)
    from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction

    # import the HaldaneMapFunction class (a concrete class)
    from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction

    # import the KosambiMapFunction class (a concrete class)
    from pybrops.popgen.gmap.KosambiMapFunction import KosambiMapFunction

Constructing Genetic Map Function Objects
=========================================

Haldane and Kosambi map function objects can be created as seen in the example below. Object construction requires no arguments.

.. code-block:: python

    # create a Haldane map function object
    haldane = HaldaneMapFunction()

    # create a Kosambi map function object
    kosambi = KosambiMapFunction()

Calculating Recombination Probabilities
=======================================

Recombination probabilities for arbitrary genetic map distances can be calculated using the ``mapfn`` method. An example of its usage is below:

.. code-block:: python

    # create an array of map distances in Morgans
    d = numpy.array([1.2, 0.5, 0.1, 0.8, 1.4])

    # calculate recombination probability using the Haldane map function
    r_haldane = haldane.mapfn(d)

    # calculate recombination probability using the Kosambi map function
    r_kosambi = kosambi.mapfn(d)

Calculating Genetic Map Distances
=================================

Genetic map distances for arbitrary recombination probabilities can be calculated by using the ``invmapfn`` method, the inverse of the ``mapfn`` function. An example of its usage is below:

.. code-block:: python

    # invert the recombination probability calculated previously
    d_haldane = haldane.invmapfn(r_haldane)

    # invert the recombination probability calculated previously
    d_kosambi = kosambi.invmapfn(r_kosambi)

Calculating Sequential and Pairwise Recombination Probabilities
===============================================================

In some instances, it may be useful to calculate sequential and/or pairwise recombination probabilities to simulate recombination. Genetic map function objects have four methods which allow for the calculation of these probabilities: ``rprob1g``, ``rprob1p``, ``rprob2g``, and ``rprob2p``.

Calculating sequential recombination probabilities
--------------------------------------------------

The ``rprob1g`` and ``rprob1p`` methods are used to calculate sequential recombination probabilities between successive markers, provided chromosome group assignments, and genetic positions or physical positions, respectively. Chromosome boundaries always start with probabilities of 0.5.

The code segment below demonstrates the use of the ``rprob1g`` and ``rprob1p`` methods.

.. code-block:: python

    # calculate from genetic positions
    xoprob = haldane.rprob1g(
        gmap = gmap,
        vrnt_chrgrp = gmap.vrnt_chrgrp,
        vrnt_genpos = gmap.vrnt_genpos
    )

    # calculate from physical positions
    xoprob = haldane.rprob1p(
        gmap = gmap,
        vrnt_chrgrp = gmap.vrnt_chrgrp,
        vrnt_phypos = gmap.vrnt_phypos
    )

Calculating pairwise recombination probabilities
------------------------------------------------

The ``rprob2g`` and ``rprob2p`` methods are used to calculate pairwise recombination probabilities between markers, provided chromosome group assignments, and genetic positions or physical positions, respectively. Inter-chromosome pairings are always assigned probabilities of 0.5.

The code segment below demonstrates the use of the ``rprob2g`` and ``rprob2p`` methods.

.. code-block:: python

    # calculate from genetic positions
    xoprob = haldane.rprob2g(
        gmap = gmap,
        vrnt_chrgrp = gmap.vrnt_chrgrp,
        vrnt_genpos = gmap.vrnt_genpos
    )

    # calculate from physical positions
    xoprob = haldane.rprob2p(
        gmap = gmap,
        vrnt_chrgrp = gmap.vrnt_chrgrp,
        vrnt_phypos = gmap.vrnt_phypos
    )
