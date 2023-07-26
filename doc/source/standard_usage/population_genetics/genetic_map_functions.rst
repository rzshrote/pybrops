Genetic Map Functions
#####################

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

Calculating Map Distances
=========================

Genetic map distances for arbitrary recombination probabilities can be calculated by using the ``invmapfn`` method, the inverse of the ``mapfn`` function. An example of its usage is below:

.. code-block:: python

    # invert the recombination probability calculated previously
    d_haldane = haldane.invmapfn(r_haldane)

    # invert the recombination probability calculated previously
    d_kosambi = kosambi.invmapfn(r_kosambi)
