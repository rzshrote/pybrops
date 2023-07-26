#!/usr/bin/env python3

import numpy

# import the GeneticMapFunction class (an abstract interface class)
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction

# import the HaldaneMapFunction class (a concrete class)
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction

# import the KosambiMapFunction class (a concrete class)
from pybrops.popgen.gmap.KosambiMapFunction import KosambiMapFunction

###
### Genetic Map Function Construction
###

# create a Haldane map function object
haldane = HaldaneMapFunction()

# create a Kosambi map function object
kosambi = KosambiMapFunction()

###
### Recombination Probability Calculation
###

# create an array of map distances in Morgans
d = numpy.array([1.2, 0.5, 0.1, 0.8, 1.4])
print("Original map distances:", d)

# calculate recombination probability using the Haldane map function
r_haldane = haldane.mapfn(d)
print("Haldane recombination probabilities:", r_haldane.round(2))

# calculate recombination probability using the Kosambi map function
r_kosambi = kosambi.mapfn(d)
print("Kosambi recombination probabilities:", r_kosambi.round(2))

###
### Genetic map distance
###

# invert the recombination probability calculated previously
d_haldane = haldane.invmapfn(r_haldane)
print("Haldane map distances:", d_haldane.round(2))

# invert the recombination probability calculated previously
d_kosambi = kosambi.invmapfn(r_kosambi)
print("Kosambi map distance:", d_kosambi.round(2))
