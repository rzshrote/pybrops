#!/usr/bin/env python3

import numpy

# import the GeneticMap class (an abstract interface class)
from pybrops.popgen.gmap.GeneticMap import GeneticMap

# import the ExtendedGeneticMap class (a concrete class)
from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap

###
### Reading/Import/Object Construction
###

# read genetic map from file
# for the purpose of this example, do not automatically group markers 
# and build an interpolation spline after reading genetic map data.
gmap = ExtendedGeneticMap.from_egmap(
    "McMullen_2009_US_NAM.egmap",
    auto_group = False,
    auto_build_spline = False
)

###
### Grouping
###

# group markers based on their chromosome/linkage group
gmap.group()

# determine whether a GeneticMap is grouped using the ``is_grouped`` method
value = gmap.is_grouped()
print("GeneticMap is grouped:", value)

###
### Congruency testing
###

# elementwise test of marker congruence
value = gmap.congruence()
print("Congruence of each marker:", value)

# whole genetic map congruence test
value = gmap.is_congruent()
print("Whole GeneticMap is congruent:", value)

###
### Building Spline
###

# construct a linear spline to interpolate genetic map positions
gmap.build_spline()

# determine whether a GeneticMap has an interpolation spline using the 
# ``has_spline`` method
value = gmap.has_spline()
print("GeneticMap has interpolation spline:", value)

###
### Interpolating New Positions
###

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

# print results
print("New Variant Linkage Group:", new_vrnt_chrgrp)
print("New Variant Physical Positions:", new_vrnt_phypos)
print("Interpolated Genetic Positions:", new_vrnt_genpos.round(3))

# export gmap as a CSV
gmap.to_csv("McMullen_2009_US_NAM.csv")

