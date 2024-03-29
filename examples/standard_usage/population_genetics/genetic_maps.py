#!/usr/bin/env python3

import copy
import numpy

# import the GeneticMap class (an abstract interface class)
from pybrops.popgen.gmap.GeneticMap import GeneticMap

# import the ExtendedGeneticMap class (a concrete class)
from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap

###
### Creating Genetic Maps
### =====================

##
## Creating genetic maps from NumPy arrays
## ---------------------------------------

# define number of variants and number of chromosomes
nvrnt = 100
nchrom = 10

# create random chromosome groups 1-10
vrnt_chrgrp = numpy.random.choice(nchrom, nvrnt, True)+1
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
    vrnt_chrgrp = vrnt_chrgrp,
    vrnt_phypos = vrnt_phypos,
    vrnt_stop   = vrnt_phypos,
    vrnt_genpos = vrnt_genpos,
    vrnt_name   = vrnt_name,
    vrnt_fncode = None, # not needed
)

##
## Reading genetic maps from a file
## --------------------------------

# read genetic map from file
# for the purpose of this example, do not automatically group markers 
# and build an interpolation spline after reading genetic map data.
gmap = ExtendedGeneticMap.from_egmap(
    "McMullen_2009_US_NAM.egmap",
    auto_group = False,
    auto_build_spline = False
)

###
### Genetic Map Properties
### ======================

##
## Marker variant properties
## -------------------------
tmp = gmap.nvrnt            # Number of variants in the Genetic Map
tmp = gmap.vrnt_chrgrp      # Marker variant chromosome group labels
tmp = gmap.vrnt_phypos      # Marker variant chromosome physical positions
tmp = gmap.vrnt_genpos      # Marker variant chromosome genetic positions
tmp = gmap.vrnt_name        # Marker variant names
tmp = gmap.vrnt_chrgrp_name # Names of chromosome groups
tmp = gmap.vrnt_chrgrp_stix # Chromosome group start indices
tmp = gmap.vrnt_chrgrp_spix # Chromosome group stop indices
tmp = gmap.vrnt_chrgrp_len  # Number of marker variants on each chromosome group

##
## Spline properties
## -----------------
tmp = gmap.spline               # Interpolation splines
tmp = gmap.spline_kind          # Interpolation spline type
tmp = gmap.spline_fill_value    # Interpolation spline default fill value

###
### Copying Genetic Maps
### ====================

##
## Shallow copying
## ---------------

# copy the genetic map
tmp = copy.copy(gmap)
tmp = gmap.copy()

##
## Deep copying
## ------------

# deep copy the genetic map
tmp = copy.deepcopy(gmap)
tmp = gmap.deepcopy()

###
### Sorting and Grouping Genetic Maps
### =================================

##
## Reordering map elements
## -----------------------

# create reordering indices
indices = numpy.arange(gmap.nvrnt)
numpy.random.shuffle(indices)
tmp = gmap.deepcopy()

# reorder values
tmp.reorder(indices)

##
## Lexsorting map elements
## -----------------------

# create lexsort keys
key1 = numpy.random.randint(0, 10, gmap.nvrnt)
key2 = numpy.random.choice(gmap.nvrnt, gmap.nvrnt, False)

# lexsort using keys
out = gmap.lexsort((key2,key1))

##
## Sorting map elements
## --------------------

# sort the genetic map
gmap.sort()

##
## Grouping map elements
## ---------------------

# group markers based on their chromosome/linkage group
gmap.group()

# determine whether a GeneticMap is grouped using the ``is_grouped`` method
value = gmap.is_grouped()
print("GeneticMap is grouped:", value)

### 
### Genetic Map Congruency
### ======================

## 
## Checking for congruency
## -----------------------

# elementwise test of marker congruence
value = gmap.congruence()
print("Congruence of each marker:", value)

# whole genetic map congruence test
value = gmap.is_congruent()
print("Whole GeneticMap is congruent:", value)

## 
## Removing map discrepancies
## --------------------------

# automatically remove discrepancies
# it may be better to manually remove these
gmap.remove_discrepancies()

###
### Building Interpolation Splines
### ==============================

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
chrgrp = numpy.array([1, 1, 1, 1, 1], dtype = int)

# construct physical position array
phypos = numpy.array([18203210,19293034,20110347,20474722,21398386], dtype = int)

# interpolate new gentic map positions
genpos = gmap.interp_genpos(
    vrnt_chrgrp = chrgrp,
    vrnt_phypos = phypos
)

# print results
print("New Variant Linkage Group:", chrgrp)
print("New Variant Physical Positions:", phypos)
print("Interpolated Genetic Positions:", genpos.round(3))

# export gmap as a CSV
gmap.to_csv("McMullen_2009_US_NAM.csv")

