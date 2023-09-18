#!/usr/bin/env python3

import numpy
from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap

### 
### Loading Genetic Map Function Modules
### ====================================

# import the GeneticMapFunction class (an abstract interface class)
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction

# import the HaldaneMapFunction class (a concrete class)
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction

# import the KosambiMapFunction class (a concrete class)
from pybrops.popgen.gmap.KosambiMapFunction import KosambiMapFunction

###
### Constructing Genetic Map Function Objects
### =========================================

# create a Haldane map function object
haldane = HaldaneMapFunction()

# create a Kosambi map function object
kosambi = KosambiMapFunction()

###
### Calculating Recombination Probabilities
### =======================================

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
### Calculating Genetic Map Distances
### =================================

# invert the recombination probability calculated previously
d_haldane = haldane.invmapfn(r_haldane)
print("Haldane map distances:", d_haldane.round(2))

# invert the recombination probability calculated previously
d_kosambi = kosambi.invmapfn(r_kosambi)
print("Kosambi map distance:", d_kosambi.round(2))

###
### Calculating Sequential and Pairwise Recombination Probabilities
### ===============================================================

# define number of variants and number of chromosomes
nvrnt = 100
nchrom = 10

# create random chromosome groups 1-10
vrnt_chrgrp = numpy.random.choice(nchrom, nvrnt, True)+1

# create variant physical positions in range [1, 2**28]
vrnt_phypos = numpy.random.randint(1, 2**20, nvrnt)

# create variant genetic positions in range [0,1]
vrnt_genpos = (1.0 / 2.0**20) * vrnt_phypos

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
## Calculating sequential recombination probabilities
## --------------------------------------------------

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

##
## Calculating pairwise recombination probabilities
## ------------------------------------------------

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
