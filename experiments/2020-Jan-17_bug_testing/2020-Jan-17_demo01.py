# import 3rd party libraries
import numpy
import time
import pandas

# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.realpath(__file__)
    )))
)
# import our libraries
from pybropt import objfn

################################################################################
# define several constants
seed = 117119
n_indiv = 200
n_loci = 10000
n_trials = 2

# number of lines to select
sel_size = 10

# generate binary marker data
numpy.random.seed(seed)
geno = numpy.empty((2,n_indiv,n_loci), dtype=numpy.uint8)
geno[0,:,:] = numpy.random.binomial(1, 0.1, (1,n_indiv,n_loci))
geno[1,:,:] = geno[0,:,:]

effect = numpy.random.normal(0, 1, n_loci)      # make marker effects
wcoeff = numpy.absolute(effect)                 # make wcoeff
tfreq = numpy.where(effect >= 0, 1.0, 0.0)      # make target frequency
tld = 1.0                                       # make target ld

# make genetic maps
gmap = numpy.random.uniform(0.0, 1.5, n_loci)
for b in range(0, n_loci, n_loci//10):
    gmap[b:b+(n_loci//10)] = numpy.sort(gmap[b:b+(n_loci//10)])
lgroup = numpy.repeat(n_loci//10, 10)

# number of random breeding generations
generations = 10

# dimension coefficients
dcoeff = numpy.array([1.0,1.0,1.0]) # even weight on all objectives

################################################################################

# calculate results
out = objfn.pa_v2_max(
    rsel = numpy.array([1,3,5,7,9,11,13,15,17,19,21]),
    geno = geno,
    wcoeff = wcoeff,
    tfreq = tfreq,
    tld = tld,
    gmap = gmap,
    lgroup = lgroup,
    generations = generations,
    dcoeff = dcoeff,
    mapfn = 'haldane',
    ldfn = 'r_sq',
    mtype = 'tril',
    mem = 1000
)

print(out)
