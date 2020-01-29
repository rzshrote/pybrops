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
from pybropt import gs

################################################################################
# define several constants
seed = 117119
n_indiv = 100
n_loci = 500
n_trials = 2

# number of lines to select
sel_size = 10

# generate binary marker data
numpy.random.seed(seed)
geno = numpy.empty((2,n_indiv,n_loci), dtype=numpy.uint8)
geno[0,:,:] = numpy.random.binomial(1, 0.1, (1,n_indiv,n_loci))
geno[1,:,:] = geno[0,:,:]

coeff = numpy.random.normal(0, 1, n_loci)      # make marker effects
wcoeff = numpy.absolute(coeff)                 # make wcoeff
tfreq = numpy.where(coeff >= 0, 1.0, 0.0)      # make target frequency
tld = 1.0                                      # make target ld

# make genetic maps
gmap = numpy.random.uniform(0.0, 1.5, n_loci)
for b in range(0, n_loci, n_loci//10):
    gmap[b:b+(n_loci//10)] = numpy.sort(gmap[b:b+(n_loci//10)])
lgroup = numpy.repeat(n_loci//10, 10)

# number of random breeding generations after crossing
generations = 5

# dimension coefficients
dcoeff = numpy.array([1.0,1.0,1.0]) # even weight on all objectives

################################################################################

# calculate results
results_df, history_df = gs.pa_v2(
    sel_size = sel_size,
    geno = geno,
    coeff = coeff,
    tfreq = tfreq,
    tld = tld,
    gmap = gmap,
    lgroup = lgroup,
    generations = generations,
    dcoeff = dcoeff,
    algorithm = 'hc_sa_set',
    mapfn = 'haldane',
    ldfn = 'r_sq',
    mtype = 'tril',
    mem = 1000,
    verbose = True,
    verbose_algo = True
)

print("Results:")
print(results_df)
print("History:")
print(history_df)
