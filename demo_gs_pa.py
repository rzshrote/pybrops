# import 3rd party libraries
import numpy
import time
import pandas

# import our libraries
from pybropt import gs
from pybropt import stpfn
from pybropt import objfn

################################################################################
# define several constants
seed = 112119
n_indiv = 50
n_loci = 200
n_trials = 2

# number of lines to select
sel_size = 10

# generate binary marker data
numpy.random.seed(seed)
geno = numpy.empty((2,n_indiv,n_loci), dtype=numpy.uint8)
geno[0,:,:] = numpy.random.binomial(1, 0.1, (1,n_indiv,n_loci))
geno[1,:,:] = geno[0,:,:]
coeff = numpy.random.normal(0, 1, n_loci)
pos_set = numpy.arange(n_indiv, dtype='uint32')
pos_state = numpy.tile(pos_set, sel_size) # select 10
dimsizes = numpy.repeat(n_indiv, sel_size)

# make wcoeff
wcoeff = numpy.absolute(coeff)

# make target frequency
tfreq = numpy.where(coeff >= 0, 1.0, 0.0)

# make target ld
tld = 1.0

# make genetic maps
gmap = numpy.random.uniform(0.0, 1.5, n_loci)
for b in range(0, n_loci, n_loci//10):
    gmap[b:b+(n_loci//10)] = numpy.sort(gmap[b:b+(n_loci//10)])
lgroup = numpy.repeat(n_loci//10, 10)

# number of random breedings
cycles = 10

################################################################################

# calculate results
results, history = gs.pa(
    sel_size = sel_size,
    geno = geno,
    coeff = coeff,
    tfreq = tfreq,
    tld = tld,
    gmap = gmap,
    lgroup = lgroup,
    cycles = cycles,
    algorithm = "hc_sa_set",
    mtype = "tril",
    verbose = True
)

print(results)
print(history)
