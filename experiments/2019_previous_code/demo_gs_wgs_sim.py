# import 3rd party libraries
import numpy
import time
import pandas

# import our libraries
from pybropt import gs
from pybropt import stpfn

################################################################################
# define several constants
seed = 112119
n_indiv = 200
n_loci = 10000
n_sel = 10
n_gen = 20
n_trials = 2

# generate binary marker data
numpy.random.seed(seed)
markers = numpy.empty((2,n_indiv,n_loci), dtype=numpy.uint8)
markers[0,:,:] = numpy.random.binomial(1, 0.1, (1,n_indiv,n_loci))
markers[1,:,:] = markers[0,:,:]
effects = numpy.random.normal(0, 1, n_loci)
pos_set = numpy.arange(n_indiv, dtype='uint32')
pos_state = numpy.tile(pos_set, n_sel) # select 10
dimsizes = numpy.repeat(n_indiv, n_sel)

# make genetic maps
d = numpy.random.uniform(0.0, 1.5, n_loci)
for b in range(0, n_loci, n_loci//10):
    d[b:b+(n_loci//10)] = numpy.sort(d[b:b+(n_loci//10)])
d_size = numpy.repeat(n_loci//10, 10)
################################################################################

# calculate results
population_df, selections_df = gs.wgs_sim(
    geno = markers,
    coeff = effects,
    pop_size = n_indiv,
    sel_size = n_sel,
    cycles = n_gen,
    d = d,
    lgroup_size = d_size,
    wgtfn = None,
    algorithm = None,
    interference = None,
    seed = None,
    nthreads = 1,
    zwidth = 3,
    verbose = True
)

print(population_df)
print(selections_df)
