# import 3rd party libraries
import numpy
import time
import pandas
import os, sys

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.realpath(__file__)
    )))
)

# import our libraries
from pybropt import gs
from pybropt import stpfn

################################################################################
# define several constants
seed = 112119
n_indiv = 200
n_loci = 10000
n_sel = 10
n_gen = 5
n_trials = 120

# generate binary marker data
numpy.random.seed(seed)
markers = numpy.empty((2,n_indiv,n_loci), dtype=numpy.uint8)
markers[0,:,:] = numpy.random.binomial(1, 0.1, (1,n_indiv,n_loci))
markers[1,:,:] = markers[0,:,:]
effects = numpy.random.normal(0, 1, n_loci)
pos_set = numpy.arange(n_indiv, dtype='int64')
pos_state = numpy.tile(pos_set, n_sel) # select 10
dim_sizes = numpy.repeat(n_indiv, n_sel)

# make genetic maps
gmap = numpy.random.uniform(0.0, 1.5, n_loci)
for b in range(0, n_loci, n_loci//10):
    gmap[b:b+(n_loci//10)] = numpy.sort(gmap[b:b+(n_loci//10)])
gmap_size = numpy.repeat(n_loci//10, 10)
################################################################################

# calculate results
icpso_history = gs.opv(
    geno = markers,
    coeff = effects,
    sel_size = n_sel,
    algorithm = "icpso",
    algorithm_varg = {
        "n": 1000,
        "inertia_wt": 0.33,
        "accel_coeff_pbest": 0.33,
        "accel_coeff_gbest": 0.34,
        "scale_factor": 0.1,
        "stpfn": lambda x: x<30,
        "states": pos_state,
        "dim_sizes": dim_sizes
    },
    seed = None,
    nthreads = 1,
    zwidth = 3,
    verbose = True
)

icpso_history.history_to_df().to_csv(
    "icpso_history.tsv",
    sep = '\t',
    index = False
)

hc_sa_set_history = gs.opv(
    geno = markers,
    coeff = effects,
    sel_size = n_sel,
    algorithm = "hc_sa_set",
    # algorithm_varg = {
    #     'n': 100,
    #     'states': numpy.tile(numpy.arange(n_indiv), (n_sel,)),
    #     'dim_sizes': numpy.repeat(n_indiv, n_sel),
    #     'inertia_wt': 0.33,
    #     'accel_coeff_pbest': 0.33,
    #     'accel_coeff_gbest': 0.34,
    #     'scale_factor': 0.1,
    #     'stpfn': lambda x: x < 20,
    # },
    seed = None,
    nthreads = 1,
    zwidth = 3,
    verbose = True
)

hc_sa_set_history.history_to_df().to_csv(
    "hc_sa_set_history.tsv",
    sep = "\t",
    index = False
)
