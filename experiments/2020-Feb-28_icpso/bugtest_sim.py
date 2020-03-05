# import 3rd party libraries
import numpy
import time
import pandas
import sys, os

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
)

# import our libraries
from pybropt import gs
from pybropt import stpfn
from pybropt import objfn

################################################################################
# define several constants
seed = 112119
n_indiv = 200
n_loci = 10000
n_sel = 10
n_gen = 5
n_trials = 1

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
gmap = numpy.random.uniform(0.0, 1.5, n_loci)
for b in range(0, n_loci, n_loci//10):
    gmap[b:b+(n_loci//10)] = numpy.sort(gmap[b:b+(n_loci//10)])
gmap_size = numpy.repeat(n_loci//10, 10)
################################################################################

for i in range(n_trials):
    # calculate results
    simout = gs.opv_sim(
        geno = markers,
        coeff = effects,
        pop_size = n_indiv,
        sel_size = n_sel,
        cycles = n_gen,
        gmap = gmap,
        lgroup = gmap_size,
        algorithm = "hc_sa_set",
        # algorithm = "icpso",
        # algorithm_varg = {
        #     "n": 1000,
        #     "inertia_wt": 0.33,
        #     "accel_coeff_pbest": 0.33,
        #     "accel_coeff_gbest": 0.34,
        #     "scale_factor": 0.1,
        #     "stpfn": lambda x: x<30,
        #     "states": numpy.tile(numpy.arange(n_indiv), (n_sel,)),
        #     "dim_sizes": numpy.repeat(n_indiv, n_sel)
        # },
        interference = None,
        seed = None,
        nthreads = 1,
        zwidth = 3,
        verbose = True
    )

    population_df, selections_df = simout.history_to_df()

    # write population to file
    pop_fname = "data/" + "opv_icpso_pop" + str(i+1).zfill(4) + ".tsv"
    population_df.to_csv(
        pop_fname,
        sep = '\t',
        index = False
    )

    # write selections to file
    sel_fname = "data/" + "opv_icpso_sel" + str(i+1).zfill(4) + ".tsv"
    selections_df.to_csv(
        sel_fname,
        sep = '\t',
        index = False
    )

    # print progress
    print("Trial:", i+1)
