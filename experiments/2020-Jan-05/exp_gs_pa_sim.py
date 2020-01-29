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
from pybropt import gs
from pybropt import stpfn
from pybropt import objfn

################################################################################
# define several constants
seed = 112119
n_indiv = 50
n_loci = 200
n_trials = 20

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

# calculate initial condition GEBVs
gebvs = numpy.dot(geno, coeff).sum(0).tolist()
pop_opv_score = objfn.opv(slice(None), geno * coeff)

# calculate initial conditions
initial_df = pandas.DataFrame(
    [["method", "algorithm", seed, 0, pop_opv_score] + gebvs],
    columns = (["method", "algorithm", "seed", "bcycle", "opv_score"] +
        ["gebv"+str(i).zfill(3) for i in range(n_indiv)])
)

# write initial_df to file
initial_df.to_csv("initial.csv")


################################################################################

# declare dataframes for storing things
# pa_population_df = None
# pa_selections_df = None
# opv_population_df = None
# opv_selections_df = None
cgs_population_df = None
cgs_selections_df = None

# simulate breeding cycles
for trial in range(n_trials):
    # simulate
    # pa_pop_df, pa_sel_df = gs.pa_sim(
    #     pop_size = n_indiv,
    #     sel_size = sel_size,
    #     geno = geno,
    #     coeff = coeff,
    #     tfreq = tfreq,
    #     tld = tld,
    #     gmap = gmap,
    #     lgroup = lgroup,
    #     cycles = cycles,
    #     algorithm = "hc_sa_set",
    #     mtype = "tril",
    #     verbose = False
    # )
    # opv_pop_df, opv_sel_df = gs.opv_sim(
    #     geno = geno,
    #     coeff = coeff,
    #     pop_size = n_indiv,
    #     sel_size = sel_size,
    #     cycles = cycles,
    #     d = gmap,
    #     lgroup_size = lgroup,
    #     algorithm = "hc_sa_set",
    #     algorithm_varg = None,
    #     interference = None,
    #     seed = None,
    #     nthreads = 1,
    #     zwidth = 3,
    #     verbose = False
    # )
    cgs_pop_df, cgs_sel_df = gs.cgs_sim(
        pop_size = n_indiv,
        sel_size = sel_size,
        geno = geno,
        coeff = coeff,
        gmap = gmap,
        lgroup = lgroup,
        bcycles = cycles
    )
    # concatenate
    # pa_population_df = pandas.concat([pa_population_df, pa_pop_df])
    # pa_selections_df = pandas.concat([pa_selections_df, pa_sel_df])
    # opv_population_df = pandas.concat([opv_population_df, opv_pop_df])
    # opv_selections_df = pandas.concat([opv_selections_df, opv_sel_df])
    cgs_population_df = pandas.concat([cgs_population_df, cgs_pop_df])
    cgs_selections_df = pandas.concat([cgs_selections_df, cgs_sel_df])
    # print cycle trial number
    print("Trial:", trial+1)

# write to files
# pa_population_df.to_csv("pa_population.csv")
# pa_selections_df.to_csv("pa_selections.csv")
# opv_population_df.to_csv("opv_population.csv")
# opv_selections_df.to_csv("opv_selections.csv")
cgs_population_df.to_csv("cgs_population.csv")
cgs_selections_df.to_csv("cgs_selections.csv")
