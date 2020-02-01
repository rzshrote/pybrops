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
from pybropt import objfn, gs, gm

################################################################################
# define several constants
seed = 118119   # RNG seed for generating genotype
pop_size = 100  # population size
sel_size = 10   # number of lines to select
n_loci = 300    # number of markers
bcycles = 10    # number of breeding cycles
n_trials = 1   # number of breeding cycle trajectories to simulate

# generate binary marker data
numpy.random.seed(seed)
geno = numpy.empty((2,pop_size,n_loci), dtype=numpy.uint8)
geno[0,:,:] = numpy.random.binomial(1, 0.1, (1,pop_size,n_loci))
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

# dimension coefficients; even weighting on all three objectives
dcoeff_t = numpy.repeat(1.0, 3*bcycles).reshape(bcycles, 3)

# calculate initial condition GEBVs
gebvs = numpy.dot(geno, coeff).sum(0).tolist()
pop_opv_score = objfn.opv(slice(None), geno * coeff)

# calculate initial conditions
initial_df = pandas.DataFrame(
    [["method", "algorithm", seed, 0, pop_opv_score] + gebvs],
    columns = (["method", "algorithm", "seed", "bcycle", "opv_score"] +
        ["gebv"+str(i).zfill(3) for i in range(pop_size)])
)

# write initial_df to file
initial_df.to_csv("initial.csv")

################################################################################
# TODO: clean up opv_sim internals
# declare dataframes for storing things
# pa_v2_population_df = None
# pa_v2_selections_df = None
# opv_population_df = None
# opv_selections_df = None
# cgs_population_df = None
# cgs_selections_df = None
# mogs_population_df = None
# mogs_selections_df = None
mogm_population_df = None
mogm_selections_df = None

# simulate breeding cycles
for trial in range(n_trials):
    # pa_v2_pop_df, pa_v2_sel_df = gs.pa_v2_sim(
    #     pop_size = pop_size,
    #     sel_size = sel_size,
    #     geno = geno,
    #     coeff = coeff,
    #     tfreq = tfreq,
    #     tld = tld,
    #     gmap = gmap,
    #     lgroup = lgroup,
    #     dcoeff_t = dcoeff_t,
    #     bcycles = bcycles,
    #     algorithm = 'hc_sa_set',
    #     mapfn = 'haldane',
    #     ldfn = 'r_sq',
    #     mtype = 'tril',
    #     mem = 1000,
    #     verbose = True,
    #     verbose_algo = False
    # )
    # opv_pop_df, opv_sel_df = gs.opv_sim(
    #     geno = geno,
    #     coeff = coeff,
    #     pop_size = pop_size,
    #     sel_size = sel_size,
    #     cycles = bcycles,
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
    # cgs_pop_df, cgs_sel_df = gs.cgs_sim(
    #     pop_size = pop_size,
    #     sel_size = sel_size,
    #     geno = geno,
    #     coeff = coeff,
    #     gmap = gmap,
    #     lgroup = lgroup,
    #     bcycles = bcycles
    # )
    # mogs_pop_df, mogs_sel_df = gs.mogs_sim(
    #     pop_size = pop_size,
    #     sel_size = sel_size,
    #     geno = geno,
    #     coeff = coeff,
    #     tfreq = tfreq,
    #     gmap = gmap,
    #     lgroup = lgroup,
    #     dcoeff_t = numpy.repeat(1.0, bcycles*2).reshape(bcycles, 2),
    #     bcycles = bcycles,
    #     algorithm = "hc_sa_set"
    # )
    mogm_pop_df, mogm_sel_df = gm.mogm_sim(
        pop_size = pop_size,
        sel_size = sel_size *2,
        geno = geno,
        coeff = coeff,
        tfreq = tfreq,
        gmap = gmap,
        lgroup = lgroup,
        dcoeff_t = numpy.tile([1.0,1.0,0.0], (bcycles,1)),
        bcycles = bcycles,
        algorithm = "hc_sa_state"
    )
    # concatenate dataframes
    # pa_v2_population_df = pandas.concat([pa_v2_population_df, pa_v2_pop_df])
    # pa_v2_selections_df = pandas.concat([pa_v2_selections_df, pa_v2_sel_df])
    # opv_population_df = pandas.concat([opv_population_df, opv_pop_df])
    # opv_selections_df = pandas.concat([opv_selections_df, opv_sel_df])
    # cgs_population_df = pandas.concat([cgs_population_df, cgs_pop_df])
    # cgs_selections_df = pandas.concat([cgs_selections_df, cgs_sel_df])
    # mogs_population_df = pandas.concat([mogs_population_df ,mogs_pop_df])
    # mogs_selections_df = pandas.concat([mogs_selections_df, mogs_sel_df])
    mogm_population_df = pandas.concat([mogm_population_df ,mogm_pop_df])
    mogm_selections_df = pandas.concat([mogm_selections_df, mogm_sel_df])
    # print cycle trial number
    print("Trial:", trial+1)

# write to files
# pa_v2_population_df.to_csv("pa_v2_population.csv")
# pa_v2_selections_df.to_csv("pa_v2_selections.csv")
# opv_population_df.to_csv("opv_population.csv")
# opv_selections_df.to_csv("opv_selections.csv")
# cgs_population_df.to_csv("cgs_population.csv")
# cgs_selections_df.to_csv("cgs_selections.csv")
# mogs_population_df.to_csv("mogs_population.csv")
# mogs_selections_df.to_csv("mogs_selections.csv")
mogm_population_df.to_csv("mogm_population_x2state_0.0.csv")
mogm_selections_df.to_csv("mogm_selections_x2state_0.0.csv")
