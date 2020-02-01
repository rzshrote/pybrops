# import 3rd party libraries
import numpy
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
seed = 1312020   # RNG seed for generating genotype
pop_size = 100  # population size
sel_size = 20   # number of lines to select (can have replicates)
n_loci = 300    # number of markers
bcycles = 10    # number of breeding cycles
n_trials = 5   # number of breeding cycle trajectories to simulate

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

# output variables
mogm_population_df = None
mogm_selections_df = None

for cross_variance in numpy.linspace(0.0,0.002,11):
    # define breeding priorities
    diversity = (1.0 - cross_variance) / 2.0
    allele = (1.0 - cross_variance) / 2.0
    dcoeff = numpy.array([diversity, allele, cross_variance])

    # simulate breeding cycles
    for trial in range(n_trials):
        # optimize
        mogm_pop_df, mogm_sel_df = gm.mogm_sim(
            pop_size = pop_size,
            sel_size = sel_size, # divide this by 2 to get number of crosses
            geno = geno,
            coeff = coeff,
            tfreq = tfreq,
            gmap = gmap,
            lgroup = lgroup,
            dcoeff_t = numpy.tile(dcoeff, (bcycles,1)),
            bcycles = bcycles,
            algorithm = "hc_sa_state"
        )

        # concatenate dataframes
        mogm_population_df = pandas.concat([mogm_population_df ,mogm_pop_df])
        mogm_selections_df = pandas.concat([mogm_selections_df, mogm_sel_df])

        # print cycle trial number
        print("stdA:", cross_variance, "Trial:", trial+1)

# write to files
mogm_population_df.to_csv("mogm_population_variable3.csv", index=False)
mogm_selections_df.to_csv("mogm_selections_variable3.csv", index=False)
