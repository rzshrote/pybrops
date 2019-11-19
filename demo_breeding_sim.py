import numpy
import time

# import optimizing functions
from pybropt.algo.hc import hc_sa_set, hc_sa_state
from pybropt.algo.icpso import icpso
from pybropt.algo.ga import ga_fps
# import meiosis simulating engine
from pybropt.mate.meiosis import meiosis
# import scoring function
from pybropt.objfn.opv import opv
from pybropt.objfn.cgs import cgs

################################################################################

# define several constants
seed = 111819
n_indiv = 200
n_loci = 10000
n_sel = 10
n_gen = 10
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

# define a stopping function for the PSO algorithm
def stop(x):
    return x < 50

# test the efficacy of each method for long term genetic gain
method1_pop = list()    # hc_sa_set population statistics
method1_sel = list()    # hc_sa_set selected statistics
for trial in range(1, n_trials+1):
    # initialize the population
    phases = markers.copy()
    varg = {"hcoeff": phases * effects}

    # append population stats to list
    method1_pop.append([
        "hc_sa_set",
        "pop",
        trial,
        0,
        cgs(phases,effects)
    ])

    # begin simulating breeding cycles
    for gen in range(1, n_gen+1):
        # make selections from founding generation
        output = hc_sa_set(opv, varg, 10, pos_set, verbose=False)

        # append selected individual statistics
        method1_sel.append([
            "hc_sa_set",
            "sel",
            trial,
            gen,
            cgs(phases[:,output["X_gbest_pos"],:],effects)
        ])

        # make copies of parents
        females = output["X_gbest_pos"].copy()
        males = output["X_gbest_pos"].copy()

        # shuffle parent matings (simulates equal contribution random matings)
        numpy.random.shuffle(females)
        numpy.random.shuffle(males)

        # make hybrid genotypes
        hybrids = numpy.empty((2, n_sel, n_loci), dtype=phases.dtype)
        hybrids[0,:,:] = phases[0,females,:]
        hybrids[1,:,:] = phases[1,males,:]

        # make gametes (in our case, make 200/10=20 gamets from one hybrid)
        sources = numpy.repeat(numpy.arange(n_sel), n_indiv//n_sel)
        gout = meiosis(hybrids, d, d_size, sources, verbose=False)

        # double these gamets to make DH; replace prev. data in 'phases' array
        phases[0,:,:] = gout
        phases[1,:,:] = gout

        # append population stats to list
        method1_pop.append([
            "hc_sa_set",
            "pop",
            trial,
            gen,
            cgs(phases,effects)
        ])

        # recalculate haplotype coefficients
        varg = {"hcoeff": phases * effects}

        # print progress
        print("Trial:", trial, "\tGen:", gen)

# write population stats
with open("demo_breed_pop.tsv", "w") as pout:
    for item in method1_pop:
        pout.write(
            "%s\t%s\t%s\t%s" %
            (item[0], item[1], item[2], item[3])
        )
        for i in item[4]:
            pout.write("\t%s" % i)
        pout.write("\n")

# write selection stats
with open("demo_breed_sel.tsv", "w") as sout:
    for item in method1_sel:
        sout.write(
            "%s\t%s\t%s\t%s" %
            (item[0], item[1], item[2], item[3])
        )
        for i in item[4]:
            sout.write("\t%s" % i)
        sout.write("\n")

# method2 = list()    # hc_sa_state
# method3 = list()    # icpso

# method_scores = list()
# for i in range(20):
#     method1 = hc_sa_set(opv, varg, 10, set_states, verbose=False)
#     method2 = hc_sa_state(opv, varg, state_states, dimsizes, verbose=False)
#     method3 = icpso(opv, varg, 1000, state_states, dimsizes,
#                     0.33, 0.33, 0.34, 0.1, stop, verbose=False)
#     # append scores
#     method_scores.append([method1["X_gbest_scr"],
#                           method2["X_gbest_scr"],
#                           method3["X_gbest_scr"]])
#     # print iteration
#     print("iter", i, "done", flush=True)
