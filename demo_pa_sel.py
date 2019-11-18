import numpy
import time

# import optimizing functions
from pybropt.algo.hc import hc_sa_set, hc_sa_state
from pybropt.algo.icpso import icpso
from pybropt.algo.ga import ga_fps

# import scoring function
from pybropt.objfn.opb import opv
from pybropt.objfn.pa import pa_sel



################################################################################
################################################################################
################################################################################


# generate binary marker data
n_lines = 20
n_markers = 1000
n_phases = 2
n_selected = 10

# seed random number
numpy.random.seed(111019)
# make markers
markers = numpy.random.binomial(1, 0.1, (n_phases,n_lines,n_markers))
# make effects
effects = numpy.random.normal(0,2,n_markers)

#rslice, geno, coeff, tfreq, tfreq_edge, tld, tld_edge,
#        cycles, d, lgroup_size

# calculate target frequencies: if effect >= 0, then tfreq = 1; else tfreq = 0
tfreq = numpy.where(effects >= 0.0, 1.0, 0.0)
# calculate the edge frequencies: if tfreq > 0.5, then edge = 0; else edge = 1
tfreq_edge = numpy.where(tfreq >= 0.5, 0.0, 1.0)
# set the target ld
tr_sq = 1.0
# calculate edge target ld
tr_sq_edge = 0.0 if tr_sq >= 0.5 else 1.0
# calculate abs(effects)
abs_effects = numpy.abs(effects)
# simulate a single chromosome of length 'n_markers'
gmap = numpy.sort(numpy.random.uniform(0.0,1.0,n_markers))
# make linkage group sizes
gmap_group_sizes = numpy.array([n_markers])
# set the number of breeding cycles
generations = 5
# select a random subset of individuals to score
selections = numpy.random.choice(numpy.arange(n_lines), n_selected, False)


# rslice, geno, coeff, tfreq, tfreq_edge, tld, tld_edge, cycles, d, lgroup_size
score = pa_sel(selections, markers, abs_effects, tfreq, tfreq_edge, tr_sq,
               tr_sq_edge, generations, gmap, gmap_group_sizes)

print("Selections:", selections)
print("Score:", score)

# varg = {"geno": markers,
#         "coeff": abs_effects,
#         "tfreq": tfreq,
#         "tfreq_edge": tfreq_edge,
#         "tld": tr_sq,
#         "tld_edge": tr_sq_edge,
#         "cycles": generations,
#         "d": gmap,
#         "lgroup_size": gmap_group_sizes}
# set_states = numpy.arange(n_lines, dtype='uint32')
# method1 = hc_sa_set(pas, varg, 10, set_states, verbose=True)

# state_states = numpy.tile(set_states, 10) # select 10
# dimsizes = numpy.repeat(n_lines, 10)

# hcoeff = markers * effects
# varg = {"hcoeff": hcoeff}

# def stop(x):
#     return x < 50
#
# # test the efficacy of each method
# method_scores = list()
# for i in range(100):
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
#
# # write results to file
# with open('method_scores.txt', 'w') as fout:
#     for item in method_scores:
#         fout.write('%s\t%s\t%s\n' % (item[0], item[1], item[2]))

# print()
# print("history =")
# for item in pso_pop["history"]:
#     print(item)
# print("X_gbest_smpl", pso_pop["X_gbest_smpl"])
# print("X_gbest_scr", pso_pop["X_gbest_scr"])
#
# print("\n\n\n")
#
# print("history =")
# for item in state_history["history"]:
#     print(item)
# print("X_gbest_pos", state_history["X_gbest_pos"])
# print("X_gbest_scr", state_history["X_gbest_scr"])
#
# print("\n\n\n")
#
# print("history =")
# for item in set_history["history"]:
#     print(item)
# print("X_gbest_pos", set_history["X_gbest_pos"])
# print("X_gbest_scr", set_history["X_gbest_scr"])
