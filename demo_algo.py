import numpy
import time

# import optimizing functions
from pybropt.algo.hc import hc_sa_set, hc_sa_state
from pybropt.algo.icpso import icpso
from pybropt.algo.ga import ga_fps

# import scoring function
from pybropt.objfn.opb import opv


################################################################################
################################################################################
################################################################################


# generate binary marker data
numpy.random.seed(111019)
markers = numpy.random.binomial(1, 0.1, (2,200,10000))
effects = numpy.random.normal(0,2,10000)
set_states = numpy.arange(200, dtype='uint32')
state_states = numpy.tile(set_states, 10) # select 10
dimsizes = numpy.repeat(200, 10)

hcoeff = markers * effects
varg = {"hcoeff": hcoeff}

def stop(x):
    return x < 50

# test the efficacy of each method
method_scores = list()
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
#
# # write results to file
# with open('method_scores.txt', 'a') as fout:
#     for item in method_scores:
#         fout.write('%s\t%s\t%s\n' % (item[0], item[1], item[2]))

# specify a custom function for calculating fps
def fpstrans(scr):
    si = numpy.argsort(scr) # get sort indices
    up = numpy.int64(len(si))
    lo = numpy.int64(numpy.floor(up * 0.5))
    tmp = numpy.empty(up, dtype=numpy.float64)
    for i in range(up):
        tmp[i] = 1.0/(up-lo) if i in si[lo:up] else 0.0
    print(tmp)
    return numpy.cumsum(tmp)

def stop2(i):
    return i < 100

# calculate ga
method4 = ga_fps(opv, varg, 50, state_states,
                 dimsizes, 0.1, 0.0001, stop2, verbose=False)
for item in method4["history"]:
    print(item)

# for i in range(120):
#     method4 = ga_fps(opv, varg, 200, state_states,
#                      dimsizes, 0.5, 0.01, stop, verbose=False)
#     method_scores.append([method4["X_gbest_scr"]])
#     print("iter", i, "done", flush=True)
# with open('ga_scores.txt', 'a') as fout:
#     for item in method_scores:
#         fout.write('%s\n' % (item))

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
