import numpy
import pandas
import os, sys

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from util.error_subroutines import *
#from util.proc_subroutines import *


class breeding_program:
    """docstring for selection."""

    ############################################################################
    ############################# Reserved methods #############################
    ############################################################################

    def __init__(self, dtype_method, dtype_algorithm, dtype_seed, dtype_cycle,
    dtype_score, dtype_gebv, dtype_phase):
        check_is_string_or_object_dtype(dtype_method, "dtype_method")
        check_is_string_or_object_dtype(dtype_algorithm, "dtype_algorithm")
        check_is_numeric_dtype(dtype_seed, "dtype_seed")
        check_is_numeric_dtype(dtype_cycle, "dtype_cycle")
        check_is_numeric_dtype(dtype_score, "dtype_score")
        check_is_numeric_dtype(dtype_gebv, "dtype_gebv")
        check_is_numeric_dtype(dtype_phase, "dtype_phase")

        self._dtype_method = dtype_method
        self._dtype_algorithm = dtype_algorithm
        self._dtype_seed = dtype_seed
        self._dtype_cycle = dtype_cycle
        self._dtype_score = dtype_score
        self._dtype_gebv = dtype_gebv
        self._dtype_phase = dtype_phase

        self._pop_method = []
        self._pop_algorithm = []
        self._pop_seed = []
        self._pop_cycle = []
        self._pop_score = []
        self._pop_gebv = []
        self._pop_phase1 = []
        self._pop_phase2 = []

        self._sel_method = []
        self._sel_algorithm = []
        self._sel_seed = []
        self._sel_cycle = []
        self._sel_score = []
        self._sel_gebv = []
        self._sel_phase1 = []
        self._sel_phase2 = []


    ############################################################################
    ################################ Properties ################################
    ############################################################################

    def dtype_method():
        doc = "The dtype_method property."
        def fget(self):
            return self._dtype_method
        def fset(self, value):
            self._dtype_method = value
        def fdel(self):
            del self._dtype_method
        return locals()
    dtype_method = property(**dtype_method())

    def dtype_algorithm():
        doc = "The dtype_algorithm property."
        def fget(self):
            return self._dtype_algorithm
        def fset(self, value):
            self._dtype_algorithm = value
        def fdel(self):
            del self._dtype_algorithm
        return locals()
    dtype_algorithm = property(**dtype_algorithm())

    def dtype_seed():
        doc = "The dtype_seed property."
        def fget(self):
            return self._dtype_seed
        def fset(self, value):
            self._dtype_seed = value
        def fdel(self):
            del self._dtype_seed
        return locals()
    dtype_seed = property(**dtype_seed())

    def dtype_cycle():
        doc = "The dtype_cycle property."
        def fget(self):
            return self._dtype_cycle
        def fset(self, value):
            self._dtype_cycle = value
        def fdel(self):
            del self._dtype_cycle
        return locals()
    dtype_cycle = property(**dtype_cycle())

    def dtype_score():
        doc = "The dtype_score property."
        def fget(self):
            return self._dtype_score
        def fset(self, value):
            self._dtype_score = value
        def fdel(self):
            del self._dtype_score
        return locals()
    dtype_score = property(**dtype_score())

    def dtype_gebv():
        doc = "The dtype_gebv property."
        def fget(self):
            return self._dtype_gebv
        def fset(self, value):
            self._dtype_gebv = value
        def fdel(self):
            del self._dtype_gebv
        return locals()
    dtype_gebv = property(**dtype_gebv())

    def dtype_phase():
        doc = "The dtype_phase property."
        def fget(self):
            return self._dtype_phase
        def fset(self, value):
            self._dtype_phase = value
        def fdel(self):
            del self._dtype_phase
        return locals()
    dtype_phase = property(**dtype_phase())

    def pop_method():
        doc = "The pop_method property."
        def fget(self):
            return self._pop_method
        def fset(self, value):
            self._pop_method = value
        def fdel(self):
            del self._pop_method
        return locals()
    pop_method = property(**pop_method())

    def sel_method():
        doc = "The sel_method property."
        def fget(self):
            return self._sel_method
        def fset(self, value):
            self._sel_method = value
        def fdel(self):
            del self._sel_method
        return locals()
    sel_method = property(**sel_method())

    def pop_algorithm():
        doc = "The pop_algorithm property."
        def fget(self):
            return self._pop_algorithm
        def fset(self, value):
            self._pop_algorithm = value
        def fdel(self):
            del self._pop_algorithm
        return locals()
    pop_algorithm = property(**pop_algorithm())

    def sel_algorithm():
        doc = "The sel_algorithm property."
        def fget(self):
            return self._sel_algorithm
        def fset(self, value):
            self._sel_algorithm = value
        def fdel(self):
            del self._sel_algorithm
        return locals()
    sel_algorithm = property(**sel_algorithm())

    def pop_seed():
        doc = "The pop_seed property."
        def fget(self):
            return self._pop_seed
        def fset(self, value):
            self._pop_seed = value
        def fdel(self):
            del self._pop_seed
        return locals()
    pop_seed = property(**pop_seed())

    def sel_seed():
        doc = "The sel_seed property."
        def fget(self):
            return self._sel_seed
        def fset(self, value):
            self._sel_seed = value
        def fdel(self):
            del self._sel_seed
        return locals()
    sel_seed = property(**sel_seed())

    def pop_cycle():
        doc = "The pop_cycle property."
        def fget(self):
            return self._pop_cycle
        def fset(self, value):
            self._pop_cycle = value
        def fdel(self):
            del self._pop_cycle
        return locals()
    pop_cycle = property(**pop_cycle())

    def sel_cycle():
        doc = "The sel_cycle property."
        def fget(self):
            return self._sel_cycle
        def fset(self, value):
            self._sel_cycle = value
        def fdel(self):
            del self._sel_cycle
        return locals()
    sel_cycle = property(**sel_cycle())

    def pop_score():
        doc = "The pop_score property."
        def fget(self):
            return self._pop_score
        def fset(self, value):
            self._pop_score = value
        def fdel(self):
            del self._pop_score
        return locals()
    pop_score = property(**pop_score())

    def sel_score():
        doc = "The sel_score property."
        def fget(self):
            return self._sel_score
        def fset(self, value):
            self._sel_score = value
        def fdel(self):
            del self._sel_score
        return locals()
    sel_score = property(**sel_score())

    def pop_phase1():
        doc = "The pop_phase1 property."
        def fget(self):
            return self._pop_phase1
        def fset(self, value):
            self._pop_phase1 = value
        def fdel(self):
            del self._pop_phase1
        return locals()
    pop_phase1 = property(**pop_phase1())

    def pop_phase2():
        doc = "The pop_phase2 property."
        def fget(self):
            return self._pop_phase2
        def fset(self, value):
            self._pop_phase2 = value
        def fdel(self):
            del self._pop_phase2
        return locals()
    pop_phase2 = property(**pop_phase2())

    def sel_phase1():
        doc = "The sel_phase1 property."
        def fget(self):
            return self._sel_phase1
        def fset(self, value):
            self._sel_phase1 = value
        def fdel(self):
            del self._sel_phase1
        return locals()
    sel_phase1 = property(**sel_phase1())

    def sel_phase2():
        doc = "The sel_phase2 property."
        def fget(self):
            return self._sel_phase2
        def fset(self, value):
            self._sel_phase2 = value
        def fdel(self):
            del self._sel_phase2
        return locals()
    sel_phase2 = property(**sel_phase2())

    def sel_gebv():
        doc = "The sel_gebv property."
        def fget(self):
            return self._sel_gebv
        def fset(self, value):
            self._sel_gebv = value
        def fdel(self):
            del self._sel_gebv
        return locals()
    sel_gebv = property(**sel_gebv())

    def pop_gebv():
        doc = "The pop_gebv property."
        def fget(self):
            return self._pop_gebv
        def fset(self, value):
            self._pop_gebv = value
        def fdel(self):
            del self._pop_gebv
        return locals()
    pop_gebv = property(**pop_gebv())

    ############################################################################
    ################################# Methods ##################################
    ############################################################################

    def history_add_population(self, method, algorithm, seed, cycle, score,
        gebv, phase1, phase2):
        ####
        method = numpy.object_(method)
        algorithm = numpy.object_(algorithm)
        seed = numpy.array(seed)
        cycle = numpy.array(cycle)
        score = numpy.array(score)
        gebv = numpy.array(gebv)
        phase1 = numpy.array(phase1)
        phase2 = numpy.array(phase2)

        if method.ndim == 1:
            method = numpy.repeat(method, len(score))
        if algorithm.ndim == 1:
            algorithm = numpy.repeat(algorithm, len(score))
        if seed.ndim == 1:
            seed = numpy.repeat(seed, len(score))
        if cycle.ndim == 1:
            cycle = numpy.repeat(cycle, len(score))
        if phase1.ndim == 1:
            phase1 = phase1[numpy.newaxis]
        if phase2.ndim == 1:
            phase2 = phase2[numpy.newaxis]

        self._pop_method.append(method)
        self._pop_algorithm.append(algorithm)
        self._pop_seed.append(seed)
        self._pop_score.append(score)
        self._pop_gebv.append(gebv)
        self._pop_cycle.append(cycle)
        self._pop_phase1.append(phase1)
        self._pop_phase2.append(phase2)

    def history_add_selection(self, method, algorithm, seed, cycle, score,
        gebv, phase1, phase2):
        ####
        method = numpy.object_(method)
        algorithm = numpy.object_(algorithm)
        seed = numpy.array(seed)
        cycle = numpy.array(cycle)
        score = numpy.array(score)
        gebv = numpy.array(gebv)
        phase1 = numpy.array(phase1)
        phase2 = numpy.array(phase2)

        if method.ndim == 1:
            method = numpy.repeat(method, len(score))
        if algorithm.ndim == 1:
            algorithm = numpy.repeat(algorithm, len(score))
        if seed.ndim == 1:
            seed = numpy.repeat(seed, len(score))
        if cycle.ndim == 1:
            cycle = numpy.repeat(cycle, len(score))
        if phase1.ndim == 1:
            phase1 = phase1[numpy.newaxis]
        if phase2.ndim == 1:
            phase2 = phase2[numpy.newaxis]

        self._sel_method.append(method)
        self._sel_algorithm.append(algorithm)
        self._sel_seed.append(seed)
        self._sel_cycle.append(cycle)
        self._sel_score.append(score)
        self._sel_gebv.append(gebv)
        self._sel_phase1.append(phase1)
        self._sel_phase2.append(phase2)

    def _concatenate(self):
        self._pop_method = [numpy.concatenate(self._pop_method, axis=0)]
        self._sel_method = [numpy.concatenate(self._sel_method, axis=0)]
        self._pop_algorithm = [numpy.concatenate(self._pop_algorithm, axis=0)]
        self._sel_algorithm = [numpy.concatenate(self._sel_algorithm, axis=0)]
        self._pop_seed = [numpy.concatenate(self._pop_seed, axis=0)]
        self._sel_seed = [numpy.concatenate(self._sel_seed, axis=0)]
        self._pop_cycle = [numpy.concatenate(self._pop_cycle, axis=0)]
        self._sel_cycle = [numpy.concatenate(self._sel_cycle, axis=0)]
        self._pop_score = [numpy.concatenate(self._pop_score, axis=0)]
        self._sel_score = [numpy.concatenate(self._sel_score, axis=0)]
        self._pop_gebv = [numpy.concatenate(self._pop_gebv, axis=0)]
        self._sel_gebv = [numpy.concatenate(self._sel_gebv, axis=0)]
        self._pop_phase1 = [numpy.concatenate(self._pop_phase1, axis=0)]
        self._sel_phase1 = [numpy.concatenate(self._sel_phase1, axis=0)]
        self._pop_phase2 = [numpy.concatenate(self._pop_phase2, axis=0)]
        self._sel_phase2 = [numpy.concatenate(self._sel_phase2, axis=0)]

    def history_to_df(self, zfill = 7):
        self._concatenate()

        pop_dict = {
            "method" : self._pop_method[0],
            "algorithm" : self._pop_algorithm[0],
            "seed" : self._pop_seed[0],
            "cycle" : self._pop_cycle[0],
            "score" : self._pop_score[0],
            "gebv" : self._pop_gebv[0]
        }
        sel_dict = {
            "method" : self._sel_method[0],
            "algorithm" : self._sel_algorithm[0],
            "seed" : self._sel_seed[0],
            "cycle" : self._sel_cycle[0],
            "score" : self._sel_score[0],
            "gebv" : self._sel_gebv[0]
        }

        pop_p1head = [
            "p1"+str(i).zfill(zfill)
            for i in range(self._pop_phase1[0].shape[1])
        ]
        pop_p2head = [
            "p2"+str(i).zfill(zfill)
            for i in range(self._pop_phase2[0].shape[1])
        ]
        sel_p1head = [
            "p1"+str(i).zfill(zfill)
            for i in range(self._sel_phase1[0].shape[1])
        ]
        sel_p2head = [
            "p2"+str(i).zfill(zfill)
            for i in range(self._sel_phase2[0].shape[1])
        ]

        for i,header in enumerate(pop_p1head):
            pop_dict[header] = self._pop_phase1[0][:,i].copy()
        for i,header in enumerate(pop_p2head):
            pop_dict[header] = self._pop_phase2[0][:,i].copy()
        for i,header in enumerate(sel_p1head):
            sel_dict[header] = self._sel_phase1[0][:,i].copy()
        for i,header in enumerate(sel_p2head):
            sel_dict[header] = self._sel_phase2[0][:,i].copy()

        pop_df = pandas.DataFrame(pop_dict)
        sel_df = pandas.DataFrame(sel_dict)

        return pop_df, sel_df
