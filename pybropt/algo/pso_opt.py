import sys, os
import pandas
import numpy

# HACK: add more directories to path.
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
# import our libraries
from util.error_subroutines import *
from algo.algo_opt import algo_opt

class pso_opt(algo_opt):
    """docstring for pso_opt."""

    ############################################################################
    ############################# Reserved methods #############################
    ############################################################################

    def __init__(self,
                 ssize, inertia_wt, pbest_comp, gbest_comp,
                 dtype_iter, dtype_score, dtype_pos, dtype_vel):
        """
        ssize : int
            Size of swarm
        inertia_wt : floating
        pbest_comp : floating
        gbest_comp : floatting
        """
        # call algo_opt constructor to initilize some variables
        super(pso_opt, self).__init__(dtype_iter, dtype_score, dtype_pos)

        # check for correct input data types.
        check_is_numeric_or_bool_dtype(dtype_vel, "dtype_vel")
        check_is_integer(ssize, "ssize")
        check_is_numeric(inertia_wt, "inertia_wt")
        check_is_numeric(pbest_comp, "pbest_comp")
        check_is_numeric(gbest_comp, "gbest_comp")

        # assign arguments (or their default) to private variables
        self._dtype_vel = dtype_vel
        self._vel = []                # assign None to vel private variable
        self._ssize = ssize
        self._inertia_wt = inertia_wt
        self._pbest_comp = pbest_comp
        self._gbest_comp = gbest_comp

    ############################################################################
    ################################ Properties ################################
    ############################################################################

    def dtype_vel():
        doc = "The dtype_vel property."
        def fget(self):
            return self._dtype_vel
        def fset(self, value):
            error_readonly("dtype_vel")
        def fdel(self):
            error_readonly("dtype_vel")
        return locals()
    dtype_vel = property(**dtype_vel())

    def vel():
        doc = "The vel property."
        def fget(self):
            return self._vel
        def fset(self, value):
            self._vel = value
        def fdel(self):
            del self._vel
        return locals()
    vel = property(**vel())

    def ssize():
        doc = "The ssize property."
        def fget(self):
            return self._ssize
        def fset(self, value):
            self._ssize = value
        def fdel(self):
            del self._ssize
        return locals()
    ssize = property(**ssize())

    def inertia_wt():
        doc = "The inertia_wt property."
        def fget(self):
            return self._inertia_wt
        def fset(self, value):
            self._inertia_wt = value
        def fdel(self):
            del self._inertia_wt
        return locals()
    inertia_wt = property(**inertia_wt())

    def pbest_comp():
        doc = "The pbest_comp property."
        def fget(self):
            return self._pbest_comp
        def fset(self, value):
            self._pbest_comp = value
        def fdel(self):
            del self._pbest_comp
        return locals()
    pbest_comp = property(**pbest_comp())

    def gbest_comp():
        doc = "The gbest_comp property."
        def fget(self):
            return self._gbest_comp
        def fset(self, value):
            self._gbest_comp = value
        def fdel(self):
            del self._gbest_comp
        return locals()
    gbest_comp = property(**gbest_comp())

    ############################################################################
    ################################# Methods ##################################
    ############################################################################

    def history_add(self, iter, score, pos, vel, check_ssize = True):
        """
        Add algorithm history to this object.

        Parameters
        ==========
        iter : int, array-like
            Iteration number.
        score : int, array-like
            Score of the position.
        pos : array-like
            1D or 2D matrix of positions.
        check_ssize : boolean, default = True
            Check the number of rows in the input data to make sure it aligns
            with ssize. Basically, the default is put all particles in at once
            for a given iteration.
        """
        # We do not call super method because we don't want to prematurely
        # append data to the internal matrices.

        ### convert everything to numpy.ndarray and do some preprocessing
        iter = numpy.array(iter)                    # convert iter to ndarray
        score = numpy.array(score)                  # convert score to ndarray
        pos = numpy.array(pos)                      # convert pos to ndarray
        vel = numpy.array(vel)                      # convert vel to ndarray
        if pos.ndim == 1:                           # if 1D array
            pos = pos[numpy.newaxis]                # add dimension
        if vel.ndim == 1:                           # if 1D array
            vel = vel[numpy.newaxis]                # add dimension
        if iter.size == 1:                          # if single value
            iter = numpy.repeat(iter, score.size)   # repeat to fill

        ### shape checking
        check_matrix_ndim(pos, "pos", 2)                    # check if pos is 2D
        check_matrix_ndim(vel, "vel", 2)                    # check is vel is 2D
        check_matrix_axis_len(pos, "pos", 0, score.size)    # check pos row count
        check_matrix_axis_len(vel, "vel", 0, score.size)    # check vel row count
        check_matrix_size(iter, "iter", score.size)         # check iter element count
        if check_ssize:
            check_matrix_size(score, "score", self._ssize)      # check score element count
            check_matrix_axis_len(pos, "pos", 0, self._ssize)   # check pos row count
            check_matrix_axis_len(vel, "vel", 0, self._ssize)   # check vel row count


        ### dtype checking
        check_matrix_dtype(iter, "iter", self._dtype_iter)      # check iter matrix dtype
        check_matrix_dtype(score, "score", self._dtype_score)   # check score matrix dtype
        check_matrix_dtype(pos, "pos", self._dtype_pos)         # check pos matrix dtype
        check_matrix_dtype(vel, "vel", self._dtype_vel)         # check vel matrix dtype

        ### finally append the matrices
        self._iter.append(iter)
        self._score.append(score)
        self._pos.append(pos)
        self._vel.append(vel)

        return

    def _concatenate(self):
        super(pso_opt, self)._concatenate()
        if len(self._vel) > 1:
            self._vel = [numpy.concatenate(self._vel, axis=0)]


    def history_to_df(self, zfill = 3):
        """
        Convert internal history to pandas.DataFrame

        Parameters
        ==========
        zfill : int
            Number of zeros to fill for X position labels

        Returns
        =======
        df : pandas.DataFrame
            A pandas DataFrame of the results.
        """
        self._concatenate()

        # make a dictionary to construct the pandas.DataFrame
        df_dict = {
            "iter" : self._iter[0],
            "score" : self._score[0]
        }

        # make labels for the X position headers
        xhead = ["x"+str(i).zfill(zfill) for i in range(self._pos[0].shape[1])]
        vhead = ["v"+str(i).zfill(zfill) for i in range(self._vel[0].shape[1])]

        # add columns + header name to df_dict
        for i,header in enumerate(xhead):
            df_dict[header] = self._pos[0][:,i].copy()
        for i,header in enumerate(vhead):
            df_dict[header] = self._vel[0][:,i].copy()

        # make DataFrame
        df = pandas.DataFrame(df_dict)

        return df

    def gbest_ix(self, iter = None):
        """
        Grab global best array index

        iter : int, default = None
            Assume that all iterations are in order. Calculate the personal best
            values for each particle up to, but not including the 'iter'th
            iteration.
            If None, grab gbest using all iterations.
            Example:
                iter = 0 -> grab nothing
                iter = 1 -> grab gbest from iter 0
                iter = n -> grab gbest from iter 0, 1, ..., n-1
        """
        self._concatenate()

        # test if length of score is divisible by ssize; raise error if not
        check_divisibility(self._score[0].size, "score.size", self._ssize, "ssize")

        # convert iter to an index on the score matrix
        iter = self._score[0].size if iter is None else iter * self._ssize

        # grab argmax index
        ix = self._score[0][:iter].argmax()

        return ix

    def gbest(self, iter = None):
        """
        Get global best iter, score, position, velocity as a tuple. If multiple
        global best scores are equivalent, return the first instance.

        Returns
        =======
        gbest : tuple
            A tuple of global best iter, score, pos, vel. The tuple order is:
                (iter, score, pos, vel)
        """
        # get best score index
        ix = self.gbest_ix(iter)

        # construct tuple
        gbest = (
            self._iter[0][ix],
            self._score[0][ix],
            self._pos[0][ix,:],
            self._vel[0][ix,:]
        )

        # return gbest
        return gbest

    def pbest_ix(self, iter = None):
        """
        Assume that internal data is divisible by

        iter : int, default = None
            Assume that all iterations are in order. Calculate the personal best
            values for each particle up to, but not including the 'iter'th
            iteration.
            If None, grab pbest using all iterations.
            Example:
                iter = 0 -> grab nothing
                iter = 1 -> grab pbest from iter 0
                iter = n -> grab pbest from iter 0, 1, ..., n-1
        """
        # test if length of score is divisible by ssize; raise error if not
        check_divisibility(self._score[0].size, "score.size", self._ssize, "ssize")

        # convert iter to an index on the score matrix
        iter = self._score[0].size if iter is None else iter * self._ssize

        # grab indices along score matrix
        ix = numpy.from_iter(
            (
                i + (self._score[0][i:iter:self._ssize].argmax() * self._ssize)
                for i in range(self._ssize)
            )
        )

        return ix

    def pbest(self, ptcl = None, iter = None):
        """
        Grab personal best data
        """
        # grab indices
        ix = self.pbest_ix(iter)

        # grab pbest iter, score, pos, vel
        pbest_iter = self._iter[ix][ptcl]
        pbest_score = self._score[ix][ptcl]
        pbest_pos = self._pos[ix][ptcl,:]
        pbest_vel = self._vel[ix][ptcl,:]

        # pack into tuple
        pbest = (pbest_iter, pbest_score, pbest_pos, pbest_vel)

        return pbest
