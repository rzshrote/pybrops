import sys, os
import numpy
import pandas

# HACK: add more directories to path.
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
# import our libraries
from util.error_subroutines import *
from algo.pso_opt import pso_opt

class icpso_opt(pso_opt):
    """docstring for icpso_opt."""

    ############################################################################
    ############################# Reserved methods #############################
    ############################################################################

    def __init__(self,
                 ssize, inertia_wt, pbest_comp, gbest_comp, scale_factor,
                 dtype_iter, dtype_score, dtype_pos, dtype_vel, dtype_smpl):
        super(icpso_opt, self).__init__(
            ssize, inertia_wt, pbest_comp, gbest_comp,
            dtype_iter, dtype_score, dtype_pos, dtype_vel
        )

        # check for correct input data types.
        check_is_numeric_or_bool_dtype(dtype_smpl, "dtype_smpl")
        check_is_numeric(scale_factor, "scale_factor")

        # assign arguments or their defaults to private variables
        self._dtype_smpl = dtype_smpl
        self._smpl = []
        self._scale_factor = scale_factor

    ############################################################################
    ################################ Properties ################################
    ############################################################################

    def dtype_smpl():
        doc = "The dtype_smpl property."
        def fget(self):
            return self._dtype_smpl
        def fset(self, value):
            error_readonly("dtype_smpl")
        def fdel(self):
            error_readonly("dtype_smpl")
        return locals()
    dtype_smpl = property(**dtype_smpl())

    def smpl():
        doc = "The smpl property."
        def fget(self):
            return self._smpl
        def fset(self, value):
            self._smpl = value
        def fdel(self):
            del self._smpl
        return locals()
    smpl = property(**smpl())

    def scale_factor():
        doc = "The scale_factor property."
        def fget(self):
            return self._scale_factor
        def fset(self, value):
            self._scale_factor = value
        def fdel(self):
            del self._scale_factor
        return locals()
    scale_factor = property(**scale_factor())

    ############################################################################
    ################################# Methods ##################################
    ############################################################################

    def history_add(self, iter, score, pos, vel, smpl, check_ssize = True):
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
        vel : array-like
            1D or 2D matrix of velocities.
        smpl : array-like
            1D or 2D matrix of distribution samples.
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
        smpl = numpy.array(smpl)                    # convert smpl to ndarray
        if pos.ndim == 1:                           # if 1D array
            pos = pos[numpy.newaxis]                # add dimension
        if vel.ndim == 1:                           # if 1D array
            vel = vel[numpy.newaxis]                # add dimension
        if smpl.ndim == 1:                          # if 1D array
            smpl = smpl[numpy.newaxis]              # add dimension
        if iter.size == 1:                          # if single value
            iter = numpy.repeat(iter, score.size)   # repeat to fill

        ### shape checking
        check_matrix_ndim(pos, "pos", 2)                    # check if pos is 2D
        check_matrix_ndim(vel, "vel", 2)                    # check is vel is 2D
        check_matrix_axis_len(pos, "pos", 0, score.size)    # check pos row count
        check_matrix_axis_len(vel, "vel", 0, score.size)    # check vel row count
        check_matrix_axis_len(smpl, "smpl", 0, score.size)  # check smpl row count
        check_matrix_size(iter, "iter", score.size)         # check iter element count
        if check_ssize:
            check_matrix_size(score, "score", self._ssize)      # check score element count
            check_matrix_axis_len(pos, "pos", 0, self._ssize)   # check pos row count
            check_matrix_axis_len(vel, "vel", 0, self._ssize)   # check vel row count
            check_matrix_axis_len(smpl, "smpl", 0, self._ssize) # check vel row count


        ### dtype checking
        check_matrix_dtype(iter, "iter", self._dtype_iter)      # check iter matrix dtype
        check_matrix_dtype(score, "score", self._dtype_score)   # check score matrix dtype
        check_matrix_dtype(pos, "pos", self._dtype_pos)         # check pos matrix dtype
        check_matrix_dtype(vel, "vel", self._dtype_vel)         # check vel matrix dtype
        check_matrix_dtype(smpl, "smpl", self._dtype_smpl)      # check smpl matrix dtype

        ### finally append the matrices
        self._iter.append(iter)
        self._score.append(score)
        self._pos.append(pos)
        self._vel.append(vel)
        self._smpl.append(smpl)

        return

    def _concatenate(self):
        super(icpso_opt, self)._concatenate()
        if len(self._smpl) > 1:
            self._smpl = [numpy.concatenate(self._smpl, axis=0)]

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
        shead = ["s"+str(i).zfill(zfill) for i in range(self._smpl[0].shape[1])]

        # add columns + header name to df_dict
        for i,header in enumerate(xhead):
            df_dict[header] = self._pos[0][:,i].copy()
        for i,header in enumerate(vhead):
            df_dict[header] = self._vel[0][:,i].copy()
        for i,header in enumerate(shead):
            df_dict[header] = self._smpl[0][:,i].copy()

        # make DataFrame
        df = pandas.DataFrame(df_dict)

        return df

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
            self._vel[0][ix,:],
            self._smpl[0][ix,:]
        )

        # return gbest
        return gbest

    def pbest(self, ptcl = None, iter = None):
        """
        Grab personal best data
        """
        # grab indices
        ix = self.pbest_ix(iter)

        # grab pbest iter, score, pos, vel
        pbest_iter = self._iter[0][ix][ptcl]
        pbest_score = self._score[0][ix][ptcl]
        pbest_pos = self._pos[0][ix][ptcl,:]
        pbest_vel = self._vel[0][ix][ptcl,:]
        pbest_smpl = self._smpl[0][ix][ptcl,:]

        # pack into tuple
        pbest = (pbest_iter, pbest_score, pbest_pos, pbest_vel, pbest_smpl)

        return pbest
