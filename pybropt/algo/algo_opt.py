import numpy
import pandas

# HACK: add more directories to path.
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
# import our libraries
from util.error_subroutines import *

# hold the algo_opt class
class algo_opt:
    """
    Holds algorithm results for a generic algorithm.

    Attributes
    ==========
    dtype_iter : numpy.dtype
        Data type of the iter vector.
    dtype_score : numpy.dtype
        Data type of the X score vector.
    dtype_pos : numpy.dtype
        Data type of the X position matrix.
    iter : numpy.ndarray, None
        Vector of iteration numbers.
    score : numpy.ndarray, None
        Vector of X scores.
    pos : numpy.ndarray, None
        2D matrix of X positions. Each row is a position.
    """

    ############################################################################
    ############################# Reserved methods #############################
    ############################################################################

    def __init__(self, dtype_iter, dtype_score, dtype_pos):
        """
        Constructor for the algo_result class.

        Parameters
        ==========
        dtype_iter : numpy.dtype
            Data type of the iter vector.
        dtype_score : numpy.dtype
            Data type of the X score vector.
        dtype_pos : numpy.dtype
            Data type of the X position vector.
        """
        # check for correct input data types.
        check_is_integer_dtype(dtype_iter, "dtype_iter")
        check_is_numeric_dtype(dtype_score, "dtype_score")
        check_is_numeric_or_bool_dtype(dtype_pos, "dtype_pos")

        # assign dtype arguments to private variables
        self._dtype_iter = dtype_iter
        self._dtype_score = dtype_score
        self._dtype_pos = dtype_pos

        # assign empty variables to iter, score, pos private variables
        self._iter = None
        self._score = None
        self._pos = None

    def __len__(self):
        """
        Get the length of the score array.
        """
        l = len(self._score)
        return l

    ############################################################################
    ################################ Properties ################################
    ############################################################################

    def dtype_iter():
        doc = "The dtype_iter property. Holds the data type of "
        def fget(self):
            return self._dtype_iter
        def fset(self, value):
            error_readonly("dtype_iter")
        def fdel(self):
            error_readonly("dtype_iter")
        return locals()
    dtype_iter = property(**dtype_iter())

    def dtype_score():
        doc = "The dtype_score property."
        def fget(self):
            return self._dtype_score
        def fset(self, value):
            error_readonly("dtype_score")
        def fdel(self):
            error_readonly("dtype_score")
        return locals()
    dtype_score = property(**dtype_score())

    def dtype_pos():
        doc = "The dtype_pos property."
        def fget(self):
            return self._dtype_pos
        def fset(self, value):
            error_readonly("dtype_pos")
        def fdel(self):
            error_readonly("dtype_pos")
        return locals()
    dtype_pos = property(**dtype_pos())

    def iter():
        doc = "The iter property."
        def fget(self):
            return self._iter
        def fset(self, value):
            self._iter = value
        def fdel(self):
            del self._iter
        return locals()
    iter = property(**iter())

    def score():
        doc = "The score property."
        def fget(self):
            return self._score
        def fset(self, value):
            self._score = value
        def fdel(self):
            del self._score
        return locals()
    score = property(**score())

    def pos():
        doc = "The pos property."
        def fget(self):
            return self._pos
        def fset(self, value):
            self._pos = value
        def fdel(self):
            del self._pos
        return locals()
    pos = property(**pos())

    ############################################################################
    ################################# Methods ##################################
    ############################################################################

    def history_add(self, iter, score, pos):
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
        """
        ### convert everything to numpy.ndarray and do some preprocessing
        iter = numpy.array(iter)                    # convert iter to ndarray
        score = numpy.array(score)                  # convert score to ndarray
        pos = numpy.array(pos)                      # convert pos to ndarray
        if pos.ndim == 1:                           # if 1D array
            pos = pos[numpy.newaxis]                # add dimension
        if iter.size == 1:                          # if single value
            iter = numpy.repeat(iter, score.size)   # repeat to fill

        ### shape checking
        check_matrix_ndim(pos, "pos", 2)                    # check that pos is a 2d matrix.
        check_matrix_axis_len(pos, "pos", 0, score.size)    # check that number of rows in 'pos' aligns with size of 'score'
        check_matrix_size(iter, "iter", score.size)         # check if matrix sizes are compatible

        ### dtype checking
        check_matrix_dtype(iter, "iter", self._dtype_iter)      # check iter matrix dtype
        check_matrix_dtype(score, "score", self._dtype_score)   # check score matrix dtype
        check_matrix_dtype(pos, "pos", self._dtype_pos)         # check pos matrix dtype

        ### finally append the matrices
        # append iter
        if self._iter is None:
            self._iter = iter.copy()
        else:
            self._iter = numpy.append(self._iter, iter)
        # append score
        if self._score is None:
            self._score = score.copy()
        else:
            self._score = numpy.append(self._score, score)
        # append pos
        if self._pos is None:
            self._pos = pos.copy()
        else:
            self._pos = numpy.append(self._pos, pos)

        return

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
        # make a dictionary to construct the pandas.DataFrame
        df_dict = {
            "iter" : self._iter.copy(),
            "score" : self._score.copy()
        }

        # make labels for the X position headers
        headers = ["x"+str(i).zfill(zfill) for i in range(self._pos.shape[1])]

        # add columns + header name to df_dict
        for i,header in enumerate(headers):
            df_dict[header] = self._pos[:,i].copy()

        # make DataFrame
        df = pandas.DataFrame(df_dict)

        return df

    def gbest_ix(self):
        """
        Get the index of the global best position. If multiple global best
        scores are equivalent, return the first instance.

        Returns
        =======
        ix : int
            An index of the global best score.
        """
        # use argmax to get the first instance of the best score.
        ix = self._score.argmax()

        return ix

    def gbest(self):
        """
        Get global best iter, score, position as a tuple. If multiple global
        best scores are equivalent, return the first instance.

        Returns
        =======
        gbest : tuple
            A tuple of global best iter, score, pos. The tuple order is:
                (iter, score, pos)
        """
        # get best score index
        ix = self.gbest_ix()

        # construct tuple
        gbest = (
            self._iter[ix],
            self._score[ix],
            self._pos[ix,:]
        )

        # return gbest
        return gbest
