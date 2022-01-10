# TODO: figure out how rpy2 works and write a good, fast converter

# from rpy2 import robjects
# from rpy2.robjects.robject import RObjectMixin

# from rpy2.robjects.vectors import Array
# from rpy2.robjects.vectors import ByteArray
# from rpy2.robjects.vectors import ByteMatrix
# from rpy2.robjects.vectors import ByteSexpVector
# from rpy2.robjects.vectors import ByteVector
# from rpy2.robjects.vectors import ComplexArray
# from rpy2.robjects.vectors import ComplexMatrix
# from rpy2.robjects.vectors import ComplexSexpVector
# from rpy2.robjects.vectors import ComplexVector
# from rpy2.robjects.vectors import DataFrame
# from rpy2.robjects.vectors import DateVector
# from rpy2.robjects.vectors import FactorVector
# from rpy2.robjects.vectors import ListSexpVector
# from rpy2.robjects.vectors import ListVector
# from rpy2.robjects.vectors import Matrix
# from rpy2.robjects.vectors import PairlistSexpVector
# from rpy2.robjects.vectors import PairlistVector
# from rpy2.robjects.vectors import StrArray
# from rpy2.robjects.vectors import StrMatrix
# from rpy2.robjects.vectors import StrSexpVector
# from rpy2.robjects.vectors import StrVector
# from rpy2.robjects.vectors import Vector


from rpy2.robjects import baseenv
from rpy2.robjects import NULL
# vectors
from rpy2.robjects.vectors import BoolVector
from rpy2.robjects.vectors import FloatVector
from rpy2.robjects.vectors import IntVector
from rpy2.robjects.vectors import StrVector
# matrices
from rpy2.robjects.vectors import BoolMatrix
from rpy2.robjects.vectors import FloatMatrix
from rpy2.robjects.vectors import IntMatrix
from rpy2.robjects.vectors import StrMatrix
# arrays
# from rpy2.robjects.vectors import BoolArray
# from rpy2.robjects.vectors import FloatArray
# from rpy2.robjects.vectors import IntArray
# from rpy2.robjects.vectors import StrArray

__all__ = [
    "generic_numpy_to_R_Vector",
    "generic_numpy_to_R_Matrix",
    "rpy2_to_R_FactorVector",
    "numpy_to_R_BoolVector",
    "numpy_to_R_BoolFactorVector",
    "numpy_to_R_BoolMatrix",
    "numpy_to_R_FloatVector",
    "numpy_to_R_FloatFactorVector",
    "numpy_to_R_FloatMatrix",
    "numpy_to_R_IntVector",
    "numpy_to_R_IntFactorVector",
    "numpy_to_R_IntMatrix",
    "numpy_to_R_StrVector",
    "numpy_to_R_StrFactorVector",
    "numpy_to_R_StrMatrix",
    "numpy_to_R"
]

################################################################################
############################## Generic Functions ###############################
################################################################################
def generic_numpy_to_R_Vector(a, dtype, cls, Rfn_name = 'as.vector'):
    """
    Generic convert a numpy array to an R Matrix.
    Does not type check for speed.

    Parameters
    ----------
    a : numpy.ndarray
        Input numpy array. Can be any shape, but will be raveled/flattened.
    dtype : str, numpy.dtype
        Numpy dtype conversion requried to construct an R Vector of class 'cls'.
    cls : class
        R Vector class. Must have class method 'from_memoryview'.
    Rfn_name : str
        Name of function in R baseenv to call to create R Vector.
        Must be:
            as.vector
            as.logical
            as.double
            as.integer
            as.complex
            as.raw
            ...
            etc.

    Returns
    -------
    out : Matrix
        R Matrix class with base Vector type of 'cls'.
    """
    a_dtype = numpy.array(          # convert to dtype array
        a,                          # input array
        dtype = dtype,              # convert to dtype
        copy = False,               # only copy if necessary
        order = 'K'                 # keep ordering of input array
    )
    a_ravel = numpy.ravel(          # ravel/flatten array (no memory allocation)
        a_dtype,                    # input array
        order = 'K'                 # keep ordering of input array
    )
    R_vec = cls.from_memoryview(    # create IntVector from memory
        a_ravel.data                # raveled array memory
    )
    Rfn = baseenv[Rfn_name]         # get vector construction function from R
    R_res = Rfn(R_vec)              # call vector construction function from R
    return R_res

def generic_numpy_to_R_Matrix(a, dtype, cls, dimnames = NULL):
    """
    Generic convert a numpy array to an R Matrix.
    Does not type check for speed.

    Parameters
    ----------
    a : numpy.ndarray
        Input numpy matrix.
    dtype : str, numpy.dtype
        Numpy dtype conversion requried to construct an R Vector of class 'cls'.
    cls : class
        R Vector class. Must have class method 'from_memoryview'.
    **kwargs
        Additional keyword arguments for calling 'matrix' function from R.

    Returns
    -------
    out : Matrix
        R Matrix class with base Vector type of 'cls'.
    """
    nrow, ncol = a.shape            # get number of rows and columns
    a_dtype = numpy.array(          # convert to int32 array
        a,                          # input array
        dtype = dtype,              # convert to dtype
        copy = False,               # only copy if necessary
        order = 'K'                 # keep ordering of input array
    )
    a_ravel = numpy.ravel(          # ravel/flatten array (no memory allocation)
        a_dtype,                    # input array
        order = 'K'                 # keep ordering of input array
    )
    R_vec = cls.from_memoryview(    # create IntVector from memory
        a_ravel.data                # raveled array memory
    )
    Rfn = baseenv['matrix']         # get 'matrix' function from R
    a_flags = a_dtype.flags         # get array flags
    if a_flags['C_CONTIGUOUS']:     # if memory is C contiguous
        R_res = Rfn(                # call 'matrix' function from R
            R_vec,                  # input R vector
            nrow = nrow,            # number of rows
            ncol = ncol,            # number of columns
            byrow = True,           # fill in C order
            dimnames = dimnames     # supply axis labels
        )
    elif a_flags['F_CONTIGUOUS']:   # if memory is Fortran contiguous
        R_res = Rfn(                # call 'matrix' function from R
            R_vec,                  # input R vector
            nrow = nrow,            # number of rows
            ncol = ncol,            # number of columns
            byrow = False,          # fill in Fortran order
            dimnames = dimnames     # supply axis labels
        )
    else:                           # otherwise raise error
        raise ValueError("input matrix is neither C nor Fortran contiguous")
    return R_res

################################################################################
############################# Factor Type Objects ##############################
################################################################################
def rpy2_to_R_FactorVector(a):
    return baseenv['as.factor'](a)

################################################################################
############################# Boolean Type Objects #############################
################################################################################
def numpy_to_R_BoolVector(a):
    return generic_numpy_to_R_Vector(a, "int32", BoolVector, 'as.logical')

def numpy_to_R_BoolFactorVector(a):
    return rpy2_to_R_FactorVector(numpy_to_R_BoolVector(a))

def numpy_to_R_BoolMatrix(a, dimnames = NULL):
    return generic_numpy_to_R_Matrix(a, "int32", BoolVector, dimnames)

################################################################################
######################### Floating-Point Type Objects ##########################
################################################################################
def numpy_to_R_FloatVector(a):
    return generic_numpy_to_R_Vector(a, "float64", FloatVector, 'as.double')

def numpy_to_R_FloatFactorVector(a):
    return rpy2_to_R_FactorVector(numpy_to_R_FloatVector(a))

def numpy_to_R_FloatMatrix(a, dimnames = NULL):
    return generic_numpy_to_R_Matrix(a, "float64", FloatVector, dimnames)

################################################################################
############################# Integer Type Objects #############################
################################################################################
def numpy_to_R_IntVector(a):
    return generic_numpy_to_R_Vector(a, "int32", IntVector, 'as.integer')

def numpy_to_R_IntFactorVector(a):
    return rpy2_to_R_FactorVector(numpy_to_R_IntVector(a))

def numpy_to_R_IntMatrix(a, dimnames = NULL):
    return generic_numpy_to_R_Matrix(a, "int32", IntVector, dimnames)

################################################################################
############################# String Type Objects ##############################
################################################################################
def numpy_to_R_StrVector(a):
    a_ravel = numpy.ravel(          # ravel/flatten array (no memory allocation)
        a,                          # input array
        order = 'K'                 # keep ordering of input array
    )
    R_res = StrVector(a_ravel)      # convert to string vector
    return R_res

def numpy_to_R_StrFactorVector(a):
    return rpy2_to_R_FactorVector(numpy_to_R_StrVector(a))

def numpy_to_R_StrMatrix(a, dimnames = NULL):
    nrow, ncol = a.shape            # get number of rows and columns
    a_ravel = numpy.ravel(          # ravel/flatten array (no memory allocation)
        a,                          # input array
        order = 'K'                 # keep ordering of input array
    )
    R_vec = StrVector(a_ravel)      # convert to string vector
    Rfn = baseenv['matrix']         # get 'matrix' function from R
    a_flags = a.flags               # get array flags
    if a_flags['C_CONTIGUOUS']:     # if memory is C contiguous
        R_res = Rfn(                # call 'matrix' function from R
            R_vec,                  # input R vector
            nrow = nrow,            # number of rows
            ncol = ncol,            # number of columns
            byrow = True,           # fill in C order
            dimnames = dimnames     # supply axis labels
        )
    elif a_flags['F_CONTIGUOUS']:   # if memory is Fortran contiguous
        R_res = Rfn(                # call 'matrix' function from R
            R_vec,                  # input R vector
            nrow = nrow,            # number of rows
            ncol = ncol,            # number of columns
            byrow = False,          # fill in Fortran order
            dimnames = dimnames     # supply axis labels
        )
    else:                           # otherwise raise error
        raise ValueError("input matrix is neither C nor Fortran contiguous")
    return R_res

################################################################################
############################# General Type Objects #############################
################################################################################

# dictionary to get vector conversion functions
_kind_numpy_to_R_Vector = {
    'b' : numpy_to_R_BoolVector,
    'i' : numpy_to_R_IntVector,
    'u' : numpy_to_R_IntVector,
    'f' : numpy_to_R_FloatVector,
    # 'c' : None, # TODO: fix conversion
    # 'm' : None, # TODO: fix conversion
    # 'M' : None, # TODO: fix conversion
    'O' : numpy_to_R_StrVector,
    # 'S' : None, # TODO: fix conversion
    # 'U' : None, # TODO: fix conversion
    # 'V' : None  # TODO: fix conversion
}

# dictionary to get matrix conversion functions
_kind_numpy_to_R_Matrix = {
    'b' : numpy_to_R_BoolMatrix,
    'i' : numpy_to_R_IntMatrix,
    'u' : numpy_to_R_IntMatrix,
    'f' : numpy_to_R_FloatMatrix,
    # 'c' : None, # TODO: fix conversion
    # 'm' : None, # TODO: fix conversion
    # 'M' : None, # TODO: fix conversion
    'O' : numpy_to_R_StrMatrix,
    # 'S' : None, # TODO: fix conversion
    # 'U' : None, # TODO: fix conversion
    # 'V' : None  # TODO: fix conversion
}

def numpy_to_R(a, Rtype = None, factor = False):
    """
    Convert a numpy.ndarray to corresponding R object (vector, matrix, array)

    Parameters
    ----------
    a : numpy.ndarray
        Input array.
    Rtype : str, default = None
        R vector data type to coerce. If None, infer from input numpy.ndarray.
        Options (correspond to numpy.dtype kinds):
        Code    | Description
        --------+----------------------------------
        'b'     | Coerce to a logical R data type
        'i'     | Coerce to an integer R data type
        'u'     | Coerce to an integer R data type
        'f'     | Coerce to a numeric R data type
        'c'     | UNDER CONSTRUCTION
        'm'     | UNDER CONSTRUCTION
        'M'     | UNDER CONSTRUCTION
        'O'     | Coerce to a character R data type
        'S'     | UNDER CONSTRUCTION
        'U'     | UNDER CONSTRUCTION
        'V'     | UNDER CONSTRUCTION
    factor : boolean, default = False
        Whether to convert vector to a factor.

    """
    ##############################
    ### argument processing
    if Rtype is None:
        Rtype = a.dtype.kind

    ##############################
    ### get conversion dictionary
    if a.ndim == 1:
        conv_dict = _kind_numpy_to_R_Vector
    elif a.ndim == 2:
        conv_dict = _kind_numpy_to_R_Matrix
    else:
        raise ValueError("numpy.ndarray of ndim == {0} not supported".format(ndim))

    ##############################
    ### get conversion function
    try:
        convfn = conv_dict[Rtype]
    except:
        raise TypeError("conversion routine for numpy.ndarray kind '{0}' not supported".format(kind))

    ##############################
    ### convert array
    R_res = convfn(a)
    if factor:
        R_res = rpy2_to_R_FactorVector(R_res)

    return R_res


# setup = """
# import numpy
# from rpy2.robjects.vectors import BoolVector
# a = numpy.random.randint(0,2,10000)
# # fastest thing developed
# def fn1(a):
#     if a.dtype != 'int32':
#         a = numpy.int32(a)
#     out = BoolVector.from_memoryview(a.data)
#     return out
# # extremely slow!!!
# def fn2(a):
#     out = BoolVector(a)
#     return out
# def fn3(a):
#     if a.dtype != 'int32':
#         a = numpy.int32(a)
#     out = BoolVector.from_iterable(a.data)
#     return out
# """
# timeit.timeit(stmt = "fn1(a)", setup = setup, number = 10000)
# timeit.timeit(stmt = "fn2(a)", setup = setup, number = 100)
# timeit.timeit(stmt = "fn3(a)", setup = setup, number = 100)
