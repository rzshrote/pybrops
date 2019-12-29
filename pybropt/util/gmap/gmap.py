import math
import numpy
from itertools import chain
import haldane
import kosambi
# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import human2bytes


def gmap2r(gmap, lgroup, mapfn = None, mtype = None, dtype = None):
    """
    Convert genetic map positions to a matrix of recombination probabilities
    using a provided mapping function.

    Parameters
    ----------
    gmap : numpy.ndarray
        A 1D array of genetic map positions to convert to recombination
        probabilities. The positions units should be in Morgans.
    lgroup : numpy.ndarray
        A 1D array of linkage group sizes. The sum of the elements should equal
        the length of 'gmap'.
    mapfn : None, callable, {'haldane', 'kosambi'}
        Mapping function to use. If None is provided, default to 'haldane'.
        Options:
            callable    : A custom callable Python function.
                          Must take the following arguments:
                            mapfn(d)
                            Parameters:
                                d : numpy.ndarray
                                    An array of genetic distances in Morgans.
                                    This can be an array of any shape.
                            Returns:
                                r : numpy.ndarray
                                    An array of recombination probabilities.
                                    The shape of the array is the same shape
                                    as that of 'd'.
            'haldane'   : Haldane (1919) mapping function
            'kosambi'   : Kosambi (1944) mapping function
    mtype : None, {'sym', 'triu', 'tril'}
        Matrix type to return. If None is provided, default to 'sym'.
        Options:
            'sym'   : Symmetric matrix
            'triu'  : Upper triangular matrix
            'tril'  : Lower triangular matrix
    dtype : None, numpy.dtype
        The type of the output array. If dtype is None, infer the data
        type from the 'gmap' argument.

    Returns
    -------
    r : numpy.ndarray
        An array of recombination probabilities between loci as defined by 'd'.
    """
    ######################## Parse 'gmap' and 'lgroup' ########################
    # check to make sure 'gmap' is the right shape
    if len(gmap.shape) != 1:
        raise ValueError(
            "Incorrect 'gmap' dimension. Needs to be 1D.\n"\
            "    len(gmap.shape) = %s\n"\
            "    gmap.shape = %s" % (len(gmap.shape), gmap.shape)
        )
    # check to make sure 'lgroup' is the right shape
    if len(lgroup.shape) != 1:
        raise ValueError(
            "Incorrect 'lgroup' dimension. Needs to be 1D.\n"\
            "    len(lgroup.shape) = %s\n"\
            "    lgroup.shape = %s" % (len(lgroup.shape), lgroup.shape)
        )
    # check to make sure gmap and lgroup align
    if len(gmap) != lgroup.sum():
        raise ValueError(
            "The 'gmap' and 'lgroup' do not align:\n"\
            "    len(gmap)    = %s\n"\
            "    lgroup.sum() = %s\n"
            % (len(gmap), lgroup.sum())
        )

    ############################## Parse 'mapfn' ##############################
    # if mapfn is None, set it to 'kosambi'
    if mapfn is None:
        mapfn = 'haldane'
    # process mapfn which may be a function or a string
    if callable(mapfn):
        pass
    elif isinstance(mapfn, str):
        if mapfn == 'haldane':
            mapfn = haldane2r
        elif mapfn == 'kosambi':
            mapfn = kosambi2r
        else:
            raise ValueError(
                "Mapping function '%s' not supported.\n"\
                "    Options:\n"\
                "        'haldane'   : Haldane (1919) mapping function\n"\
                "        'kosambi'   : Kosambi (1944) mapping function\n"
                % mapfn
            )
    else:
        raise TypeError(
            "'mapfn' not correct type."\
            "'mapfn' must be a callable() function or a string.\n"\
            "    Received: type(mapfn) = %s\n"\
            "    Options:\n"\
            "        callable    : A custom callable Python function.\n"\
            "        'haldane'   : Haldane (1919) mapping function\n"\
            "        'kosambi'   : Kosambi (1944) mapping function\n"
            % type(mapfn)
        )

    ############################## Parse 'dtype' ##############################
    # if there is no provided dtype, set it to the genetic map dtype
    if dtype is None:
        dtype = gmap.dtype
    # if dtype is not of type numpy.dtype, raise an error
    if not isinstance(dtype, numpy.dtype):
        raise TypeError(
            "'dtype' must be of type 'numpy.dtype'.\n"\
            "    Received: type(dtype) = %s" % type(dtype)
        )

    ############################## Parse 'mtype' ##############################
    # if 'mtype' is None, set to compute entire matrix in blocks
    if mtype is None:
        mtype = 'sym'

    # if the mtype is not in the list of options, raise an error
    if mtype not in ['sym', 'triu', 'tril']:
        raise ValueError(
            "Invalid 'mtype'. Options are:\n"\
            "    'sym'   : Symmetric matrix\n"\
            "    'triu'  : Upper triangular matrix\n"\
            "    'tril'  : Lower triangular matrix\n"
        )

    ########################### Compute whole matrix ###########################
    # allocate an empty matrix for recombination probabilities
    r = numpy.empty(            # make empty matrix
        (len(gmap), len(gmap)), # make a square matrix of width len(gmap)
        dtype=dtype             # set the data type as dtype
    )

    # calculate start, stop indices for matrix navigation
    dsp = numpy.cumsum(lgroup)     # calculate stop indices
    dst = dsp - lgroup[0]          # calculate start indices

    # for each start, stop index
    for st,sp in zip(dst,dsp):
        # make mesh for that array region
        tmpmesh = numpy.meshgrid(d[st:sp], d[st:sp])

        # calculate recombination for that region
        r[st:sp,st:sp] = mapfn(numpy.abs(tmpmesh[0] - tmpmesh[1]))

        # fill everything else assuming independent assortment
        r[:st,st:sp].fill(0.5)  # fill matrix above with 0.5
        r[sp:,st:sp].fill(0.5)  # fill matrix below with 0.5

    # finally, modify our 'r' matrix if we need a specific mtype
    if mtype == 'sym':              # if 'sym', do nothing
        pass
    elif mtype == 'triu':           # if 'triu', multiply by a triu mask
        r *= numpy.logical_not(     # invert triangle mask for triu (fast)
            numpy.tri(              # make inverse of triu mask
                *r.shape[-2:],      # get first two dimensions
                k=-1,               # offset to get lower triangle
                dtype=numpy.bool    # use bool (1 byte) data type (fast)
            )
        )
    elif mtype == 'tril':           # if 'tril', multiply by a tril mask
        r *= numpy.tri(             # make lower triangle mask
            *r.shape[-2:],          # get first two dimensions
            k=0,                    # offset of zero to get lower triangle mask
            dtype=numpy.bool        # use bool (1 byte) data type (fast)
        )

    # return probability matrix
    return r


def gmap2r_gen(gmap, lgroup, mapfn = None, mem = None, mtype = None,
               dtype = None):
    """
    Convert genetic map positions in Morgans to recombination probabilities
    using the Haldane mapping function (Haldane, 1919).

    This is a generator function, so it will iteratively yield different
    chunks of the matrix based on memory chunk size. Most chunks will be
    square if the length of gmap is larger than the width of the chunk.

    Parameters
    ==========
    gmap : numpy.ndarray
        A 1D array of genetic map positions to convert to recombination
        probabilities. The positions units should be in Morgans.
    lgroup : numpy.ndarray
        A 1D array of linkage group sizes. The sum of the elements should equal
        the length of 'gmap'.
    mapfn : None, callable, {'haldane', 'kosambi'}
        Mapping function to use. If None is provided, default to 'haldane'.
        Options:
            callable    : A custom callable Python function.
                          Must take the following arguments:
                            mapfn(d)
                            Parameters:
                                d : numpy.ndarray
                                    An array of genetic distances in Morgans.
                                    This can be an array of any shape.
                            Returns:
                                r : numpy.ndarray
                                    An array of recombination probabilities.
                                    The shape of the array is the same shape
                                    as that of 'd'.
            'haldane'   : Haldane (1919) mapping function
            'kosambi'   : Kosambi (1944) mapping function
    mem : None, {int, numpy.integer, str}
        Size of memory chunk to compute at a time. This has several options.
        Options:
            None:
                Compute everything at once; no memory restrictions.
            int, numpy.int:
                Specify chunks of maximum width 'mem'. Memory is restricted,
                but is proportional to dtype.itemsize * (mem**2)
            str:
                String specifying amount of memory to use for each chunk.
                This is a string that is parsed. Can use SI or IEC format.
                Examples:
                    mem = "1M"  # specifies 1 megabyte chunks
                    mem = "1Mi" # specifies 1 mebibits chunks
    mtype : None, {'sym', 'triu', 'tril'}
        Matrix type to return. If None is provided, default to 'sym'
        Options:
            'sym'   : Symmetric matrix
            'triu'  : Upper triangular matrix
            'tril'  : Lower triangular matrix
    dtype : numpy.dtype
        The type of the output array. If dtype is not given, infer the data
        type from the 'gmap' argument.

    Returns
    -------
    r : numpy.ndarray
        An array of recombination probabilities between loci as defined by 'd'.
    """
    ######################## Parse 'gmap' and 'lgroup' ########################
    # check to make sure 'gmap' is the right shape
    if len(gmap.shape) != 1:
        raise ValueError(
            "Incorrect 'gmap' dimension. Needs to be 1D.\n"\
            "    len(gmap.shape) = %s\n"\
            "    gmap.shape = %s" % (len(gmap.shape), gmap.shape)
        )
    # check to make sure 'lgroup' is the right shape
    if len(lgroup.shape) != 1:
        raise ValueError(
            "Incorrect 'lgroup' dimension. Needs to be 1D.\n"\
            "    len(lgroup.shape) = %s\n"\
            "    lgroup.shape = %s" % (len(lgroup.shape), lgroup.shape)
        )
    # check to make sure gmap and lgroup align
    if len(gmap) != lgroup.sum():
        raise ValueError(
            "The 'gmap' and 'lgroup' do not align:\n"\
            "    len(gmap)    = %s\n"\
            "    lgroup.sum() = %s\n"
            % (len(gmap), lgroup.sum())
        )

    ############################## Parse 'dtype' ##############################
    # if there is no provided dtype, set it to the genetic map dtype
    if dtype is None:
        dtype = gmap.dtype
    # if dtype is not of type numpy.dtype, raise an error
    if not isinstance(dtype, numpy.dtype):
        raise TypeError(
            "'dtype' must be of type 'numpy.dtype'.\n"\
            "    Received: type(dtype) = %s" % type(dtype)
        )

    ############################## Parse 'mapfn' ##############################
    # if mapfn is None, set it to 'kosambi'
    if mapfn is None:
        mapfn = 'haldane'
    # process mapfn which may be a function or a string
    if callable(mapfn):
        pass
    elif isinstance(mapfn, str):
        if mapfn == 'haldane':
            mapfn = haldane2r
        elif mapfn == 'kosambi':
            mapfn = kosambi2r
        else:
            raise ValueError(
                "Mapping function '%s' not supported.\n"\
                "    Options:\n"\
                "        'haldane'   : Haldane (1919) mapping function\n"\
                "        'kosambi'   : Kosambi (1944) mapping function\n"
                % mapfn
            )
    else:
        raise TypeError(
            "'mapfn' not correct type."\
            "'mapfn' must be a callable() function or a string.\n"\
            "    Received: type(mapfn) = %s\n"\
            "    Options:\n"\
            "        callable    : A custom callable Python function.\n"\
            "        'haldane'   : Haldane (1919) mapping function\n"\
            "        'kosambi'   : Kosambi (1944) mapping function\n"
            % type(mapfn)
        )

    ############################### Parse 'mem' ###############################
    # if 'mem' is None, set to compute whole matrix
    if mem is None:
        mem = len(gmap)
    # if 'mem' is an integer type, do nothing
    elif isinstance(mem, (int, numpy.int)):
        pass
    # if 'mem' is a string, parse it and calculate chunk size
    elif isinstance(mem, str):
        mem = int(                              # cast as integer
            math.floor(                         # get floor
                math.sqrt(                      # take square root
                    human2bytes(mem) /          # parse number of required bytes
                    numpy.dtype(dtype).itemsize # divide by dtype size in bytes
                )
            )
        )
    else:
        raise TypeError(
            "Incorrect 'mem' data type.\n"\
            "    Options: None, int, numpy.int, str\n"\
            "    type(mem) = %s" % type(mem)
        )
    # if mem is =< 0, raise an error
    if mem =< 0:
        raise ValueError(
            "'mem' must be greater than zero:\n"\
            "    Received: %s" % mem
        )

    ############################## Parse 'mtype' ##############################
    # if 'mtype' is None, set to compute entire matrix in blocks
    if mtype is None:
        mtype = 'sym'
    # invalid 'mtype' options handled in compute matrix chunks if-else chain

    ########################## Compute matrix chunks ##########################
    # get length of gmap
    l = len(gmap)

    # take cumulative sum to calculate stop indices
    lsp_arr = numpy.cumsum(lgroup)

    # if 'mtype' equals 'sym', we're doing the entire matrix
    if mtype == 'sym':
        # row loop
        for rst,rsp in zip(             # zip start, stop indices
            range(0, l, mem),           # rst) range start indices
            chain(                      # rsp) chain stop indices
                range(mem, l, mem),     # 1) first several indices
                (l,)                    # 2) last index
            )
        ):
            # column loop
            for cst,csp in zip(             # zip start, stop indices
                range(0, l, mem),           # cst) start indices
                chain(                      # csp) chain stop indices
                    range(mem, l, mem),     # 1) first several indices
                    (l,)                    # 2) last index
                )
            ):
                # Allocate an empty matrix containing random bit garbage.
                out = numpy.empty(
                    (rsp-rst, csp-cst),     # (rows, columns) dimensions
                    dtype = dtype           # set data type
                )

                # variables to track columns we haven't filled yet
                lcols = csp  # cols to the left not filled (index)
                rcols = cst  # cols to the right not filled (index)

                # Algorithm:
                #     For each linkage group, calculate the overlap with the
                #     current computational chunk. If there is an overlap,
                #     compute the overlap. Otherwise, fill with 0.5.

                # start coordinate variable for the linkage groups.
                lst = 0

                # cycle through all the linkage group sectors trying to identify
                # regions of overlap between the computational chunk and the
                # linkage group block
                for lsp in lsp_arr:
                    # calculate an overlap region defining a box: (x1,y1,x2,y2)
                    # where x1 < x2, y1 < y2
                    x1 = max(lst, cst)   # x1; matrix start column
                    y1 = max(lst, rst)   # y1; matrix start row
                    x2 = min(lsp, csp)   # x2; matrix stop column
                    y2 = min(lsp, rsp)   # y2; matrix start row

                    # calclulate the differences between x2-x1, y2-y1
                    dx = x2 - x1  # x overlap; if negative, no overlap
                    dy = y2 - y1  # y overlap; if negative, no overlap

                    # if we've found a linkage group that overlaps
                    if (dx >= 0) and (dy >= 0):
                        # get the cols overlap (x range) as a slice
                        cols = slice(x1-cst, x2-cst)

                        # make distance mesh for the overlap region
                        # first argument specifies number of columns
                        # second argument specifies number of rows
                        dmesh = numpy.meshgrid(gmap[x1:x2], gmap[y1:y2])

                        # fill the recombination frequencies in the matrix
                        #   rst-y1:rst-y2     # rows (y)
                        #   cols              # columns (x)
                        out[y1-rst:y2-rst,cols] = mapfn(
                            numpy.abs(dmesh[0] - dmesh[1])
                        )

                        # fill rows above the linkage block with 0.5
                        out[:y1-rst,cols].fill(0.5)

                        # fill rows below the linkage block with 0.5
                        out[y2-rst:,cols].fill(0.5)

                        # alter column tracking variables
                        if x1 < lcols:  # if x1 is before lcols index
                            lcols = x1  # set the index to x1
                        if x2 > rcols:  # if x2 is after rcols index
                            rcols = x2  # set the index to x2

                    # advance the linkage group start index to the current stop
                    lst = lsp

                # fill the remaining untouched columns with 0.5
                if lcols > rcols:               # if lcols and rcols overlap
                    out.fill(0.5)               # fill the whole thing in one go
                else:                           # else fill each side
                    out[:,:lcols-cst].fill(0.5) # fill left
                    out[:,rcols-cst:].fill(0.5) # fill right

                # yield the output matrix
                yield out
    elif mtype == 'triu':
        # row loop
        for rst,rsp in zip(             # zip start, stop indices
            range(0, l, mem),           # rst) range start indices
            chain(                      # rsp) chain stop indices
                range(mem, l, mem),     # 1) first several indices
                (l,)                    # 2) last index
            )
        ):
            # column loop
            for cst,csp in zip(             # zip start, stop indices
                range(rst, l, mem),         # cst) start indices
                chain(                      # csp) chain stop indices
                    range(rst+mem, l, mem), # 1) first several indices
                    (l,)                    # 2) last index
                )
            ):
                # Allocate an empty matrix containing random bit garbage.
                out = numpy.empty(
                    (rsp-rst, csp-cst),     # (rows, columns) dimensions
                    dtype = dtype           # set data type
                )

                # variables to track columns we haven't filled yet
                lcols = csp  # cols to the left not filled (index)
                rcols = cst  # cols to the right not filled (index)

                # Algorithm:
                #     For each linkage group, calculate the overlap with the
                #     current computational chunk. If there is an overlap,
                #     compute the overlap. Otherwise, fill with 0.5.

                # start coordinate variable for the linkage groups.
                lst = 0

                # cycle through all the linkage group sectors trying to identify
                # regions of overlap between the computational chunk and the
                # linkage group block
                for lsp in lsp_arr:
                    # calculate an overlap region defining a box: (x1,y1,x2,y2)
                    # where x1 < x2, y1 < y2
                    x1 = max(lst, cst)   # x1; matrix start column
                    y1 = max(lst, rst)   # y1; matrix start row
                    x2 = min(lsp, csp)   # x2; matrix stop column
                    y2 = min(lsp, rsp)   # y2; matrix start row

                    # calclulate the differences between x2-x1, y2-y1
                    dx = x2 - x1  # x overlap; if negative, no overlap
                    dy = y2 - y1  # y overlap; if negative, no overlap

                    # if we've found a linkage group that overlaps
                    if (dx > 0) and (dy > 0):
                        # get the cols overlap (x range) as a slice
                        cols = slice(x1-cst, x2-cst)

                        # make distance mesh for the overlap region
                        # first argument specifies number of columns
                        # second argument specifies number of rows
                        dmesh = numpy.meshgrid(gmap[x1:x2], gmap[y1:y2])

                        # fill the recombination frequencies in the matrix
                        out[y1-rst:y2-rst,cols] = mapfn(
                            numpy.abs(dmesh[0] - dmesh[1])
                        )

                        # fill rows above the linkage block with 0.5
                        out[:y1-rst,cols].fill(0.5)

                        # fill rows below the linkage block with 0.5
                        out[y2-rst:,cols].fill(0.5)

                        # alter column tracking variables
                        if x1 < lcols:  # if x1 is before lcols index
                            lcols = x1  # set the index to x1
                        if x2 > rcols:  # if x2 is after rcols index
                            rcols = x2  # set the index to x2

                    # advance the linkage group start index to the current stop
                    lst = lsp

                # fill the remaining untouched columns with 0.5
                if lcols > rcols:               # if lcols and rcols overlap
                    out.fill(0.5)               # fill the whole thing in one go
                else:                           # else fill each side
                    out[:,:lcols-cst].fill(0.5) # fill left
                    out[:,rcols-cst:].fill(0.5) # fill right

                # if we are along the diagonal
                if rst == cst:
                    # multiply in place by an upper triangle mask
                    out *= numpy.logical_not(   # invert triangle mask for triu
                        numpy.tri(              # make inverse of triu mask
                            *out.shape[-2:],    # get first two dimensions
                            k=-1,               # offset to get lower triangle
                            dtype=numpy.bool    # use bool data type
                        )
                    )

                # yield the output matrix
                yield out
    elif mtype == 'tril':
        # row loop
        for rst,rsp in zip(             # zip start, stop indices
            range(0, l, mem),           # rst) range start indices
            chain(                      # rsp) chain stop indices
                range(mem, l, mem),     # 1) first several indices
                (l,)                    # 2) last index
            )
        ):
            # column loop
            for cst,csp in zip(             # zip start, stop indices
                range(0, rsp, mem),         # cst) start indices
                chain(                      # csp) chain stop indices
                    range(mem, rsp, mem),   # 1) first several indices
                    (rsp,)                  # 2) last index
                )
            ):
                # Allocate an empty matrix containing random bit garbage.
                out = numpy.empty(
                    (rsp-rst, csp-cst),     # (rows, columns) dimensions
                    dtype = dtype           # set data type
                )

                # variables to track columns we haven't filled yet
                lcols = csp  # cols to the left not filled (index)
                rcols = cst  # cols to the right not filled (index)

                # Algorithm:
                #     For each linkage group, calculate the overlap with the
                #     current computational chunk. If there is an overlap,
                #     compute the overlap. Otherwise, fill with 0.5.

                # start coordinate variable for the linkage groups.
                lst = 0

                # cycle through all the linkage group sectors trying to identify
                # regions of overlap between the computational chunk and the
                # linkage group block
                for lsp in lsp_arr:
                    # calculate an overlap region defining a box: (x1,y1,x2,y2)
                    # where x1 < x2, y1 < y2
                    x1 = max(lst, cst)   # x1; matrix start column
                    y1 = max(lst, rst)   # y1; matrix start row
                    x2 = min(lsp, csp)   # x2; matrix stop column
                    y2 = min(lsp, rsp)   # y2; matrix start row

                    # calclulate the differences between x2-x1, y2-y1
                    dx = x2 - x1  # x overlap; if negative, no overlap
                    dy = y2 - y1  # y overlap; if negative, no overlap

                    # if we've found a linkage group that overlaps
                    if (dx > 0) and (dy > 0):
                        # get the cols overlap (x range) as a slice
                        cols = slice(x1-cst, x2-cst)

                        # make distance mesh for the overlap region
                        # first argument specifies number of columns
                        # second argument specifies number of rows
                        dmesh = numpy.meshgrid(gmap[x1:x2], gmap[y1:y2])

                        # fill the recombination frequencies in the matrix
                        out[y1-rst:y2-rst,cols] = mapfn(
                            numpy.abs(dmesh[0] - dmesh[1])
                        )

                        # fill rows above the linkage block with 0.5
                        out[:y1-rst,cols].fill(0.5)

                        # fill rows below the linkage block with 0.5
                        out[y2-rst:,cols].fill(0.5)

                        # alter column tracking variables
                        if x1 < lcols:  # if x1 is before lcols index
                            lcols = x1  # set the index to x1
                        if x2 > rcols:  # if x2 is after rcols index
                            rcols = x2  # set the index to x2

                    # advance the linkage group start index to the current stop
                    lst = lsp

                # fill the remaining untouched columns with 0.5
                if lcols > rcols:               # if lcols and rcols overlap
                    out.fill(0.5)               # fill the whole thing in one go
                else:                           # else fill each side
                    out[:,:lcols-cst].fill(0.5) # fill left
                    out[:,rcols-cst:].fill(0.5) # fill right

                # if we are along the diagonal
                if rst == cst:
                    # multiply in place by an lower triangle mask
                    out *= numpy.tri(       # make tril mask
                        *out.shape[-2:],    # get first two dimensions
                        k=0,                # offset to get lower triangle
                        dtype=numpy.bool    # use bool data type
                    )

                # yield the output matrix
                yield out
    else:
        raise ValueError(
            "Invalid 'mtype' value '%s'.\n"\
            "    Options: None, {'sym', 'triu', 'tril'}" % mtype
        )
