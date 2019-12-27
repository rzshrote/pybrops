import math
import numpy
# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import human2bytes

def haldane2r_gen(gmap, lgroup, mem=None, mtype=None, dtype=None):
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
        the length of 'd'.
    mem : None, int, numpy.int, str
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
        type from the 'd' argument.

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
        dtype = d.dtype

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

    ############################## Parse 'mtype' ##############################
    # if 'mtype' is None, set to compute entire matrix in blocks
    if mtype is None:
        mtype = 'sym'
    # if 'mtype' is not in valid options, throw an error
    if mtype not in ['sym', 'triu', 'tril']:
        raise ValueError(
            "Invalid 'mtype' value '%s'.\n"\
            "    Options: None, {'sym', 'triu', 'tril'}" % mtype
        )

    ########################## Compute matrix chunks ##########################
    # get length of gmap
    l = len(gmap)

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
            for cst,csp in zip(         # zip start, stop indices
                range(0, l, mem),       # cst) start indices
                chain(                  # csp) chain stop indices
                    range(mem, l, mem), # 1) first several indices
                    (l,)                # 2) last index
                )
            ):
                # Step 1) Allocate an empty matrix containing random bit garbage.
                #         The dimensions correspond to start, stop positions
                out = numpy.empty(
                    shape = (rsp-rst, csp-cst),     # (rows, columns) dimensions
                    dtype = dtype                   # set data type
                )

                # Step 2) make a counter to track current column we've filled
                #         This starts at 'cst' and advances with each processed
                #         linkage group.
                curcol = cst

                # Pseudocode:
                # declare start variable
                # st = 0
                # for sp in lsp_arr:
                #     overlap lgroup_box{ (st,st), (sp,sp) } and chunk_box{ (cst, rst), (csp, rsp) }
                #     if any are greater than or equal to 0 there is an overlap

                # Step 3) declare start coordinate variable for the linkage groups.
                # note: this is independent of the starting sector, so all linkage
                #       groups are tested.
                lst = 0

                # Step 4) cycle through all the linkage group sectors trying to
                # identify regions of overlap between the computational chunk and
                # the linkage group block
                for lsp in lsp_arr:
                    # calculate an overlap region defining a box: (x1,y1,x2,y2)
                    # where x1 < x2, y1 < y2
                    x1 = max(lst, cst),  # x1
                    y1 = max(lst, rst),  # y1
                    x2 = min(lsp, csp),  # x2
                    y2 = min(lsp, rsp)   # y2

                    # calclulate the differences between x2-x1, y2-y1
                    dx = x2 - x1  # x overlap; if negative, no overlap
                    dy = y2 - y1  # y overlap; if negative, no overlap

                    # if we've found a linkage group that overlaps
                    if (dx >= 0) and (dy >= 0):
                        # get the cols overlap (x range) as a slice; this slice will
                        # be reused multiple times, so we want to save time for
                        # allocation
                        cslice = slice(cst - x1, cst - x2)

                        # make mesh for the overlap region
                        # first argument specifies number of columns
                        # second argument specifies number of rows
                        tmpmesh = numpy.meshgrid(d[x1:x2], d[y1:y2])

                        # fill the recombination frequencies in the matrix
                        out[
                            slice(rst - y1, rst - y2),  # rows (y)
                            cslice                      # cols (x)
                        ] = 0.5 * (
                            1 - numpy.exp(-2 * numpy.abs(tmpmesh[0] - tmpmesh[1]))
                        )

                        # fill rows above the linkage block with 0.5
                        out[
                            slice(0, rst - y1),     # rows above linkage group block
                            cslice                  # cols (x)
                        ].fill(0.5)

                        # fill rows below the linkage block with 0.5
                        out[
                            slice(rst - y2, rsp),   # rows above linkage group block
                            cslice                  # cols (x)
                        ].fill(0.5)

                        # advance the curcol variable to the next location
                        curcol = x2

                    # advance the linkage group start index to the current stop
                    lst = lsp

                # Step 5) fill the rest of the remaining matrix with 0.5, since it
                # does not intersect with a linkage block.
                # linkage group blocks should always be touching at their
                # vertices, so advancing curcol will not skip empty sectors.
                out[:,curcol:].fill(0.5)

                # Step 6) if the computational chunk is along the diagonal
                # (rst == cst since we are dealing with square matrices), get the
                # upper triangle.
                if rst == cst:
                    # multiply in place by an upper triangle mask
                    # performance notes: data conversion from bool to float will
                    # reduce speed. May need to do some testing to determine if
                    # doing 1 - tri(dtype=dtype) is faster
                    out *= numpy.logical_not(   # invert triangle mask for triu (fast)
                        numpy.tri(              # make inverse of triu mask
                            *out.shape[-2:],    # get first two dimensions
                            k=-1,               # offset to get lower triangle
                            dtype=numpy.bool    # use bool (1 byte) data type (fast)
                        )
                    )

                # yield the output matrix
                yield out
    elif

    # take cumulative sum to calculate stop indices
    # calculate start indices from stop indices
    lsp_arr = numpy.cumsum(lgroup)

    yield 0

def haldane_to_r_triu(d, lgroup, sector = None, chunk = None, dtype = None):
    """
    Convert genetic map distances in Morgans to recombination probabilities
    using the Haldane mapping function (Haldane, 1919).

    This will only calculate an upper triangle matrix for the provided sector.

    Parameters
    ==========
    d : numpy.ndarray
        A 1D array of genetic map positions. Data type should be floating point.
    lgroup : numpy.ndarray
        A 1D array of linkage group sizes. The sum of the elements should equal
        the length of 'd'. Data type should be integral.
    sector : None, tuple, array-like
        A 2-mer tuple or 2-mer array-like representing the region of the
        recombination matrix to calculate. The format of this tuple is:
            (x1, x2)
            Where:
                x1 is an integer index *Inclusive*
                x2 is an integer index *Exclusive*
                x1 value corresponds to index start in a 2D numpy.ndarray.
                x2 value corresponds to index stop in a 2D numpy.ndarray.
                This will always specify along the diagonal

          +---------------------------------> x (index increases to right)
          |   rest of matrix
          |
          |   (x1,x1)+-------------------+
          |          |                   |
          |          |      sector       |
          |          |                   |
          |          +-------------------+(x2,x2)
          |
          |
          v
          y (index increases downwards)

    chunk : None, integer, numpy.int
        Length of elements in 'sector' to process at once. This will cut up a
        large sector into a number of mostly square, but also rectangular
        computational chunks. This can be used to limit memory usage.

    Yields
    ======
    yield : numpy.ndarray
        A numpy array of recombination probabilities.
    """
    # if sector is None, compute the entire matrix
    # format: (x1, x2)
    if sector is None:
        sector = (0, len(d))

    # if the chunk size is None, set chunk to compute whole matrix
    if chunk is None:
        chunk = sector[1]-sector[0]

    # if there is no provided dtype, set it to the genetic map dtype
    if dtype is None:
        dtype = d.dtype

    # take cumulative sum to calculate stop indices
    # calculate start indices from stop indices
    lsp_arr = numpy.cumsum(lgroup)


    # row loop
    for rst,rsp in zip(             # enumerate, zip start, stop indices
        range(                      # generator for start indices
            sector[0],              # start in sector.x1 (inclusive)
            sector[1],              # stop in x2 (exclusive)
            chunk[1]                # advance by chunk.y
        ),
        chain(                      # calculate stop indices
            range(                  # calculate first several indices
                sector[0]+chunk,    # start in sector.x1 + chunk (inclusive)
                sector[1],          # stop in sector.x2 (exclusive)
                chunk[1]            # advance by chunk
            ),
            (sector[3],)            # tack on last index
        )
    ):
        # column loop
        for cst,csp in zip(                     # enumerate, zip start, stop indices
            range(rst,stop,chunk),              # calculate start indices
            chain(                              # calculate stop indices
                range(rst+chunk,stop,chunk),    # calculate first several
                (stop,)                         # tack on last index
            )
        )):
            # Step 1) Allocate an empty matrix containing random bit garbage.
            #         The dimensions correspond to start, stop positions
            out = numpy.empty(
                shape = (rsp-rst, csp-cst),     # (rows, columns) dimensions
                dtype = dtype                   # set data type
            )

            # Step 2) make a counter to track current column we've filled
            #         This starts at 'cst' and advances with each processed
            #         linkage group.
            curcol = cst

            # Pseudocode:
            # declare start variable
            # st = 0
            # for sp in lsp_arr:
            #     overlap lgroup_box{ (st,st), (sp,sp) } and chunk_box{ (cst, rst), (csp, rsp) }
            #     if any are greater than or equal to 0 there is an overlap

            # Step 3) declare start coordinate variable for the linkage groups.
            # note: this is independent of the starting sector, so all linkage
            #       groups are tested.
            lst = 0

            # Step 4) cycle through all the linkage group sectors trying to
            # identify regions of overlap between the computational chunk and
            # the linkage group block
            for lsp in lsp_arr:
                # calculate an overlap region defining a box: (x1,y1,x2,y2)
                # where x1 < x2, y1 < y2
                x1 = max(lst, cst),  # x1
                y1 = max(lst, rst),  # y1
                x2 = min(lsp, csp),  # x2
                y2 = min(lsp, rsp)   # y2

                # calclulate the differences between x2-x1, y2-y1
                dx = x2 - x1  # x overlap; if negative, no overlap
                dy = y2 - y1  # y overlap; if negative, no overlap

                # if we've found a linkage group that overlaps
                if (dx >= 0) and (dy >= 0):
                    # get the cols overlap (x range) as a slice; this slice will
                    # be reused multiple times, so we want to save time for
                    # allocation
                    cslice = slice(cst - x1, cst - x2)

                    # make mesh for the overlap region
                    # first argument specifies number of columns
                    # second argument specifies number of rows
                    tmpmesh = numpy.meshgrid(d[x1:x2], d[y1:y2])

                    # fill the recombination frequencies in the matrix
                    out[
                        slice(rst - y1, rst - y2),  # rows (y)
                        cslice                      # cols (x)
                    ] = 0.5 * (
                        1 - numpy.exp(-2 * numpy.abs(tmpmesh[0] - tmpmesh[1]))
                    )

                    # fill rows above the linkage block with 0.5
                    out[
                        slice(0, rst - y1),     # rows above linkage group block
                        cslice                  # cols (x)
                    ].fill(0.5)

                    # fill rows below the linkage block with 0.5
                    out[
                        slice(rst - y2, rsp),   # rows above linkage group block
                        cslice                  # cols (x)
                    ].fill(0.5)

                    # advance the curcol variable to the next location
                    curcol = x2

                # advance the linkage group start index to the current stop
                lst = lsp

            # Step 5) fill the rest of the remaining matrix with 0.5, since it
            # does not intersect with a linkage block.
            # linkage group blocks should always be touching at their
            # vertices, so advancing curcol will not skip empty sectors.
            out[:,curcol:].fill(0.5)

            # Step 6) if the computational chunk is along the diagonal
            # (rst == cst since we are dealing with square matrices), get the
            # upper triangle.
            if rst == cst:
                # multiply in place by an upper triangle mask
                # performance notes: data conversion from bool to float will
                # reduce speed. May need to do some testing to determine if
                # doing 1 - tri(dtype=dtype) is faster
                out *= numpy.logical_not(   # invert triangle mask for triu (fast)
                    numpy.tri(              # make inverse of triu mask
                        *out.shape[-2:],    # get first two dimensions
                        k=-1,               # offset to get lower triangle
                        dtype=numpy.bool    # use bool (1 byte) data type (fast)
                    )
                )

            # yield the output matrix
            yield out


def haldane2r_fast(d):
    """
    Convert genetic distances in Morgans to recombination probabilities.

    This is a bare-bones function and does no parameter checking for errors.

    Parameters
    ==========
    d : numpy.ndarray
        An array of genetic distances in Morgans. This can be an array of any
        shape.

    Returns
    =======
    out : numpy.ndarray
        An array of recombination probabilities. The shape of the array is
        the same shape as that of 'd'.
    """
    return 0.5 * (1 - numpy.exp(-2 * d))

def haldane2r(gmap, lgroup, mtype = None, dtype = None):
    """
    Convert genetic map positions to a matrix of recombination probabilities
    using the Haldane mapping function (Haldane, 1919).

    Parameters
    ----------
    gmap : numpy.ndarray
        A 1D array of genetic map positions to convert to recombination
        probabilities. The positions units should be in Morgans.
    lgroup : numpy.ndarray
        A 1D array of linkage group sizes. The sum of the elements should equal
        the length of 'd'.
    mtype : None, {'sym', 'triu', 'tril'}
        Matrix type to return. If None is provided, default to 'sym'
        Options:
            'sym'   : Symmetric matrix
            'triu'  : Upper triangular matrix
            'tril'  : Lower triangular matrix
    dtype : dtype
        The type of the output array. If dtype is not given, infer the data
        type from the 'd' argument.

    Returns
    -------
    r : numpy.ndarray
        An array of recombination probabilities between loci as defined by 'd'.
    """
    # check if there is
    # if there is no provided mtype, set it to 'sym'
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

    # if there is no provided dtype, set it to the genetic map dtype
    if dtype is None:
        dtype = d.dtype

    # allocate an empty matrix for recombination probabilities
    r = numpy.empty(        # make empty matrix
        (len(d), len(d)),   # make a square matrix of size (r=len(d), c=len(d))
        dtype=dtype         # set the data type as dtype
    )

    # calculate start, stop indices for matrix navigation
    dsp = numpy.cumsum(lgroup)     # calculate stop indices
    dst = dsp - lgroup[0]          # calculate start indices

    # for each start, stop index
    for st,sp in zip(dst,dsp):
        # make mesh for that array region
        tmpmesh = numpy.meshgrid(d[st:sp], d[st:sp])

        # calculate Haldane distance for that region
        r[st:sp,st:sp] = 0.5 * (
            1 - numpy.exp(-2 * numpy.abs(tmpmesh[0] - tmpmesh[1]))
        )

        # fill everything else assuming independent assortment
        # NOTE: use of '=' vs 'fill': 'fill' looks to be slightly faster
        r[:st,st:sp].fill(0.5)  # fill matrix above with 0.5
        r[sp:,st:sp].fill(0.5)  # fill matrix below with 0.5

    # finally, modify our 'r' matrix if we need a specific mtype
    if mtype == 'sym':              # if 'sym', do nothing
        pass
    if mtype == 'triu':             # if 'triu', multiply by a triu mask
        r *= numpy.logical_not(     # invert triangle mask for triu (fast)
            numpy.tri(              # make inverse of triu mask
                *r.shape[-2:],      # get first two dimensions
                k=-1,               # offset to get lower triangle
                dtype=numpy.bool    # use bool (1 byte) data type (fast)
            )
        )
    if mtype == 'tril':             # if 'tril', multiply by a tril mask
        r *= numpy.tri(             # make lower triangle mask
            *r.shape[-2:],          # get first two dimensions
            k=0,                    # offset of zero to get lower triangle mask
            dtype=numpy.bool        # use bool (1 byte) data type (fast)
        )

    # return probability matrix
    return r



def r2haldane_fast(r):
    """
    Convert recombination probabilities between two loci to genetic map
    distances using the Haldane mapping function (Haldane, 1919).

    Parameters
    ----------
    r : numpy.ndarray
        An array of recombination probabilities between two loci.

    Returns
    -------
    d : numpy.ndarray
        An array of genetic map distances as defined by the Haldane mapping
        function.
    """
    return -0.5 * numpy.log(1 - (2 * r))
