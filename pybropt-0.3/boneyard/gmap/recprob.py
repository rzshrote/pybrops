import numpy

# import a function lookup table
from .const import GMAP_MAP2R_DICT


def recprob(rst, rsp, cst, csp, gmap, lgroup, mapfn = 'haldane'):
    """
    Determine recombination probabilities of a region of a recombination matrix
    given a genetic map and a mapping function.

    Parameters
    ==========
    rst : int
        If the full recombination matrix is (n x n), rst is the starting row
        index (inclusive) to calculate.
    rsp : int
        If the full recombination matrix is (n x n), rsp is the stopping row
        index (exclusive) to calculate.
    cst : int
        If the full recombination matrix is (n x n), cst is the starting column
        index (inclusive) to calculate.
    csp : int
        If the full recombination matrix is (n x n), csp is the stopping column
        index (exclusive) to calculate.
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

    Returns
    =======
    r : numpy.ndarray
        A numpy.ndarray of size (rsp-rst, csp-cst) containing recombination
        probabilities.
    """
    # test if mapfn is a callable function. If not, lookup in lookup table.
    if not callable(mapfn):
        mapfn = GMAP_MAP2R_DICT[mapfn]   # assume that mapfn is a valid key

    # allocate recombination probability matrix
    r = numpy.empty(            # allocate empty matrix
        (rsp-rst, csp-cst),     # (rows, columns) dimensions
        dtype=numpy.float64     # use float64 dtype
    )

    # variables to track columns we haven't filled yet
    lcols = csp  # cols to the left not filled (index)
    rcols = cst  # cols to the right not filled (index)

    # Algorithm:
    #     For each linkage group, calculate the overlap with the
    #     current computational chunk. If there is an overlap,
    #     compute the overlap. Otherwise, fill with 0.5.

    # start and stop coordinate variable for the linkage groups.
    lst = 0
    lsp = 0

    # cycle through all the linkage group sectors trying to identify
    # regions of overlap between the computational chunk and the
    # linkage group block
    for g in lgroup:
        # advance the lsp index by the size of the linkage group.
        lsp += g

        # calculate an overlap region defining a box: (c1,r1,c2,r2)
        # where c1 < c2, r1 < r2
        c1 = max(lst, cst)   # c1; matrix start column
        r1 = max(lst, rst)   # r1; matrix start row
        c2 = min(lsp, csp)   # c2; matrix stop column
        r2 = min(lsp, rsp)   # r2; matrix start row

        # if we've found a linkage group that overlaps
        if (c2 >= c1) and (r2 >= r1):
            # get the cols overlap (x range) as a slice
            cols = slice(c1-cst, c2-cst)

            # make distance mesh for the overlap region
            dmesh = numpy.meshgrid(
                gmap[r1:r2],    # number of rows
                gmap[c1:c2],    # number of columns
                indexing='ij'   # use ij indexing for speed
            )

            # fill the recombination frequencies in the matrix
            #   rst-r1:rst-r2     # rows (y)
            #   cols              # columns (x)
            r[r1-rst:r2-rst,cols] = mapfn(
                numpy.absolute(dmesh[0] - dmesh[1])
            )

            # fill rows above the linkage block with 0.5
            r[:r1-rst,cols].fill(0.5)

            # fill rows below the linkage block with 0.5
            r[r2-rst:,cols].fill(0.5)

            # alter column tracking variables
            if c1 < lcols:  # if c1 is before lcols index
                lcols = c1  # set the index to c1
            if c2 > rcols:  # if c2 is after rcols index
                rcols = c2  # set the index to c2

        # advance the linkage group start index to the current stop
        lst = lsp

    # fill the remaining untouched columns with 0.5
    if lcols > rcols:             # if lcols and rcols overlap
        r.fill(0.5)               # fill the whole thing in one go
    else:                         # else fill each side
        r[:,:lcols-cst].fill(0.5) # fill left
        r[:,rcols-cst:].fill(0.5) # fill right

    return r
