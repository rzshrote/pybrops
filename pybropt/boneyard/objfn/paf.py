import numpy

def paf(rsel, geno, wcoeff, dtype=numpy.dtype('float64')):
    """
    Population allele frequency (PAF) score objective function.

    Bare bones objective function; does not check input types, etc.

    Computes all internal calculations in 'float64' to reduce floating point
    errors.
    """
    # generate a view of the geno matrix that only contains 'rsel' rows.
    sgeno = geno[:,rsel,:]

    # calculate population frequencies in 'float64' precision
    pfreq = sgeno.sum((0,1)) / numpy.float64(sgeno.shape[0] * sgeno.shape[1])

    # calculate dot product
    # use numpy.dot over pfreq[None,:] @ wcoeff[:,None] for speed
    paf = numpy.dot(pfreq, wcoeff)

    # return the dot product
    return dtype.type(paf)
