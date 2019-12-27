import numpy
from itertools import chain

################################################################################
# Create a couple of helper functions that will assist the main 'pa' function
################################################################################


# population architect selection
def pa(rslice,
       geno,
       coeff,
       tfreq,
       tfreq_edge,
       tld,
       tld_edge,
       cycles,
       d,
       lgroup_size,
       ldfn = None,
       chunk = None):
    """
    Population Architect (PA) Selection

    Parameters
    ==========
    rslice : numpy.ndarray, tuple, list, None
        A 1D array of indices to use for slicing and summing the matrix by row.
        Each index in 'rslice' represents a single individual's row. If 'rslice'
        is None, use all individuals.
        # TODO: performance testing to see which array type is better.
    geno : numpy.ndarray
        An array of allele states. The dtype of 'geno' should be 'uint8'. Array
        shape should be (depth, row, column) = (M, N, L) where 'M' represents
        number of chromosome phases, 'N' represents number of individuals, 'L'
        represents number of markers. Array format should be the 'C' format.
        # TODO: performance testing to see which array format is better 'F' or
                'C'.
    coeff : numpy.ndarray
        An array of coefficients for allele effects. The dtype of 'coeff'
        should be either 'float32' or 'float64'. This array should be single
        dimensional (column,) = (N,) where 'N' represents number of markers.
    tfreq : numpy.ndarray
        An array of target allele frequencies.
        Example:
            tfreq = numpy.array([0.2, 0.6, 0.7])
    tfreq_edge : numpy.ndarray
        An array of allele frequencies maximally distant from the target
        frequency.
        Example:
            tfreq      = numpy.array([0.2, 0.6, 0.7])
            tfreq_edge = numpy.array([1.0, 0.0, 0.0])
    tld : numpy.ndarray
        An array of target LD values. Current measure is in r squared.
        Example:
            tld = numpy.array([1.0, 1.0, 1.0])
    tld_edge : numpy.ndarray
        An array of LD values maximally distant from the target LD. Current
        measure is in r squared.
        Example:
            tld      = numpy.array([1.0, 1.0, 1.0])
            tld_edge = numpy.array([0.0, 0.0, 0.0])
    cycles : int, numpy.integer
        Number of breeding cycles left to the deadline.
    d : numpy.ndarray
        Genetic map in Morgan units.
    lgroup_size : numpy.ndarray
        Array of linkage group sizes. The sum of the elements should equal the
        length of 'd'.
    ldfn : None, function
        Linkage disequilibrium function. Customize this to get custom LD
        measures or crossing schemes.
    chunk : None, {int, numpy.integer}
        Matrix chunk size to compute. This is useful in situations with large
        numbers of markers (>5000). Large matrix calculations are exceedingly
        memory hungry (memory needed = n^2). Set this to a lower number to
        reduce memory consumption during matrix operations.

    Returns
    =======
    numpy.float64
        A floating point number representing the Population Architect score for
        the selected subset.
    """
    if chunk is None:
        # calculate pairwise marker weights
        # make triangle matrix
        marker_score = numpy.triu(coeff * coeff[:,None])
        # scale to sum to 1
        marker_score /= marker_score.sum()

        # calculate probabilities for recombination
        # fill with 0.5 (for independent assortment)
        prec = numpy.empty((len(d), len(d)), dtype=numpy.float64)
        prec.fill(0.5)

        # stop indices
        dsp = numpy.cumsum(lgroup_size)

        # start indices
        dst = dsp - lgroup_size[0]

        # for each start, stop index
        for st,sp in zip(dst,dsp):
            # make mesh for that array region
            tmpmesh = numpy.meshgrid(d[st:sp], d[st:sp])
            # calculate Haldane distance for that region
            prec[st:sp,st:sp] = 0.5 * (1 - numpy.exp(
                            -2 * numpy.abs(tmpmesh[0] - tmpmesh[1])))

        # select rows, THIS MAKES A COPY OF THE ARRAY (PERFORMANCE ISSUE)
        phase = geno[:,rslice,:].reshape(len(rslice)*geno.shape[0], geno.shape[2])

        # get the dimensions of the matrix
        rows, columns = phase.shape

        # convert the number of sequences (rows) to float64
        nseq_f64 = numpy.float64(rows)

        # calculate allele sums (axis=0), divide by nseq_f64 to get probabilities
        pfreq = phase.sum(0) / nseq_f64

        # calculate D; multiply by decay over 'cycles' based on recombination rate
        D = ((numpy.matmul(phase.transpose(), phase) / nseq_f64) - \
              numpy.matmul(pfreq[:,None], pfreq[None,:])) * \
              (1 - prec)**cycles

        # calculate A = p(1-p)
        A = pfreq * (1 - pfreq)

        # calculate r^2 matrix
        r_sq = (D**2) / numpy.matmul(A, A.transpose())

        # determine allele availability
        allele_avail = numpy.where(
            tfreq >= 1.0,
            pfreq > 0.0,
            numpy.where(
                tfreq > 0.0,
                numpy.logical_and(pfreq > 0.0, pfreq < 1.0),
                pfreq < 1.0
            )
        )

        # calculate pairwise allele availability
        pair_avail = allele_avail * allele_avail[:,None]

        # calculate difference between target and population
        diff_tp = tfreq - pfreq

        # calculate difference between target and edge
        diff_te = tfreq - tfreq_edge

        # make a mesh for differences between target and population; this is a list
        mesh_tp = numpy.meshgrid(diff_tp, diff_tp)

        # make a mesh for differences between edge and population; this is a list
        mesh_te = numpy.meshgrid(diff_te, diff_te)

        # calculate pairwise distances between target and population freq. loci
        # calculate pairwise distances between the edge and population freq. loci
        # finally divide population distance from target by distance to edge
        # do 1 - dist_prime
        allele_score = 1 - (numpy.sqrt(mesh_tp[0]**2, mesh_tp[1]**2) /
                            numpy.sqrt(mesh_te[0]**2, mesh_te[1]**2))

        # calculate the difference between target LD and population LD
        # calculate the difference between target LD and edge LD
        # calculate LD score component
        # modify the ld_score based on ambiguous LD states (r_sq == 0)
        # TODO: replace this with an approx. equal due to floating point errors
        ld_score = numpy.where(
            r_sq == 0.0,
            allele_avail,
            1 - ((tld - r_sq) / (tld - tld_edge))
        )

        tmpdiag = (marker_score * allele_score * ld_score)

        # marker_score is a triangle matrix, so it should zero out.
        # calculate product and return the sum of the triangle matrix.
        return tmpdiag.sum()
    else:
        return 0
