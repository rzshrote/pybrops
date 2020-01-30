import numpy

def mogm(rsel, geno, wcoeff, tfreq, varA, varAfn, dcoeff,
         dtype = numpy.dtype('float64')):
    """
    Multi-Objective Genomic Mating (MOGM) objective function.
        The goal is to minimize this function. Lower is better.
        This is a bare bones function. Minimal error checking is done.

    Given a weight vector 'dcoeff', calculate the dot product of F(x) and the
    weight vector.
        score = dot( dcoeff, F(x) )
        Where:
            F(x) is a vector of objective functions:
                F(x) = < f_PAU(x), f_PAFD(x), f_stdA(x) >
            wcoeff is a vector of objective weights:
                wcoeff = < f_PAU_weight, f_PAFD_weight, f_stdA_weight >

    f_PAU(x):

    Given the provided genotype matrix 'geno' and row selections from it 'rsel',
    calculate the selection allele freq. From the selection allele frequencies
    and the target allele frequencies, determine if the target frequencies
    cannot be attained after unlimited generations and selection rounds.
    Multiply this vector by a weight coefficients vector 'wcoeff'.

    f_PAFD(x):

    Given a genotype matrix, a target allele frequency vector, and a vector of
    weights, calculate the distance between the selection frequency and the
    target frequency.

    f_stdA(x)

    Given a progeny variance matrix and a crossing structure function, take
    the sum of standard deviations for each cross.

    """
    # if varAfn is not callable, we assume it's a key for STDA_DICT
    if not callable(varAfn):
        varAfn = STDA_DICT[varAfn]

    ###############################################################
    ## create several variables for in matrix calculations below ##
    # calculate number of phases
    phases = geno[:,rsel,:].shape[0] * geno[:,rsel,:].shape[1]

    # calculate allele frequencies
    pfreq = geno[:,rsel,:].sum((0,1)) / numpy.float64(phases)

    # calculate allele availability for the selected population
    pfreq_gt_0 = (pfreq > 0.0)  # make mask for pfreq > 0
    pfreq_lt_1 = (pfreq < 1.0)  # make mask for pfreq < 1
    allele_avail = numpy.where( # assess allele availability
        tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
        pfreq_gt_0,             # then set TRUE if population has >0 allele freq
        numpy.where(            # else
            tfreq > 0.0,        # if target freq > 0.0
            numpy.logical_and(  # then set TRUE if pop freq is between (0.0,1.0)
                pfreq_gt_0,
                pfreq_lt_1
            ),
            pfreq_lt_1          # else set TRUE if pop freq is less than 1.0
        )
    )

    # make a vector of doubles to store scores for F(x)
    F_x = numpy.zeros(3, dtype='float64')
    ###############################################################

    ############################# Compute f_PAU(x) #############################
    F_x[0] = numpy.logical_not(allele_avail).dot(wcoeff)

    ############################# Compute f_PAU(x) #############################
    F_x[1] = numpy.absolute(tfreq - pfreq).dot(wcoeff)

    ############################ Compute f_stdA(x) #############################
    F_x[2] = varAfn(rsel, varA)

    ###################### Compute F(x) and return score #######################
    # calculate weight sum method
    score = F_x.dot(dcoeff)

    # return the score casted as the output data type
    return dtype.type(score)


def mogm_max(rsel, geno, wcoeff, tfreq, varA, varAfn, dcoeff,
         dtype = numpy.dtype('float64')):
    """
    Negate Multi-Objective Genomic Mating (MOGM) function 'mogm' such that the
    objective becomes to maximize the function instead.
    """
    return -1 * mogs(**locals())
