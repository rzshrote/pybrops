# 3rd party libraries
import numpy

# our libraries
from . import GenomicMating
import pybropt.util

class MOGM(GenomicMating):
    """docstring for MOGM."""

    ############################################################################
    ######################### Reserved Object Methods ##########################
    ############################################################################
    def __init__(self, population, cross, wcoeff = None, tfreq = None, method = "MOGM"):
        super(MOGM, self).__init__(population, cross, method)

        # check that we have marker coefficients
        pybropt.util.check_is_ParametricGenomicModel(
            self._population.genomic_model,
            "population.genomic_model"
        )

        # set wcoeff if needed
        if wcoeff is None:
            wcoeff = self.calc_wcoeff()
        self._wcoeff = wcoeff

        # set tfreq if needed
        if tfreq is None:
            tfreq = self.calc_tfreq()
        self._tfreq = tfreq

    ############################################################################
    ################################ Properties ################################
    ############################################################################
    def wcoeff():
        doc = "The wcoeff property."
        def fget(self):
            return self._wcoeff
        def fset(self, value):
            self._wcoeff = value
        def fdel(self):
            del self._wcoeff
        return locals()
    wcoeff = property(**wcoeff())

    def tfreq():
        doc = "The tfreq property."
        def fget(self):
            return self._tfreq
        def fset(self, value):
            self._tfreq = value
        def fdel(self):
            del self._tfreq
        return locals()
    tfreq = property(**tfreq())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    def calc_wcoeff(self):
        """
        Calculate weight coefficients.
        """
        # calculate wcoeff
        wcoeff = MOGM.wcoeff_mat(self._population.genomic_model.coeff)

        return wcoeff

    def calc_tfreq(self):
        """
        Calculate target frequencies
        """
        # calculate tfreq
        tfreq = MOGM.tfreq_mat(self._population.genomic_model.coeff)

        return tfreq

    def objfn(self, sel, objcoeff = None, minimizing = True):
        """
        MOGM objective function.

        Parameters
        ----------
        sel : numpy.ndarray
            A selection indices matrix of shape (k,)
            Where:
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, use all individuals.
        objcoeff : numpy.ndarray, None
            An objective coefficients matrix of shape (t,).
            Where:
                't' is the number of objectives.
            These are used to weigh objectives in the weight sum method.
            If None, do not multiply score by a weight sum vector.
        minimizing : bool, default = True
            If True, MOGM scores are adjusted such that OPV becomes a
            minimizing function: lower is better.
            If False, MOGM scores are adjusted such that OPV becomes a
            maximizing function: higher is better.
            Adjusted MOGM scores are used in the dot product with 'objcoeff'.

        Returns
        -------
        mogm : numpy.ndarray
            An MOGM score matrix of shape (k,) or (k,t) depending on whether
            'objcoeff' was specified or not, respectively.
        """
        # calculate MOGM values
        mogm = MOGM.objfn_mat(
            sel,
            self._population.geno,
            self._population.genomic_model.coeff,
            self._cross.varA,
            wcoeff = self._wcoeff,
            tfreq = self._tfreq
        )

        # MOGM is vanilla minimizing, but if we want pos
        if not minimizing:
            mogm = -mogm

        # if we have objective weights, take dot product for weight sum method
        if objcoeff is not None:
            if objcoeff.ndim == 1:
                mogm = mogm.dot(objcoeff)
            elif objcoeff.ndim == 2:
                # get dimensions
                ma, mb = mogm.shape
                oa, ob = objcoeff.shape

                # get dot products
                if (oa, ob) == (ma, mb):
                    mogm = (mogm * objcoeff).sum()
                elif (oa, ob) == (1, mb):
                    mogm = (mogm * objcoeff).sum(0)
                elif (oa, ob) == (ma, 1):
                    mogm = (mogm * objcoeff).sum(1)
                else:
                    raise ValueError("'mogm' and 'objcoeff' shapes not aligned")
            else:
                raise ValueError("'mogm' and 'objcoeff' shapes not aligned")

        return mogm

    def objfn_vec(self, sel, objcoeff = None, minimizing = True):
        # calculate MOGM values
        mogm = MOGM.objfn_vec_mat(
            sel,
            self._population.geno,
            self._population.genomic_model.coeff,
            wcoeff = self._wcoeff,
            tfreq = self._tfreq
        )

        # negate MOGM scores if necessary.
        if not minimizing:
            mogm = -mogm

        # if we have objective weights, take dot product for weight sum method
        if objcoeff is not None:
            mogm = mogm.dot(objcoeff)

        return mogm

    # optimize() function does not need to be overridden

    def simulate(self):
        raise NotImplementedError

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def objfn_mat(sel, geno, coeff, varAfn, wcoeff = None, tfreq = None):
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

        Given the provided genotype matrix 'geno' and row selections from it 'sel',
        calculate the selection allele freq. From the selection allele frequencies
        and the target allele frequencies, determine if the target frequencies
        cannot be attained after unlimited generations and selection rounds.
        Multiply this vector by a weight coefficients vector 'wcoeff'.

        f_PAFD(x):

        Given a genotype matrix, a target allele frequency vector, and a vector of
        weights, calculate the distance between the selection frequency and the
        target frequency.

        f_SPstdA(x)

        Given a progeny variance matrix and a crossing structure function, take
        the sum of standard deviations for each cross.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape (k,)
            Where:
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, use all individuals.
        geno : numpy.ndarray, None
            A int8 binary genotype matrix of shape (m, n, p).
            Where:
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.
            Remarks:
                Shape of the matrix is most critical. Underlying matrix
                operations will support other numeric data types.
        coeff : numpy.ndarray, None
            A trait prediction coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        varAfn : callable
            A function for getting variance components.
            Function signature is as follows:
            varAfn(sel, geno, coeff)
                Where:
                    sel : numpy.ndarray, None
                        A selection indices matrix of shape (k,).
                    geno : numpy.ndarray, None
                        A int8 binary genotype matrix of shape (m, n, p).
                    coeff : numpy.ndarray, None
                        A trait prediction coefficients matrix of shape (p, t).
                Remark:
                    varAfn does not have to utilize geno, coeff, but must
                    accept them.
                Returns:
                    varA : numpy.ndarray
                        A *variance* matrix of shape (a, t).
                        Where:
                            'a' is the number of mate configurations.
                            't' is the number of traits.
        wcoeff : numpy.ndarray, None
            A marker weight coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            Remarks: Values in 'wcoeff' have an assumption:
                All values must be non-negative.
        tfreq : numpy.ndarray, None
            A target allele frequency matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.

        Returns
        -------
        mogm : numpy.ndarray
            A MOGM score matrix of shape (3,t).
            Where:
                mogm[0] are PAU scores.
                mogm[1] are PAFD scores.
                mogs[2] are -SPstdA scores.
        """
        # if 'sel' is None, set to all individuals
        if sel is None:
            sel = slice(None)

        # if weight coefficient is not provided, calculate it.
        if wcoeff is None:
            wcoeff = MOGM.wcoeff_mat(coeff)

        # if tfreq is not provided, calculate it.
        if tfreq is None:
            tfreq = MOGM.tfreq_mat(coeff)

        # generate a view of the geno matrix that only contains 'sel' rows.
        sgeno = geno[:,sel,:]

        # calculate num chromosomes in the selection (m * n) as double
        phases = numpy.float64(sgeno.shape[0] * sgeno.shape[1])

        # calculate population frequencies; add axis for correct broadcast
        pfreq = (sgeno.sum((0,1)) / phases)[:,None]

        # calculate some inequalities for use multiple times
        pfreq_lteq_0 = (pfreq <= 0.0)   # is population freq <= 0.0
        pfreq_gteq_1 = (pfreq >= 1.0)   # is population freq >= 1.0

        # calculate allele unavailability: True == unavailable
        allele_unavail = numpy.where(
            tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
            pfreq_lteq_0,           # then set True if sel has allele freq == 0
            numpy.where(            # else
                tfreq > 0.0,        # if 0.0 < target freq < 1.0
                numpy.logical_or(   # then set TRUE if pop freq is outside (0.0,1.0)
                    pfreq_lteq_0,
                    pfreq_gteq_1
                ),
                pfreq_gteq_1        # else set TRUE if pop freq is >= 1.0
            )
        )

        # calculate distance between target and population
        dist = numpy.absolute(tfreq - pfreq) # (p,t)-(p,1) -> (p,t)


        ####################################################
        # compute f_PAU(x)
        pau = (wcoeff * allele_unavail).sum(0) # (p,t)[0] -> (t,)

        # compute f_PAFD(x)
        pafd = (wcoeff * dist).sum(0) # (p,t)[0] -> (t,)

        # compute f_SPstdA(x)
        spstda = numpy.sqrt(varAfn(sel)).sum(1) # (t,a)[0] -> (t,)

        # allocate an empty vector of float64 to hold pau, pafd, stdA
        mogm = numpy.empty((3, wcoeff.shape[1]), dtype='float64')

        # copy scores into rows of allocated matrix
        mogm[0] = pau
        mogm[1] = pafd
        mogm[2] = -spstda # negate before copy

        return mogm

    @staticmethod
    def objfn_vec_mat(sel, geno, coeff, varAfn, wcoeff = None, tfreq = None):
        """

        """
        # if 'sel' is None, set to all individuals
        if sel is None:
            sel = numpy.arange(geno.shape[1])[None,:]

        # if weight coefficient is not provided, calculate it.
        if wcoeff is None:
            wcoeff = MOGM.wcoeff_mat(coeff)

        # if tfreq is not provided, calculate it.
        if tfreq is None:
            tfreq = MOGM.tfreq_mat(coeff)

        # generate a view of the geno matrix that only contains 'sel' rows.
        sgeno = geno[:,sel,:] # (m,n,p)[1] -> (m,j,k,p)

        # calculate num chromosomes in the selection (m * n) as double
        phases = numpy.float64(sgeno.shape[0] * sgeno.shape[2]) # (m*k)

        # calculate population frequencies; add axis for correct broadcast
        pfreq = (sgeno.sum((0,2)) / phases)[:,None] # (m,j,k,p)[0,2] -> (j,p,1)

        # calculate some inequalities for use multiple times
        pfreq_lteq_0 = (pfreq <= 0.0)   # is population freq <= 0.0
        pfreq_gteq_1 = (pfreq >= 1.0)   # is population freq >= 1.0

        # calculate allele unavailability: True == unavailable
        allele_unavail = numpy.where(
            tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
            pfreq_lteq_0,           # then set True if sel has allele freq == 0
            numpy.where(            # else
                tfreq > 0.0,        # if 0.0 < target freq < 1.0
                numpy.logical_or(   # then set TRUE if pop freq is outside (0.0,1.0)
                    pfreq_lteq_0,
                    pfreq_gteq_1
                ),
                pfreq_gteq_1        # else set TRUE if pop freq is >= 1.0
            )
        )

        # calculate distance between target and population
        dist = numpy.absolute(tfreq - pfreq) # (p,t)-(j,p,1) -> (j,p,t)


        ####################################################
        # compute f_PAU(x)
        pau = (wcoeff * allele_unavail).sum(1) # (j,p,t)[1] -> (j,t)

        # compute f_PAFD(x)
        pafd = (wcoeff * dist).sum(1) # (j,p,t)[1] -> (j,t)

        # compute f_SPstdA(x)
        spstda = numpy.sqrt(varAfn(sel, geno, coeff)).sum(1) # (j,a,t)[1] -> (j,t)

        # allocate an empty vector of float64 to hold pau, pafd, stdA
        mogm = numpy.empty((3, pau.shape[0], pau.shape[1]), dtype='float64')

        # copy scores into rows of allocated matrix
        mogm[0] = pau
        mogm[1] = pafd
        mogm[2] = -spstda # negate before copy

        return mogm # (3,j,t)

    @staticmethod
    def wcoeff_mat(coeff):
        wcoeff = numpy.absolute(coeff)
        return wcoeff

    @staticmethod
    def tfreq_mat(coeff):
        tfreq = (coeff >= 0).astype('float64')
        return tfreq
