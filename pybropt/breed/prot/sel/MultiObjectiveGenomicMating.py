import numpy
import math
import types

from deap import algorithms
from deap import base
from deap import creator
from deap import tools
from deap import benchmarks

import pybropt.core.random

from pybropt.breed.prot.sel.SelectionProtocol import SelectionProtocol

from pybropt.core.error import check_is_int
from pybropt.core.error import check_is_str

class MultiObjectiveGenomicMating(SelectionProtocol):
    """docstring for MultiObjectiveGenomicMating."""

    def __init__(self, arg):
        super(MultiObjectiveGenomicMating, self).__init__()
        self.arg = arg

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn(sel, mat, tfreq, mkrwt, trans, kwargs):
        """
        Multi-objective genomic selection objective function.
            The goal is to minimize this function. Lower is better.
            This is a bare bones function. Minimal error checking is done.

        Given a 2D weight vector 'dcoeff', calculate the Euclidian distance from the
        origin according to:

        dist = dot( dcoeff, F(x) )

        Where :math:`F(\\textbf{x})` is a vector of objective functions:

        .. math::
            F(\\textbf{x}) = \begin{bmatrix}
            f^{\\textup{PAU}}(\\textbf{x})  \\\\
            f^{\\textup{PAFD}}(\\textbf{x})
            \end{bmatrix}

        Objectives:

        Population Allele Unavailability (PAU): :math:`f^{\\textup{PAU}}(\\textbf{x})`

        Given the provided genotype matrix 'geno' and row selections from it 'sel',
        calculate the selection allele freq. From the selection allele frequencies
        and the target allele frequencies, determine if the target frequencies
        cannot be attained after unlimited generations and selection rounds.
        Multiply this vector by a weight coefficients vector 'wcoeff'.

        Population Allele Frequency Distance (PAFD): :math:`f^{\\textup{PAFD}}(\\textbf{x})`

        Given a genotype matrix, a target allele frequency vector, and a vector of
        weights, calculate the distance between the selection frequency and the
        target frequency.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape (k,)
            Where:
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, use all individuals.
        mat : numpy.ndarray, None
            A int8 binary genotype matrix of shape (m, n, p).
            Where:
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.
            Remarks:
                Shape of the matrix is most critical. Underlying matrix
                operations will support other numeric data types.
        tfreq : floating, numpy.ndarray
            A target allele frequency matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            Example:
                tfreq = numpy.array([0.2, 0.6, 0.7])
        mkrwt : numpy.ndarray
            A marker weight coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            Remarks: Values in 'mkrwt' have an assumption:
                All values must be non-negative.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:
                Must accept a single numpy.ndarray argument.
                Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to 'trans' function.

        Returns
        -------
        mogs : numpy.ndarray
            A MOGS score matrix of shape (2,t) or other.
        """
        # if no selection, select all
        if sel is None:
            sel = slice(None)

        # generate a view of the genotype matrix that only contains 'sel' rows.
        # (m,(k,),p) -> (m,k,p)
        sgeno = mat[:,sel,:]

        # calculate reciprocal number of phases
        rphase = 1.0 / (sgeno.shape[0] * sgeno.shape[1])

        # calculate population frequencies; add axis for correct broadcast
        # (m,k,p).sum[0,1] -> (p,)
        # (p,) * scalar -> (p,)
        # (p,None) -> (p,1) # need (p,1) for broadcasting with (p,t) arrays
        pfreq = (sgeno.sum((0,1)) * rphase)[:,None]

        # calculate some inequalities for use multiple times
        pfreq_lteq_0 = (pfreq <= 0.0)   # is population freq <= 0.0
        pfreq_gteq_1 = (pfreq >= 1.0)   # is population freq >= 1.0

        # calculate allele unavailability
        allele_unavail = numpy.where(
            tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
            pfreq_lteq_0,           # then set True if sel has allele freq == 0
            numpy.where(            # else
                tfreq > 0.0,        # if 0.0 < target freq < 1.0
                numpy.logical_or(   # then set True if pop freq is outside (0.0,1.0)
                    pfreq_lteq_0,
                    pfreq_gteq_1
                ),
                pfreq_gteq_1        # else set True if pop freq is >= 1.0
            )
        )

        # calculate distance between target and population
        # (p,t)-(p,1) -> (p,t)
        dist = numpy.absolute(tfreq - pfreq)

        # compute f_PAU(x)
        # (p,t) * (p,t) -> (p,t)
        # (p,t).sum[0] -> (t,)
        pau = (mkrwt * allele_unavail).sum(0)

        # compute f_PAFD(x)
        # (p,t) * (p,t) -> (p,t)
        # (p,t).sum[0] -> (t,)
        pafd = (mkrwt * dist).sum(0)

        # stack to make MOGS matrix
        # (2,t)
        mogs = numpy.stack([pau, pafd])

        # transform and return
        return trans(mogs, **kwargs)
