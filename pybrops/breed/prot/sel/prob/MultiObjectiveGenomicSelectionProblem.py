"""
Module implementing Multi-Objective Genomic Selection (MOGS) optimization problems.
"""

__all__ = [
    "MultiObjectiveGenomicSubsetSelectionProblem",
]

from abc import ABCMeta
from abc import abstractmethod
from numbers import Integral
from numbers import Real
from typing import Callable
from typing import Optional
from typing import Union

import numpy
# from pybrops.breed.prot.sel.prob.IntegerSelectionProblem import IntegerSelectionProblem
# from pybrops.breed.prot.sel.prob.RealSelectionProblem import RealSelectionProblem
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix


class MultiObjectiveGenomicSelectionProblemMixin(
        metaclass = ABCMeta,
    ):
    """Helper class to implement properties common to MOGS."""

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in haplotype matrix
        return 2 * self._mkrwt.shape[1]

    ################# genotype information #################
    @property
    def geno(self) -> numpy.ndarray:
        """Genotype matrix of shape (n,p) in {0,1,2} format."""
        return self._geno
    @geno.setter
    def geno(self, value: numpy.ndarray) -> None:
        """Set genotype matrix."""
        check_is_ndarray(value, "geno")
        check_ndarray_ndim(value, "geno", 2)
        self._geno = value
    
    ##################### Ploidy level #####################
    @property
    def ploidy(self) -> Integral:
        """ploidy."""
        return self._ploidy
    @ploidy.setter
    def ploidy(self, value: Integral) -> None:
        """Set ploidy."""
        check_is_Integral(value, "ploidy")
        check_is_gt(value, "ploidy", 0)
        self._ploidy = value
    
    @property
    def mkrwt(self) -> numpy.ndarray:
        """Marker weights."""
        return self._mkrwt
    @mkrwt.setter
    def mkrwt(self, value: numpy.ndarray) -> None:
        """Set marker weights."""
        check_is_ndarray(value, "mkrwt")
        check_ndarray_ndim(value, "mkrwt", 2)
        self._mkrwt = value

    ############ target allele frequency values ############
    @property
    def tfreq(self) -> numpy.ndarray:
        """Target allele frequency."""
        return self._tfreq
    @tfreq.setter
    def tfreq(self, value: numpy.ndarray) -> None:
        """Set target allele frequency."""
        check_is_ndarray(value, "tfreq")
        check_ndarray_ndim(value, "tfreq", 2)
        self._tfreq = value
        self._tminor = self._calc_tminor(self._tfreq)
        self._thet = self._calc_thet(self._tfreq)
        self._tmajor = self._calc_tminor(self._tfreq)
    
    @property
    def tminor(self) -> numpy.ndarray:
        """Whether the target allele frequency is fixation of a minor allele."""
        return self._tminor
    
    @property
    def thet(self) -> numpy.ndarray:
        """Whether the target allele frequency is heterozygous."""
        return self._thet
    
    @property
    def tmajor(self) -> numpy.ndarray:
        """Whether the target allele frequency is fixation of a major allele."""
        return self._tmajor
    
    ######################### Private Object Methods ###########################
    @staticmethod
    def _calc_mkrwt(weight: Union[numpy.ndarray,Callable], u_a: numpy.ndarray):
        if callable(weight):
            return weight(u_a)
        elif isinstance(weight, numpy.ndarray):
            return weight
        else:
            raise TypeError("variable 'weight' must be a callable function or numpy.ndarray")
    
    @staticmethod
    def _calc_tfreq(target: Union[numpy.ndarray,Callable], u_a: numpy.ndarray):
        if callable(target):
            return target(u_a)
        elif isinstance(target, numpy.ndarray):
            return target
        else:
            raise TypeError("variable 'target' must be a callable function or numpy.ndarray")

    @staticmethod
    def _calc_tminor(tfreq: numpy.ndarray):
        return (tfreq == 0.0)
    
    @staticmethod
    def _calc_tmajor(tfreq: numpy.ndarray):
        return (tfreq == 1.0)
    
    @staticmethod
    def _calc_thet(tfreq: numpy.ndarray):
        return (tfreq > 0.0) & (tfreq < 1.0)
    
    ############################## Class Methods ###############################
    @classmethod
    @abstractmethod
    def from_gmat_gpmod(
            cls,
            gmat: GenotypeMatrix,
            weight: Union[numpy.ndarray,Callable],
            target: Union[numpy.ndarray,Callable],
            gpmod: AdditiveLinearGenomicModel,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            **kwargs: dict
        ) -> "MultiObjectiveGenomicSelectionProblemMixin":
        raise NotImplementedError("class method is abstract")

class MultiObjectiveGenomicSubsetSelectionProblem(
        MultiObjectiveGenomicSelectionProblemMixin,
        SubsetSelectionProblem,
    ):
    """
    docstring for SubsetMultiObjectiveSelectionProblem.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            geno: numpy.ndarray,
            ploidy: Integral,
            mkrwt: numpy.ndarray,
            tfreq: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetMultiObjectiveSelectionProblem.
        
        Parameters
        ----------
        geno : numpy.ndarray
            A genotype matrix of shape ``(n,p)`` representing only biallelic
            loci. One of the two alleles at a locus is coded using a ``1``. The
            other allele is coded as a ``0``. ``mat`` holds the counts of the
            allele coded by ``1``.

            Where:

            - ``n`` is the number of individuals.
            - ``p`` is the number of markers.

            Example::

                # matrix of shape (n = 3, p = 4)
                geno = numpy.array([[0,2,1,0],
                                    [2,2,1,1],
                                    [0,1,0,2]])
        ploidy : Integral
            Number of phases that the genotype matrix ``mat`` represents.
        tfreq : numpy.ndarray
            A target allele frequency matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Example::

                tfreq = numpy.array([0.2, 0.6, 0.7, 0.5])
        mkrwt : numpy.ndarray
            A marker weight coefficients matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Remarks:

            - All values in ``mkrwt`` must be non-negative.
        ndecn : Integral
            Number of decision variables.
        decn_space: numpy.ndarray, None
            An array of shape ``(2,ndecn)`` defining the decision space.
            If None, do not set a decision space.
        decn_space_lower: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing lower limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a lower limit for the decision variables.
        decn_space_upper: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing upper limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a upper limit for the decision variables.
        nobj: Integral
            Number of objectives.
        obj_wt: numpy.ndarray
            Objective function weights.
        obj_trans: Callable, None
            A transformation function transforming a latent space vector to an objective space vector.
            The transformation function must be of the form: ``obj_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the identity transformation function: copy the latent space vector to the objective space vector.
        obj_trans_kwargs: dict, None
            Keyword arguments for the latent space to objective space transformation function.
            If None, an empty dictionary is used.
        nineqcv: Integral,
            Number of inequality constraints.
        ineqcv_wt: numpy.ndarray,
            Inequality constraint violation weights.
        ineqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an inequality constraint violation vector.
            The transformation function must be of the form: ``ineqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        ineqcv_trans_kwargs: Optional[dict],
            Keyword arguments for the latent space to inequality constraint violation space transformation function.
            If None, an empty dictionary is used.
        neqcv: Integral
            Number of equality constraints.
        eqcv_wt: numpy.ndarray
            Equality constraint violation weights.
        eqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an equality constraint violation vector.
            The transformation function must be of the form: ``eqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        eqcv_trans_kwargs: dict, None
            Keyword arguments for the latent space to equality constraint violation space transformation function.
            If None, an empty dictionary is used.
        kwargs : dict
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        super(MultiObjectiveGenomicSubsetSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            **kwargs
        )
        # assignments
        self.geno = geno
        self.ploidy = ploidy
        self.mkrwt = mkrwt
        self.tfreq = tfreq

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Multi-objective genomic selection objective function.

        - The goal is to minimize all objectives for this function.
        - This is a bare bones function. Minimal error checking is done.

        Objectives: :math:`F(\\textbf{x})`

        .. math::

            F(\\textbf{x}) = {[f^{\\textup{PAU}}(\\textbf{x}), f^{\\textup{PAFD}}(\\textbf{x})]}'

        Population Allele Unavailability (PAU): :math:`f^{\\textup{PAU}}(\\textbf{x})`

        .. math::

            f^{\\textup{PAU}}(\\textbf{x}) = \\textbf{w} \\cdot \\textbf{u}

        Given a genotype matrix ``mat`` and a selection indices vector
        :math:`\\textbf{x} =` ``sel``, calculate the selection allele frequency.
        From the selection allele frequencies and the target allele frequencies
        ``tfreq``, determine if the target frequencies can be attained after
        unlimited generations of selection. If the target allele frequency at a
        locus cannot be attained, score locus as ``1``, otherwise score as
        ``0``. Store this into a binary score vector :math:`\\textbf{u}`.
        Take the dot product between the binary score vector and the marker
        weight vector :math:`\\textbf{w} =` ``mkrwt`` to calculate
        :math:`f^{\\textup{PAU}}(\\textbf{x})` and return the result.

        Population Allele Frequency Distance (PAFD): :math:`f^{\\textup{PAFD}}(\\textbf{x})`

        .. math::
            f^{\\textup{PAFD}}(\\textbf{x}) = \\textbf{w} \\cdot \\left | \\textbf{p}_{x} - \\textbf{p}_{t} \\right |

        Given a genotype matrix ``mat`` and a selection indices vector
        :math:`\\textbf{x} =` ``sel``, calculate the selection allele frequency
        :math:`\\textbf{p}_{x}`. From the selection allele frequencies and the
        target allele frequencies :math:`\\textbf{p}_{t} =` ``tfreq``,
        calculate the absolute value of the difference between the two vectors.
        Finally, take the dot product between the difference vector and the marker
        weight vector :math:`\\textbf{w} =` ``mkrwt`` to calculate
        :math:`f^{\\textup{PAFD}}(\\textbf{x})` and return the result.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,)``.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            An MOGS matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # calculate the allele frequency of the selected subset
        # (n,p)[(k,),:,None] -> (p,1)
        pfreq = (1.0 / (self.ploidy * len(x))) * self.geno[x,:,None].sum(0)

        # determine where allele frequencies are < 1.0
        # (p,1)
        p_ltmajor = (pfreq < 1.0)

        # determine where allele frequencies are > 0.0
        # (p,1)
        p_gtminor = (pfreq > 0.0)

        # determine where allele frequencies are < 1.0 and > 0.0
        # (p,1)
        p_het = numpy.logical_and(p_ltmajor, p_gtminor)

        # determine where alleles are unavailable using precomputed arrays
        # (p,t)
        allele_unavail = numpy.logical_not(
            numpy.logical_or(
                numpy.logical_and(p_ltmajor, self.tminor), 
                numpy.logical_or(
                    numpy.logical_and(p_het, self.thet), 
                    numpy.logical_and(p_gtminor, self.tmajor)
                )
            )
        )

        # calculate the manhattan distance and PAFD
        # (p,t) -> (t,)
        pafd = (self.mkrwt * numpy.absolute(self.tfreq - pfreq)).sum(0)
        
        # calculate the allele unavailability
        # (p,t) -> (t,)
        pau = (self.mkrwt * allele_unavail).sum(0)

        # concatenate to make MOGS vector
        # (t,) and (t,) -> (t + t,)
        out = numpy.concatenate([pau, pafd])

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_gmat_gpmod(
            cls,
            gmat: GenotypeMatrix,
            weight: Union[numpy.ndarray,Callable],
            target: Union[numpy.ndarray,Callable],
            gpmod: AdditiveLinearGenomicModel,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            **kwargs: dict
        ) -> "MultiObjectiveGenomicSubsetSelectionProblem":
        # extract genotype matrix
        geno = gmat.mat_asformat("{0,1,2}")
        ploidy = gmat.ploidy
        mkrwt = cls._calc_mkrwt(weight, gpmod.u_a)
        tfreq = cls._calc_tfreq(target, gpmod.u_a)

        # construct class
        out = cls(
            geno = geno,
            ploidy = ploidy,
            mkrwt = mkrwt,
            tfreq = tfreq,
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            **kwargs
        )

        return out
