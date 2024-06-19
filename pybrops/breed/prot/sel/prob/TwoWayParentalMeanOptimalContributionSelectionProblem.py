"""
Module defining optimization problems for parental mean optimal contribution selection.
"""

__all__ = [
    "TwoWayParentalMeanOptimalContributionSelectionProblemMixin",
    "TwoWayParentalMeanOptimalContributionSubsetSelectionProblem",
]

from abc import ABCMeta, abstractmethod
from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy

from pybrops.breed.prot.sel.prob.SubsetMateSelectionProblem import SubsetMateSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len, check_ndarray_is_square, check_ndarray_is_triu, check_ndarray_ndim
from pybrops.core.util.crossix import twowayix, twowayix_len
from pybrops.model.pmebvmat.DenseTwoWayParentalMeanEstimatedBreedingValueMatrix import DenseTwoWayParentalMeanEstimatedBreedingValueMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix


class TwoWayParentalMeanOptimalContributionSelectionProblemMixin(
        metaclass = ABCMeta,
    ):
    """Helper class containing common properties for Parental Mean Optimal Contribution Selection Problems."""

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        return 1 + self._pmebv.shape[1]
    
    @property
    def pmebv(self) -> numpy.ndarray:
        """Parental mean estimated breeding value matrix."""
        return self._pmebv
    @pmebv.setter
    def pmebv(self, value: numpy.ndarray) -> None:
        """Set parental mean estimated breeding value matrix."""
        check_is_ndarray(value, "pmebv")
        check_ndarray_ndim(value, "pmebv", 2)
        self._pmebv = value
    
    @property
    def xmap(self) -> numpy.ndarray:
        """Cross map."""
        return self._xmap
    @xmap.setter
    def xmap(self, value: numpy.ndarray) -> None:
        """Set cross map."""
        check_is_ndarray(value, "xmap")
        check_ndarray_ndim(value, "xmap", 2)
        check_ndarray_axis_len(value, "xmap", 0, len(self._pmebv))
        check_ndarray_axis_len(value, "xmap", 1, 2)
        self._xmap = value
    
    @property
    def xcontrib(self) -> numpy.ndarray:
        """Cross contribution."""
        return self._xcontrib
    @xcontrib.setter
    def xcontrib(self, value: numpy.ndarray) -> None:
        """Set cross contribution."""
        check_is_ndarray(value, "xcontrib")
        check_ndarray_ndim(value, "xcontrib", 2)
        check_ndarray_axis_len(value, "xcontrib", 0, len(self._pmebv))
        check_ndarray_axis_len(value, "xcontrib", 1, 2)
        self._xcontrib = value
    
    @property
    def C(self) -> numpy.ndarray:
        """Cholesky decomposition of the kinship matrix."""
        return self._C
    @C.setter
    def C(self, value: numpy.ndarray) -> None:
        """Set Cholesky decomposition of the kinship matrix."""
        check_is_ndarray(value, "C")
        check_ndarray_ndim(value, "C", 2)
        check_ndarray_is_square(value, "C")
        check_ndarray_is_triu(value, "C")
        self._C = value

    ######################### Private Object Methods ###########################
    @staticmethod
    def _calc_pmebv(
            bvmat: BreedingValueMatrix, 
            xmap: numpy.ndarray,
            unscale: bool = True,
        ) -> numpy.ndarray:
        """
        Construct a parental mean breeding value matrix.

        Parameters
        ----------
        bvmat : BreedingValueMatrix
            Input breeding value matrix.
        
        Returns
        -------
        out : numpy.ndarray
            A breeding value matrix of shape ``(n,n,t)``. May be scaled or unscaled.
        """
        # calculate pmEBV matrix
        pmebvmat = DenseTwoWayParentalMeanEstimatedBreedingValueMatrix.from_bvmat(bvmat = bvmat)

        # unscale
        if unscale:
            pmebvmat.unscale(inplace=True)
        
        # get pmEBV matrix
        tmp = pmebvmat.mat
        
        # get female indices
        female = xmap[:,0]
        male = xmap[:,1]

        # (s,d)
        out = tmp[female,male,:]

        return out

    @staticmethod
    def _calc_xmap(
            ntaxa: Integral, 
            symab: bool = True,
            mateab: str = "uniq",
        ) -> numpy.ndarray:
        """
        Calculate the cross map.

        Parameters
        ----------
        ntaxa : int
            Number of taxa.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(s,d)`` containing cross map indices.

            Where:

            - ``s`` is the number of elements in the upper triangle, including
              or not including the diagonal (depending on ``unique_parents``).
            - ``d`` is the number of parents in the cross.
        """
        # construct the output from an iterator
        out = numpy.fromiter(
            iter = twowayix(ntaxa, symab, mateab),
            dtype = numpy.dtype((int, 2)),
            count = twowayix_len(ntaxa, symab, mateab)
        )
        return out

    @staticmethod
    def _calc_xcontrib(
            ntaxa: Integral, 
            symab: bool = True,
            mateab: str = "uniq",
        ) -> numpy.ndarray:
        """
        Calculate the cross contribution map.
        """
        out = numpy.tile((0.5,0.5), (twowayix_len(ntaxa, symab, mateab),1))
        return out

    @staticmethod
    def _calc_C(gmat: GenotypeMatrix, cmatfcty: CoancestryMatrixFactory) -> numpy.ndarray:
        """
        Construct a Cholesky decomposition of a coancestry matrix.

        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotype matrix from which to calculate a coancestry matrix.
        
        Returns
        -------
        out : numpy.ndarray
            A Cholesky decomposition of a coancestry matrix of shape ``(n,n)``.
        """
        # get genomic relationship matrix: (n,n)
        G = cmatfcty.from_gmat(gmat)

        # to ensure we're able to perform cholesky decomposition, apply jitter if needed.
        # if we are unable to fix, then raise value error
        if not G.apply_jitter():
            raise ValueError(
                "Unable to construct objective function: Kinship matrix is not positive definite.\n"+
                "    This could be caused by lack of genetic diversity.\n"
            )

        # convert G to (1/2)G (kinship analogue): (n,n)
        K = G.mat_asformat("kinship")

        # cholesky decomposition of K matrix: (n,n)
        out = numpy.linalg.cholesky(K).T

        return out

    ############################## Class Methods ###############################
    @classmethod
    @abstractmethod
    def from_bvmat_gmat(
            cls,
            bvmat: BreedingValueMatrix,
            gmat: GenotypeMatrix, 
            cmatfcty: CoancestryMatrixFactory,
            unscale: bool,
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
        ) -> "TwoWayParentalMeanOptimalContributionSelectionProblemMixin":
        raise NotImplementedError("class method is abstract")

class TwoWayParentalMeanOptimalContributionSubsetSelectionProblem(
        TwoWayParentalMeanOptimalContributionSelectionProblemMixin,
        SubsetMateSelectionProblem,
    ):
    """
    Class representing a Parental Mean Optimal Contribution Selection Problem for subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            pmebv: numpy.ndarray,
            xmap: numpy.ndarray,
            xcontrib: numpy.ndarray,
            C: numpy.ndarray,
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
        Constructor for TwoWayParentalMeanOptimalContributionSubsetMateSelectionProblem.
        
        Parameters
        ----------
        pmebv : numpy.ndarray
            A breeding value matrix of shape ``(s,t)``. 
            If you are using a penalization transformation function, preferably
            these breeding values are centered and scaled to make the penalies 
            less extreme.

            Where:

            - ``s`` is the number of crosses.
            - ``t`` is the number of traits.
        xmap : numpy.ndarray
            A cross map of shape ``(s,d)``.
        xcontrib : numpy.ndarray
            A cross parental contribution matrix of shape ``(s,2)``.
        C : numpy.ndarray
            An upper triangle matrix of shape ``(n,n)`` resulting from a Cholesky 
            decomposition of a kinship matrix: K = C'C.

            Where:

            - ``n`` is the number of individuals.
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
            Additional keyword arguments passed to the parent class (SubsetMateSelectionProblem) constructor.
        """
        super(TwoWayParentalMeanOptimalContributionSubsetSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            decn_space_xmap = xmap,
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
        # order dependent assignments
        self.pmebv = pmebv
        self.xmap = xmap
        self.xcontrib = xcontrib
        self.C = C

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Encode a candidate solution for the given Problem into an ``l`` 
        dimensional latent evaluation space.
        
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
            A matrix of shape (1+t,).

            The first index in the array is the mean genomic relationship 
            (a minimizing objective):

            .. math::
                MGR = || \\textbf{C} \\textbf{(sel)} ||_2

            The next `t` indices in the array are the sum of breeding values for 
            each of ``t`` traits for the selection (all maximizing objectives).

            Where:

            - ``t`` is the number of traits.
        """
        # FIXME
        # calculate individual contribution
        # scalar
        indcontrib = 1.0 / len(x)

        # get cross contributions
        # (s,d)[ndecn,:] -> (ndecn,d)
        xcontrib = self._xcontrib[x,:]

        # get cross map
        # (s,d)[ndecn,:] -> (ndecn,d)
        xmap = self._xmap[x,:]

        # get cross weights
        # (s,1) * (s,d) -> (s,d)
        # scalar * (s,d)[ndecn,:] -> (ndecn,d)
        xweight = indcontrib * xcontrib

        # bin based on their parents
        # (n,)
        xbins = numpy.bincount(
            xmap.flatten(), 
            weights = xweight.flatten(), 
            minlength = len(self._C),
        )

        # (n,n) @ (n,) -> (n,)
        CQx = self._C.dot(xbins)

        # calculate mean genomic relationship
        # norm2( (n,), keepdims=True ) -> (1,)
        mgr = numpy.linalg.norm(CQx, ord = 2, keepdims = True)

        # calculate negative mean breeding value of the selection
        # (n,t)[(ndecn,),:] -> (ndecn,t)
        # (ndecn,t).sum(0) -> (t,)
        gain = -indcontrib * self._pmebv[x,:].sum(0)

        # concatenate everything
        # (1,) concat (t,) -> (1+t,)
        out = numpy.concatenate([mgr,gain])

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_bvmat_gmat_xmap_xcontrib(
            cls,
            bvmat: BreedingValueMatrix,
            gmat: GenotypeMatrix, 
            cmatfcty: CoancestryMatrixFactory,
            xmap: numpy.ndarray,
            xcontrib: numpy.ndarray,
            unscale: bool,
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
        ) -> "TwoWayParentalMeanOptimalContributionSubsetSelectionProblem":
        # calculate estimated breeding values and relationships
        pmebv = cls._calc_pmebv(bvmat, xmap, unscale)
        C = cls._calc_C(gmat, cmatfcty)

        # construct class
        out = cls(
            pmebv = pmebv,
            xmap = xmap,
            xcontrib = xcontrib,
            C = C,
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

    @classmethod
    def from_bvmat_gmat(
            cls,
            bvmat: BreedingValueMatrix,
            gmat: GenotypeMatrix, 
            cmatfcty: CoancestryMatrixFactory,
            symab: bool,
            mateab: str,
            unscale: bool,
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
        ) -> "TwoWayParentalMeanOptimalContributionSubsetSelectionProblem":
        # calculate estimated breeding values and relationships
        ntaxa = bvmat.ntaxa
        xmap = cls._calc_xmap(ntaxa, symab, mateab)
        xcontrib = cls._calc_xcontrib(ntaxa, symab, mateab)
        pmebv = cls._calc_pmebv(bvmat, xmap, unscale)
        C = cls._calc_C(gmat, cmatfcty)

        # construct class
        out = cls(
            pmebv = pmebv,
            xmap = xmap,
            xcontrib = xcontrib,
            C = C,
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
