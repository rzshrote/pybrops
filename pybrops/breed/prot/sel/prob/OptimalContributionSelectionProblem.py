"""
Module defining optimization problems for binary optimal constribution selection.
"""

__all__ = [
    "OptimalContributionSubsetSelectionProblem",
    "OptimalContributionRealSelectionProblem",
    "OptimalContributionIntegerSelectionProblem",
    "OptimalContributionBinarySelectionProblem",
]

from abc import ABCMeta
from abc import abstractmethod
from numbers import Integral
from numbers import Real
from typing import Callable
from typing import Optional
from typing import Union

import numpy
from pybrops.breed.prot.sel.prob.BinarySelectionProblem import BinarySelectionProblem
from pybrops.breed.prot.sel.prob.IntegerSelectionProblem import IntegerSelectionProblem
from pybrops.breed.prot.sel.prob.RealSelectionProblem import RealSelectionProblem
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_is_square
from pybrops.core.error.error_value_numpy import check_ndarray_is_triu
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix


class OptimalContributionSelectionProblemMixin(
        metaclass = ABCMeta,
    ):
    """Helper class containing common properties for Optimal Contribution Selection Problems."""

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in BV matrix plus 1
        return 1 + self._ebv.shape[1]

    @property
    def ebv(self) -> numpy.ndarray:
        """Breeding value matrix."""
        return self._ebv
    @ebv.setter
    def ebv(self, value: numpy.ndarray) -> None:
        """Set breeding value matrix."""
        check_is_ndarray(value, "ebv")
        check_ndarray_ndim(value, "ebv", 2)
        self._ebv = value

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
    def _calc_ebv(bvmat: BreedingValueMatrix, unscale = True) -> numpy.ndarray:
        """
        Construct a breeding value matrix for use by a Problem specification.

        Parameters
        ----------
        bvmat : BreedingValueMatrix
            Input breeding value matrix.
        
        Returns
        -------
        out : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``. May be scaled or unscaled.
        """
        # get breeding value matrix (n,t)
        out = bvmat.unscale() if unscale else bvmat.mat

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
        ) -> "OptimalContributionSelectionProblemMixin":
        raise NotImplementedError("class method is abstract")

class OptimalContributionSubsetSelectionProblem(
        OptimalContributionSelectionProblemMixin,
        SubsetSelectionProblem,
    ):
    """
    Class representing an Optimal Contribution Selection Problem for subset
    search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ebv: numpy.ndarray,
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
        Constructor for SubsetOptimalContributionSelectionProblem.
        
        Parameters
        ----------
        ebv : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``. 
            If you are using a penalization transformation function, preferably
            these breeding values are centered and scaled to make the penalies 
            less extreme.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
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
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        # call SubsetSelectionProblem constructor
        super(OptimalContributionSubsetSelectionProblem, self).__init__(
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
        # order dependent assignments
        self.ebv = ebv
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
        # calculate individual contribution
        # scalar
        indcontrib = 1.0 / len(x)

        # calculate MEH
        # (n,n)[:,(k,)] -> (n,k)
        # scalar * (n,k).sum(1) -> (n,)
        Cx = indcontrib * self.C[:,x].sum(1)

        # calculate mean genomic relationship
        # norm2( (n,), keepdims=True ) -> (1,)
        mgr = numpy.linalg.norm(Cx, ord = 2, keepdims = True)

        # calculate negative mean breeding value of the selection
        # (n,t)[(k,),:] -> (k,t)
        # (k,t).sum(0) -> (t,)
        gain = -indcontrib * self.ebv[x,:].sum(0)
        
        # concatenate everything
        # (1,) concat (t,) -> (1+t,)
        out = numpy.concatenate([mgr,gain])

        # return (1+t,)
        return out

    ############################## Class Methods ###############################
    @classmethod
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
        ) -> "OptimalContributionSubsetSelectionProblem":
        # calculate estimated breeding values and relationships
        ebv = cls._calc_ebv(bvmat, unscale)
        C = cls._calc_C(gmat, cmatfcty)

        # construct class
        out = cls(
            ebv = ebv,
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

class OptimalContributionRealSelectionProblem(
        OptimalContributionSelectionProblemMixin,
        RealSelectionProblem,
    ):
    """
    Class representing an Optimal Contribution Selection Problem for real
    search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ebv: numpy.ndarray,
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
        Constructor for RealOptimalContributionSelectionProblem.
        
        Parameters
        ----------
        ebv : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``. 
            If you are using a penalization transformation function, preferably
            these breeding values are centered and scaled to make the penalies 
            less extreme.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
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
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        # call SubsetSelectionProblem constructor
        super(OptimalContributionRealSelectionProblem, self).__init__(
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
        # order dependent assignments
        self.ebv = ebv
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
            A candidate solution vector of shape ``(ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(kx)`` where ``k`` is any number.
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
        # calculate sum(x)
        xsum = x.sum()

        # if sum(x) ~== 0, then set to 1
        xsum = xsum if abs(xsum) >= 1e-10 else 1.0

        # scale x to have a sum of 1 (contribution)
        # (n,) -> (n,)
        contrib = (1.0 / xsum) * x

        # calculate mean genomic contribution
        # (n,n) . (n,) -> (n,)
        # scalar * (n,) -> (n,)
        # norm2( (n,), keepdims=True ) -> (1,)
        mgc = numpy.linalg.norm(self.C.dot(contrib), ord = 2, keepdims = True)

        # calculate negative mean breeding value of the selection
        # (n,) . (n,t) -> (t,)
        gain = -contrib.dot(self._ebv)
        
        # concatenate everything
        # (1,) concat (t,) -> (1+t,)
        out = numpy.concatenate([mgc,gain])

        # return (1+t,)
        return out

    ############################## Class Methods ###############################
    @classmethod
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
        ) -> "OptimalContributionRealSelectionProblem":
        # calculate estimated breeding values and relationships
        ebv = cls._calc_ebv(bvmat, unscale)
        C = cls._calc_C(gmat, cmatfcty)

        # construct class
        out = cls(
            ebv = ebv,
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

class OptimalContributionIntegerSelectionProblem(
        OptimalContributionSelectionProblemMixin,
        IntegerSelectionProblem,
    ):
    """
    Class representing an Optimal Contribution Selection Problem for integer
    search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ebv: numpy.ndarray,
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
        Constructor for IntegerOptimalContributionSelectionProblem.
        
        Parameters
        ----------
        ebv : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``. 
            If you are using a penalization transformation function, preferably
            these breeding values are centered and scaled to make the penalies 
            less extreme.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
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
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        # call SubsetSelectionProblem constructor
        super(OptimalContributionIntegerSelectionProblem, self).__init__(
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
        # order dependent assignments
        self.ebv = ebv
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
            A candidate solution vector of shape ``(ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(kx)`` where ``k`` is any number.
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
        # calculate sum(x)
        xsum = x.sum()

        # if sum(x) ~== 0, then set to 1
        xsum = xsum if abs(xsum) >= 1e-10 else 1.0

        # scale x to have a sum of 1 (contribution)
        # (n,) -> (n,)
        contrib = (1.0 / xsum) * x

        # calculate mean genomic contribution
        # (n,n) . (n,) -> (n,)
        # scalar * (n,) -> (n,)
        # norm2( (n,), keepdims=True ) -> (1,)
        mgc = numpy.linalg.norm(self.C.dot(contrib), ord = 2, keepdims = True)

        # calculate negative mean breeding value of the selection
        # (n,) . (n,t) -> (t,)
        gain = -contrib.dot(self._ebv)
        
        # concatenate everything
        # (1,) concat (t,) -> (1+t,)
        out = numpy.concatenate([mgc,gain])

        # return (1+t,)
        return out

    ############################## Class Methods ###############################
    @classmethod
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
        ) -> "OptimalContributionIntegerSelectionProblem":
        # calculate estimated breeding values and relationships
        ebv = cls._calc_ebv(bvmat, unscale)
        C = cls._calc_C(gmat, cmatfcty)

        # construct class
        out = cls(
            ebv = ebv,
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

class OptimalContributionBinarySelectionProblem(
        OptimalContributionSelectionProblemMixin,
        BinarySelectionProblem,
    ):
    """
    Class representing an Optimal Contribution Selection Problem for integer
    search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ebv: numpy.ndarray,
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
        Constructor for BinaryOptimalContributionSelectionProblem.
        
        Parameters
        ----------
        ebv : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``. 
            If you are using a penalization transformation function, preferably
            these breeding values are centered and scaled to make the penalies 
            less extreme.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
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
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        # call SubsetSelectionProblem constructor
        super(OptimalContributionBinarySelectionProblem, self).__init__(
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
        # order dependent assignments
        self.ebv = ebv
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
            A candidate solution vector of shape ``(ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(kx)`` where ``k`` is any number.
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
        # calculate sum(x)
        xsum = x.sum()

        # if sum(x) ~== 0, then set to 1
        xsum = xsum if abs(xsum) >= 1e-10 else 1.0

        # scale x to have a sum of 1 (contribution)
        # (n,) -> (n,)
        contrib = (1.0 / xsum) * x

        # calculate mean genomic contribution
        # (n,n) . (n,) -> (n,)
        # scalar * (n,) -> (n,)
        # norm2( (n,), keepdims=True ) -> (1,)
        mgc = numpy.linalg.norm(self.C.dot(contrib), ord = 2, keepdims = True)

        # calculate negative mean breeding value of the selection
        # (n,) . (n,t) -> (t,)
        gain = -contrib.dot(self._ebv)
        
        # concatenate everything
        # (1,) concat (t,) -> (1+t,)
        out = numpy.concatenate([mgc,gain])

        # return (1+t,)
        return out

    ############################## Class Methods ###############################
    @classmethod
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
        ) -> "OptimalContributionBinarySelectionProblem":
        # calculate estimated breeding values and relationships
        ebv = cls._calc_ebv(bvmat, unscale)
        C = cls._calc_C(gmat, cmatfcty)

        # construct class
        out = cls(
            ebv = ebv,
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
