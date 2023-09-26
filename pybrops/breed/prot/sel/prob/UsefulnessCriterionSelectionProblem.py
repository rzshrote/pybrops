"""
Module implementing Usefulness Criterion (UC) Selection problems.
"""

from abc import ABCMeta, abstractmethod
from numbers import Integral, Real
from typing import Callable, Optional, Union
import numpy
import scipy.stats
from pybrops.breed.prot.sel.prob.BinaryMateSelectionProblem import BinaryMateSelectionProblem
from pybrops.breed.prot.sel.prob.IntegerMateSelectionProblem import IntegerMateSelectionProblem
from pybrops.breed.prot.sel.prob.RealMateSelectionProblem import RealMateSelectionProblem
from pybrops.breed.prot.sel.prob.SubsetMateSelectionProblem import SubsetMateSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral, check_is_Real, check_is_bool
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len_eq, check_ndarray_axis_len_gteq, check_ndarray_ndim
from pybrops.core.error.error_value_python import check_is_in_interval_inclusive
from pybrops.core.util.arrayix import triudix, triuix
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory, check_is_GeneticVarianceMatrixFactory
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction, check_is_GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix


class UsefulnessCriterionSelectionProblemMixin(metaclass=ABCMeta):
    """Mixin class containing properties common to UC selection problems."""

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################

    ############## Number of latent variables ##############
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in EMBV matrix
        return self._ucmat.shape[1]

    ##################### EMBV matrix ######################
    @property
    def ucmat(self) -> numpy.ndarray:
        """Usefulness criterion matrix of shape ``(s,t)``."""
        return self._ucmat
    @ucmat.setter
    def ucmat(self, value: numpy.ndarray) -> None:
        """Set usefulness criterion matrix."""
        check_is_ndarray(value, "ucmat")
        check_ndarray_ndim(value, "ucmat", 2)
        # most (binary, real, integer) problems require decisons for each cross
        check_ndarray_axis_len_eq(value, "ucmat", 0, self.ndecn)
        self._ucmat = value

    ######################### Private Object Methods ###########################
    @staticmethod
    def _calc_xmap(
            ntaxa: Integral, 
            nparent: Integral, 
            unique_parents: bool = True
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
        if unique_parents:
            return numpy.array(list(triudix(ntaxa,nparent)))
        else:
            return numpy.array(list(triuix(ntaxa,nparent)))

    @staticmethod
    def _calc_uc(
            vmatfcty: GeneticVarianceMatrixFactory, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral, 
            gmapfn: GeneticMapFunction, 
            selection_intensity: Real,
            pgmat: PhasedGenotypeMatrix, 
            gmod: GenomicModel, 
            xmap: numpy.ndarray
        ) -> numpy.ndarray:
        # calculate breeding values
        bvmat_obj = gmod.gebv(pgmat)
        
        # calculate variance matrix
        vmat_obj = vmatfcty.from_gmod(
            gmod = gmod, 
            pgmat = pgmat, 
            ncross = ncross, 
            nprogeny = nprogeny, 
            nself = nself,
            gmapfn = gmapfn
        )

        # get expected genome contributions
        # (p,)
        epgc = numpy.array(vmat_obj.epgc)
        
        # extract decentered and unscaled
        bvmat = bvmat_obj.unscale()     # (n,t)

        # extract variances
        vmat = vmat_obj.mat             # (n,...,n,t)

        # allocate memory for usefulness criterion
        uc = numpy.empty((len(xmap),bvmat_obj.ntrait), dtype = float)

        # for each cross configuration
        for i,cconfig in enumerate(xmap):
            # take dot product with expected genome contributions to get progeny mean
            # (p,) . (p,t) -> (t,)
            pmean = epgc.dot(bvmat[cconfig,:])

            # extract the variance at (i,...,i,:)
            # (t,)
            pvar = vmat[tuple(cconfig) + (slice(None),)]

            # calculate the usefulness criterion
            # (t,) + (t,) -> (t,)
            uc[i,:] = pmean + selection_intensity * numpy.sqrt(pvar)
        
        return uc

    ############################## Class Methods ###############################
    @classmethod
    @abstractmethod
    def from_pgmat_gpmod_xmap(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool, 
            pgmat: PhasedGenotypeMatrix,
            gpmod: GenomicModel,
            xmap: numpy.ndarray,
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
        ) -> "UsefulnessCriterionSelectionProblemMixin":
        raise NotImplementedError("class method is abstract")

    @classmethod
    @abstractmethod
    def from_pgmat_gpmod(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool, 
            pgmat: PhasedGenotypeMatrix,
            gpmod: GenomicModel,
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
        ) -> "UsefulnessCriterionSelectionProblemMixin":
        raise NotImplementedError("class method is abstract")

class UsefulnessCriterionBinaryMateSelectionProblem(UsefulnessCriterionSelectionProblemMixin,BinaryMateSelectionProblem):
    """
    Class representing Usefulness Criterion (UC) selection problems in subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ucmat: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
            decn_space_xmap: numpy.ndarray,
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
        Constructor for UsefulnessCriterionBinarySelectionProblem

        Parameters
        ----------
        ucmat : numpy.ndarray
            An usefulness criterion matrix of shape ``(s,t)``.

            Where:

            - ``s`` is the number of cross candidates.
            - ``t`` is the number of traits.
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
            Additional keyword arguments passed to the parent class (BinaryMateSelectionProblem) constructor.
        """
        super(UsefulnessCriterionBinaryMateSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            decn_space_xmap = decn_space_xmap,
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
        self.ucmat = ucmat

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a subset of individuals based on Usefulness Criterion (UC) 
        selection. 
        
        Scoring for UC selection is defined as the negative mean of UC values 
        for a selected subset. A smaller value indicates a better UC score.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(s,) == (ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(ax)`` where ``a`` is any number.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            An EMBV selection matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # scale x to have a sum of 1 (contribution)
        # (s,) -> (s,)
        contrib = (1.0 / x.sum()) * x

        # select individuals and take the negative mean of their EMBVs
        # CGS calculation explanation
        # Step 1: (s,) . (s,t) -> (t,)  # take dot product with contributions
        out = -contrib.dot(self._ucmat)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_pgmat_gpmod_xmap(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool, 
            pgmat: PhasedGenotypeMatrix,
            gpmod: GenomicModel,
            xmap: numpy.ndarray,
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
        ) -> "UsefulnessCriterionBinaryMateSelectionProblem":
        # type checks
        check_is_Integral(nparent, "nparent")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nself, "nself")
        check_is_Real(upper_percentile, "upper_percentile")
        check_is_in_interval_inclusive(upper_percentile, "upper_percentile", 0.0, 1.0)
        check_is_GeneticVarianceMatrixFactory(vmatfcty, "vmatfcty")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        check_is_ndarray(xmap, "xmap")
        
        # convert percentile to selection intensity
        selection_intensity = scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - upper_percentile)) / upper_percentile
        
        # calculate usefulness criterion
        ucmat = cls._calc_uc(
            vmatfcty, 
            ncross, 
            nprogeny, 
            nself, 
            gmapfn, 
            selection_intensity, 
            pgmat, 
            gpmod, 
            xmap
        )

        # construct object
        out = cls(
            ucmat = ucmat,
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

        return out

    @classmethod
    def from_pgmat_gpmod(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool, 
            pgmat: PhasedGenotypeMatrix,
            gpmod: GenomicModel,
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
        ) -> "UsefulnessCriterionBinaryMateSelectionProblem":
        # type checks
        check_is_Integral(nparent, "nparent")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nself, "nself")
        check_is_Real(upper_percentile, "upper_percentile")
        check_is_in_interval_inclusive(upper_percentile, "upper_percentile", 0.0, 1.0)
        check_is_GeneticVarianceMatrixFactory(vmatfcty, "vmatfcty")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # convert percentile to selection intensity
        selection_intensity = scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - upper_percentile)) / upper_percentile
        
        # calculate cross map
        xmap = cls._calc_xmap(
            pgmat.ntaxa, 
            nparent, 
            unique_parents
        )
        
        # calculate usefulness criterion
        ucmat = cls._calc_uc(
            vmatfcty, 
            ncross, 
            nprogeny, 
            nself, 
            gmapfn, 
            selection_intensity, 
            pgmat, 
            gpmod, 
            xmap
        )

        # construct object
        out = cls(
            ucmat = ucmat,
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

        return out

class UsefulnessCriterionIntegerMateSelectionProblem(UsefulnessCriterionSelectionProblemMixin,IntegerMateSelectionProblem):
    """
    Class representing Usefulness Criterion (UC) selection problems in integer search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ucmat: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
            decn_space_xmap: numpy.ndarray,
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
        Constructor for UsefulnessCriterionIntegerSelectionProblem

        Parameters
        ----------
        ucmat : numpy.ndarray
            An usefulness criterion matrix of shape ``(s,t)``.

            Where:

            - ``s`` is the number of cross candidates.
            - ``t`` is the number of traits.
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
            Additional keyword arguments passed to the parent class (IntegerMateSelectionProblem) constructor.
        """
        super(UsefulnessCriterionIntegerMateSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            decn_space_xmap = decn_space_xmap,
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
        self.ucmat = ucmat

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a subset of individuals based on Usefulness Criterion (UC) 
        selection. 
        
        Scoring for UC selection is defined as the negative mean of UC values 
        for a selected subset. A smaller value indicates a better UC score.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(s,) == (ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(ax)`` where ``a`` is any number.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            An EMBV selection matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # scale x to have a sum of 1 (contribution)
        # (s,) -> (s,)
        contrib = (1.0 / x.sum()) * x

        # select individuals and take the negative mean of their EMBVs
        # CGS calculation explanation
        # Step 1: (s,) . (s,t) -> (t,)  # take dot product with contributions
        out = -contrib.dot(self._ucmat)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_pgmat_gpmod_xmap(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool, 
            pgmat: PhasedGenotypeMatrix,
            gpmod: GenomicModel,
            xmap: numpy.ndarray,
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
        ) -> "UsefulnessCriterionIntegerMateSelectionProblem":
        # type checks
        check_is_Integral(nparent, "nparent")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nself, "nself")
        check_is_Real(upper_percentile, "upper_percentile")
        check_is_in_interval_inclusive(upper_percentile, "upper_percentile", 0.0, 1.0)
        check_is_GeneticVarianceMatrixFactory(vmatfcty, "vmatfcty")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        check_is_ndarray(xmap, "xmap")
        
        # convert percentile to selection intensity
        selection_intensity = scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - upper_percentile)) / upper_percentile
        
        # calculate usefulness criterion
        ucmat = cls._calc_uc(
            vmatfcty, 
            ncross, 
            nprogeny, 
            nself, 
            gmapfn, 
            selection_intensity, 
            pgmat, 
            gpmod, 
            xmap
        )

        # construct object
        out = cls(
            ucmat = ucmat,
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

        return out

    @classmethod
    def from_pgmat_gpmod(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool, 
            pgmat: PhasedGenotypeMatrix,
            gpmod: GenomicModel,
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
        ) -> "UsefulnessCriterionIntegerMateSelectionProblem":
        # type checks
        check_is_Integral(nparent, "nparent")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nself, "nself")
        check_is_Real(upper_percentile, "upper_percentile")
        check_is_in_interval_inclusive(upper_percentile, "upper_percentile", 0.0, 1.0)
        check_is_GeneticVarianceMatrixFactory(vmatfcty, "vmatfcty")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # convert percentile to selection intensity
        selection_intensity = scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - upper_percentile)) / upper_percentile
        
        # calculate cross map
        xmap = cls._calc_xmap(
            pgmat.ntaxa, 
            nparent, 
            unique_parents
        )
        
        # calculate usefulness criterion
        ucmat = cls._calc_uc(
            vmatfcty, 
            ncross, 
            nprogeny, 
            nself, 
            gmapfn, 
            selection_intensity, 
            pgmat, 
            gpmod, 
            xmap
        )

        # construct object
        out = cls(
            ucmat = ucmat,
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

        return out

class UsefulnessCriterionRealMateSelectionProblem(UsefulnessCriterionSelectionProblemMixin,RealMateSelectionProblem):
    """
    Class representing Usefulness Criterion (UC) selection problems in real search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ucmat: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
            decn_space_xmap: numpy.ndarray,
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
        Constructor for UsefulnessCriterionRealSelectionProblem

        Parameters
        ----------
        ucmat : numpy.ndarray
            An usefulness criterion matrix of shape ``(s,t)``.

            Where:

            - ``s`` is the number of cross candidates.
            - ``t`` is the number of traits.
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
            Additional keyword arguments passed to the parent class (RealMateSelectionProblem) constructor.
        """
        super(UsefulnessCriterionRealMateSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            decn_space_xmap = decn_space_xmap,
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
        self.ucmat = ucmat

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a subset of individuals based on Usefulness Criterion (UC) 
        selection. 
        
        Scoring for UC selection is defined as the negative mean of UC values 
        for a selected subset. A smaller value indicates a better UC score.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(s,) == (ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(ax)`` where ``a`` is any number.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            An EMBV selection matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # scale x to have a sum of 1 (contribution)
        # (s,) -> (s,)
        contrib = (1.0 / x.sum()) * x

        # select individuals and take the negative mean of their EMBVs
        # CGS calculation explanation
        # Step 1: (s,) . (s,t) -> (t,)  # take dot product with contributions
        out = -contrib.dot(self._ucmat)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_pgmat_gpmod_xmap(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool, 
            pgmat: PhasedGenotypeMatrix,
            gpmod: GenomicModel,
            xmap: numpy.ndarray,
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
        ) -> "UsefulnessCriterionRealMateSelectionProblem":
        # type checks
        check_is_Integral(nparent, "nparent")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nself, "nself")
        check_is_Real(upper_percentile, "upper_percentile")
        check_is_in_interval_inclusive(upper_percentile, "upper_percentile", 0.0, 1.0)
        check_is_GeneticVarianceMatrixFactory(vmatfcty, "vmatfcty")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        check_is_ndarray(xmap, "xmap")
        
        # convert percentile to selection intensity
        selection_intensity = scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - upper_percentile)) / upper_percentile
        
        # calculate usefulness criterion
        ucmat = cls._calc_uc(
            vmatfcty, 
            ncross, 
            nprogeny, 
            nself, 
            gmapfn, 
            selection_intensity, 
            pgmat, 
            gpmod, 
            xmap
        )

        # construct object
        out = cls(
            ucmat = ucmat,
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

        return out

    @classmethod
    def from_pgmat_gpmod(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool, 
            pgmat: PhasedGenotypeMatrix,
            gpmod: GenomicModel,
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
        ) -> "UsefulnessCriterionRealMateSelectionProblem":
        # type checks
        check_is_Integral(nparent, "nparent")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nself, "nself")
        check_is_Real(upper_percentile, "upper_percentile")
        check_is_in_interval_inclusive(upper_percentile, "upper_percentile", 0.0, 1.0)
        check_is_GeneticVarianceMatrixFactory(vmatfcty, "vmatfcty")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # convert percentile to selection intensity
        selection_intensity = scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - upper_percentile)) / upper_percentile
        
        # calculate cross map
        xmap = cls._calc_xmap(
            pgmat.ntaxa, 
            nparent, 
            unique_parents
        )
        
        # calculate usefulness criterion
        ucmat = cls._calc_uc(
            vmatfcty, 
            ncross, 
            nprogeny, 
            nself, 
            gmapfn, 
            selection_intensity, 
            pgmat, 
            gpmod, 
            xmap
        )

        # construct object
        out = cls(
            ucmat = ucmat,
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

        return out

class UsefulnessCriterionSubsetMateSelectionProblem(UsefulnessCriterionSelectionProblemMixin,SubsetMateSelectionProblem):
    """
    Class representing Usefulness Criterion (UC) selection problems in subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ucmat: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
            decn_space_xmap: numpy.ndarray,
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
        Constructor for UsefulnessCriterionSubsetSelectionProblem

        Parameters
        ----------
        ucmat : numpy.ndarray
            An usefulness criterion matrix of shape ``(s,t)``.

            Where:

            - ``s`` is the number of cross candidates.
            - ``t`` is the number of traits.
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
        super(UsefulnessCriterionSubsetMateSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            decn_space_xmap = decn_space_xmap,
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
        self.ucmat = ucmat

    ############################ Object Properties #############################

    ##################### EMBV matrix ######################
    @UsefulnessCriterionSelectionProblemMixin.ucmat.setter
    def ucmat(self, value: numpy.ndarray) -> None:
        """Set expected maximum breeding value matrix."""
        check_is_ndarray(value, "ucmat")
        check_ndarray_ndim(value, "ucmat", 2)
        # for subset problems, must have more crosses than decision variables
        check_ndarray_axis_len_gteq(value, "ucmat", 0, self.ndecn)
        self._ucmat = value

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a subset of individuals based on Usefulness Criterion (UC) 
        selection. 
        
        Scoring for UC selection is defined as the negative mean of UC values 
        for a selected subset. A smaller value indicates a better UC score.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(k,) == (ndecn,)``.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            An EMBV selection matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # select individuals and take the negative mean of their EMBVs
        # Step 1: (s,t)[(k,),:] -> (k,t)    # select individuals
        # Step 2: (k,t).sum(0)  -> (t,)     # sum across all individuals
        # Step 3: scalar * (t,) -> (t,)     # take mean across selection
        out = -(1.0 / len(x)) * (self._ucmat[x,:].sum(0))

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_pgmat_gpmod_xmap(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool, 
            pgmat: PhasedGenotypeMatrix,
            gpmod: GenomicModel,
            xmap: numpy.ndarray,
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
        ) -> "UsefulnessCriterionSubsetMateSelectionProblem":
        # type checks
        check_is_Integral(nparent, "nparent")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nself, "nself")
        check_is_Real(upper_percentile, "upper_percentile")
        check_is_in_interval_inclusive(upper_percentile, "upper_percentile", 0.0, 1.0)
        check_is_GeneticVarianceMatrixFactory(vmatfcty, "vmatfcty")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        check_is_ndarray(xmap, "xmap")
        
        # convert percentile to selection intensity
        selection_intensity = scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - upper_percentile)) / upper_percentile
        
        # calculate usefulness criterion
        ucmat = cls._calc_uc(
            vmatfcty, 
            ncross, 
            nprogeny, 
            nself, 
            gmapfn, 
            selection_intensity, 
            pgmat, 
            gpmod, 
            xmap
        )

        # construct object
        out = cls(
            ucmat = ucmat,
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

        return out

    @classmethod
    def from_pgmat_gpmod(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool, 
            pgmat: PhasedGenotypeMatrix,
            gpmod: GenomicModel,
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
        ) -> "UsefulnessCriterionSubsetMateSelectionProblem":
        # type checks
        check_is_Integral(nparent, "nparent")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nself, "nself")
        check_is_Real(upper_percentile, "upper_percentile")
        check_is_in_interval_inclusive(upper_percentile, "upper_percentile", 0.0, 1.0)
        check_is_GeneticVarianceMatrixFactory(vmatfcty, "vmatfcty")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # convert percentile to selection intensity
        selection_intensity = scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - upper_percentile)) / upper_percentile
        
        # calculate cross map
        xmap = cls._calc_xmap(
            pgmat.ntaxa, 
            nparent, 
            unique_parents
        )
        
        # calculate usefulness criterion
        ucmat = cls._calc_uc(
            vmatfcty, 
            ncross, 
            nprogeny, 
            nself, 
            gmapfn, 
            selection_intensity, 
            pgmat, 
            gpmod, 
            xmap
        )

        # construct object
        out = cls(
            ucmat = ucmat,
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

        return out

