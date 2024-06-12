from abc import ABCMeta, abstractmethod
from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
from pybrops.breed.prot.sel.prob.SubsetMateSelectionProblem import SubsetMateSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_bool, check_is_str
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len_gteq, check_ndarray_ndim
from pybrops.core.util.crossix import twowayix, twowayix_len
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.model.pmgebvmat.DenseTwoWayParentalMeanGenomicEstimatedBreedingValueMatrix import DenseTwoWayParentalMeanGenomicEstimatedBreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix, check_is_GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix


class TwoWayParentalMeanGenomicEstimatedBreedingValueProblemMixin(
        metaclass = ABCMeta,
    ):
    """Mixin class containing properties common to pmGEBV selection problems."""

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################

    ############## Number of latent variables ##############
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in EMBV matrix
        return self._pmgebvmat.shape[1]

    ##################### EMBV matrix ######################
    @property
    def pmgebvmat(self) -> numpy.ndarray:
        """Usefulness criterion matrix of shape ``(s,t)``."""
        return self._pmgebvmat
    @pmgebvmat.setter
    def pmgebvmat(self, value: numpy.ndarray) -> None:
        """Set usefulness criterion matrix."""
        check_is_ndarray(value, "pmgebvmat")
        check_ndarray_ndim(value, "pmgebvmat", 2)
        # most (binary, real, integer) problems require decisons for each cross
        check_ndarray_axis_len_gteq(value, "pmgebvmat", 0, self.ndecn)
        self._pmgebvmat = value

    ######################### Private Object Methods ###########################
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

    ############################## Class Methods ###############################
    @classmethod
    @abstractmethod
    def from_gmat_gpmod_xmap(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            gmat: PhasedGenotypeMatrix,
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
        ) -> "TwoWayParentalMeanGenomicEstimatedBreedingValueProblemMixin":
        raise NotImplementedError("class method is abstract")

    @classmethod
    @abstractmethod
    def from_gmat_gpmod(
            cls,
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            gmat: PhasedGenotypeMatrix,
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
        ) -> "TwoWayParentalMeanGenomicEstimatedBreedingValueProblemMixin":
        raise NotImplementedError("class method is abstract")

class TwoWayParentalMeanGenomicEstimatedBreedingValueSubsetMateSelectionProblem(
        TwoWayParentalMeanGenomicEstimatedBreedingValueProblemMixin,
        SubsetMateSelectionProblem,
    ):
    ########################## Special Object Methods ##########################
    def __init__(
            self,
            pmgebvmat: numpy.ndarray,
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
        pmgebvmat : numpy.ndarray
            An usefulness criterion matrix of shape ``(s,t)``.

            Where:

            - ``s`` is the number of cross candidates.
            - ``t`` is the number of traits.
        ndecn : Integral
            Number of decision variables.
        decn_space : numpy.ndarray, None
            An array of shape ``(2,ndecn)`` defining the decision space.
            If None, do not set a decision space.
        decn_space_lower : numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing lower limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a lower limit for the decision variables.
        decn_space_upper : numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing upper limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a upper limit for the decision variables.
        decn_space_xmap : numpy.ndarray
            Cross map corresponding to the decision space.
        nobj : Integral
            Number of objectives.
        obj_wt : numpy.ndarray
            Objective function weights.
        obj_trans : Callable, None
            A transformation function transforming a latent space vector to an objective space vector.
            The transformation function must be of the form: ``obj_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the identity transformation function: copy the latent space vector to the objective space vector.
        obj_trans_kwargs : dict, None
            Keyword arguments for the latent space to objective space transformation function.
            If None, an empty dictionary is used.
        nineqcv : Integral,
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
        super(TwoWayParentalMeanGenomicEstimatedBreedingValueSubsetMateSelectionProblem, self).__init__(
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
        self.pmgebvmat = pmgebvmat

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
        out = -(1.0 / len(x)) * (self._pmgebvmat[x,:].sum(0))

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_gmat_gpmod_xmap(
            cls,
            gmat: GenotypeMatrix,
            gpmod: GenomicModel,
            xmap: numpy.ndarray,
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
        ) -> "TwoWayParentalMeanGenomicEstimatedBreedingValueProblemMixin":
        # type checks
        check_is_GenotypeMatrix(gmat, "gmat")
        check_is_GenomicModel(gpmod, "gpmod")
        check_is_ndarray(xmap, "xmap")

        # calculate pmGEBV matrix
        tmp = DenseTwoWayParentalMeanGenomicEstimatedBreedingValueMatrix.from_gmod(
            gmod = gpmod,
            gmat = gmat,
        )

        # get matrix
        # (n,n,t)
        mat = tmp.unscale(False) if unscale else tmp.mat

        # get female, male indicies: (s,2) -> ((s,), (s,))
        female = xmap[:,0]
        male = xmap[:,1]

        # extract and flatten selected crosses
        # (n,n,t)[(s,),(s,),:] -> (s,t)
        pmgebvmat = mat[female,male,:]

        # construct object
        out = cls(
            pmgebvmat = pmgebvmat,
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
    def from_gmat_gpmod(
            cls,
            gmat: GenotypeMatrix,
            gpmod: GenomicModel,
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
        ) -> "TwoWayParentalMeanGenomicEstimatedBreedingValueProblemMixin":
        # type checks
        check_is_bool(symab, "symab")
        check_is_str(mateab, "mateab")
        check_is_PhasedGenotypeMatrix(gmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # calculate cross map
        xmap = cls._calc_xmap(
            ntaxa = gmat.ntaxa,
            symab = symab,
            mateab = mateab,
        )
        
        # construct object
        out = cls.from_gmat_gpmod_xmap(
            gmat = gmat,
            gpmod = gpmod,
            xmap = xmap,
            unscale = unscale,
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

