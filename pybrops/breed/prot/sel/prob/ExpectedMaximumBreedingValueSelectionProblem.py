"""
Module implementing Expected Maximum Breeding Value (EMBV) Selection problems.
"""

__all__ = [
    "ExpectedMaximumBreedingValueSubsetSelectionProblem",
    "ExpectedMaximumBreedingValueRealSelectionProblem",
    "ExpectedMaximumBreedingValueIntegerSelectionProblem",
    "ExpectedMaximumBreedingValueBinarySelectionProblem"
]

from abc import ABCMeta
from numbers import Integral, Real
from typing import Callable, Optional, Union
import numpy
from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol, check_is_MatingProtocol
from pybrops.breed.prot.sel.prob.BinarySelectionProblem import BinarySelectionProblem
from pybrops.breed.prot.sel.prob.IntegerSelectionProblem import IntegerSelectionProblem
from pybrops.breed.prot.sel.prob.RealSelectionProblem import RealSelectionProblem
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral, check_is_bool
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len_eq, check_ndarray_axis_len_gteq, check_ndarray_ndim
from pybrops.core.util.arrayix import triudix, triuix
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix


class ExpectedMaximumBreedingValueSelectionProblemMixin(metaclass=ABCMeta):
    """Helper class to implement properties common to EMBV selection problems."""

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################

    ############## Number of latent variables ##############
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in EMBV matrix
        return self._embv.shape[1]

    ##################### EMBV matrix ######################
    @property
    def embv(self) -> numpy.ndarray:
        """Expected maximum breeding value matrix of shape ``(s,t)``."""
        return self._embv
    @embv.setter
    def embv(self, value: numpy.ndarray) -> None:
        """Set expected maximum breeding value matrix."""
        check_is_ndarray(value, "embv")
        check_ndarray_ndim(value, "embv", 2)
        # most (binary, real, integer) problems require decisons for each cross
        check_ndarray_axis_len_eq(value, "embv", 0, self.ndecn)
        self._embv = value

    ########################## Private Object Methods ##########################
    @staticmethod
    def _calc_xmap(ntaxa, nparent, unique_parents = True):
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
    def _calc_embv(
            nparent: int,
            nmating: int,
            nprogeny: int,
            nrep: int,
            unique_parents: bool,
            pgmat: PhasedGenotypeMatrix, 
            gpmod: GenomicModel, 
            mateprot: MatingProtocol
        ) -> numpy.ndarray:
        # calculate cross map for our genotype matrix
        # (s,d)
        xmap = ExpectedMaximumBreedingValueSelectionProblemMixin._calc_xmap(
            pgmat.ntaxa,
            nparent,
            unique_parents
        )

        # allocate matrix for output EMBVs
        # (s,t)
        embv = numpy.empty((xmap.shape[0],gpmod.ntrait), dtype = float)

        # for each cross configuration
        # (d,)
        for i,xconfig in enumerate(xmap):
            # variable for tracking the average
            # (t,)
            avg = 0

            # run progeny simulations
            for i in range(nrep):
                # create progeny
                progeny = mateprot.mate(
                    pgmat = pgmat,
                    xconfig = xconfig,
                    nmating = nmating,
                    nprogeny = nprogeny,
                    miscout = None
                )

                # predict progeny breeding values
                # (nprogeny,t)
                bvmat = gpmod.gebv(progeny)

                # find max trait values and add to avg
                # (nprogeny,t).max(0) -> (t,)
                # (nprogeny,).max(0) -> scalar
                avg = avg + bvmat.tmax(True)

            # divide by the number of replicates
            # (t,)
            avg = avg / nrep

            embv[i,:] = avg

        return embv

class ExpectedMaximumBreedingValueSubsetSelectionProblem(ExpectedMaximumBreedingValueSelectionProblemMixin,SubsetSelectionProblem):
    """
    Class representing Expected Maximum Breeding Value (EMBV) selection problems in subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            embv: numpy.ndarray,
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
        Constructor for ExpectedMaximumBreedingValueSubsetSelectionProblem

        Parameters
        ----------
        embv : numpy.ndarray
            An expected maximum breeding value matrix of shape ``(s,t)``.

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
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        super(ExpectedMaximumBreedingValueSubsetSelectionProblem, self).__init__(
            embv = embv,
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
        self.embv = embv

    ############################ Object Properties #############################

    ##################### EMBV matrix ######################
    @ExpectedMaximumBreedingValueSelectionProblemMixin.embv.setter
    def embv(self, value: numpy.ndarray) -> None:
        """Set expected maximum breeding value matrix."""
        check_is_ndarray(value, "embv")
        check_ndarray_ndim(value, "embv", 2)
        # for subset problems, must have more crosses than decision variables
        check_ndarray_axis_len_gteq(value, "embv", 0, self.ndecn)
        self._embv = value

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a subset of individuals based on Expected Maximum Breeding Value
        (EMBV) selection (Müller et al., 2018). 
        
        Scoring for EMBV selection is defined as the negative mean of EMBVs for 
        a selected subset. A smaller value indicates a better EMBV score.

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
        out = -(1.0 / len(x)) * (self._embv[x,:].sum(0))

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_pgmat_gpmod(
            cls,
            nparent: Integral,
            nmating: Integral,
            nprogeny: Integral,
            nrep: Integral,
            unique_parents: bool,
            pgmat: PhasedGenotypeMatrix, 
            gpmod: GenomicModel, 
            mateprot: MatingProtocol,
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
        ) -> "ExpectedMaximumBreedingValueSubsetSelectionProblem":
        # check input types
        check_is_Integral(nparent, "nparent")
        check_is_Integral(nmating, "nmating")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nrep, "nrep")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        check_is_MatingProtocol(mateprot, "mateprot")

        # calculate estimated maximum breeding values
        embv = cls._calc_embv(nparent, nmating, nprogeny, nrep, unique_parents, pgmat, gpmod, mateprot)

        # construct class
        out = cls(
            embv = embv,
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

class ExpectedMaximumBreedingValueRealSelectionProblem(ExpectedMaximumBreedingValueSelectionProblemMixin,RealSelectionProblem):
    """
    Class representing Expected Maximum Breeding Value (EMBV) selection problems in real search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            embv: numpy.ndarray,
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
        Constructor for RealExpectedMaximumBreedingValueSelectionProblem

        Parameters
        ----------
        embv : numpy.ndarray
            An expected maximum breeding value matrix of shape ``(s,t)``.

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
            Additional keyword arguments passed to the parent class (RealSelectionProblem) constructor.
        """
        super(ExpectedMaximumBreedingValueRealSelectionProblem, self).__init__(
            embv = embv,
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
        self.embv = embv

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a subset of individuals based on Expected Maximum Breeding Value
        (EMBV) selection (Müller et al., 2018). 
        
        Scoring for EMBV selection is defined as the negative mean of EMBVs for 
        a selected subset. A smaller value indicates a better EMBV score.

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
        out = -contrib.dot(self._embv)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_pgmat_gpmod(
            cls,
            nparent: Integral,
            nmating: Integral,
            nprogeny: Integral,
            nrep: Integral,
            unique_parents: bool,
            pgmat: PhasedGenotypeMatrix, 
            gpmod: GenomicModel, 
            mateprot: MatingProtocol,
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
        ) -> "ExpectedMaximumBreedingValueRealSelectionProblem":
        # check input types
        check_is_Integral(nparent, "nparent")
        check_is_Integral(nmating, "nmating")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nrep, "nrep")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        check_is_MatingProtocol(mateprot, "mateprot")

        # calculate estimated maximum breeding values
        embv = cls._calc_embv(nparent, nmating, nprogeny, nrep, unique_parents, pgmat, gpmod, mateprot)

        # construct class
        out = cls(
            embv = embv,
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

class ExpectedMaximumBreedingValueIntegerSelectionProblem(ExpectedMaximumBreedingValueSelectionProblemMixin,IntegerSelectionProblem):
    """
    Class representing Expected Maximum Breeding Value (EMBV) selection problems in integer search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            embv: numpy.ndarray,
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
        Constructor for IntegerExpectedMaximumBreedingValueSelectionProblem

        Parameters
        ----------
        embv : numpy.ndarray
            An expected maximum breeding value matrix of shape ``(s,t)``.

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
            Additional keyword arguments passed to the parent class (IntegerSelectionProblem) constructor.
        """
        super(ExpectedMaximumBreedingValueIntegerSelectionProblem, self).__init__(
            embv = embv,
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
        self.embv = embv

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a subset of individuals based on Expected Maximum Breeding Value
        (EMBV) selection (Müller et al., 2018). 
        
        Scoring for EMBV selection is defined as the negative mean of EMBVs for 
        a selected subset. A smaller value indicates a better EMBV score.

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
        out = -contrib.dot(self._embv)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_pgmat_gpmod(
            cls,
            nparent: Integral,
            nmating: Integral,
            nprogeny: Integral,
            nrep: Integral,
            unique_parents: bool,
            pgmat: PhasedGenotypeMatrix, 
            gpmod: GenomicModel, 
            mateprot: MatingProtocol,
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
        ) -> "ExpectedMaximumBreedingValueIntegerSelectionProblem":
        # check input types
        check_is_Integral(nparent, "nparent")
        check_is_Integral(nmating, "nmating")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nrep, "nrep")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        check_is_MatingProtocol(mateprot, "mateprot")

        # calculate estimated maximum breeding values
        embv = cls._calc_embv(nparent, nmating, nprogeny, nrep, unique_parents, pgmat, gpmod, mateprot)

        # construct class
        out = cls(
            embv = embv,
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

class ExpectedMaximumBreedingValueBinarySelectionProblem(ExpectedMaximumBreedingValueSelectionProblemMixin,BinarySelectionProblem):
    """
    Class representing Expected Maximum Breeding Value (EMBV) selection problems in binary search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            embv: numpy.ndarray,
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
        Constructor for ExpectedMaximumBreedingValueBinarySelectionProblem

        Parameters
        ----------
        embv : numpy.ndarray
            An expected maximum breeding value matrix of shape ``(s,t)``.

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
            Additional keyword arguments passed to the parent class (IntegerSelectionProblem) constructor.
        """
        super(ExpectedMaximumBreedingValueBinarySelectionProblem, self).__init__(
            embv = embv,
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
        self.embv = embv

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a subset of individuals based on Expected Maximum Breeding Value
        (EMBV) selection (Müller et al., 2018). 
        
        Scoring for EMBV selection is defined as the negative mean of EMBVs for 
        a selected subset. A smaller value indicates a better EMBV score.

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
        out = -contrib.dot(self._embv)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_pgmat_gpmod(
            cls,
            nparent: Integral,
            nmating: Integral,
            nprogeny: Integral,
            nrep: Integral,
            unique_parents: bool,
            pgmat: PhasedGenotypeMatrix, 
            gpmod: GenomicModel, 
            mateprot: MatingProtocol,
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
        ) -> "ExpectedMaximumBreedingValueBinarySelectionProblem":
        # check input types
        check_is_Integral(nparent, "nparent")
        check_is_Integral(nmating, "nmating")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral(nrep, "nrep")
        check_is_bool(unique_parents, "unique_parents")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        check_is_MatingProtocol(mateprot, "mateprot")

        # calculate estimated maximum breeding values
        embv = cls._calc_embv(nparent, nmating, nprogeny, nrep, unique_parents, pgmat, gpmod, mateprot)

        # construct class
        out = cls(
            embv = embv,
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
