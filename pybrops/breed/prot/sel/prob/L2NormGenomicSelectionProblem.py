"""
Module implementing L2-norm Genomic Selection (L2GS) problems for multiple search space types.
"""

__all__ = [
    "L2NormGenomicBinarySelectionProblem",
    "L2NormGenomicIntegerSelectionProblem",
    "L2NormGenomicRealSelectionProblem",
    "L2NormGenomicSubsetSelectionProblem",
]

from abc import ABCMeta
from abc import abstractmethod
from numbers import Integral
from numbers import Real
from typing import Callable
from typing import Optional
from typing import Union
import numpy
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pybrops.breed.prot.sel.prob.BinarySelectionProblem import BinarySelectionProblem
from pybrops.breed.prot.sel.prob.IntegerSelectionProblem import IntegerSelectionProblem
from pybrops.breed.prot.sel.prob.RealSelectionProblem import RealSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_is_square
from pybrops.core.error.error_value_numpy import check_ndarray_is_triu
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import check_is_CoancestryMatrixFactory
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.GenotypeMatrix import check_is_GenotypeMatrix


class L2NormGenomicSelectionProblemMixin(
        metaclass = ABCMeta,
    ):
    """Helper class containing common properties for L2GS problems."""

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################

    ############## Number of latent variables ##############
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in BV matrix plus 1
        return self._C.shape[0]

    ########## Cholesky decomposition of kinship ###########
    @property
    def C(self) -> numpy.ndarray:
        """Cholesky decomposition of the kinship matrix."""
        return self._C
    @C.setter
    def C(self, value: numpy.ndarray) -> None:
        """Set Cholesky decomposition of the kinship matrix."""
        check_is_ndarray(value, "C")
        check_ndarray_ndim(value, "C", 3)
        check_ndarray_is_square(value, "C")
        check_ndarray_is_triu(value, "C")
        self._C = value

    ############################## Class Methods ###############################
    @classmethod
    @abstractmethod
    def from_gmat(
            cls,
            gmat: GenotypeMatrix,
            cmatfcty: CoancestryMatrixFactory,
            mkrwt: numpy.ndarray,
            afreq: numpy.ndarray,
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
        ) -> "L2NormGenomicSelectionProblemMixin":
        """
        Construct an L2 Norm Genomic Selection Problem from a Genotype Matrix.

        Parameters
        ----------
        mkrwt : numpy.ndarray
            A marker weight matrix of shape ``(p,t)``.
        afreq : numpy.ndarray
            A marker target allele frequency matrix of shape ``(p,t)``.
        """
        raise NotImplementedError("class method is abstract")

class L2NormGenomicSubsetSelectionProblem(
        L2NormGenomicSelectionProblemMixin,
        SubsetSelectionProblem,
    ):
    """
    Class representing L2-norm Genomic Selection (L2GS) in subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
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
        Constructor for SubsetConventionalGenomicSelectionProblem.
        
        Parameters
        ----------
        C : numpy.ndarray
            An upper triangle matrix of shape ``(t,n,n)`` resulting from a Cholesky 
            decomposition of a distance relationship matrix: K = C'C.

            Where:

            - ``t`` is the number of traits.
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
        super(L2NormGenomicSubsetSelectionProblem, self).__init__(
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
        # assignments
        self.C = C

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on the L2-norm from a utopian 
        allele frequency, specified by a relationship matrix.

        This is a minimizing objective. A lower score means a smaller distance 
        to the utopian point.

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
            A distance metric matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # calculate vector
        # (t,n,n)[:,:,(k,)] -> (t,n,k)
        # (1/k) * (t,n,k).sum(2) -> (t,n)
        Cx = (1.0 / len(x)) * self._C[:,:,x].sum(2)

        # calculate distance to utopian point for each trait
        # norm2( (t,n) ) -> (t,)
        out = numpy.linalg.norm(Cx, ord = 2, axis = 1)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_gmat(
            cls,
            gmat: GenotypeMatrix,
            cmatfcty: CoancestryMatrixFactory,
            mkrwt: numpy.ndarray,
            afreq: numpy.ndarray,
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
        ) -> "L2NormGenomicSubsetSelectionProblem":
        """
        Construct an L2 Norm Genomic Selection Problem from a Genotype Matrix.

        Parameters
        ----------
        mkrwt : numpy.ndarray
            A marker weight matrix of shape ``(p,t)``.
        afreq : numpy.ndarray
            A marker target allele frequency matrix of shape ``(p,t)``.
        """
        # type checks
        check_is_GenotypeMatrix(gmat, "gmat")
        check_is_CoancestryMatrixFactory(cmatfcty, "cmatfcty")
        check_is_ndarray(mkrwt, "mkrwt")
        check_is_ndarray(afreq, "afreq")

        # get shapes
        ntrait = mkrwt.shape[1]
        ntaxa = gmat.ntaxa

        # allocate memory for cholesky decomposition
        Ctensor = numpy.empty((ntrait,ntaxa,ntaxa), dtype=float)

        # for each trait, calculate cholesky decomposition and store
        for i in range(ntrait):
            # calculate coancestry
            G = cmatfcty.from_gmat(gmat, mkrwt = mkrwt, afreq = afreq)

            # to ensure we're able to perform cholesky decomposition, apply jitter if needed.
            # if we are unable to fix, then raise value error
            if not G.apply_jitter():
                raise ValueError(
                    "Unable to construct objective function: Kinship matrix is not positive definite.\n"+
                    "    This could be caused by lack of genetic diversity.\n"
                )

            K = G.mat_asformat("kinship")       # convert G to (1/2)G (kinship analogue): (n,n)
            C = numpy.linalg.cholesky(K).T      # cholesky decomposition of K matrix: (n,n)

            # store cholesky decomposition in tensor
            Ctensor[i,:,:] = C

        # construct class
        out = cls(
            C = Ctensor,
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

class L2NormGenomicRealSelectionProblem(
        L2NormGenomicSelectionProblemMixin,
        RealSelectionProblem,
    ):
    """
    Class representing L2-norm Genomic Selection (L2GS) in real search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
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
        Constructor for RealConventionalGenomicSelectionProblem.
        
        Parameters
        ----------
        C : numpy.ndarray
            An upper triangle matrix of shape ``(t,n,n)`` resulting from a Cholesky 
            decomposition of a distance relationship matrix: K = C'C.

            Where:

            - ``t`` is the number of traits.
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
            Additional keyword arguments passed to the parent class (RealSelectionProblem) constructor.
        """
        super(L2NormGenomicRealSelectionProblem, self).__init__(
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
        # assignments
        self.C = C

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on the L2-norm from a utopian 
        allele frequency, specified by a relationship matrix.

        This is a minimizing objective. A lower score means a smaller distance 
        to the utopian point.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(k,) == (ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(ax)`` where ``a`` is any number.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A distance metric matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # scale x to have a sum of 1 (contribution)
        contrib = (1.0 / x.sum()) * x

        # calculate distance
        # (t,n,n) . (n,) -> (t,n)
        # norm2( (t,n), axis=1 ) -> (t,)
        out = numpy.linalg.norm(self.C.dot(contrib), ord = 2, axis = 1)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_gmat(
            cls,
            gmat: GenotypeMatrix,
            cmatfcty: CoancestryMatrixFactory,
            mkrwt: numpy.ndarray,
            afreq: numpy.ndarray,
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
        ) -> "L2NormGenomicRealSelectionProblem":
        """
        Construct an L2 Norm Genomic Selection Problem from a Genotype Matrix.

        Parameters
        ----------
        mkrwt : numpy.ndarray
            A marker weight matrix of shape ``(p,t)``.
        afreq : numpy.ndarray
            A marker target allele frequency matrix of shape ``(p,t)``.
        """
        # type checks
        check_is_GenotypeMatrix(gmat, "gmat")
        check_is_CoancestryMatrixFactory(cmatfcty, "cmatfcty")
        check_is_ndarray(mkrwt, "mkrwt")
        check_is_ndarray(afreq, "afreq")

        # get shapes
        ntrait = mkrwt.shape[1]
        ntaxa = gmat.ntaxa

        # allocate memory for cholesky decomposition
        Ctensor = numpy.empty((ntrait,ntaxa,ntaxa), dtype=float)

        # for each trait, calculate cholesky decomposition and store
        for i in range(ntrait):
            # calculate coancestry
            G = cmatfcty.from_gmat(gmat, mkrwt = mkrwt, afreq = afreq)

            # to ensure we're able to perform cholesky decomposition, apply jitter if needed.
            # if we are unable to fix, then raise value error
            if not G.apply_jitter():
                raise ValueError(
                    "Unable to construct objective function: Kinship matrix is not positive definite.\n"+
                    "    This could be caused by lack of genetic diversity.\n"
                )

            K = G.mat_asformat("kinship")       # convert G to (1/2)G (kinship analogue): (n,n)
            C = numpy.linalg.cholesky(K).T      # cholesky decomposition of K matrix: (n,n)

            # store cholesky decomposition in tensor
            Ctensor[i,:,:] = C

        # construct class
        out = cls(
            C = Ctensor,
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

class L2NormGenomicIntegerSelectionProblem(
        L2NormGenomicSelectionProblemMixin,
        IntegerSelectionProblem,
    ):
    """
    Class representing L2-norm Genomic Selection (L2GS) in integer search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
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
        Constructor for IntegerConventionalGenomicSelectionProblem.
        
        Parameters
        ----------
        C : numpy.ndarray
            An upper triangle matrix of shape ``(t,n,n)`` resulting from a Cholesky 
            decomposition of a distance relationship matrix: K = C'C.

            Where:

            - ``t`` is the number of traits.
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
            Additional keyword arguments passed to the parent class (IntegerSelectionProblem) constructor.
        """
        super(L2NormGenomicIntegerSelectionProblem, self).__init__(
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
        # assignments
        self.C = C

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on the L2-norm from a utopian 
        allele frequency, specified by a relationship matrix.

        This is a minimizing objective. A lower score means a smaller distance 
        to the utopian point.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(k,) == (ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(ax)`` where ``a`` is any number.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A distance metric matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # scale x to have a sum of 1 (contribution)
        contrib = (1.0 / x.sum()) * x

        # calculate distance
        # (t,n,n) . (n,) -> (t,n)
        # norm2( (t,n), axis=1 ) -> (t,)
        out = numpy.linalg.norm(self.C.dot(contrib), ord = 2, axis = 1)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_gmat(
            cls,
            gmat: GenotypeMatrix,
            cmatfcty: CoancestryMatrixFactory,
            mkrwt: numpy.ndarray,
            afreq: numpy.ndarray,
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
        ) -> "L2NormGenomicIntegerSelectionProblem":
        """
        Construct an L2 Norm Genomic Selection Problem from a Genotype Matrix.

        Parameters
        ----------
        mkrwt : numpy.ndarray
            A marker weight matrix of shape ``(p,t)``.
        afreq : numpy.ndarray
            A marker target allele frequency matrix of shape ``(p,t)``.
        """
        # type checks
        check_is_GenotypeMatrix(gmat, "gmat")
        check_is_CoancestryMatrixFactory(cmatfcty, "cmatfcty")
        check_is_ndarray(mkrwt, "mkrwt")
        check_is_ndarray(afreq, "afreq")

        # get shapes
        ntrait = mkrwt.shape[1]
        ntaxa = gmat.ntaxa

        # allocate memory for cholesky decomposition
        Ctensor = numpy.empty((ntrait,ntaxa,ntaxa), dtype=float)

        # for each trait, calculate cholesky decomposition and store
        for i in range(ntrait):
            # calculate coancestry
            G = cmatfcty.from_gmat(gmat, mkrwt = mkrwt, afreq = afreq)

            # to ensure we're able to perform cholesky decomposition, apply jitter if needed.
            # if we are unable to fix, then raise value error
            if not G.apply_jitter():
                raise ValueError(
                    "Unable to construct objective function: Kinship matrix is not positive definite.\n"+
                    "    This could be caused by lack of genetic diversity.\n"
                )

            K = G.mat_asformat("kinship")       # convert G to (1/2)G (kinship analogue): (n,n)
            C = numpy.linalg.cholesky(K).T      # cholesky decomposition of K matrix: (n,n)

            # store cholesky decomposition in tensor
            Ctensor[i,:,:] = C

        # construct class
        out = cls(
            C = Ctensor,
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

class L2NormGenomicBinarySelectionProblem(
        L2NormGenomicSelectionProblemMixin,
        BinarySelectionProblem,
    ):
    """
    Class representing L2-norm Genomic Selection (L2GS) in binary search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
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
        Constructor for BinaryConventionalGenomicSelectionProblem.
        
        Parameters
        ----------
        C : numpy.ndarray
            An upper triangle matrix of shape ``(t,n,n)`` resulting from a Cholesky 
            decomposition of a distance relationship matrix: K = C'C.

            Where:

            - ``t`` is the number of traits.
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
            Additional keyword arguments passed to the parent class (BinarySelectionProblem) constructor.
        """
        super(L2NormGenomicBinarySelectionProblem, self).__init__(
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
        # assignments
        self.C = C

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on the L2-norm from a utopian 
        allele frequency, specified by a relationship matrix.

        This is a minimizing objective. A lower score means a smaller distance 
        to the utopian point.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(k,) == (ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(ax)`` where ``a`` is any number.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A distance metric matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # scale x to have a sum of 1 (contribution)
        contrib = (1.0 / x.sum()) * x

        # calculate distance
        # (t,n,n) . (n,) -> (t,n)
        # norm2( (t,n), axis=1 ) -> (t,)
        out = numpy.linalg.norm(self.C.dot(contrib), ord = 2, axis = 1)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_gmat(
            cls,
            gmat: GenotypeMatrix,
            cmatfcty: CoancestryMatrixFactory,
            mkrwt: numpy.ndarray,
            afreq: numpy.ndarray,
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
        ) -> "L2NormGenomicBinarySelectionProblem":
        """
        Construct an L2 Norm Genomic Selection Problem from a Genotype Matrix.

        Parameters
        ----------
        mkrwt : numpy.ndarray
            A marker weight matrix of shape ``(p,t)``.
        afreq : numpy.ndarray
            A marker target allele frequency matrix of shape ``(p,t)``.
        """
        # type checks
        check_is_GenotypeMatrix(gmat, "gmat")
        check_is_CoancestryMatrixFactory(cmatfcty, "cmatfcty")
        check_is_ndarray(mkrwt, "mkrwt")
        check_is_ndarray(afreq, "afreq")

        # get shapes
        ntrait = mkrwt.shape[1]
        ntaxa = gmat.ntaxa

        # allocate memory for cholesky decomposition
        Ctensor = numpy.empty((ntrait,ntaxa,ntaxa), dtype=float)

        # for each trait, calculate cholesky decomposition and store
        for i in range(ntrait):
            # calculate coancestry
            G = cmatfcty.from_gmat(gmat, mkrwt = mkrwt, afreq = afreq)

            # to ensure we're able to perform cholesky decomposition, apply jitter if needed.
            # if we are unable to fix, then raise value error
            if not G.apply_jitter():
                raise ValueError(
                    "Unable to construct objective function: Kinship matrix is not positive definite.\n"+
                    "    This could be caused by lack of genetic diversity.\n"
                )

            K = G.mat_asformat("kinship")       # convert G to (1/2)G (kinship analogue): (n,n)
            C = numpy.linalg.cholesky(K).T      # cholesky decomposition of K matrix: (n,n)

            # store cholesky decomposition in tensor
            Ctensor[i,:,:] = C

        # construct class
        out = cls(
            C = Ctensor,
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
