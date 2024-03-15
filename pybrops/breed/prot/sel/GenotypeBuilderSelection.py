"""
Module defining protocols for Genotype Builder (GB) selection.
"""

__all__ = [
    "GenotypeBuilderSelectionMixin",
    "GenotypeBuilderSubsetSelection",
]

from abc import ABCMeta
from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
from numpy.random import Generator, RandomState
import pandas

from pybrops.breed.prot.sel.SubsetSelectionProtocol import SubsetSelectionProtocol
from pybrops.breed.prot.sel.prob.GenotypeBuilderSelectionProblem import GenotypeBuilderSubsetSelectionProblem
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_python import check_is_gt, check_is_in_interval_inclusive
from pybrops.model.gmod.AdditiveLinearGenomicModel import check_is_AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.opt.algo.SubsetOptimizationAlgorithm import SubsetOptimizationAlgorithm
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix


class GenotypeBuilderSelectionMixin(
        metaclass = ABCMeta,
    ):
    """
    docstring for GenotypeBuilderBaseSelection.
    """

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################
    @property
    def ntrait(self) -> Integral:
        """Number of traits to expect from matrix inputs."""
        return self._ntrait
    @ntrait.setter
    def ntrait(self, value: Integral) -> None:
        """Set number of traits to expect."""
        check_is_Integral(value, "ntrait")
        check_is_gt(value, "ntrait", 0)
        self._ntrait = value

    @property
    def nhaploblk(self) -> Integral:
        """Number of haplotype blocks to consider."""
        return self._nhaploblk
    @nhaploblk.setter
    def nhaploblk(self, value: Integral) -> None:
        """Set number of haplotype blocks to consider."""
        check_is_Integral(value, "nhaploblk")
        check_is_gt(value, "nhaploblk", 0)
        self._nhaploblk = value
    
    @property
    def nbestfndr(self) -> Integral:
        """Number of best founders to consider."""
        return self._nbestfndr
    @nbestfndr.setter
    def nbestfndr(self, value: Integral) -> None:
        """Set number of best founders to consider."""
        check_is_Integral(value, "nbestfndr")           # must be int
        check_is_gt(value, "nbestfndr", 0)
        self._nbestfndr = value

class GenotypeBuilderSubsetSelection(GenotypeBuilderSelectionMixin,SubsetSelectionProtocol):
    """
    Class defining Genotype Builder (GB) Selection for subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nhaploblk: Integral,
            nbestfndr: Integral,
            ncross: Integral,
            nparent: Integral,
            nmating: Union[Integral,numpy.ndarray],
            nprogeny: Union[Integral,numpy.ndarray],
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
            ndset_wt: Optional[Real] = None,
            ndset_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            rng: Optional[Union[Generator,RandomState]] = None, 
            soalgo: Optional[SubsetOptimizationAlgorithm] = None,
            moalgo: Optional[SubsetOptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class GenotypeBuilderSubsetSelection.

        Parameters
        ----------
        ntrait : Integral
            Number of traits to expect from matrix inputs.
        
        nhaploblk : Integral
            Number of haplotype blocks to consider.

        nbestfndr : Integral
            Number of best founders to consider.

        ncross : Integral
            Number of cross configurations to consider.
        
        nparent : Integral
            Number of parents per cross configuration.
        
        nmating : Integral, numpy.ndarray
            Number of matings per configuration.

            If ``nmating`` is ``Integral``, then broadcast to a ``numpy.ndarray`` 
            of shape ``(ncross,)``.
            
            If ``nmating`` is ``numpy.ndarray``, then the array must be of type 
            ``Integral`` and of shape ``(ncross,)``.
        
        nprogeny : Integral, numpy.ndarray
            Number of progeny to derive from each mating event.

            If ``nprogeny`` is ``Integral``, then broadcast to a ``numpy.ndarray`` 
            of shape ``(ncross,)``.
            
            If ``nprogeny`` is ``numpy.ndarray``, then the array must be of type 
            ``Integral`` and of shape ``(ncross,)``.

        nobj : Integral
            Number of optimization objectives when constructing a 
            ``SelectionProblem``. This is equivalent to the vector length 
            returned by the ``obj_trans`` function. Must be ``Integral`` greater 
            than 0.
        
        obj_wt : numpy.ndarray, Real, None
            Objective function weights. Weights from this vector are applied 
            to objective function values via the Hadamard product. If values 
            are ``1.0`` or ``-1.0``, this can be used to specify minimizing 
            and maximizing objectives, respectively.

            If ``obj_wt`` is ``numpy.ndarray``, then the array must be of shape 
            ``(nobj,)``.

            If ``obj_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(nobj,)``.

            If ``obj_wt`` is ``None``, then the value ``1.0`` is broadcast to a 
            ``numpy.ndarray`` of shape ``(nobj,)``. This assumes that all 
            objectives are to be minimized.
        
        obj_trans : Callable, None
            A function which transforms values from a latent objective space to 
            the objective space. This transformation function must have the 
            following signature::

                def obj_trans(
                        decnvec: numpy.ndarray,
                        latentvec: numpy.ndarray, 
                        **kwargs: dict
                    ) -> numpy.ndarray:
                    # do stuff
                    return output
            
            Where:

            - ``decnvec`` is a ``numpy.ndarray`` containing the decision vector.
            - ``latentvec`` is a ``numpy.ndarray`` containing the latent space 
                objective function values which are to be transformed.
            - ``kwargs`` is a ``dict`` containing additional keyword arguments.

            If ``obj_trans`` is ``None``, then default to an identity objective 
            transformation function.

        obj_trans_kwargs : dict
            Keyword arguments for the latent space to objective space 
            transformation function. 

            If `obj_trans_kwargs`` is ``None``, then default to an empty 
            dictionary.

        nineqcv : Integral, None
            Number of inequality constraint violation functions. This is 
            equivalent to the vector length returned by the ``ineqcv_trans`` 
            function. Must be ``Integral`` greater than or equal to zero.

            If ``nineqcv`` is ``None``, then set to zero.

        ineqcv_wt : numpy.ndarray, None
            Inequality constraint violation function weights. Weights from this 
            vector are applied to inequality constraint violation function 
            values via the Hadamard product. If values are ``1.0`` or ``-1.0``, 
            this can be used to specify minimizing and maximizing constraints, 
            respectively.

            If ``ineqcv_wt`` is ``numpy.ndarray``, then the array must be of 
            shape ``(nineqcv,)``.

            If ``ineqcv_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(nineqcv,)``.

            If ``ineqcv_wt`` is ``None``, then the value ``1.0`` is broadcast 
            to a ``numpy.ndarray`` of shape ``(nineqcv,)``. This assumes that 
            all constraints are to be minimized.

        ineqcv_trans : Callable, None
            A function which transforms values from a latent objective space to 
            the inequality constraint violation space. This transformation 
            function must have the following signature::

                def ineqcv_trans(
                        decnvec: numpy.ndarray,
                        latentvec: numpy.ndarray, 
                        **kwargs: dict
                    ) -> numpy.ndarray:
                    # do stuff
                    return output
            
            Where:

            - ``decnvec`` is a ``numpy.ndarray`` containing the decision vector.
            - ``latentvec`` is a ``numpy.ndarray`` containing the latent space 
                objective function values which are to be transformed.
            - ``kwargs`` is a ``dict`` containing additional keyword arguments.

            If ``ineqcv_trans`` is ``None``, then default to a transformation 
            function returning an empty vector.
        
        ineqcv_trans_kwargs : dict, None
            Keyword arguments for the latent space to inequality constraint 
            violation transformation function.
        
            If `ineqcv_trans_kwargs`` is ``None``, then default to an empty 
            dictionary.

        neqcv : Integral, None
            Number of equality constraint violations. This is equivalent to the 
            vector length returned by the ``eqcv_trans`` function. Must be 
            ``Integral`` greater than or equal to zero.
        
            If ``neqcv`` is ``None``, then set to zero.

        eqcv_wt : numpy.ndarray, None
            Equality constraint violation function weights. Weights from this 
            vector are applied to equality constraint violation function 
            values via the Hadamard product. If values are ``1.0`` or ``-1.0``, 
            this can be used to specify minimizing and maximizing constraints, 
            respectively.

            If ``eqcv_wt`` is ``numpy.ndarray``, then the array must be of 
            shape ``(neqcv,)``.

            If ``eqcv_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(neqcv,)``.

            If ``eqcv_wt`` is ``None``, then the value ``1.0`` is broadcast 
            to a ``numpy.ndarray`` of shape ``(neqcv,)``. This assumes that 
            all constraints are to be minimized.

        eqcv_trans : Callable, None
            A function which transforms values from a latent objective space to 
            the equality constraint violation space. This transformation 
            function must have the following signature::

                def eqcv_trans(
                        decnvec: numpy.ndarray,
                        latentvec: numpy.ndarray, 
                        **kwargs: dict
                    ) -> numpy.ndarray:
                    # do stuff
                    return output
            
            Where:

            - ``decnvec`` is a ``numpy.ndarray`` containing the decision vector.
            - ``latentvec`` is a ``numpy.ndarray`` containing the latent space 
                objective function values which are to be transformed.
            - ``kwargs`` is a ``dict`` containing additional keyword arguments.

            If ``eqcv_trans`` is ``None``, then default to a transformation 
            function returning an empty vector.

        eqcv_trans_kwargs : dict, None
            Keyword arguments for the latent space to equality constraint 
            violation transformation function.

            If `eqcv_trans_kwargs`` is ``None``, then default to an empty 
            dictionary.

        ndset_wt : Real, None
            Nondominated set weight. The weight from this function is applied 
            to outputs from ``ndset_trans``. If values are ``1.0`` or ``-1.0``, 
            this can be used to specify minimizing and maximizing objectives, 
            respectively.

            If ``ndset_wt`` is ``None``, then it is set to the default value of ``1.0``.
            This assumes that the objective is to be minimized.

        ndset_trans : Callable, None
            A function which transforms values from the non-dominated set 
            objective space to the single-objective space. This transformation 
            function must have the following signature::

                def ndset_trans(
                        mat: numpy.ndarray, 
                        **kwargs: dict
                    ) -> numpy.ndarray:
                    # do stuff
                    return output
            
            Where:

            - ``mat`` is a ``numpy.ndarray`` containing a point coordinate array 
                of shape ``(npt, nobj)`` where ``npt`` is the number of points 
                and ``nobj`` is the number of objectives (dimensions). This 
                array contains input points for calculating the distance between 
                a point to the vector ``vec_wt``.
            - ``kwargs`` is a ``dict`` containing additional keyword arguments.

            If ``ndset_trans`` is ``None``, then default to a transformation 
            function calculating the distance between a weight vector and 
            provided points

        ndset_trans_kwargs : dict, None
            Nondominated set transformation function keyword arguments.

            If ``ndset_trans_kwargs`` is ``None``, then default to defaults for 
            the default ``ndset_trans`` function::

                ndset_trans_kwargs = {
                    "obj_wt": numpy.repeat(1.0, nobj),
                    "vec_wt": numpy.repeat(1.0, nobj)
                }

        rng : numpy.random.Generator, numpy.random.RandomState, None
            Random number source.

            If ``rng`` is ``None``, default to the global random number 
            generator.

        soalgo : SubsetOptimizationAlgorithm, None
            Single-objective optimization algorithm.

            If ``soalgo`` is ``None``, then use a default single-objective 
            optimization algorithm.

        moalgo : SubsetOptimizationAlgorithm, None
            Multi-objective opimization algorithm.

            If ``moalgo`` is ``None``, then use a default multi-objective 
            optimization algorithm.

        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        # make assignments from unmodified Mixin class first
        self.ntrait = ntrait
        self.nhaploblk = nhaploblk
        # make assignments from SubsetSelectionProtocol second
        super(GenotypeBuilderSubsetSelection, self).__init__(
            ncross = ncross,
            nparent = nparent,
            nmating = nmating,
            nprogeny = nprogeny,
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
            ndset_wt = ndset_wt,
            ndset_trans = ndset_trans,
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )
        # make assignments from modified Mixin class properties third
        self.nbestfndr = nbestfndr

    ############################ Object Properties #############################
    @GenotypeBuilderSelectionMixin.nbestfndr.setter
    def nbestfndr(self, value: Integral) -> None:
        """Set number of best founders to consider."""
        check_is_Integral(value, "nbestfndr")           # must be int
        check_is_in_interval_inclusive(value, "nbestfndr", 0, self.nselindiv) # must be in [0,nparent]
        self._nbestfndr = value

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: pandas.DataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            **kwargs: dict
        ) -> SelectionProblem:
        """
        Create an optimization problem definition using provided inputs.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : pandas.DataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : SelectionProblem
            An optimization problem definition.
        """
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_AdditiveLinearGenomicModel(gpmod, "gpmod")
        
        # get number of taxa
        ntaxa = pgmat.ntaxa

        # get decision space parameters
        decn_space = numpy.arange(ntaxa)
        decn_space_lower = numpy.repeat(0, self.nparent)
        decn_space_upper = numpy.repeat(ntaxa-1, self.nparent)

        # construct problem
        prob = GenotypeBuilderSubsetSelectionProblem.from_pgmat_gpmod(
            pgmat = pgmat,
            gpmod = gpmod,
            nhaploblk = self.nhaploblk,
            nbestfndr = self.nbestfndr,
            ndecn = self.nparent,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = self.nobj,
            obj_wt = self.obj_wt,
            obj_trans = self.obj_trans,
            obj_trans_kwargs = self.obj_trans_kwargs,
            nineqcv = self.nineqcv,
            ineqcv_wt = self.ineqcv_wt,
            ineqcv_trans = self.ineqcv_trans,
            ineqcv_trans_kwargs = self.ineqcv_trans_kwargs,
            neqcv = self.neqcv,
            eqcv_wt = self.eqcv_wt,
            eqcv_trans = self.eqcv_trans,
            eqcv_trans_kwargs = self.eqcv_trans_kwargs
        )

        return prob

    ################ Single Objective Solve ################
    # inherit sosolve() from SubsetSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from SubsetSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from SubsetSelectionProtocol
