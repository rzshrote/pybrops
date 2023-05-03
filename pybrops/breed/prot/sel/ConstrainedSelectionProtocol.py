"""
Module defining a general class for selection protocols.
"""

# list of all public objects in this module
__all__ = [
    "ConstrainedSelectionProtocol"
]

# imports
from numbers import Integral, Number, Real
from typing import Callable, Optional, Union

import numpy
from pybrops.breed.prot.sel.prob.trans import trans_empty, trans_identity, trans_ndpt_to_vec_dist
from pybrops.core.error.error_type_python import check_is_Callable, check_is_Integral, check_is_Real, check_is_dict, check_is_str
from pybrops.core.error.error_value_numpy import check_ndarray_len_eq, check_ndarray_ndim
from pybrops.core.error.error_value_python import check_is_gt, check_is_gteq
from pybrops.breed.prot.sel.ConstrainedSelectionProtocolType import ConstrainedSelectionProtocolType

class ConstrainedSelectionProtocol(ConstrainedSelectionProtocolType):
    """
    A semi-abstract class implementing several key properties common to most, 
    if not all, constrained selection protocols.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            method: str,
            nobj: Integral,
            obj_wt: Optional[numpy.ndarray] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[numpy.ndarray] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[numpy.ndarray] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            ndset_wt: Optional[Real] = None,
            ndset_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedSelectionProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(ConstrainedSelectionProtocol, self).__init__(**kwargs)
        # order dependent assignments
        self.method = method
        self.nobj = nobj
        self.obj_wt = obj_wt
        self.obj_trans = obj_trans
        self.obj_trans_kwargs = obj_trans_kwargs
        self.nineqcv = nineqcv
        self.ineqcv_wt = ineqcv_wt
        self.ineqcv_trans = ineqcv_trans
        self.ineqcv_trans_kwargs = ineqcv_trans_kwargs
        self.neqcv = neqcv
        self.eqcv_wt = eqcv_wt
        self.eqcv_trans = eqcv_trans
        self.eqcv_trans_kwargs = eqcv_trans_kwargs
        self.ndset_wt = ndset_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def method(self) -> str:
        """Selection method."""
        return self._method
    @method.setter
    def method(self, value: str) -> None:
        """Set selection method."""
        check_is_str(value, "method")       # must be string
        value = value.lower()               # convert to lowercase
        options = ("single", "pareto")      # method options
        # if not method supported raise ValueError
        if value not in options:
            raise ValueError("Unsupported 'method'. Options are: " + ", ".join(map(str, options)))
        self._method = value

    @property
    def nobj(self) -> Integral:
        """Number of objectives."""
        return self._nobj
    @nobj.setter
    def nobj(self, value: Integral) -> None:
        """Set number of objectives."""
        check_is_Integral(value, "nobj")
        check_is_gt(value, "nobj", 0)     # cannot have 0 objectives
        self._nobj = value
    
    @property
    def obj_wt(self) -> numpy.ndarray:
        """Objective function weights."""
        return self._obj_wt
    @obj_wt.setter
    def obj_wt(self, value: Union[numpy.ndarray,Number,None]) -> None:
        """Set objective function weights. If None, set to 1.0."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "obj_wt", 1)
            check_ndarray_len_eq(value, "obj_wt", self.nobj)
        elif isinstance(value, Number):
            value = numpy.repeat(value, self.nobj)
        elif value is None:
            value = numpy.repeat(1.0, self.nobj)
        else:
            raise TypeError("'obj_wt' must be of type numpy.ndarray, a numeric type, or None")
        self._obj_wt = value

    @property
    def obj_trans(self) -> Callable[[numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to objective function values."""
        return self._obj_trans
    @obj_trans.setter
    def obj_trans(self, value: Union[Callable[[numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to objective space transformation function. If None, set to identity function."""
        if value is None:
            value = trans_identity
        check_is_Callable(value, "obj_trans")
        self._obj_trans = value
    
    @property
    def obj_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to objective space transformation function."""
        return self._obj_trans_kwargs
    @obj_trans_kwargs.setter
    def obj_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to objective space transformation function. If None, set to empty dict."""
        if value is None:
            value = {}
        check_is_dict(value, "obj_trans_kwargs")
        self._obj_trans_kwargs = value
    
    @property
    def nineqcv(self) -> Integral:
        """Number of inequality constraint violation functions."""
        return self._nineqcv
    @nineqcv.setter
    def nineqcv(self, value: Union[Integral,None]) -> None:
        """Set number of inequality constraint violation functions. If None, set to 0."""
        if value is None:
            value = 0
        check_is_Integral(value, "nineqcv")
        check_is_gteq(value, "nineqcv", 0)  # possible to have 0 inequality constraints
        self._nineqcv = value

    @property
    def ineqcv_wt(self) -> numpy.ndarray:
        """Inequality constraint violation function weights."""
        return self._ineqcv_wt
    @ineqcv_wt.setter
    def ineqcv_wt(self, value: Union[numpy.ndarray,Number,None]) -> None:
        """Set inequality constraint violation function weights. If None, set to 1.0."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "ineqcv_wt", 1)
            check_ndarray_len_eq(value, "ineqcv_wt", self.nineqcv)
        elif isinstance(value, Number):
            value = numpy.repeat(value, self.nineqcv)
        elif value is None:
            value = numpy.repeat(1.0, self.nineqcv)
        else:
            raise TypeError("'ineqcv_wt' must be of type numpy.ndarray, a numeric type, or None")
        self._ineqcv_wt = value

    @property
    def ineqcv_trans(self) -> Callable[[numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to inequality constraint violation values."""
        return self._ineqcv_trans
    @ineqcv_trans.setter
    def ineqcv_trans(self, value: Union[Callable[[numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to inequality constraint violation transformation function. If None, set to the empty function."""
        if value is None:
            value = trans_empty
        check_is_Callable(value, "ineqcv_trans")
        self._ineqcv_trans = value
    
    @property
    def ineqcv_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to inequality constraint violation transformation function."""
        return self._ineqcv_trans_kwargs
    @ineqcv_trans_kwargs.setter
    def ineqcv_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to inequality constraint violation transformation function. If None, set to empty dict."""
        if value is None:
            value = {}
        check_is_dict(value, "ineqcv_trans_kwargs")
        self._ineqcv_trans_kwargs = value
    
    @property
    def neqcv(self) -> Integral:
        """Number of equality constraint violations."""
        return self._neqcv
    @neqcv.setter
    def neqcv(self, value: Union[Integral,None]) -> None:
        """Set number of equality constraint violations. If None, set to 0."""
        if value is None:
            value = 0
        check_is_Integral(value, "neqcv")
        check_is_gteq(value, "neqcv", 0)    # possible to have 0 equality constraints
        self._neqcv = value
    
    @property
    def eqcv_wt(self) -> numpy.ndarray:
        """Equality constraint violation function weights."""
        return self._eqcv_wt
    @eqcv_wt.setter
    def eqcv_wt(self, value: Union[numpy.ndarray,Number,None]) -> None:
        """Set equality constraint violation function weights. If None, set to 1.0."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "eqcv_wt", 1)
            check_ndarray_len_eq(value, "eqcv_wt", self.neqcv)
        elif isinstance(value, Number):
            value = numpy.repeat(value, self.neqcv)
        elif value is None:
            value = numpy.repeat(1.0, self.neqcv)
        else:
            raise TypeError("'eqcv_wt' must be of type numpy.ndarray or a numeric type")
        self._eqcv_wt = value

    @property
    def eqcv_trans(self) -> Callable[[numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to equality constraint violation values."""
        return self._eqcv_trans
    @eqcv_trans.setter
    def eqcv_trans(self, value: Union[Callable[[numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to equality constraint violation transformation function. If None, set to the empty function."""
        if value is None:
            value = trans_empty
        check_is_Callable(value, "eqcv_trans")
        self._eqcv_trans = value 
    
    @property
    def eqcv_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to equality constraint violation transformation function."""
        return self._eqcv_trans_kwargs
    @eqcv_trans_kwargs.setter
    def eqcv_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to equality constraint violation transformation function. If None, set to empty dict."""
        if value is None:
            value = {}
        check_is_dict(value, "eqcv_trans_kwargs")
        self._eqcv_trans_kwargs = value

    @property
    def ndset_wt(self) -> Real:
        """Nondominated set weights."""
        return self._ndset_wt
    @ndset_wt.setter
    def ndset_wt(self, value: Union[Real,None]) -> None:
        """Set nondominated set weights. If None, set to 1.0."""
        if value is None:
            value = 1.0
        check_is_Real(value, "ndset_wt")
        self._ndset_wt = value

    @property
    def ndset_trans(self) -> Callable[[numpy.ndarray,dict],numpy.ndarray]:
        """Nondominated set transformation function."""
        return self._ndset_trans
    @ndset_trans.setter
    def ndset_trans(self, value: Union[Callable[[numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set nondominated set transformation function. If None, set to closest point to vector function."""
        if value is None:
            value = trans_ndpt_to_vec_dist
        check_is_Callable(value, "ndset_trans")
        self._ndset_trans = value

    @property
    def ndset_trans_kwargs(self) -> dict:
        """Nondominated set transformation function keyword arguments."""
        return self._ndset_trans_kwargs
    @ndset_trans_kwargs.setter
    def ndset_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set nondominated set transformation function keyword arguments. If None, set to dict with one vectors."""
        if value is None:                           # if given None
            value = {                               # set default to empty dict
                "obj_wt": numpy.repeat(1.0, self.nobj),
                "vec_wt": numpy.repeat(1.0, self.nobj)
            }
        check_is_dict(value, "ndset_trans_kwargs")  # check is dict
        self._ndset_trans_kwargs = value

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ########## Optimization Problem Construction ###########
    # leave `problem` abstract

    ############## Pareto Frontier Functions ###############
    # leave `pareto` abstract
    # def pareto(
    #         self, 
    #         pgmat: PhasedGenotypeMatrix, 
    #         gmat: GenotypeMatrix, 
    #         ptdf: PhenotypeDataFrame, 
    #         bvmat: BreedingValueMatrix, 
    #         gpmod: GenomicModel, 
    #         t_cur: Integral, 
    #         t_max: Integral, 
    #         miscout: Optional[dict] = None, 
    #         **kwargs: dict
    #     ) -> tuple:
    #     """
    #     Calculate a Pareto frontier for objectives.

    #     Parameters
    #     ----------
    #     pgmat : PhasedGenotypeMatrix
    #         Genomes
    #     gmat : GenotypeMatrix
    #         Genotypes
    #     ptdf : PhenotypeDataFrame
    #         Phenotype dataframe
    #     bvmat : BreedingValueMatrix
    #         Breeding value matrix
    #     gpmod : GenomicModel
    #         Genomic prediction model
    #     t_cur : int
    #         Current generation number.
    #     t_max : int
    #         Maximum (deadline) generation number.
    #     miscout : dict, None
    #         Pointer to a dictionary for miscellaneous user defined output.
    #         If ``dict``, write to dict (may overwrite previously defined fields).
    #         If ``None``, user defined output is not calculated or stored.
    #     kwargs : dict
    #         Additional keyword arguments.

    #     Returns
    #     -------
    #     out : tuple
    #         A tuple containing two objects ``(frontier, sel_config)``.

    #         Where:

    #         - ``frontier`` is a ``numpy.ndarray`` of shape ``(q,v)`` containing
    #           Pareto frontier points.
    #         - ``sel_config`` is a ``numpy.ndarray`` of shape ``(q,k)`` containing
    #           parent selection decisions for each corresponding point in the
    #           Pareto frontier.

    #         Where:

    #         - ``q`` is the number of points in the frontier.
    #         - ``v`` is the number of objectives for the frontier.
    #         - ``k`` is the number of search space decision variables.
    #     """
    #     # construct the problem
    #     prob = self.problem(
    #         pgmat = pgmat,
    #         gmat = gmat,
    #         ptdf = ptdf,
    #         bvmat = bvmat,
    #         gpmod = gpmod,
    #         t_cur = t_cur,
    #         t_max = t_max
    #     )

    #     # optimize the problem
    #     soln = self.moalgo.minimize(
    #         prob = prob,
    #         miscout = miscout
    #     )

    #     if miscout is not None:
    #         miscout["soln"] = soln

    #     frontier = soln.soln_obj
    #     sel_config = soln.soln_decn

    #     return frontier, sel_config

    ################# Selection Functions ##################
    # leave `select` abstract
    # def select(
    #         self, 
    #         pgmat: PhasedGenotypeMatrix, 
    #         gmat: GenotypeMatrix, 
    #         ptdf: PhenotypeDataFrame, 
    #         bvmat: BreedingValueMatrix, 
    #         gpmod: GenomicModel, 
    #         t_cur: Integral, 
    #         t_max: Integral, 
    #         miscout: Optional[dict] = None, 
    #         **kwargs: dict
    #     ) -> tuple:
    #     """
    #     Select individuals for breeding.

    #     Parameters
    #     ----------
    #     pgmat : PhasedGenotypeMatrix
    #         Genomes
    #     gmat : GenotypeMatrix
    #         Genotypes
    #     ptdf : PhenotypeDataFrame
    #         Phenotype dataframe
    #     bvmat : BreedingValueMatrix
    #         Breeding value matrix
    #     gpmod : GenomicModel
    #         Genomic prediction model
    #     t_cur : int
    #         Current generation number.
    #     t_max : int
    #         Maximum (deadline) generation number.
    #     miscout : dict, None
    #         Pointer to a dictionary for miscellaneous user defined output.
    #         If ``dict``, write to dict (may overwrite previously defined fields).
    #         If ``None``, user defined output is not calculated or stored.
    #     kwargs : dict
    #         Additional keyword arguments.

    #     Returns
    #     -------
    #     out : tuple
    #         A tuple containing four objects: ``(pgmat, sel, ncross, nprogeny)``.

    #         Where:

    #         - ``pgmat`` is a PhasedGenotypeMatrix of parental candidates.
    #         - ``sel`` is a ``numpy.ndarray`` of indices specifying a cross
    #           pattern. Each index corresponds to an individual in ``pgmat``.
    #         - ``ncross`` is a ``numpy.ndarray`` specifying the number of
    #           crosses to perform per cross pattern.
    #         - ``nprogeny`` is a ``numpy.ndarray`` specifying the number of
    #           progeny to generate per cross.
    #     """
    #     # single-objective method: objfn_trans returns a single value for each
    #     # selection configuration
    #     if self.method == "single":
    #         # construct the problem
    #         prob = self.problem(
    #             pgmat = pgmat,
    #             gmat = gmat,
    #             ptdf = ptdf,
    #             bvmat = bvmat,
    #             gpmod = gpmod,
    #             t_cur = t_cur,
    #             t_max = t_max
    #         )

    #         # optimize the problem
    #         soln = self.soalgo.minimize(
    #             prob = prob,
    #             miscout = miscout
    #         )

    #         if miscout is not None:
    #             miscout["soln"] = soln

    #         # extract decision variables
    #         sel = soln.soln_decn[0]

    #         return pgmat, sel, self.ncross, self.nprogeny

    #     # multi-objective method: objfn_trans returns a multiple values for each
    #     # selection configuration
    #     elif self.method == "pareto":
    #         # get the pareto frontier
    #         frontier, sel_config = self.pareto(
    #             pgmat = pgmat,
    #             gmat = gmat,
    #             ptdf = ptdf,
    #             bvmat = bvmat,
    #             gpmod = gpmod,
    #             t_cur = t_cur,
    #             t_max = t_max,
    #             miscout = miscout,
    #             **kwargs
    #         )

    #         # get scores for each of the points along the pareto frontier
    #         score = self.ndset_wt * self.ndset_trans(frontier, **self.ndset_trans_kwargs)

    #         # get index of maximum score
    #         ix = score.argmax()

    #         # add fields to miscout
    #         if miscout is not None:
    #             miscout["frontier"] = frontier
    #             miscout["sel_config"] = sel_config

    #         return pgmat, sel_config[ix], self.ncross, self.nprogeny
    #     else:
    #         raise ValueError("argument 'method' must be either 'single' or 'pareto'")
