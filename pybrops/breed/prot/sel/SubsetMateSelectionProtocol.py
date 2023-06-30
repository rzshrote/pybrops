"""
Module defining a general class for subset selection protocols.
"""

from abc import ABCMeta, abstractmethod
from numbers import Integral
from typing import Optional

from pybrops.breed.prot.sel.MateSelectionProtocol import MateSelectionProtocol
from pybrops.breed.prot.sel.SubsetSelectionProtocol import SubsetSelectionProtocol
from pybrops.breed.prot.sel.cfg.SubsetMateSelectionConfiguration import SubsetMateSelectionConfiguration
from pybrops.breed.prot.sel.soln.SubsetMateSelectionSolution import SubsetMateSelectionSolution
from pybrops.breed.prot.sel.prob.SubsetMateSelectionProblem import SubsetMateSelectionProblem
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame


class SubsetMateSelectionProtocol(SubsetSelectionProtocol,MateSelectionProtocol,metaclass=ABCMeta):
    """
    Semi-abstract class for creating subset selection protocols.
    """
    ########################## Special Object Methods ##########################
    # inherit __init__ from SubsetSelectionProtocol

    ############################ Object Properties #############################
    # inherit properties from SubsetSelectionProtocol

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    @abstractmethod
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            **kwargs: dict
        ) -> SubsetMateSelectionProblem:
        """
        Create an optimization problem definition using provided inputs.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
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
        out : SubsetMateSelectionProblem
            An optimization problem definition.
        """
        raise NotImplementedError("method is abstract")

    ################ Single Objective Solve ################
    def sosolve(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> SubsetMateSelectionSolution:
        """
        Solve the selection problem using a single-objective optimization algorithm.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : SubsetMateSelectionSolution
            A single-objective solution to the posed selection problem.
        """
        # check the number of objectives and raise error if needed
        if self.nobj != 1:
            raise RuntimeError("{0} instance is not single-objective in nature: expected nobj == 1 but received nobj == {1}".format(type(self).__name__,self.nobj))

        # construct the problem
        prob = self.problem(
            pgmat = pgmat, 
            gmat = gmat, 
            ptdf = ptdf, 
            bvmat = bvmat, 
            gpmod = gpmod, 
            t_cur = t_cur, 
            t_max = t_max, 
            **kwargs
        )

        # optimize the problem
        soln = self.soalgo.minimize(
            prob = prob,
            miscout = miscout
        )

        # convert subset solution to subset selection solution
        # add cross map metadata from problem to metadata from SubsetSolution
        out = SubsetMateSelectionSolution(
            ndecn = soln.ndecn,
            decn_space = soln.decn_space,
            decn_space_lower = soln.decn_space_lower,
            decn_space_upper = soln.decn_space_upper,
            decn_space_xmap = prob.decn_space_xmap,
            nobj = soln.nobj,
            obj_wt = soln.obj_wt,
            nineqcv = soln.nineqcv,
            ineqcv_wt = soln.ineqcv_wt,
            neqcv = soln.neqcv,
            eqcv_wt = soln.eqcv_wt,
            nsoln = soln.nsoln,
            soln_decn = soln.soln_decn,
            soln_obj = soln.soln_obj,
            soln_ineqcv = soln.soln_ineqcv,
            soln_eqcv = soln.soln_eqcv
        )

        return out

    ############## Pareto Frontier Functions ###############
    def mosolve(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> SubsetMateSelectionSolution:
        """
        Calculate a Pareto frontier for objectives.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : SubsetMateSelectionSolution
            A multi-objective solution to the posed selection problem.
        """
        # check the number of objectives and raise error if needed
        if self.nobj <= 1:
            raise RuntimeError("{0} instance is not multi-objective in nature: expected nobj > 1 but received nobj == {1}".format(type(self).__name__,self.nobj))

        # construct the problem
        prob = self.problem(
            pgmat = pgmat,
            gmat = gmat,
            ptdf = ptdf,
            bvmat = bvmat,
            gpmod = gpmod,
            t_cur = t_cur,
            t_max = t_max
        )

        # optimize the problem
        soln = self.moalgo.minimize(
            prob = prob,
            miscout = miscout
        )

        # convert subset solution to subset selection solution
        # add cross map metadata from problem to metadata from SubsetSolution
        out = SubsetMateSelectionSolution(
            ndecn = soln.ndecn,
            decn_space = soln.decn_space,
            decn_space_lower = soln.decn_space_lower,
            decn_space_upper = soln.decn_space_upper,
            decn_space_xmap = prob.decn_space_xmap,
            nobj = soln.nobj,
            obj_wt = soln.obj_wt,
            nineqcv = soln.nineqcv,
            ineqcv_wt = soln.ineqcv_wt,
            neqcv = soln.neqcv,
            eqcv_wt = soln.eqcv_wt,
            nsoln = soln.nsoln,
            soln_decn = soln.soln_decn,
            soln_obj = soln.soln_obj,
            soln_ineqcv = soln.soln_ineqcv,
            soln_eqcv = soln.soln_eqcv
        )

        return out

    ################# Selection Functions ##################
    def select(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> SubsetMateSelectionConfiguration:
        """
        Select individuals for breeding.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : SubsetMateSelectionConfiguration
            A selection configuration object, requiring all necessary information to mate individuals.
        """
        # if the number of objectives is 1, then we use a single objective algorithm
        if self.nobj == 1:
            # solve the single-objective problem
            sosoln = self.sosolve(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                miscout = miscout,
                **kwargs
            )

            # add solution to miscout if provided.
            if miscout is not None:
                miscout["sosoln"] = sosoln

            # construct a SubsetMateSelectionConfiguration
            selcfg = SubsetMateSelectionConfiguration(
                ncross = self.ncross,
                nparent = self.nparent,
                nmating = self.nmating,
                nprogeny = self.nprogeny,
                pgmat = pgmat,
                xconfig_decn = sosoln.soln_decn[0],
                xconfig_xmap = sosoln.decn_space_xmap,
                rng = None
            )

            return selcfg

        # else, we use a multi-objective algorithm and apply a transformation on the non-dominated points to identify a 
        # multi-objective method: objfn_trans returns a multiple values for each
        # selection configuration
        elif self.nobj > 1:
            # solve the single-objective problem
            mosoln = self.mosolve(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                miscout = miscout,
                **kwargs
            )

            # get scores for each of the points along the non-dominated set
            # (nndpts,)
            score = self.ndset_wt * self.ndset_trans(
                mosoln.soln_obj, 
                **self.ndset_trans_kwargs
            )

            # get index of maximum score
            ix = score.argmax()

            # add solution to miscout if provided.
            if miscout is not None:
                miscout["mosoln"] = mosoln

            # construct a SubsetMateSelectionConfiguration
            selcfg = SubsetMateSelectionConfiguration(
                ncross = self.ncross,
                nparent = self.nparent,
                nmating = self.nmating,
                nprogeny = self.nprogeny,
                pgmat = pgmat,
                xconfig_decn = mosoln.soln_decn[ix],
                xconfig_xmap = mosoln.decn_space_xmap,
                rng = None
            )

            return selcfg

        # else raise an error as the number of objectives is an illegal value
        else:
            raise ValueError("number of objectives must be greater than zero")



################################## Utilities ###################################
def check_is_SubsetMateSelectionProtocol(v: object, vname: str) -> None:
    """
    Check if object is of type SubsetMateSelectionProtocol, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SubsetMateSelectionProtocol):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,SubsetMateSelectionProtocol.__name__,type(v).__name__))
