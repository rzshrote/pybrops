"""
Module defining a general class for subset selection protocols.
"""

from numbers import Integral
from typing import Optional, Union
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.opt.algo.NSGA2SubsetGeneticAlgorithm import NSGA2SubsetGeneticAlgorithm
from pybrops.opt.algo.SubsetGeneticAlgorithm import SubsetGeneticAlgorithm
from pybrops.opt.algo.SubsetOptimizationAlgorithm import SubsetOptimizationAlgorithm, check_is_SubsetOptimizationAlgorithm
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame


class SubsetSelectionProtocol(SelectionProtocol):
    """
    """
    ########################## Special Object Methods ##########################
    # __init__ abstract method from SelectionProtocol

    ############################ Object Properties #############################
    @property
    def soalgo(self) -> SubsetOptimizationAlgorithm:
        """Single-objective optimization algorithm."""
        return self._soalgo
    @soalgo.setter
    def soalgo(self, value: Union[SubsetOptimizationAlgorithm,None]) -> None:
        """Set single-objective optimization algorithm."""
        if value is None:
            # construct default multi-objective algorithm
            value = SubsetGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                pop_size = 100  # number of parents in population
            )
            # construct default hillclimber
            # value = SteepestDescentSubsetHillClimber(self.rng)
        check_is_SubsetOptimizationAlgorithm(value, "soalgo")
        self._soalgo = value

    @property
    def moalgo(self) -> SubsetOptimizationAlgorithm:
        """Multi-objective opimization algorithm."""
        return self._moalgo
    @moalgo.setter
    def moalgo(self, value: Union[SubsetOptimizationAlgorithm,None]) -> None:
        """Set multi-objective opimization algorithm."""
        if value is None:
            # construct default multi-objective algorithm
            value = NSGA2SubsetGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                pop_size = 100  # number of parents in population
            )
        check_is_SubsetOptimizationAlgorithm(value, "moalgo")
        self._moalgo = value

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    # abstract `problem` method from SelectionProtocol

    ############## Pareto Frontier Functions ###############
    def pareto(
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
        ) -> tuple:
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
        out : tuple
            A tuple containing two objects ``(frontier, sel_config)``.

            Where:

            - ``frontier`` is a ``numpy.ndarray`` of shape ``(q,v)`` containing
              Pareto frontier points.
            - ``sel_config`` is a ``numpy.ndarray`` of shape ``(q,k)`` containing
              parent selection decisions for each corresponding point in the
              Pareto frontier.

            Where:

            - ``q`` is the number of points in the frontier.
            - ``v`` is the number of objectives for the frontier.
            - ``k`` is the number of search space decision variables.
        """
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

        if miscout is not None:
            miscout["soln"] = soln

        frontier = soln.soln_obj
        sel_config = soln.soln_decn

        return frontier, sel_config

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
        ) -> tuple:
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
        out : tuple
            A tuple containing four objects: ``(pgmat, sel, ncross, nprogeny)``.

            Where:

            - ``pgmat`` is a PhasedGenotypeMatrix of parental candidates.
            - ``sel`` is a ``numpy.ndarray`` of indices specifying a cross
              pattern. Each index corresponds to an individual in ``pgmat``.
            - ``ncross`` is a ``numpy.ndarray`` specifying the number of
              crosses to perform per cross pattern.
            - ``nprogeny`` is a ``numpy.ndarray`` specifying the number of
              progeny to generate per cross.
        """
        # if the number of objectives is 1, then we use a single objective algorithm
        if self.nobj == 1:
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
            soln = self.soalgo.minimize(
                prob = prob,
                miscout = miscout
            )

            if miscout is not None:
                miscout["soln"] = soln

            # extract decision variables as a vector (ndecn,)
            sel = soln.soln_decn[0]

            return pgmat, sel, self.nmating, self.nprogeny
        # else, we use a multi-objective algorithm and apply a transformation on the non-dominated points to identify a 
        # multi-objective method: objfn_trans returns a multiple values for each
        # selection configuration
        elif self.nobj > 1:
            # get the pareto frontier
            frontier, sel_config = self.pareto(
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

            # get scores for each of the points along the pareto frontier
            score = self.ndset_wt * self.ndset_trans(frontier, **self.ndset_trans_kwargs)

            # get index of maximum score
            ix = score.argmax()

            # add fields to miscout
            if miscout is not None:
                miscout["frontier"] = frontier
                miscout["sel_config"] = sel_config

            return pgmat, sel_config[ix], self.nmating, self.nprogeny
        else:
            raise ValueError("number of objectives must be greater than zero")
