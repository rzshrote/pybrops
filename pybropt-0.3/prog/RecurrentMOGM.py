# 3rd party
import numpy
import pandas

# our libraries
import pybropt.util
import pybropt.breed
from . import BreedingProgram

class RecurrentMOGM(BreedingProgram):
    """docstring for RecurrentMOGM."""

    def __init__(self):
        super(RecurrentMOGM, self).__init__(method = "Recurrent MOGM")


    def simulate(self, population, cross, algorithm, bcycle = 0, objfn_kwargs = None, algo_kwargs = None, savegeno = False, seed = None):
        """
        Run breeding simulations.
        """
        # set seed if needed
        pybropt.util.cond_seed_rng(seed)

        # create MOGM object
        mogm = pybropt.breed.MOGM(population, cross)

        # get initial population score not applicable
        # pop_score = mogm.objfn([i for i in range(mogm.population.ntaxa)], **objfn_kwargs)

        # get initial population GEBVs
        pop_gebv = mogm.population.gebv(None)

        # record history
        self.add_population_history(
            method = mogm.method,
            algorithm = algorithm.name,
            seed = seed,
            cycle = 0,
            score = 0.0,
            gebv = pop_gebv,
            geno = population.geno if savegeno else None
        )

        # simulate the breeding cycles
        for cycle in range(1, bcycle+1):
            # create a new Breeding object
            mogm = mogm.from_self(
                population = population,
                cross = cross
            )

            # get number of taxa
            ntaxa = mogm.population.ntaxa

            # make states list
            sstates = [i for i in range(ntaxa)]

            # make search space dimensions; duplicate 20 times for 10 two-way crosses
            sdims = [sstates for i in range(20)]

            # make search space
            sspace = pybropt.algo.CategoricalSearchSpace(*sdims)

            # reset algorithm
            algorithm.reset()

            # update algorithm search space
            algorithm.search_space = sspace

            # make selections
            print("optimizing...")
            sel = mogm.optimize(
                algorithm = algorithm,
                **algo_kwargs,
                **objfn_kwargs
            )

            print("Cycle:", cycle, "Selected:", sel)

            # get score of the selection
            sel_score = algorithm.gbest_score()

            # get selection GEBVs
            sel_gebv = mogm.population.gebv(sel)

            # add history selection
            self.add_selection_history(
                method = mogm.method,
                algorithm = algorithm.name,
                seed = seed,
                cycle = cycle,
                score = sel_score,
                gebv = sel_gebv,
                geno = population.geno[:,sel,:] if savegeno else None
            )

            # mate and generate new population
            population = breeding.cross.mate(sel)

            # get score of the entire population
            pop_score = mogm.objfn(None, **objfn_kwargs)

            # get entire population GEBVs
            pop_gebv = mogm.population.gebv(None)

            # record history
            self.add_population_history(
                method = mogm.method,
                algorithm = algorithm.name,
                seed = seed,
                cycle = cycle,
                score = pop_score,
                gebv = pop_gebv,
                geno = population.geno if savegeno else None
            )

            # make new cross object
            cross = cross.from_self(
                population = population,
            )
