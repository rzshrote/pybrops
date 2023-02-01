import numpy

# https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
def is_pareto_efficient(fitnesses, return_mask = True):
    """
    Find the pareto-efficient points (maximizing function)
    :param fitnesses: An (n_points, n_fitnesses) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    is_efficient = numpy.arange(fitnesses.shape[0])
    n_points = fitnesses.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index<len(fitnesses):
        nondominated_point_mask = numpy.any(fitnesses>fitnesses[next_point_index], axis=1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        fitnesses = fitnesses[nondominated_point_mask]
        next_point_index = numpy.sum(nondominated_point_mask[:next_point_index])+1
    if return_mask:
        is_efficient_mask = numpy.zeros(n_points, dtype = bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient

def dcn(indiv, pop):
    closest = numpy.in1d(indiv, pop[0]).sum()
    for i in range(1,len(pop)):
        tmp = numpy.in1d(indiv, pop[i]).sum()
        if tmp < closest:
            closest = tmp
    return closest

class MULTISelectionOperator(SelectionOperator):
    """docstring for MULTISelectionOperator."""

    def __init__(self, objwt, **kwargs: dict):
        super(MULTISelectionOperator, self).__init__(**kwargs)
        self.objwt = objwt

    def select(self, t_cur, t_max, k, ppop, pscore, opop, oscore, objwt = None, **kwargs: dict):
        if objwt is None:
            objwt = self.objwt

        # concatenate populations
        pop = ppop + opop
        score = numpy.concatenate([pscore, oscore])
        wscore = numpy.dot(score, objwt)

        #
        new_pop = []
        new_score = []

        # find best individual
        ix = wscore.argmax()

        # append best individual and delete from main pool
        new_pop.append(pop.pop(ix))
        new_score.append(score[ix])
        score = numpy.delete(score, ix)
        wscore = numpy.delete(wscore, ix)

        while len(new_pop) < k:
            fitnesses = numpy.empty(
                (len(wscore), 2),
                dtype = wscore.dtype
            )
            fitnesses[:,0] = wscore
            for i in range(len(wscore)):
                fitnesses[i,1] = dcn(pop[i], new_pop)

            mask = is_pareto_efficient(fitnesses)
            ix = self.rng.choice(
                numpy.flatnonzero(mask),
                1
            )

            # append best individual and delete from main pool
            new_pop.append(pop.pop(ix))
            new_score.append(score[ix])
            score = numpy.delete(score, ix)
            wscore = numpy.delete(wscore, ix)

        new_score = numpy.array(new_score)

        return new_pop, new_score
