import copy

from . import MemeticOperator

class NoMemeticOperator(MemeticOperator):
    """docstring for NoMemeticOperator."""

    def __init__(self, **kwargs: dict):
        super(NoMemeticOperator, self).__init__(**kwargs)

    def evolve(self, objfn, pop, score, **kwargs: dict):
        new_pop = copy.copy(pop)        # shallow copy list
        new_score = copy.copy(score)    # copy array

        return new_pop, new_score
