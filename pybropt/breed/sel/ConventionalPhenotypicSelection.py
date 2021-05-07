from pybropt.breed.psel import ConventionalPhenotypicParentSelection
from pybropt.breed.ssel import ConventionalPhenotypicSurvivorSelection

class ConventionalPhenotypicSelection(ConventionalPhenotypicParentSelection,ConventionalPhenotypicSurvivorSelection):
    """docstring for ConventionalPhenotypicSelection."""

    def __init__(self, k_p, traitwt_p, ncross, nprogeny, k_s, traitwt_s, rng, **kwargs):
        super(ConventionalPhenotypicSelection, self).__init__(
            k_p = k_p,
            traitwt_p = traitwt_p,
            ncross = ncross,
            nprogeny = nprogeny,
            k_s = k_s,
            traitwt_s = traitwt_s,
            rng = rng
        )
