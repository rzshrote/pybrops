from pybropt.breed.psel import OptimalContributionParentSelection
from pybropt.breed.ssel import OptimalContributionSurvivorSelection

class OptimalContributionSelection(OptimalContributionParentSelection,OptimalContributionSurvivorSelection):
    """docstring for OptimalContributionSelection."""

    def __init__(self, k_p, traitwt_p, inbfn_p, ncross, nprogeny, k_s, traitwt_s, inbfn_s, rng, cmatcls, bvtype = "gebv", **kwargs):
        super(OptimalContributionSelection, self).__init__(
            k_p = k_p,
            traitwt_p = traitwt_p,
            inbfn_p = inbfn_p,
            ncross = ncross,
            nprogeny = nprogeny,
            k_s = k_s,
            traitwt_s = traitwt_s,
            inbfn_s = inbfn_s,
            rng = rng,
            cmatcls = cmatcls,
            bvtype = bvtype
            **kwargs
        )
