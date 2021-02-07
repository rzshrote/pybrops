from . import WeightedGenomicParentSelection
from . import WeightedGenomicSurvivorSelection

class WeightedGenomicSelection(WeightedGenomicParentSelection,WeightedGenomicSurvivorSelection):
    """docstring for WeightedGenomicSelection."""

    def __init__(self, k_p, traitwt_p, ncross, nprogeny, k_s, traitwt_s, rng, **kwargs):
        super(WeightedGenomicSelection, self).__init__(
            k_p = k_p,
            traitwt_p = traitwt_p,
            ncross = ncross,
            nprogeny = nprogeny,
            k_s = k_s,
            traitwt_s = traitwt_s,
            rng = rng
        )
