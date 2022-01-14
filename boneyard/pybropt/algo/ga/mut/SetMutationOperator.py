import copy
import numpy

from . import MutationOperator

class SetMutationOperator(MutationOperator):
    """docstring for SetMutationOperator."""

    def __init__(self, p, setspace, rng, **kwargs):
        super(SetMutationOperator, self).__init__(**kwargs)
        self.p = p
        self.rng = rng

    def mutate(self, chrom, **kwargs):
        # initialize
        out = copy.copy(chrom)                      # copy chromosome
        clen = len(chrom)                           # get chromosome length
        mut = self.rng.uniform(0, 1, clen) < self.p # determine mutation points

        # mutate algorithm
        for i in range(clen):                       # for each element in chromosome
            if mut[i]:                              # if is has been selected for mutation
                mask = ~numpy.in1d(                 # determine available mutations
                    self.setspace,                  # set space
                    out                             # working copy of chromosome
                )
                out[i] = self.rng.choice(           # randomly select one value
                    self.setspace[mask],            # from reduced set
                    1,
                    replace = False
                )

        return out
