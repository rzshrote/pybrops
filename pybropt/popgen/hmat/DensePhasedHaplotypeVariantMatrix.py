# TODO: class write-up

from . import HaplotypeVariantMatrix, PhasedHaplotypeMatrix

class PhasedHaplotypeVariantMatrix(HaplotypeVariantMatrix,PhasedHaplotypeMatrix):
    """docstring for PhasedHaplotypeVariantMatrix."""

    def __init__(self, **kwargs):
        super(PhasedHaplotypeVariantMatrix, self).__init__(**kwargs)
