from . import TraitMatrix
from . import DenseMutableMatrix

# TODO: implement me
class DenseTraitMatrix(DenseMutableMatrix,TraitMatrix):
    """docstring for DenseTraitMatrix."""

    def __init__(self, **kwargs):
        super(DenseTraitMatrix, self).__init__(**kwargs)
