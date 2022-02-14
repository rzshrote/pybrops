class MutationOperator:
    """docstring for MutationOperator."""

    def __init__(self, **kwargs):
        super(MutationOperator, self).__init__()

    def mutate(self, chrom, **kwargs):
        raise NotImplementedError("method is abstract")
