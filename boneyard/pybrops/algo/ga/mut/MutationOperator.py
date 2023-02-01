class MutationOperator:
    """docstring for MutationOperator."""

    def __init__(self, **kwargs: dict):
        super(MutationOperator, self).__init__()

    def mutate(self, chrom, **kwargs: dict):
        raise NotImplementedError("method is abstract")
