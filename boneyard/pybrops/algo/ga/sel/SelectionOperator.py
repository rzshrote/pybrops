class SelectionOperator:
    """docstring for SelectionOperator."""

    def __init__(self, **kwargs: dict):
        super(SelectionOperator, self).__init__()

    def select(self, t_cur, t_max, k, ppop, pscore, opop, oscore, **kwargs: dict):
        raise NotImplementedError("method is abstract")
