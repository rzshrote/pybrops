class SelectionOperator:
    """docstring for SelectionOperator."""

    def __init__(self, **kwargs):
        super(SelectionOperator, self).__init__()

    def select(self, t_cur, t_max, k, ppop, pscore, opop, oscore, **kwargs):
        raise NotImplementedError("method is abstract")
