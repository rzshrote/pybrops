class Logbook:
    """docstring for Logbook."""

    def __init__(self, **kwargs: dict):
        super(Logbook, self).__init__()

    def book():
        doc = "The book property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    book = property(**book())

    def log(self, gen, pop, score):
        raise NotImplementedError("method is abstract")

    def reset(self):
        raise NotImplementedError("method is abstract")

    def write(self, fname):
        raise NotImplementedError("method is abstract")
