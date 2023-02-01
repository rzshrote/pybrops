class Logbook:
    """docstring for Logbook."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict):
        super(Logbook, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
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

    def rep():
        doc = "The rep property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    rep = property(**rep())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def log_initialize(self, t_cur, t_max, geno, bval, gmod, **kwargs: dict):
        raise NotImplementedError("method is abstract")

    def log_pselect(self, t_cur, t_max, pgvmat, sel, ncross, nprogeny, misc, **kwargs: dict):
        raise NotImplementedError("method is abstract")

    def log_mate(self, t_cur, t_max, pgvmat, misc, **kwargs: dict):
        raise NotImplementedError("method is abstract")

    def log_gintegrate(self, t_cur, t_max, geno, misc, **kwargs: dict):
        raise NotImplementedError("method is abstract")

    def log_evaluate(self, t_cur, t_max, bvmat, bvmat_true, misc, **kwargs: dict):
        raise NotImplementedError("method is abstract")

    def log_bvintegrate(self, t_cur, t_max, bval, misc, **kwargs: dict):
        raise NotImplementedError("method is abstract")

    def log_calibrate(self, t_cur, t_max, gmod, misc, **kwargs: dict):
        raise NotImplementedError("method is abstract")

    def log_sselect(self, t_cur, t_max, geno, bval, gmod, misc, **kwargs: dict):
        raise NotImplementedError("method is abstract")

    def reset(self):
        """
        Reset Logbook internals.
        """
        raise NotImplementedError("method is abstract")

    def write(self, fname):
        """
        Write Logbook to file

        Parameters
        ----------
        fname : str
            File name to which to write file.
        """
        raise NotImplementedError("method is abstract")
