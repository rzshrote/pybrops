class LogBook:
    """docstring for LogBook."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(LogBook, self).__init__()

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

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def log_initialize(self, t_cur, t_max, geno, bval, gmod, **kwargs):
        raise NotImplementedError("method is abstract")

    def log_pselect(self, t_cur, t_max, pgvmat, sel, ncross, nprogeny, misc, **kwargs):
        raise NotImplementedError("method is abstract")

    def log_mate(self, t_cur, t_max, pgvmat, misc, **kwargs):
        raise NotImplementedError("method is abstract")

    def log_evaluate(self, t_cur, t_max, bvmat, misc, **kwargs):
        raise NotImplementedError("method is abstract")

    def log_calibrate(self, t_cur, t_max, gmod, misc, **kwargs):
        raise NotImplementedError("method is abstract")

    def log_sselect(self, t_cur, t_max, geno, bval, gmod, misc, **kwargs):
        raise NotImplementedError("method is abstract")

    def reset(self):
        """
        Reset LogBook internals.
        """
        raise NotImplementedError("method is abstract")

    def write(self, fname):
        """
        Write LogBook to file

        Parameters
        ----------
        fname : str
            File name to which to write file.
        """
        raise NotImplementedError("method is abstract")
