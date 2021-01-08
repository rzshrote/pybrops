from . import LogBook

class BlankLogBook(LogBook):
    """docstring for BlankLogBook."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(BlankLogBook, self).__init__(**kwargs)
        self.reset()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def book():
        doc = "The book property."
        def fget(self):
            return self._book
        def fset(self, value):
            self._book = value
        def fdel(self):
            del self._book
        return locals()
    book = property(**book())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def log_initialize(self, t_cur, t_max, geno, bval, gmod, **kwargs):
        pass

    def log_pselect(self, t_cur, t_max, pgvmat, sel, ncross, nprogeny, misc, **kwargs):
        pass

    def log_mate(self, t_cur, t_max, pgvmat, misc, **kwargs):
        pass

    def log_evaluate(self, t_cur, t_max, bvmat, misc, **kwargs):
        pass

    def log_calibrate(self, t_cur, t_max, gmod, misc, **kwargs):
        pass

    def log_sselect(self, t_cur, t_max, geno, bval, gmod, misc, **kwargs):
        pass

    def reset(self):
        """
        Reset LogBook internals.
        """
        self.book = None

    def write(self, fname):
        """
        Write LogBook to file

        Parameters
        ----------
        fname : str
            File name to which to write file.
        """
        pass
