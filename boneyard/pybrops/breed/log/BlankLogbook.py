from . import Logbook

class BlankLogbook(Logbook):
    """docstring for BlankLogbook."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict):
        super(BlankLogbook, self).__init__(**kwargs)
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
    def log_initialize(self, t_cur, t_max, geno, bval, gmod, **kwargs: dict):
        pass

    def log_pselect(self, t_cur, t_max, pgvmat, sel, ncross, nprogeny, misc, **kwargs: dict):
        pass

    def log_mate(self, t_cur, t_max, pgvmat, misc, **kwargs: dict):
        pass

    def log_evaluate(self, t_cur, t_max, bvmat, misc, **kwargs: dict):
        pass

    def log_calibrate(self, t_cur, t_max, gmod, misc, **kwargs: dict):
        pass

    def log_sselect(self, t_cur, t_max, geno, bval, gmod, misc, **kwargs: dict):
        pass

    def reset(self):
        """
        Reset Logbook internals.
        """
        self.book = None

    def write(self, fname):
        """
        Write Logbook to file

        Parameters
        ----------
        fname : str
            File name to which to write file.
        """
        pass
