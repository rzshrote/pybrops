import pandas

from . import LogBook

class BasicLogBook(LogBook):
    """docstring for BasicLogBook."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(BasicLogBook, self).__init__(**kwargs)
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

    def rep():
        doc = "The rep property."
        def fget(self):
            return self._rep
        def fset(self, value):
            self._rep = value
        def fdel(self):
            del self._rep
        return locals()
    rep = property(**rep())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def log_initialize(self, t_cur, t_max, geno, bval, gmod, **kwargs):
        self.book["rep"].append(self.rep)
        self.book["t_cur"].append(t_cur)
        self.book["t_max"].append(t_max)
        # candidate pool
        self.book["cand_mehe"].append(geno["cand"].mehe())
        self.book["cand_usl"].append(gmod["cand"].usl(geno["cand"]).sum())
        self.book["cand_lsl"].append(gmod["cand"].lsl(geno["cand"]).sum())
        self.book["cand_mean_perf"].append(bval["cand"].mat.mean(0).sum())
        self.book["cand_max_perf"].append(bval["cand"].mat.max(0).sum())
        self.book["cand_var_A"].append(gmod["cand"].var_A(geno["cand"]).sum())
        self.book["cand_var_a"].append(gmod["cand"].var_a(geno["cand"]).sum())
        self.book["cand_bulmer"].append(gmod["cand"].bulmer(geno["cand"]).sum())
        # main pool
        self.book["main_mehe"].append(geno["main"].mehe())
        self.book["main_usl"].append(gmod["main"].usl(geno["main"]).sum())
        self.book["main_lsl"].append(gmod["main"].lsl(geno["main"]).sum())
        self.book["main_mean_perf"].append(bval["main"].mat.mean(0).sum())
        self.book["main_max_perf"].append(bval["main"].mat.max(0).sum())
        self.book["main_var_A"].append(gmod["main"].var_A(geno["main"]).sum())
        self.book["main_var_a"].append(gmod["main"].var_a(geno["main"]).sum())
        self.book["main_bulmer"].append(gmod["main"].bulmer(geno["main"]).sum())
        # candidate pool (true genomic model)
        self.book["cand_true_usl"].append(gmod["true"].usl(geno["cand"]).sum())
        self.book["cand_true_lsl"].append(gmod["true"].lsl(geno["cand"]).sum())
        self.book["cand_true_mean_perf"].append(bval["cand_true"].mat.mean(0).sum())
        self.book["cand_true_max_perf"].append(bval["cand_true"].mat.max(0).sum())
        self.book["cand_true_var_A"].append(gmod["true"].var_A(geno["cand"]).sum())
        self.book["cand_true_var_a"].append(gmod["true"].var_a(geno["cand"]).sum())
        self.book["cand_true_bulmer"].append(gmod["true"].bulmer(geno["cand"]).sum())
        # main pool (true genomic model)
        self.book["main_true_usl"].append(gmod["true"].usl(geno["main"]).sum())
        self.book["main_true_lsl"].append(gmod["true"].lsl(geno["main"]).sum())
        self.book["main_true_mean_perf"].append(bval["main_true"].mat.mean(0).sum())
        self.book["main_true_max_perf"].append(bval["main_true"].mat.max(0).sum())
        self.book["main_true_var_A"].append(gmod["true"].var_A(geno["main"]).sum())
        self.book["main_true_var_a"].append(gmod["true"].var_a(geno["main"]).sum())
        self.book["main_true_bulmer"].append(gmod["true"].bulmer(geno["main"]).sum())

    def log_pselect(self, t_cur, t_max, pgvmat, sel, ncross, nprogeny, misc, **kwargs):
        pass

    def log_mate(self, t_cur, t_max, pgvmat, misc, **kwargs):
        pass

    def log_gintegrate(self, t_cur, t_max, geno, misc, **kwargs):
        pass

    def log_evaluate(self, t_cur, t_max, bvmat, bvmat_true, misc, **kwargs):
        pass

    def log_bvintegrate(self, t_cur, t_max, bval, misc, **kwargs):
        pass

    def log_calibrate(self, t_cur, t_max, gmod, misc, **kwargs):
        pass

    def log_sselect(self, t_cur, t_max, geno, bval, gmod, misc, **kwargs):
        self.book["rep"].append(self.rep)
        self.book["t_cur"].append(t_cur)
        self.book["t_max"].append(t_max)
        # candidate pool
        self.book["cand_mehe"].append(geno["cand"].mehe())
        self.book["cand_usl"].append(gmod["cand"].usl(geno["cand"]).sum())
        self.book["cand_lsl"].append(gmod["cand"].lsl(geno["cand"]).sum())
        self.book["cand_mean_perf"].append(bval["cand"].mat.mean(0).sum())
        self.book["cand_max_perf"].append(bval["cand"].mat.max(0).sum())
        self.book["cand_var_A"].append(gmod["cand"].var_A(geno["cand"]).sum())
        self.book["cand_var_a"].append(gmod["cand"].var_a(geno["cand"]).sum())
        self.book["cand_bulmer"].append(gmod["cand"].bulmer(geno["cand"]).sum())
        # main pool
        self.book["main_mehe"].append(geno["main"].mehe())
        self.book["main_usl"].append(gmod["main"].usl(geno["main"]).sum())
        self.book["main_lsl"].append(gmod["main"].lsl(geno["main"]).sum())
        self.book["main_mean_perf"].append(bval["main"].mat.mean(0).sum())
        self.book["main_max_perf"].append(bval["main"].mat.max(0).sum())
        self.book["main_var_A"].append(gmod["main"].var_A(geno["main"]).sum())
        self.book["main_var_a"].append(gmod["main"].var_a(geno["main"]).sum())
        self.book["main_bulmer"].append(gmod["main"].bulmer(geno["main"]).sum())
        # candidate pool (true genomic model)
        self.book["cand_true_usl"].append(gmod["true"].usl(geno["cand"]).sum())
        self.book["cand_true_lsl"].append(gmod["true"].lsl(geno["cand"]).sum())
        self.book["cand_true_mean_perf"].append(bval["cand_true"].mat.mean(0).sum())
        self.book["cand_true_max_perf"].append(bval["cand_true"].mat.max(0).sum())
        self.book["cand_true_var_A"].append(gmod["true"].var_A(geno["cand"]).sum())
        self.book["cand_true_var_a"].append(gmod["true"].var_a(geno["cand"]).sum())
        self.book["cand_true_bulmer"].append(gmod["true"].bulmer(geno["cand"]).sum())
        # main pool (true genomic model)
        self.book["main_true_usl"].append(gmod["true"].usl(geno["main"]).sum())
        self.book["main_true_lsl"].append(gmod["true"].lsl(geno["main"]).sum())
        self.book["main_true_mean_perf"].append(bval["main_true"].mat.mean(0).sum())
        self.book["main_true_max_perf"].append(bval["main_true"].mat.max(0).sum())
        self.book["main_true_var_A"].append(gmod["true"].var_A(geno["main"]).sum())
        self.book["main_true_var_a"].append(gmod["true"].var_a(geno["main"]).sum())
        self.book["main_true_bulmer"].append(gmod["true"].bulmer(geno["main"]).sum())


    def reset(self):
        """
        Reset LogBook internals.
        """
        self.rep = 1
        self.book = {
            "rep" : [],
            "t_cur" : [],
            "t_max" : [],
            # candidate pool
            "cand_mehe" : [],
            "cand_usl" : [],
            "cand_lsl" : [],
            "cand_mean_perf" : [],
            "cand_max_perf" : [],
            "cand_var_A" : [],
            "cand_var_a" : [],
            "cand_bulmer" : [],
            # main pool
            "main_mehe" : [],
            "main_usl" : [],
            "main_lsl" : [],
            "main_mean_perf" : [],
            "main_max_perf" : [],
            "main_var_A" : [],
            "main_var_a" : [],
            "main_bulmer" : [],
            # candidate pool (true genomic model)
            "cand_true_usl" : [],
            "cand_true_lsl" : [],
            "cand_true_mean_perf" : [],
            "cand_true_max_perf" : [],
            "cand_true_var_A" : [],
            "cand_true_var_a" : [],
            "cand_true_bulmer" : [],
            # main pool (true genomic model)
            "main_true_usl" : [],
            "main_true_lsl" : [],
            "main_true_mean_perf" : [],
            "main_true_max_perf" : [],
            "main_true_var_A" : [],
            "main_true_var_a" : [],
            "main_true_bulmer" : [],
        }

    def write(self, fname):
        """
        Write LogBook to file

        Parameters
        ----------
        fname : str
            File name to which to write file.
        """
        df = pandas.DataFrame(self.book)
        df.to_csv(fname, index = False)
