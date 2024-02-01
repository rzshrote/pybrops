from pybrops.breed.op.psel.ParentSelectionOperator import ParentSelectionOperator


class DummyParentSelectionOperator(ParentSelectionOperator):
    def pselect(self, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, miscout: dict, **kwargs: dict) -> tuple:
        return super().pselect(genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs)

