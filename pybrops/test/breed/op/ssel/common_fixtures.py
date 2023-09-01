from pybrops.breed.op.ssel.SurvivorSelectionOperator import SurvivorSelectionOperator


class DummySurvivorSelectionOperator(SurvivorSelectionOperator):
    def sselect(self, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, miscout: dict, **kwargs: dict) -> tuple:
        return super().sselect(genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs)

