from pybrops.breed.op.mate.MatingOperator import MatingOperator


class DummyMatingOperator(MatingOperator):
    def mate(self, mcfg: dict, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, miscout: dict, **kwargs: dict) -> tuple:
        return super().mate(mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs)

