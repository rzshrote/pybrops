from pybrops.breed.op.eval.EvaluationOperator import EvaluationOperator


class DummyEvaluationOperator(EvaluationOperator):
    def evaluate(self, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, miscout: dict, **kwargs: dict) -> tuple:
        return super().evaluate(genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs)

