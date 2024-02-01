from pybrops.breed.op.log.Logbook import Logbook


class DummyLogbook(Logbook):
    @property
    def data(self) -> object:
        return super().data
    @data.setter
    def data(self, value: object) -> None:
        super().data = value
    @property
    def rep(self) -> object:
        return super().rep
    @rep.setter
    def rep(self, value: object) -> None:
        super().rep = value
    def log_initialize(self, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, **kwargs: dict) -> None:
        return super().log_initialize(genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs)
    def log_pselect(self, mcfg: dict, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, **kwargs: dict) -> None:
        return super().log_pselect(mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs)
    def log_mate(self, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, **kwargs: dict) -> None:
        return super().log_mate(genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs)
    def log_evaluate(self, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, **kwargs: dict) -> None:
        return super().log_evaluate(genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs)
    def log_sselect(self, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, **kwargs: dict) -> None:
        return super().log_sselect(genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs)
    def reset(self) -> None:
        return super().reset()
    def write(self, filename: str) -> None:
        return super().write(filename)
    