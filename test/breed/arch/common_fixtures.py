from typing import Optional
from pybrops.breed.arch.BreedingEdge import BreedingEdge
from pybrops.breed.arch.BreedingGraph import BreedingGraph
from pybrops.breed.arch.BreedingNode import BreedingNode
from pybrops.breed.arch.BreedingProgram import BreedingProgram
from pybrops.breed.arch.EmigrationOperator import EmigrationOperator
from pybrops.breed.arch.GermplasmBank import GermplasmBank
from pybrops.breed.arch.ImmigrationOperator import ImmigrationOperator
from pybrops.breed.op.eval.EvaluationOperator import EvaluationOperator
from pybrops.breed.op.init.InitializationOperator import InitializationOperator
from pybrops.breed.op.mate.MatingOperator import MatingOperator
from pybrops.breed.op.psel.ParentSelectionOperator import ParentSelectionOperator
from pybrops.breed.op.ssel.SurvivorSelectionOperator import SurvivorSelectionOperator


class DummyBreedingEdge(BreedingEdge):
    pass

class DummyBreedingGraph(BreedingGraph):
    @property
    def graph(self) -> object:
        return super().graph
    @graph.setter
    def graph(self, value: object) -> None:
        super().graph = value

class DummyBreedingNode(BreedingNode):
    @property
    def genome(self) -> object:
        return super().genome
    @genome.setter
    def genome(self, value: object) -> None:
        super().genome = value
    @property
    def geno(self) -> object:
        return super().geno
    @geno.setter
    def geno(self, value: object) -> None:
        super().geno = value
    @property
    def pheno(self) -> object:
        return super().pheno
    @pheno.setter
    def pheno(self, value: object) -> None:
        super().pheno = value
    @property
    def bval(self) -> object:
        return super().bval
    @bval.setter
    def bval(self, value: object) -> None:
        super().bval = value
    @property
    def gmod(self) -> object:
        return super().gmod
    @gmod.setter
    def gmod(self, value: object) -> None:
        super().gmod = value
    @property
    def t_cur(self) -> object:
        return super().t_cur
    @t_cur.setter
    def t_cur(self, value: object) -> None:
        super().t_cur = value
    @property
    def t_max(self) -> object:
        return super().t_max
    @t_max.setter
    def t_max(self, value: object) -> None:
        super().t_max = value

class DummyBreedingProgram(DummyBreedingNode,BreedingProgram):
    @property
    def start_genome(self) -> object:
        return super().start_genome
    @start_genome.setter
    def start_genome(self, value: object) -> None:
        super().start_genome = value
    @property
    def start_geno(self) -> object:
        return super().start_geno
    @start_geno.setter
    def start_geno(self, value: object) -> None:
        super().start_geno = value
    @property
    def start_pheno(self) -> object:
        return super().start_pheno
    @start_pheno.setter
    def start_pheno(self, value: object) -> None:
        super().start_pheno = value
    @property
    def start_bval(self) -> object:
        return super().start_bval
    @start_bval.setter
    def start_bval(self, value: object) -> None:
        super().start_bval = value
    @property
    def start_gmod(self) -> object:
        return super().start_gmod
    @start_gmod.setter
    def start_gmod(self, value: object) -> None:
        super().start_gmod = value
    @property
    def initop(self) -> object:
        return super().initop
    @initop.setter
    def initop(self, value: object) -> None:
        super().initop = value
    @property
    def pselop(self) -> object:
        return super().pselop
    @pselop.setter
    def pselop(self, value: object) -> None:
        super().pselop = value
    @property
    def mateop(self) -> object:
        return super().mateop
    @mateop.setter
    def mateop(self, value: object) -> None:
        super().mateop = value
    @property
    def evalop(self) -> object:
        return super().evalop
    @evalop.setter
    def evalop(self, value: object) -> None:
        super().evalop = value
    @property
    def sselop(self) -> object:
        return super().sselop
    @sselop.setter
    def sselop(self, value: object) -> None:
        super().sselop = value
    def initialize(self, **kwargs: dict):
        return super().initialize(**kwargs)
    def is_initialized(self, **kwargs: dict):
        return super().is_initialized(**kwargs)
    def reset(self, **kwargs: dict):
        return super().reset(**kwargs)
    def advance(self, ngen, lbook, **kwargs: dict):
        return super().advance(ngen, lbook, **kwargs)
    def evolve(self, nrep, ngen, lbook, **kwargs: dict):
        return super().evolve(nrep, ngen, lbook, **kwargs)

class DummyEmigrationOperator(DummyBreedingEdge,EmigrationOperator):
    def emigrate(self, bnode, **kwargs: dict):
        return super().emigrate(bnode, **kwargs)

class DummyGermplasmBank(DummyBreedingNode,GermplasmBank):
    pass

class DummyImmigrationOperator(DummyBreedingEdge,ImmigrationOperator):
    def immigrate(self, bnode, **kwargs: dict):
        return super().immigrate(bnode, **kwargs)

class DummyInitializationOperator(InitializationOperator):
    def initialize(self, miscout: dict | None, **kwargs: dict) -> tuple:
        return super().initialize(miscout, **kwargs)

class DummyParentSelectionOperator(ParentSelectionOperator):
    def pselect(self, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, miscout: dict, **kwargs: dict) -> tuple:
        return super().pselect(genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs)

class DummyMatingOperator(MatingOperator):
    def mate(self, mcfg: dict, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, miscout: dict, **kwargs: dict) -> tuple:
        return super().mate(mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs)

class DummyEvaluationOperator(EvaluationOperator):
    def evaluate(self, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, miscout: dict, **kwargs: dict) -> tuple:
        return super().evaluate(genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs)

class DummySurvivorSelectionOperator(SurvivorSelectionOperator):
    def sselect(self, genome: dict, geno: dict, pheno: dict, bval: dict, gmod: dict, t_cur: int, t_max: int, miscout: dict, **kwargs: dict) -> tuple:
        return super().sselect(genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs)

