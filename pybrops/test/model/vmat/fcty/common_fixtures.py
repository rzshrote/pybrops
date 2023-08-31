from numbers import Integral
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix
from pybrops.model.vmat.AdditiveGenicVarianceMatrix import AdditiveGenicVarianceMatrix
from pybrops.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix
from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix
from pybrops.model.vmat.fcty.AdditiveGeneticVarianceMatrixFactory import AdditiveGeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty.AdditiveGenicVarianceMatrixFactory import AdditiveGenicVarianceMatrixFactory
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty.GenicVarianceMatrixFactory import GenicVarianceMatrixFactory
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class DummyGenicVarianceMatrixFactory(GenicVarianceMatrixFactory):
    def from_gmod(self, gmod: GenomicModel, pgmat: PhasedGenotypeMatrix, nprogeny: int, **kwargs: dict) -> GenicVarianceMatrix:
        return super().from_gmod(gmod, pgmat, nprogeny, **kwargs)

class DummyAdditiveGenicVarianceMatrixFactory(DummyGenicVarianceMatrixFactory,AdditiveGenicVarianceMatrixFactory):
    def from_algmod(self, algmod: AdditiveLinearGenomicModel, pgmat: PhasedGenotypeMatrix, nprogeny: int, mem: int, **kwargs: dict) -> AdditiveGenicVarianceMatrix:
        return super().from_algmod(algmod, pgmat, nprogeny, mem, **kwargs)

class DummyGeneticVarianceMatrixFactory(GeneticVarianceMatrixFactory):
    def from_gmod(self, gmod: GenomicModel, pgmat: PhasedGenotypeMatrix, ncross: Integral, nprogeny: Integral, nself: Integral, gmapfn: GeneticMapFunction, **kwargs: dict) -> GeneticVarianceMatrix:
        return super().from_gmod(gmod, pgmat, ncross, nprogeny, nself, gmapfn, **kwargs)

class DummyAdditiveGeneticVarianceMatrixFactory(DummyGeneticVarianceMatrixFactory,AdditiveGeneticVarianceMatrixFactory):
    def from_algmod(self, algmod: AdditiveLinearGenomicModel, pgmat: PhasedGenotypeMatrix, ncross: int, nprogeny: int, nself: int, gmapfn: GeneticMapFunction, mem: int, **kwargs: dict) -> AdditiveGeneticVarianceMatrix:
        return super().from_algmod(algmod, pgmat, ncross, nprogeny, nself, gmapfn, mem, **kwargs)
    