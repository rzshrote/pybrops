from numbers import Integral
import pandas
import pytest
from pybrops.breed.prot.sel.MateSelectionProtocol import MateSelectionProtocol
from pybrops.breed.prot.sel.MateSelectionProtocol import check_is_MateSelectionProtocol
from pybrops.breed.prot.sel.cfg.SelectionConfiguration import SelectionConfiguration
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.test.assert_python import assert_abstract_method, assert_docstring, assert_semiabstract_class, not_raises
from .common_fixtures_large import *

class DummyMateSelectionProtocol(MateSelectionProtocol):
    @property
    def soalgo(self) -> object:
        """soalgo."""
        return self._soalgo
    @soalgo.setter
    def soalgo(self, value: object) -> None:
        """Set soalgo."""
        self._soalgo = value
    @property
    def moalgo(self) -> object:
        """moalgo."""
        return self._moalgo
    @moalgo.setter
    def moalgo(self, value: object) -> None:
        """Set moalgo."""
        self._moalgo = value
    def problem(self, pgmat: PhasedGenotypeMatrix, gmat: GenotypeMatrix, ptdf: pandas.DataFrame, bvmat: BreedingValueMatrix, gpmod: GenomicModel, t_cur: Integral, t_max: Integral, **kwargs: dict) -> SelectionProblem:
        return super().problem(pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs)
    def sosolve(self, pgmat: PhasedGenotypeMatrix, gmat: GenotypeMatrix, ptdf: pandas.DataFrame, bvmat: BreedingValueMatrix, gpmod: GenomicModel, t_cur: Integral, t_max: Integral, miscout: dict | None, **kwargs: dict) -> SelectionSolution:
        return super().sosolve(pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout, **kwargs)
    def mosolve(self, pgmat: PhasedGenotypeMatrix, gmat: GenotypeMatrix, ptdf: pandas.DataFrame, bvmat: BreedingValueMatrix, gpmod: GenomicModel, t_cur: Integral, t_max: Integral, miscout: dict | None, **kwargs: dict) -> SelectionSolution:
        return super().mosolve(pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout, **kwargs)
    def select(self, pgmat: PhasedGenotypeMatrix, gmat: GenotypeMatrix, ptdf: pandas.DataFrame, bvmat: BreedingValueMatrix, gpmod: GenomicModel, t_cur: Integral, t_max: Integral, miscout: dict | None, **kwargs: dict) -> SelectionConfiguration:
        return super().select(pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout, **kwargs)

@pytest.fixture
def selprot(
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_nobj_single
    ):
    out = DummyMateSelectionProtocol(
        ncross=common_ncross,
        nparent=common_nparent,
        nmating=common_nmating,
        nprogeny=common_nprogeny,
        nobj=common_nobj_single
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_MateSelectionProtocol_is_semiabstract():
    assert_semiabstract_class(MateSelectionProtocol)

############################## Test class docstring ############################
def test_MateSelectionProtocol_docstring():
    assert_docstring(MateSelectionProtocol)

############################ Test class properties #############################

############################## Test class methods ##############################
def test_MateSelectionProtocol_problem_is_abstract():
    assert_abstract_method(MateSelectionProtocol, "problem")

def test_MateSelectionProtocol_sosolve_is_abstract():
    assert_abstract_method(MateSelectionProtocol, "sosolve")

def test_MateSelectionProtocol_mosolve_is_abstract():
    assert_abstract_method(MateSelectionProtocol, "mosolve")

def test_MateSelectionProtocol_select_is_abstract():
    assert_abstract_method(MateSelectionProtocol, "select")
    
############################# Test class utilities #############################
def test_check_is_MateSelectionProtocol(selprot):
    with not_raises(Exception):
        check_is_MateSelectionProtocol(selprot, "selprot")
    with pytest.raises(TypeError):
        check_is_MateSelectionProtocol(object(), "selprot")
    with pytest.raises(TypeError):
        check_is_MateSelectionProtocol(None, "selprot")
