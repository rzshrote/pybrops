from numbers import Integral, Real
from typing import Callable
import pytest
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.breed.prot.sel.SelectionProtocol import check_is_SelectionProtocol
from pybrops.breed.prot.sel.cfg.SelectionConfiguration import SelectionConfiguration
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame
from pybrops.test.assert_python import assert_concrete_property, assert_docstring, assert_semiabstract_class, not_raises
from pybrops.test.breed.prot.sel.common_fixtures_large import *

class DummySelectionProtocol(SelectionProtocol):
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
    def problem(self, pgmat: PhasedGenotypeMatrix, gmat: GenotypeMatrix, ptdf: PhenotypeDataFrame, bvmat: BreedingValueMatrix, gpmod: GenomicModel, t_cur: Integral, t_max: Integral, **kwargs: dict) -> SelectionProblem:
        return super().problem(pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs)
    def sosolve(self, pgmat: PhasedGenotypeMatrix, gmat: GenotypeMatrix, ptdf: PhenotypeDataFrame, bvmat: BreedingValueMatrix, gpmod: GenomicModel, t_cur: Integral, t_max: Integral, miscout: dict | None, **kwargs: dict) -> SelectionSolution:
        return super().sosolve(pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout, **kwargs)
    def mosolve(self, pgmat: PhasedGenotypeMatrix, gmat: GenotypeMatrix, ptdf: PhenotypeDataFrame, bvmat: BreedingValueMatrix, gpmod: GenomicModel, t_cur: Integral, t_max: Integral, miscout: dict | None, **kwargs: dict) -> SelectionSolution:
        return super().mosolve(pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout, **kwargs)
    def select(self, pgmat: PhasedGenotypeMatrix, gmat: GenotypeMatrix, ptdf: PhenotypeDataFrame, bvmat: BreedingValueMatrix, gpmod: GenomicModel, t_cur: Integral, t_max: Integral, miscout: dict | None, **kwargs: dict) -> SelectionConfiguration:
        return super().select(pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout, **kwargs)

@pytest.fixture
def selprot(
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_nobj_single
    ):
    out = DummySelectionProtocol(
        ncross=common_ncross,
        nparent=common_nparent,
        nmating=common_nmating,
        nprogeny=common_nprogeny,
        nobj=common_nobj_single
    )
    yield out


################### Test class abstract/concrete properties ####################
def test_SelectionProtocol_is_semiabstract():
    assert_semiabstract_class(SelectionProtocol)

############################## Test class docstring ############################
def test_SelectionProtocol_docstring():
    assert_docstring(SelectionProtocol)

############################ Test class properties #############################

### nselindiv ###
def test_SelectionProtocol_nselindiv_is_concrete():
    assert_concrete_property(SelectionProtocol, "nselindiv")

def test_nselindiv_fget(selprot, common_ncross, common_nparent):
    assert isinstance(selprot.nselindiv, Integral)
    assert selprot.nselindiv == (common_ncross * common_nparent)

def test_nselindiv_fset_AttributeError(selprot):
    with pytest.raises(AttributeError):
        selprot.nselindiv = int(7)

def test_nselindiv_fdel_AttributeError(selprot):
    with pytest.raises(AttributeError):
        del selprot.nselindiv

### ncross ###
def test_SelectionProtocol_ncross_is_concrete():
    assert_concrete_property(SelectionProtocol, "ncross")

def test_ncross_fget(selprot, common_ncross):
    assert isinstance(selprot.ncross, Integral)
    assert selprot.ncross == common_ncross

def test_ncross_fset(selprot):
    with not_raises(Exception):
        selprot.ncross = int(1)
    with not_raises(Exception):
        selprot.ncross = numpy.int8(1)
    with not_raises(Exception):
        selprot.ncross = numpy.int16(1)
    with not_raises(Exception):
        selprot.ncross = numpy.int32(1)
    with not_raises(Exception):
        selprot.ncross = numpy.int64(1)

def test_ncross_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.ncross = object()
    with pytest.raises(TypeError):
        selprot.ncross = None
    with pytest.raises(TypeError):
        selprot.ncross = float(1)

def test_ncross_fset_ValueError(selprot):
    with pytest.raises(ValueError):
        selprot.ncross = int(-1)
    with pytest.raises(ValueError):
        selprot.ncross = int(0)

def test_ncross_fdel_AttributeError(selprot):
    with pytest.raises(AttributeError):
        del selprot.ncross

### nparent ###
def test_SelectionProtocol_nparent_is_concrete():
    assert_concrete_property(SelectionProtocol, "nparent")

def test_nparent_fget(selprot, common_nparent):
    assert isinstance(selprot.nparent, Integral)
    assert selprot.nparent == common_nparent

def test_nparent_fset(selprot):
    with not_raises(Exception):
        selprot.nparent = int(1)
    with not_raises(Exception):
        selprot.nparent = numpy.int8(1)
    with not_raises(Exception):
        selprot.nparent = numpy.int16(1)
    with not_raises(Exception):
        selprot.nparent = numpy.int32(1)
    with not_raises(Exception):
        selprot.nparent = numpy.int64(1)

def test_nparent_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.nparent = object()
    with pytest.raises(TypeError):
        selprot.nparent = None
    with pytest.raises(TypeError):
        selprot.nparent = float(1)

def test_nparent_fset_ValueError(selprot):
    with pytest.raises(ValueError):
        selprot.nparent = int(-1)
    with pytest.raises(ValueError):
        selprot.nparent = int(0)

def test_nparent_fdel_AttributeError(selprot):
    with pytest.raises(AttributeError):
        del selprot.nparent

### nmating ###
def test_SelectionProtocol_nmating_is_concrete():
    assert_concrete_property(SelectionProtocol, "nmating")

def test_nmating_fget(selprot, common_nmating):
    assert isinstance(selprot.nmating, numpy.ndarray)
    assert numpy.all(selprot.nmating == common_nmating)

def test_nmating_fset(selprot, common_ncross):
    with not_raises(Exception):
        selprot.nmating = int(0)
    with not_raises(Exception):
        selprot.nmating = numpy.repeat(int(0), common_ncross)
    with not_raises(Exception):
        selprot.nmating = int(1)
    with not_raises(Exception):
        selprot.nmating = numpy.repeat(int(1), common_ncross)
    with not_raises(Exception):
        selprot.nmating = numpy.int8(1)
    with not_raises(Exception):
        selprot.nmating = numpy.repeat(numpy.int8(1), common_ncross)
    with not_raises(Exception):
        selprot.nmating = numpy.int16(1)
    with not_raises(Exception):
        selprot.nmating = numpy.repeat(numpy.int16(1), common_ncross)
    with not_raises(Exception):
        selprot.nmating = numpy.int32(1)
    with not_raises(Exception):
        selprot.nmating = numpy.repeat(numpy.int32(1), common_ncross)
    with not_raises(Exception):
        selprot.nmating = numpy.int64(1)
    with not_raises(Exception):
        selprot.nmating = numpy.repeat(numpy.int64(1), common_ncross)

def test_nmating_fset_TypeError(selprot, common_ncross):
    with pytest.raises(TypeError):
        selprot.nmating = object()
    with pytest.raises(TypeError):
        selprot.nmating = numpy.repeat(object(), common_ncross)
    with pytest.raises(TypeError):
        selprot.nmating = None
    with pytest.raises(TypeError):
        selprot.nmating = numpy.repeat(None, common_ncross)
    with pytest.raises(TypeError):
        selprot.nmating = float(1)
    with pytest.raises(TypeError):
        selprot.nmating = numpy.repeat(float(1), common_ncross)

def test_nmating_fset_ValueError(selprot, common_ncross):
    with pytest.raises(ValueError):
        selprot.nmating = int(-1)
    with pytest.raises(ValueError):
        selprot.nmating = numpy.repeat(int(-1), common_ncross)

def test_nmating_fdel_AttributeError(selprot):
    with pytest.raises(AttributeError):
        del selprot.nmating

### nprogeny ###
def test_SelectionProtocol_nprogeny_is_concrete():
    assert_concrete_property(SelectionProtocol, "nprogeny")

def test_nprogeny_fget(selprot, common_nprogeny):
    assert isinstance(selprot.nprogeny, numpy.ndarray)
    assert numpy.all(selprot.nprogeny == common_nprogeny)

def test_nprogeny_fset(selprot, common_ncross):
    with not_raises(Exception):
        selprot.nprogeny = int(0)
    with not_raises(Exception):
        selprot.nprogeny = numpy.repeat(int(0), common_ncross)
    with not_raises(Exception):
        selprot.nprogeny = int(1)
    with not_raises(Exception):
        selprot.nprogeny = numpy.repeat(int(1), common_ncross)
    with not_raises(Exception):
        selprot.nprogeny = numpy.int8(1)
    with not_raises(Exception):
        selprot.nprogeny = numpy.repeat(numpy.int8(1), common_ncross)
    with not_raises(Exception):
        selprot.nprogeny = numpy.int16(1)
    with not_raises(Exception):
        selprot.nprogeny = numpy.repeat(numpy.int16(1), common_ncross)
    with not_raises(Exception):
        selprot.nprogeny = numpy.int32(1)
    with not_raises(Exception):
        selprot.nprogeny = numpy.repeat(numpy.int32(1), common_ncross)
    with not_raises(Exception):
        selprot.nprogeny = numpy.int64(1)
    with not_raises(Exception):
        selprot.nprogeny = numpy.repeat(numpy.int64(1), common_ncross)

def test_nprogeny_fset_TypeError(selprot, common_ncross):
    with pytest.raises(TypeError):
        selprot.nprogeny = object()
    with pytest.raises(TypeError):
        selprot.nprogeny = numpy.repeat(object(), common_ncross)
    with pytest.raises(TypeError):
        selprot.nprogeny = None
    with pytest.raises(TypeError):
        selprot.nprogeny = numpy.repeat(None, common_ncross)
    with pytest.raises(TypeError):
        selprot.nprogeny = float(1)
    with pytest.raises(TypeError):
        selprot.nprogeny = numpy.repeat(float(1), common_ncross)

def test_nprogeny_fset_ValueError(selprot, common_ncross):
    with pytest.raises(ValueError):
        selprot.nprogeny = int(-1)
    with pytest.raises(ValueError):
        selprot.nprogeny = numpy.repeat(int(-1), common_ncross)

def test_nprogeny_fdel_AttributeError(selprot):
    with pytest.raises(AttributeError):
        del selprot.nprogeny

### nobj ###
def test_SelectionProtocol_nobj_is_concrete():
    assert_concrete_property(SelectionProtocol, "nobj")

def test_nobj_fget(selprot):
    assert isinstance(selprot.nobj, Integral)

def test_nobj_fset(selprot):
    with not_raises(TypeError):
        selprot.nobj = numpy.int8(1)
    with not_raises(TypeError):
        selprot.nobj = numpy.int16(1)
    with not_raises(TypeError):
        selprot.nobj = numpy.int32(1)
    with not_raises(TypeError):
        selprot.nobj = numpy.int64(1)
    with not_raises(TypeError):
        selprot.nobj = int(1)

def test_nobj_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.nobj = None
    with pytest.raises(TypeError):
        selprot.nobj = "string"
    with pytest.raises(TypeError):
        selprot.nobj = float(1.0)

def test_nobj_fset_ValueError(selprot):
    with not_raises(ValueError):
        selprot.nobj = int(1)
    with pytest.raises(ValueError):
        selprot.nobj = int(0)
    with pytest.raises(ValueError):
        selprot.nobj = int(-1)

### obj_wt ###
def test_SelectionProtocol_obj_wt_is_concrete():
    assert_concrete_property(SelectionProtocol, "obj_wt")

def test_obj_wt_fget(selprot):
    assert isinstance(selprot.obj_wt, numpy.ndarray)

def test_obj_wt_fset(selprot, common_obj_wt_single):
    with not_raises(TypeError):
        selprot.obj_wt = common_obj_wt_single
    with not_raises(TypeError):
        selprot.obj_wt = float(1.0)
    with not_raises(TypeError):
        selprot.obj_wt = None

def test_obj_wt_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.obj_wt = object()
    with pytest.raises(TypeError):
        selprot.obj_wt = "string"

def test_obj_wt_fset_ValueError(selprot, common_obj_wt_single):
    with pytest.raises(ValueError):
        selprot.obj_wt = numpy.zeros(len(common_obj_wt_single)+1, dtype=common_obj_wt_single.dtype)

### obj_trans ###
def test_SelectionProtocol_obj_trans_is_concrete():
    assert_concrete_property(SelectionProtocol, "obj_trans")

def test_obj_trans_fget(selprot):
    assert isinstance(selprot.obj_trans, Callable)

def test_obj_trans_fset(selprot, common_obj_trans_single):
    with not_raises(Exception):
        selprot.obj_trans = common_obj_trans_single

def test_obj_trans_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.obj_trans = "string"
    with pytest.raises(TypeError):
        selprot.obj_trans = []
    with pytest.raises(TypeError):
        selprot.obj_trans = {}

def test_obj_trans_fdel(selprot):
    with pytest.raises(AttributeError):
        del selprot.obj_trans

### obj_trans_kwargs ###
def test_SelectionProtocol_obj_trans_kwargs_is_concrete():
    assert_concrete_property(SelectionProtocol, "obj_trans_kwargs")

def test_obj_trans_kwargs_fget(selprot):
    assert isinstance(selprot.obj_trans_kwargs, dict)

def test_obj_trans_kwargs_fset(selprot):
    with not_raises(Exception):
        selprot.obj_trans_kwargs = {"a":1}
    with not_raises(Exception):
        selprot.obj_trans_kwargs = None

def test_obj_trans_kwargs_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.obj_trans_kwargs = "string"
    with pytest.raises(TypeError):
        selprot.obj_trans_kwargs = []
    with pytest.raises(TypeError):
        selprot.obj_trans_kwargs = ()

def test_obj_trans_kwargs_fdel(selprot):
    with pytest.raises(AttributeError):
        del selprot.obj_trans_kwargs

### nineqcv ###
def test_SelectionProtocol_nineqcv_is_concrete():
    assert_concrete_property(SelectionProtocol, "nineqcv")

def test_nineqcv_fget(selprot):
    assert isinstance(selprot.nineqcv, Integral)

def test_nineqcv_fset(selprot):
    with not_raises(TypeError):
        selprot.nineqcv = numpy.int8(1)
    with not_raises(TypeError):
        selprot.nineqcv = numpy.int16(1)
    with not_raises(TypeError):
        selprot.nineqcv = numpy.int32(1)
    with not_raises(TypeError):
        selprot.nineqcv = numpy.int64(1)
    with not_raises(TypeError):
        selprot.nineqcv = int(1)
    with not_raises(TypeError):
        selprot.nineqcv = None

def test_nineqcv_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.nineqcv = object()
    with pytest.raises(TypeError):
        selprot.nineqcv = "string"
    with pytest.raises(TypeError):
        selprot.nineqcv = float(1.0)

def test_nineqcv_fset_ValueError(selprot):
    with not_raises(ValueError):
        selprot.nineqcv = int(1)
    with not_raises(ValueError):
        selprot.nineqcv = int(0)
    with pytest.raises(ValueError):
        selprot.nineqcv = int(-1)

### ineqcv_wt ###
def test_SelectionProtocol_ineqcv_wt_is_concrete():
    assert_concrete_property(SelectionProtocol, "ineqcv_wt")

def test_ineqcv_wt_fget(selprot):
    assert isinstance(selprot.ineqcv_wt, numpy.ndarray)

def test_ineqcv_wt_fset(selprot, common_ineqcv_wt_single):
    with not_raises(TypeError):
        selprot.ineqcv_wt = common_ineqcv_wt_single
    with not_raises(TypeError):
        selprot.ineqcv_wt = float(1.0)
    with not_raises(TypeError):
        selprot.ineqcv_wt = None

def test_ineqcv_wt_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.ineqcv_wt = "string"

def test_ineqcv_wt_fset_ValueError(selprot, common_ineqcv_wt_single):
    with pytest.raises(ValueError):
        selprot.ineqcv_wt = numpy.zeros(len(common_ineqcv_wt_single)+1, dtype=common_ineqcv_wt_single.dtype)

### ineqcv_trans ###
def test_SelectionProtocol_ineqcv_trans_is_concrete():
    assert_concrete_property(SelectionProtocol, "ineqcv_trans")

def test_ineqcv_trans_fget(selprot):
    assert isinstance(selprot.ineqcv_trans, Callable)

def test_ineqcv_trans_fset(selprot, common_ineqcv_trans_single):
    with not_raises(Exception):
        selprot.ineqcv_trans = common_ineqcv_trans_single

def test_ineqcv_trans_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.ineqcv_trans = "string"
    with pytest.raises(TypeError):
        selprot.ineqcv_trans = []
    with pytest.raises(TypeError):
        selprot.ineqcv_trans = {}

def test_ineqcv_trans_fdel(selprot):
    with pytest.raises(AttributeError):
        del selprot.ineqcv_trans

### ineqcv_trans_kwargs ###
def test_SelectionProtocol_ineqcv_trans_kwargs_is_concrete():
    assert_concrete_property(SelectionProtocol, "ineqcv_trans_kwargs")

def test_ineqcv_trans_kwargs_fget(selprot):
    assert isinstance(selprot.ineqcv_trans_kwargs, dict)

def test_ineqcv_trans_kwargs_fset(selprot):
    with not_raises(Exception):
        selprot.ineqcv_trans_kwargs = {"a":1}
    with not_raises(Exception):
        selprot.ineqcv_trans_kwargs = None

def test_ineqcv_trans_kwargs_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.ineqcv_trans_kwargs = "string"
    with pytest.raises(TypeError):
        selprot.ineqcv_trans_kwargs = []
    with pytest.raises(TypeError):
        selprot.ineqcv_trans_kwargs = ()

def test_ineqcv_trans_kwargs_fdel(selprot):
    with pytest.raises(AttributeError):
        del selprot.ineqcv_trans_kwargs

### neqcv ###
def test_SelectionProtocol_neqcv_is_concrete():
    assert_concrete_property(SelectionProtocol, "neqcv")

def test_neqcv_fget(selprot):
    assert isinstance(selprot.neqcv, Integral)

def test_neqcv_fset(selprot):
    with not_raises(TypeError):
        selprot.neqcv = numpy.int8(1)
    with not_raises(TypeError):
        selprot.neqcv = numpy.int16(1)
    with not_raises(TypeError):
        selprot.neqcv = numpy.int32(1)
    with not_raises(TypeError):
        selprot.neqcv = numpy.int64(1)
    with not_raises(TypeError):
        selprot.neqcv = int(1)
    with not_raises(TypeError):
        selprot.neqcv = None

def test_neqcv_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.neqcv = "string"
    with pytest.raises(TypeError):
        selprot.neqcv = float(1.0)

def test_neqcv_fset_ValueError(selprot):
    with not_raises(ValueError):
        selprot.neqcv = int(1)
    with not_raises(ValueError):
        selprot.neqcv = int(0)
    with pytest.raises(ValueError):
        selprot.neqcv = int(-1)

### eqcv_wt ###
def test_SelectionProtocol_eqcv_wt_is_concrete():
    assert_concrete_property(SelectionProtocol, "eqcv_wt")

def test_eqcv_wt_fget(selprot):
    assert isinstance(selprot.eqcv_wt, numpy.ndarray)

def test_eqcv_wt_fset(selprot, common_eqcv_wt_single):
    with not_raises(TypeError):
        selprot.eqcv_wt = common_eqcv_wt_single
    with not_raises(TypeError):
        selprot.eqcv_wt = float(1.0)
    with not_raises(TypeError):
        selprot.eqcv_wt = None

def test_eqcv_wt_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.eqcv_wt = object()
    with pytest.raises(TypeError):
        selprot.eqcv_wt = "string"

def test_eqcv_wt_fset_ValueError(selprot, common_eqcv_wt_single):
    with pytest.raises(ValueError):
        selprot.eqcv_wt = numpy.zeros(len(common_eqcv_wt_single)+1, dtype=common_eqcv_wt_single.dtype)

### eqcv_trans ###
def test_SelectionProtocol_eqcv_trans_is_concrete():
    assert_concrete_property(SelectionProtocol, "eqcv_trans")

def test_eqcv_trans_fget(selprot):
    assert isinstance(selprot.eqcv_trans, Callable)

def test_eqcv_trans_fset(selprot, common_eqcv_trans_single):
    with not_raises(Exception):
        selprot.eqcv_trans = common_eqcv_trans_single

def test_eqcv_trans_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.eqcv_trans = "string"
    with pytest.raises(TypeError):
        selprot.eqcv_trans = []
    with pytest.raises(TypeError):
        selprot.eqcv_trans = {}

def test_eqcv_trans_fdel(selprot):
    with pytest.raises(AttributeError):
        del selprot.eqcv_trans

### eqcv_trans_kwargs ###
def test_SelectionProtocol_eqcv_trans_kwargs_is_concrete():
    assert_concrete_property(SelectionProtocol, "eqcv_trans_kwargs")

def test_eqcv_trans_kwargs_fget(selprot):
    assert isinstance(selprot.eqcv_trans_kwargs, dict)

def test_eqcv_trans_kwargs_fset(selprot):
    with not_raises(Exception):
        selprot.eqcv_trans_kwargs = {"a":1}
    with not_raises(Exception):
        selprot.eqcv_trans_kwargs = None

def test_eqcv_trans_kwargs_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.eqcv_trans_kwargs = "string"
    with pytest.raises(TypeError):
        selprot.eqcv_trans_kwargs = []
    with pytest.raises(TypeError):
        selprot.eqcv_trans_kwargs = ()

def test_eqcv_trans_kwargs_fdel(selprot):
    with pytest.raises(AttributeError):
        del selprot.eqcv_trans_kwargs

### ndset_wt ###
def test_SelectionProtocol_ndset_wt_is_concrete():
    assert_concrete_property(SelectionProtocol, "ndset_wt")

def test_ndset_wt_fget(selprot):
    assert isinstance(selprot.ndset_wt, Real)

def test_ndset_wt_fset(selprot, common_ndset_wt_single):
    with not_raises(TypeError):
        selprot.ndset_wt = common_ndset_wt_single
    with not_raises(TypeError):
        selprot.ndset_wt = float(1.0)
    with not_raises(TypeError):
        selprot.ndset_wt = None

def test_ndset_wt_fset_TypeError(selprot, common_ndset_wt_single):
    with pytest.raises(TypeError):
        selprot.ndset_wt = object()
    with pytest.raises(TypeError):
        selprot.ndset_wt = "string"
    with pytest.raises(TypeError):
        selprot.ndset_wt = numpy.array([common_ndset_wt_single])

def test_ndset_wt_fset_ValueError(selprot):
    with pytest.raises(ValueError):
        selprot.ndset_wt = 0

### ndset_trans ###
def test_SelectionProtocol_ndset_trans_is_concrete():
    assert_concrete_property(SelectionProtocol, "ndset_trans")

def test_ndset_trans_fget(selprot):
    assert isinstance(selprot.ndset_trans, Callable)

def test_ndset_trans_fset(selprot, common_ndset_trans_single):
    with not_raises(Exception):
        selprot.ndset_trans = common_ndset_trans_single

def test_ndset_trans_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.ndset_trans = object()
    with pytest.raises(TypeError):
        selprot.ndset_trans = "string"
    with pytest.raises(TypeError):
        selprot.ndset_trans = []
    with pytest.raises(TypeError):
        selprot.ndset_trans = {}

def test_ndset_trans_fdel(selprot):
    with pytest.raises(AttributeError):
        del selprot.ndset_trans

### ndset_trans_kwargs ###
def test_SelectionProtocol_ndset_trans_kwargs_is_concrete():
    assert_concrete_property(SelectionProtocol, "ndset_trans_kwargs")

def test_ndset_trans_kwargs_fget(selprot):
    assert isinstance(selprot.ndset_trans_kwargs, dict)

def test_ndset_trans_kwargs_fset(selprot):
    with not_raises(Exception):
        selprot.ndset_trans_kwargs = {"a":1}
    with not_raises(Exception):
        selprot.ndset_trans_kwargs = None

def test_ndset_trans_kwargs_fset_TypeError(selprot):
    with pytest.raises(TypeError):
        selprot.ndset_trans_kwargs = object()
    with pytest.raises(TypeError):
        selprot.ndset_trans_kwargs = "string"
    with pytest.raises(TypeError):
        selprot.ndset_trans_kwargs = []
    with pytest.raises(TypeError):
        selprot.ndset_trans_kwargs = ()

def test_ndset_trans_kwargs_fdel(selprot):
    with pytest.raises(AttributeError):
        del selprot.ndset_trans_kwargs

############################# Test class utilities #############################
def test_check_is_SelectionProtocol(selprot):
    with not_raises(Exception):
        check_is_SelectionProtocol(selprot, "selprot")
    with pytest.raises(TypeError):
        check_is_SelectionProtocol(object(), "selprot")
    with pytest.raises(TypeError):
        check_is_SelectionProtocol(None, "selprot")
