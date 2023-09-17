import inspect
import pytest
import numpy

from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap
from pybrops.popgen.gmap.ExtendedGeneticMap import check_is_ExtendedGeneticMap

################################################################################
################################### Fixtures ###################################
################################################################################

@pytest.fixture
def egmap(shared_datadir):
    yield ExtendedGeneticMap.from_egmap(shared_datadir / "McMullen_2009_US_NAM.M.egmap")

@pytest.fixture
def egmap_sample(shared_datadir):
    yield ExtendedGeneticMap.from_egmap(shared_datadir / "sample.egmap")

################################################################################
#################################### Tests #####################################
################################################################################
def test_check_is_ExtendedGeneticMap():
    with pytest.raises(TypeError):
        check_is_ExtendedGeneticMap(None, "None")

def test_len(egmap):
    assert len(egmap) == 1122

def test_lexsort_None(egmap):
    a = egmap.lexsort()
    b = numpy.arange(1122)
    assert numpy.all(a == b)

def test_reorder(egmap):
    chrgrp = egmap.vrnt_chrgrp.copy()
    phypos = egmap.vrnt_phypos.copy()
    stop = egmap.vrnt_stop.copy()
    genpos = egmap.vrnt_genpos.copy()

    a = numpy.arange(1122)
    numpy.random.shuffle(a)
    egmap.reorder(a)

    assert numpy.all(egmap.vrnt_chrgrp == chrgrp[a])
    assert numpy.all(egmap.vrnt_phypos == phypos[a])
    assert numpy.all(egmap.vrnt_stop == stop[a])
    assert numpy.all(egmap.vrnt_genpos == genpos[a])

def test_sort_None(egmap):
    chrgrp = egmap.vrnt_chrgrp.copy()
    phypos = egmap.vrnt_phypos.copy()
    stop = egmap.vrnt_stop.copy()
    genpos = egmap.vrnt_genpos.copy()

    egmap.sort()

    assert numpy.all(egmap.vrnt_chrgrp == chrgrp)
    assert numpy.all(egmap.vrnt_phypos == phypos)
    assert numpy.all(egmap.vrnt_stop == stop)
    assert numpy.all(egmap.vrnt_genpos == genpos)

def test_group(egmap):
    egmap.group()
    assert egmap.is_grouped()

def test_build_spline(egmap_sample):
    egmap_sample.group()
    egmap_sample.build_spline()
    assert egmap_sample.has_spline()

def test_interp_genpos(egmap_sample):
    egmap_sample.group()
    egmap_sample.build_spline(kind = 'linear', fill_value = 'extrapolate')
    chrgrp = numpy.arange(5) + 1
    phypos = numpy.repeat(2, 5)
    out = numpy.repeat(0.75, 5)
    a = egmap_sample.interp_genpos(chrgrp, phypos)
    assert numpy.all(a == out)

# TODO: test for raises warning on interpolation
# def test_gmap_inversion(shared_datadir):
#     data_path = shared_datadir / "McMullen_2009_US_NAM.M.egmap"
#     with pytest.raises(ValueError):
#         m = ExtendedGeneticMap.from_egmap(data_path)
