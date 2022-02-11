import pytest
import numpy

from pybropt.core.util.haplo import calc_nhaploblk_chrom
from pybropt.core.util.haplo import calc_haplobin
from pybropt.core.util.haplo import calc_haplobin_bounds

################################################################################
################################ Test fixtures #################################
################################################################################
# nhaploblk, genpos, chrgrp_stix, chrgrp_spix
# nhaploblk_chrom, genpos, chrgrp_stix, chrgrp_spix
# haplobin

@pytest.fixture
def nhaploblk():
    yield 5

@pytest.fixture
def genpos():
    a = numpy.array([
         0.10,  1.35,  1.56,  2.10,  2.15,  2.72,  3.04,    # chromosome 1
        -0.49, -0.06,  0.59,  0.81,                         # chromosome 2
        -0.18, -0.04,  0.24,  0.25,  1.04,  1.63            # chromosome 3
    ])
    yield a

@pytest.fixture
def chrgrp_stix():
    a = numpy.array([0, 7, 11])
    yield a

@pytest.fixture
def chrgrp_spix():
    a = numpy.array([7, 11, 17])
    yield a

@pytest.fixture
def nhaploblk_chrom():
    a = numpy.array([2, 1, 2])
    yield a

@pytest.fixture
def haplobin():
    a = numpy.array([
        0, 0, 0, 1, 1, 1, 1,    # chromosome 1
        2, 2, 2, 2,             # chromosome 2
        3, 3, 3, 3, 4, 4        # chromosome 3
    ])
    yield a

@pytest.fixture
def hstix():
    a = numpy.array([0, 3, 7, 11, 15])
    yield a

@pytest.fixture
def hspix():
    a = numpy.array([3, 7, 11, 15, 17])
    yield a

@pytest.fixture
def hlen():
    a = numpy.array([3, 4, 4, 4, 2])
    yield a

@pytest.fixture
def haplobin_long():
    a = numpy.array([
        0, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
    ])
    yield a

@pytest.fixture
def hstix_long():
    a = numpy.array([0, 5, 12, 18])
    yield a

@pytest.fixture
def hspix_long():
    a = numpy.array([5, 12, 18, 29])
    yield a

@pytest.fixture
def hlen_long():
    a = numpy.array([5, 7, 6, 11])
    yield a

################################################################################
############################ Test module functions #############################
################################################################################
def test_calc_nhaploblk_chrom(nhaploblk, genpos, chrgrp_stix, chrgrp_spix, nhaploblk_chrom):
    out = calc_nhaploblk_chrom(nhaploblk, genpos, chrgrp_stix, chrgrp_spix)
    assert numpy.all(out == nhaploblk_chrom)

def test_calc_haplobin(nhaploblk_chrom, genpos, chrgrp_stix, chrgrp_spix, haplobin):
    out = calc_haplobin(nhaploblk_chrom, genpos, chrgrp_stix, chrgrp_spix)
    assert numpy.all(out == haplobin)

def test_calc_haplobin_bounds(haplobin, hstix, hspix, hlen):
    out = calc_haplobin_bounds(haplobin)
    assert numpy.all(out[0] == hstix)
    assert numpy.all(out[1] == hspix)
    assert numpy.all(out[2] == hlen)

def test_calc_haplobin_bounds_long(haplobin_long, hstix_long, hspix_long, hlen_long):
    out = calc_haplobin_bounds(haplobin_long)
    assert numpy.all(out[0] == hstix_long)
    assert numpy.all(out[1] == hspix_long)
    assert numpy.all(out[2] == hlen_long)
