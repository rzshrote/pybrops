from pybrops.core.util.crossix import threewayix
from pybrops.core.util.crossix import threewayix_asymab_anyab_anyc
from pybrops.core.util.crossix import threewayix_asymab_anyab_backc
from pybrops.core.util.crossix import threewayix_asymab_anyab_uniqc
from pybrops.core.util.crossix import threewayix_asymab_selfab_anyc
from pybrops.core.util.crossix import threewayix_asymab_selfab_backc
from pybrops.core.util.crossix import threewayix_asymab_selfab_uniqc
from pybrops.core.util.crossix import threewayix_asymab_uniqab_anyc
from pybrops.core.util.crossix import threewayix_asymab_uniqab_backc
from pybrops.core.util.crossix import threewayix_asymab_uniqab_uniqc
from pybrops.core.util.crossix import threewayix_symab_anyab_anyc
from pybrops.core.util.crossix import threewayix_symab_anyab_backc
from pybrops.core.util.crossix import threewayix_symab_anyab_uniqc
from pybrops.core.util.crossix import threewayix_symab_selfab_anyc
from pybrops.core.util.crossix import threewayix_symab_selfab_backc
from pybrops.core.util.crossix import threewayix_symab_selfab_uniqc
from pybrops.core.util.crossix import threewayix_symab_uniqab_anyc
from pybrops.core.util.crossix import threewayix_symab_uniqab_backc
from pybrops.core.util.crossix import threewayix_symab_uniqab_uniqc
from pybrops.core.util.crossix import twowayix
from pybrops.core.util.crossix import twowayix_asymab_selfab
from pybrops.core.util.crossix import twowayix_symab_selfab
from pybrops.core.util.crossix import twowayix_symab_anyab
from pybrops.core.util.crossix import twowayix_symab_uniqab
from pybrops.core.util.crossix import twowayix_asymab_anyab
from pybrops.core.util.crossix import twowayix_asymab_uniqab
from pybrops.test.assert_python import assert_generator_isconcrete
from pybrops.test.assert_python import assert_module_documentation
from pybrops.test.assert_python import assert_module_public_api

################################################################################
################################ Test fixtures #################################
################################################################################

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_error_value_numpy_module_documentation():
    import pybrops.core.util.crossix
    assert_module_documentation(pybrops.core.util.crossix)

def test_error_value_numpy_module_public_api():
    import pybrops.core.util.crossix
    assert_module_public_api(pybrops.core.util.crossix)

################################################################################
############################ Test module functions #############################
################################################################################

############################################################
############# Two-way cross index generators ###############
############################################################

### twowayix_symab_anyab
def test_twowayix_symab_anyab_is_concrete():
    assert_generator_isconcrete(twowayix_asymab_anyab)

def test_twowayix_symab_anyab():
    assert sum(1 for _ in twowayix_asymab_anyab(0)) == (0**2)
    assert sum(1 for _ in twowayix_asymab_anyab(0)) == (0**2)
    assert sum(1 for _ in twowayix_asymab_anyab(1)) == (1**2)
    assert sum(1 for _ in twowayix_asymab_anyab(2)) == (2**2)
    assert sum(1 for _ in twowayix_asymab_anyab(3)) == (3**2)
    assert sum(1 for _ in twowayix_asymab_anyab(4)) == (4**2)
    assert sum(1 for _ in twowayix_asymab_anyab(5)) == (5**2)
    assert sum(1 for _ in twowayix_asymab_anyab(6)) == (6**2)
    assert sum(1 for _ in twowayix_asymab_anyab(7)) == (7**2)
    assert sum(1 for _ in twowayix_asymab_anyab(8)) == (8**2)
    assert sum(1 for _ in twowayix_asymab_anyab(9)) == (9**2)

### twowayix_symab_selfab
def test_twowayix_symab_selfab_is_concrete():
    assert_generator_isconcrete(twowayix_symab_selfab)

### twowayix_symab_uniqab
def test_twowayix_symab_uniqab_is_concrete():
    assert_generator_isconcrete(twowayix_asymab_uniqab)

def test_twowayix_symab_uniqab():
    assert sum(1 for _ in twowayix_asymab_uniqab(0)) == (0**2 - 0)
    assert sum(1 for _ in twowayix_asymab_uniqab(1)) == (1**2 - 1)
    assert sum(1 for _ in twowayix_asymab_uniqab(2)) == (2**2 - 2)
    assert sum(1 for _ in twowayix_asymab_uniqab(3)) == (3**2 - 3)
    assert sum(1 for _ in twowayix_asymab_uniqab(4)) == (4**2 - 4)
    assert sum(1 for _ in twowayix_asymab_uniqab(5)) == (5**2 - 5)
    assert sum(1 for _ in twowayix_asymab_uniqab(6)) == (6**2 - 6)
    assert sum(1 for _ in twowayix_asymab_uniqab(7)) == (7**2 - 7)
    assert sum(1 for _ in twowayix_asymab_uniqab(8)) == (8**2 - 8)
    assert sum(1 for _ in twowayix_asymab_uniqab(9)) == (9**2 - 9)

### twowayix_asymab_anyab
def test_twowayix_asymab_anyab_is_concrete():
    assert_generator_isconcrete(twowayix_symab_anyab)

def test_twowayix_asymab_anyab():
    assert sum(1 for _ in twowayix_symab_anyab(0)) == ((0**2 - 0)/2 + 0)
    assert sum(1 for _ in twowayix_symab_anyab(1)) == ((1**2 - 1)/2 + 1)
    assert sum(1 for _ in twowayix_symab_anyab(2)) == ((2**2 - 2)/2 + 2)
    assert sum(1 for _ in twowayix_symab_anyab(3)) == ((3**2 - 3)/2 + 3)
    assert sum(1 for _ in twowayix_symab_anyab(4)) == ((4**2 - 4)/2 + 4)
    assert sum(1 for _ in twowayix_symab_anyab(5)) == ((5**2 - 5)/2 + 5)
    assert sum(1 for _ in twowayix_symab_anyab(6)) == ((6**2 - 6)/2 + 6)
    assert sum(1 for _ in twowayix_symab_anyab(7)) == ((7**2 - 7)/2 + 7)
    assert sum(1 for _ in twowayix_symab_anyab(8)) == ((8**2 - 8)/2 + 8)
    assert sum(1 for _ in twowayix_symab_anyab(9)) == ((9**2 - 9)/2 + 9)

### twowayix_asymab_selfab
def test_twowayix_asymab_selfab_is_concrete():
    assert_generator_isconcrete(twowayix_asymab_selfab)

### twowayix_asymab_uniqab
def test_twowayix_asymab_uniqab_is_concrete():
    assert_generator_isconcrete(twowayix_symab_uniqab)

def test_twowayix_asymab_uniqab():
    assert sum(1 for _ in twowayix_symab_uniqab(0)) == ((0**2 - 0)/2)
    assert sum(1 for _ in twowayix_symab_uniqab(1)) == ((1**2 - 1)/2)
    assert sum(1 for _ in twowayix_symab_uniqab(2)) == ((2**2 - 2)/2)
    assert sum(1 for _ in twowayix_symab_uniqab(3)) == ((3**2 - 3)/2)
    assert sum(1 for _ in twowayix_symab_uniqab(4)) == ((4**2 - 4)/2)
    assert sum(1 for _ in twowayix_symab_uniqab(5)) == ((5**2 - 5)/2)
    assert sum(1 for _ in twowayix_symab_uniqab(6)) == ((6**2 - 6)/2)
    assert sum(1 for _ in twowayix_symab_uniqab(7)) == ((7**2 - 7)/2)
    assert sum(1 for _ in twowayix_symab_uniqab(8)) == ((8**2 - 8)/2)
    assert sum(1 for _ in twowayix_symab_uniqab(9)) == ((9**2 - 9)/2)

### twowayix
def test_twowayix_is_concrete():
    assert_generator_isconcrete(twowayix)

############################################################
############# Three-way cross index generators #############
############################################################

### threewayix_symab_anyab_anyc
def test_threewayix_symab_anyab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_anyab_anyc)

### threewayix_symab_anyab_backc
def test_threewayix_symab_anyab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_anyab_backc)

### threewayix_symab_anyab_uniqc
def test_threewayix_symab_anyab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_anyab_uniqc)

### threewayix_symab_selfab_anyc
def test_threewayix_symab_selfab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_selfab_anyc)

### threewayix_symab_selfab_backc
def test_threewayix_symab_selfab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_selfab_backc)

### threewayix_symab_selfab_uniqc
def test_threewayix_symab_selfab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_selfab_uniqc)

### threewayix_symab_uniqab_anyc
def test_threewayix_symab_uniqab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_uniqab_anyc)

### threewayix_symab_uniqab_backc
def test_threewayix_symab_uniqab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_uniqab_backc)

### threewayix_symab_uniqab_uniqc
def test_threewayix_symab_uniqab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_uniqab_uniqc)

### threewayix_asymab_anyab_anyc
def test_threewayix_asymab_anyab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_anyab_anyc)

def test_threewayix_asymab_anyab_anyc():
    assert sum(1 for _ in threewayix_asymab_anyab_anyc(0)) == (0**3)
    assert sum(1 for _ in threewayix_asymab_anyab_anyc(1)) == (1**3)
    assert sum(1 for _ in threewayix_asymab_anyab_anyc(2)) == (2**3)
    assert sum(1 for _ in threewayix_asymab_anyab_anyc(3)) == (3**3)
    assert sum(1 for _ in threewayix_asymab_anyab_anyc(4)) == (4**3)
    assert sum(1 for _ in threewayix_asymab_anyab_anyc(5)) == (5**3)
    assert sum(1 for _ in threewayix_asymab_anyab_anyc(6)) == (6**3)
    assert sum(1 for _ in threewayix_asymab_anyab_anyc(7)) == (7**3)
    assert sum(1 for _ in threewayix_asymab_anyab_anyc(8)) == (8**3)
    assert sum(1 for _ in threewayix_asymab_anyab_anyc(9)) == (9**3)

### threewayix_asymab_anyab_backc
def test_threewayix_asymab_anyab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_anyab_backc)

def test_threewayix_asymab_anyab_backc():
    assert sum(1 for _ in threewayix_asymab_anyab_backc(0)) == (2 * 0**2)
    assert sum(1 for _ in threewayix_asymab_anyab_backc(1)) == (2 * 1**2)
    assert sum(1 for _ in threewayix_asymab_anyab_backc(2)) == (2 * 2**2)
    assert sum(1 for _ in threewayix_asymab_anyab_backc(3)) == (2 * 3**2)
    assert sum(1 for _ in threewayix_asymab_anyab_backc(4)) == (2 * 4**2)
    assert sum(1 for _ in threewayix_asymab_anyab_backc(5)) == (2 * 5**2)
    assert sum(1 for _ in threewayix_asymab_anyab_backc(6)) == (2 * 6**2)
    assert sum(1 for _ in threewayix_asymab_anyab_backc(7)) == (2 * 7**2)
    assert sum(1 for _ in threewayix_asymab_anyab_backc(8)) == (2 * 8**2)
    assert sum(1 for _ in threewayix_asymab_anyab_backc(9)) == (2 * 9**2)

### threewayix_asymab_anyab_uniqc
def test_threewayix_asymab_anyab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_anyab_uniqc)

def test_threewayix_asymab_anyab_uniqc():
    assert sum(1 for _ in threewayix_asymab_anyab_uniqc(0)) == 0
    assert sum(1 for _ in threewayix_asymab_anyab_uniqc(1)) == 0
    assert sum(1 for _ in threewayix_asymab_anyab_uniqc(2)) == 0
    assert sum(1 for _ in threewayix_asymab_anyab_uniqc(3)) == ((3 - 2) * 3**2)
    assert sum(1 for _ in threewayix_asymab_anyab_uniqc(4)) == ((4 - 2) * 4**2)
    assert sum(1 for _ in threewayix_asymab_anyab_uniqc(5)) == ((5 - 2) * 5**2)
    assert sum(1 for _ in threewayix_asymab_anyab_uniqc(6)) == ((6 - 2) * 6**2)
    assert sum(1 for _ in threewayix_asymab_anyab_uniqc(7)) == ((7 - 2) * 7**2)
    assert sum(1 for _ in threewayix_asymab_anyab_uniqc(8)) == ((8 - 2) * 8**2)
    assert sum(1 for _ in threewayix_asymab_anyab_uniqc(9)) == ((9 - 2) * 9**2)

### threewayix_asymab_selfab_anyc
def test_threewayix_asymab_selfab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_selfab_anyc)

### threewayix_asymab_selfab_backc
def test_threewayix_asymab_selfab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_selfab_backc)

### threewayix_asymab_selfab_uniqc
def test_threewayix_asymab_selfab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_selfab_uniqc)

### threewayix_asymab_uniqab_anyc
def test_threewayix_asymab_uniqab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_uniqab_anyc)

### threewayix_asymab_uniqab_backc
def test_threewayix_asymab_uniqab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_uniqab_backc)

### threewayix_asymab_uniqab_uniqc
def test_threewayix_asymab_uniqab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_uniqab_uniqc)

### threewayix
def test_threewayix_is_concrete():
    assert_generator_isconcrete(threewayix)

