### two-way crosses
# generators
from pybrops.core.util.crossix import twowayix_symab_anyab
from pybrops.core.util.crossix import twowayix_symab_selfab
from pybrops.core.util.crossix import twowayix_symab_uniqab
from pybrops.core.util.crossix import twowayix_asymab_anyab
from pybrops.core.util.crossix import twowayix_asymab_selfab
from pybrops.core.util.crossix import twowayix_asymab_uniqab
from pybrops.core.util.crossix import twowayix
# lengths
from pybrops.core.util.crossix import twowayix_symab_anyab_len
from pybrops.core.util.crossix import twowayix_symab_selfab_len
from pybrops.core.util.crossix import twowayix_symab_uniqab_len
from pybrops.core.util.crossix import twowayix_asymab_anyab_len
from pybrops.core.util.crossix import twowayix_asymab_selfab_len
from pybrops.core.util.crossix import twowayix_asymab_uniqab_len
from pybrops.core.util.crossix import twowayix_len
### three-way crosses
# generators
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
from pybrops.core.util.crossix import threewayix
# lengths
from pybrops.core.util.crossix import threewayix_asymab_anyab_anyc_len
from pybrops.core.util.crossix import threewayix_asymab_anyab_backc_len
from pybrops.core.util.crossix import threewayix_asymab_anyab_uniqc_len
from pybrops.core.util.crossix import threewayix_asymab_selfab_anyc_len
from pybrops.core.util.crossix import threewayix_asymab_selfab_backc_len
from pybrops.core.util.crossix import threewayix_asymab_selfab_uniqc_len
from pybrops.core.util.crossix import threewayix_asymab_uniqab_anyc_len
from pybrops.core.util.crossix import threewayix_asymab_uniqab_backc_len
from pybrops.core.util.crossix import threewayix_asymab_uniqab_uniqc_len
from pybrops.core.util.crossix import threewayix_symab_anyab_anyc_len
from pybrops.core.util.crossix import threewayix_symab_anyab_backc_len
from pybrops.core.util.crossix import threewayix_symab_anyab_uniqc_len
from pybrops.core.util.crossix import threewayix_symab_selfab_anyc_len
from pybrops.core.util.crossix import threewayix_symab_selfab_backc_len
from pybrops.core.util.crossix import threewayix_symab_selfab_uniqc_len
from pybrops.core.util.crossix import threewayix_symab_uniqab_anyc_len
from pybrops.core.util.crossix import threewayix_symab_uniqab_backc_len
from pybrops.core.util.crossix import threewayix_symab_uniqab_uniqc_len
from pybrops.core.util.crossix import threewayix_len
# test functions
from pybrops.test.assert_python import assert_function_isconcrete
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

######################## Generators ########################

### twowayix_symab_anyab
def test_twowayix_symab_anyab_is_concrete():
    assert_generator_isconcrete(twowayix_symab_anyab)

def test_twowayix_symab_anyab():
    for i in range(10):
        n = sum(1 for _ in twowayix_symab_anyab(i))
        m = twowayix_symab_anyab_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### twowayix_symab_selfab
def test_twowayix_symab_selfab_is_concrete():
    assert_generator_isconcrete(twowayix_symab_selfab)

def test_twowayix_symab_selfab():
    for i in range(10):
        n = sum(1 for _ in twowayix_symab_selfab(i))
        m = twowayix_symab_selfab_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### twowayix_symab_uniqab
def test_twowayix_symab_uniqab_is_concrete():
    assert_generator_isconcrete(twowayix_symab_uniqab)

def test_twowayix_symab_uniqab():
    for i in range(10):
        n = sum(1 for _ in twowayix_symab_uniqab(i))
        m = twowayix_symab_uniqab_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### twowayix_asymab_anyab
def test_twowayix_asymab_anyab_is_concrete():
    assert_generator_isconcrete(twowayix_asymab_anyab)

def test_twowayix_asymab_anyab():
    for i in range(10):
        n = sum(1 for _ in twowayix_asymab_anyab(i))
        m = twowayix_asymab_anyab_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### twowayix_asymab_selfab
def test_twowayix_asymab_selfab_is_concrete():
    assert_generator_isconcrete(twowayix_asymab_selfab)

def test_twowayix_asymab_selfab():
    for i in range(10):
        n = sum(1 for _ in twowayix_asymab_selfab(i))
        m = twowayix_asymab_selfab_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### twowayix_asymab_uniqab
def test_twowayix_asymab_uniqab_is_concrete():
    assert_generator_isconcrete(twowayix_asymab_uniqab)

def test_twowayix_asymab_uniqab():
    for i in range(10):
        n = sum(1 for _ in twowayix_asymab_uniqab(i))
        m = twowayix_asymab_uniqab_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### twowayix
def test_twowayix_is_concrete():
    assert_generator_isconcrete(twowayix)

################### Length of Generators ###################

### twowayix_symab_anyab_len
def test_twowayix_symab_anyab_len_is_concrete():
    assert_function_isconcrete(twowayix_symab_anyab_len)

def test_twowayix_symab_anyab_len():
    assert twowayix_symab_anyab_len(-1) == 0
    assert twowayix_symab_anyab_len(0) == 0
    assert isinstance(twowayix_symab_anyab_len(4), int)

### twowayix_symab_selfab_len
def test_twowayix_symab_selfab_len_is_concrete():
    assert_function_isconcrete(twowayix_symab_selfab_len)

def test_twowayix_symab_selfab_len():
    assert twowayix_symab_selfab_len(-1) == 0
    assert twowayix_symab_selfab_len(0) == 0
    assert isinstance(twowayix_symab_selfab_len(4), int)

### twowayix_symab_uniqab
def test_twowayix_symab_uniqab_len_is_concrete():
    assert_function_isconcrete(twowayix_symab_uniqab_len)

def test_twowayix_symab_uniqab_len():
    assert twowayix_symab_uniqab_len(-1) == 0
    assert twowayix_symab_uniqab_len(0) == 0
    assert isinstance(twowayix_symab_uniqab_len(4), int)

### twowayix_asymab_anyab_len
def test_twowayix_asymab_anyab_len_is_concrete():
    assert_function_isconcrete(twowayix_asymab_anyab_len)

def test_twowayix_asymab_anyab_len():
    assert twowayix_asymab_anyab_len(-1) == 0
    assert twowayix_asymab_anyab_len(0) == 0
    assert isinstance(twowayix_asymab_anyab_len(4), int)

### twowayix_asymab_selfab_len
def test_twowayix_asymab_selfab_len_is_concrete():
    assert_function_isconcrete(twowayix_asymab_selfab_len)

def test_twowayix_asymab_selfab_len():
    assert twowayix_asymab_selfab_len(-1) == 0
    assert twowayix_asymab_selfab_len(0) == 0
    assert isinstance(twowayix_asymab_selfab_len(4), int)

### twowayix_asymab_uniqab_len
def test_twowayix_asymab_uniqab_len_is_concrete():
    assert_function_isconcrete(twowayix_asymab_uniqab_len)

def test_twowayix_asymab_uniqab_len():
    assert twowayix_asymab_uniqab_len(-1) == 0
    assert twowayix_asymab_uniqab_len(0) == 0
    assert isinstance(twowayix_asymab_uniqab_len(4), int)

### twowayix
def test_twowayix_len_is_concrete():
    assert_function_isconcrete(twowayix_len)

############################################################
############# Three-way cross index generators #############
############################################################

######################## Generators ########################

### threewayix_symab_anyab_anyc
def test_threewayix_symab_anyab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_anyab_anyc)

def test_threewayix_symab_anyab_anyc():
    for i in range(10):
        n = sum(1 for _ in threewayix_symab_anyab_anyc(i))
        m = threewayix_symab_anyab_anyc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_symab_anyab_backc
def test_threewayix_symab_anyab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_anyab_backc)

def test_threewayix_symab_anyab_backc():
    for i in range(10):
        n = sum(1 for _ in threewayix_symab_anyab_backc(i))
        m = threewayix_symab_anyab_backc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_symab_anyab_uniqc
def test_threewayix_symab_anyab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_anyab_uniqc)

def test_threewayix_symab_anyab_uniqc():
    for i in range(10):
        n = sum(1 for _ in threewayix_symab_anyab_uniqc(i))
        m = threewayix_symab_anyab_uniqc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_symab_selfab_anyc
def test_threewayix_symab_selfab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_selfab_anyc)

def test_threewayix_symab_selfab_anyc():
    for i in range(10):
        n = sum(1 for _ in threewayix_symab_selfab_anyc(i))
        m = threewayix_symab_selfab_anyc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_symab_selfab_backc
def test_threewayix_symab_selfab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_selfab_backc)

def test_threewayix_symab_selfab_backc():
    for i in range(10):
        n = sum(1 for _ in threewayix_symab_selfab_backc(i))
        m = threewayix_symab_selfab_backc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_symab_selfab_uniqc
def test_threewayix_symab_selfab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_selfab_uniqc)

def test_threewayix_symab_selfab_uniqc():
    for i in range(10):
        n = sum(1 for _ in threewayix_symab_selfab_uniqc(i))
        m = threewayix_symab_selfab_uniqc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_symab_uniqab_anyc
def test_threewayix_symab_uniqab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_uniqab_anyc)

def test_threewayix_symab_uniqab_anyc():
    for i in range(10):
        n = sum(1 for _ in threewayix_symab_uniqab_anyc(i))
        m = threewayix_symab_uniqab_anyc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_symab_uniqab_backc
def test_threewayix_symab_uniqab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_uniqab_backc)

def test_threewayix_symab_uniqab_backc():
    for i in range(10):
        n = sum(1 for _ in threewayix_symab_uniqab_backc(i))
        m = threewayix_symab_uniqab_backc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_symab_uniqab_uniqc
def test_threewayix_symab_uniqab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_symab_uniqab_uniqc)

def test_threewayix_symab_uniqab_uniqc():
    for i in range(10):
        n = sum(1 for _ in threewayix_symab_uniqab_uniqc(i))
        m = threewayix_symab_uniqab_uniqc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_asymab_anyab_anyc
def test_threewayix_asymab_anyab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_anyab_anyc)

def test_threewayix_asymab_anyab_anyc():
    for i in range(10):
        n = sum(1 for _ in threewayix_asymab_anyab_anyc(i))
        m = threewayix_asymab_anyab_anyc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_asymab_anyab_backc
def test_threewayix_asymab_anyab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_anyab_backc)

def test_threewayix_asymab_anyab_backc():
    for i in range(10):
        n = sum(1 for _ in threewayix_asymab_anyab_backc(i))
        m = threewayix_asymab_anyab_backc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_asymab_anyab_uniqc
def test_threewayix_asymab_anyab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_anyab_uniqc)

def test_threewayix_asymab_anyab_uniqc():
    for i in range(10):
        n = sum(1 for _ in threewayix_asymab_anyab_uniqc(i))
        m = threewayix_asymab_anyab_uniqc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_asymab_selfab_anyc
def test_threewayix_asymab_selfab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_selfab_anyc)

def test_threewayix_asymab_selfab_anyc():
    for i in range(10):
        n = sum(1 for _ in threewayix_asymab_selfab_anyc(i))
        m = threewayix_asymab_selfab_anyc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_asymab_selfab_backc
def test_threewayix_asymab_selfab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_selfab_backc)

def test_threewayix_asymab_selfab_backc():
    for i in range(10):
        n = sum(1 for _ in threewayix_asymab_selfab_backc(i))
        m = threewayix_asymab_selfab_backc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_asymab_selfab_uniqc
def test_threewayix_asymab_selfab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_selfab_uniqc)

def test_threewayix_asymab_selfab_uniqc():
    for i in range(10):
        n = sum(1 for _ in threewayix_asymab_selfab_uniqc(i))
        m = threewayix_asymab_selfab_uniqc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_asymab_uniqab_anyc
def test_threewayix_asymab_uniqab_anyc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_uniqab_anyc)

def test_threewayix_asymab_uniqab_anyc():
    for i in range(10):
        n = sum(1 for _ in threewayix_asymab_uniqab_anyc(i))
        m = threewayix_asymab_uniqab_anyc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_asymab_uniqab_backc
def test_threewayix_asymab_uniqab_backc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_uniqab_backc)

def test_threewayix_asymab_uniqab_backc():
    for i in range(10):
        n = sum(1 for _ in threewayix_asymab_uniqab_backc(i))
        m = threewayix_asymab_uniqab_backc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix_asymab_uniqab_uniqc
def test_threewayix_asymab_uniqab_uniqc_is_concrete():
    assert_generator_isconcrete(threewayix_asymab_uniqab_uniqc)

def test_threewayix_asymab_uniqab_uniqc():
    for i in range(10):
        n = sum(1 for _ in threewayix_asymab_uniqab_uniqc(i))
        m = threewayix_asymab_uniqab_uniqc_len(i)
        if n != m:
            print("Received:", n, "Expected:", m)
        assert n == m

### threewayix
def test_threewayix_is_concrete():
    assert_generator_isconcrete(threewayix)

################### Length of Generators ###################

### threewayix_symab_anyab_anyc_len
def test_threewayix_symab_anyab_anyc_len_is_concrete():
    assert_function_isconcrete(threewayix_symab_anyab_anyc_len)

def test_threewayix_symab_anyab_anyc_len():
    assert threewayix_symab_anyab_anyc_len(-1) == 0
    assert threewayix_symab_anyab_anyc_len(0) == 0
    assert isinstance(threewayix_symab_anyab_anyc_len(4), int)

### threewayix_symab_anyab_backc_len
def test_threewayix_symab_anyab_backc_len_is_concrete():
    assert_function_isconcrete(threewayix_symab_anyab_backc_len)

def test_threewayix_symab_anyab_backc_len():
    assert threewayix_symab_anyab_backc_len(-1) == 0
    assert threewayix_symab_anyab_backc_len(0) == 0
    assert isinstance(threewayix_symab_anyab_backc_len(4), int)

### threewayix_symab_anyab_uniqc_len
def test_threewayix_symab_anyab_uniqc_len_is_concrete():
    assert_function_isconcrete(threewayix_symab_anyab_uniqc_len)

def test_threewayix_symab_anyab_uniqc_len():
    assert threewayix_symab_anyab_uniqc_len(-1) == 0
    assert threewayix_symab_anyab_uniqc_len(0) == 0
    assert isinstance(threewayix_symab_anyab_uniqc_len(4), int)

### threewayix_symab_selfab_anyc_len
def test_threewayix_symab_selfab_anyc_len_is_concrete():
    assert_function_isconcrete(threewayix_symab_selfab_anyc_len)

def test_threewayix_symab_selfab_anyc_len():
    assert threewayix_symab_selfab_anyc_len(-1) == 0
    assert threewayix_symab_selfab_anyc_len(0) == 0
    assert isinstance(threewayix_symab_selfab_anyc_len(4), int)

### threewayix_symab_selfab_backc_len
def test_threewayix_symab_selfab_backc_len_is_concrete():
    assert_function_isconcrete(threewayix_symab_selfab_backc_len)

def test_threewayix_symab_selfab_backc_len():
    assert threewayix_symab_selfab_backc_len(-1) == 0
    assert threewayix_symab_selfab_backc_len(0) == 0
    assert isinstance(threewayix_symab_selfab_backc_len(4), int)

### threewayix_symab_selfab_uniqc_len
def test_threewayix_symab_selfab_uniqc_len_is_concrete():
    assert_function_isconcrete(threewayix_symab_selfab_uniqc_len)

def test_threewayix_symab_selfab_uniqc_len():
    assert threewayix_symab_selfab_uniqc_len(-1) == 0
    assert threewayix_symab_selfab_uniqc_len(0) == 0
    assert isinstance(threewayix_symab_selfab_uniqc_len(4), int)

### threewayix_symab_uniqab_anyc_len
def test_threewayix_symab_uniqab_anyc_len_is_concrete():
    assert_function_isconcrete(threewayix_symab_uniqab_anyc_len)

def test_threewayix_symab_uniqab_anyc_len():
    assert threewayix_symab_uniqab_anyc_len(-1) == 0
    assert threewayix_symab_uniqab_anyc_len(0) == 0
    assert isinstance(threewayix_symab_uniqab_anyc_len(4), int)

### threewayix_symab_uniqab_backc_len
def test_threewayix_symab_uniqab_backc_len_is_concrete():
    assert_function_isconcrete(threewayix_symab_uniqab_backc_len)

def test_threewayix_symab_uniqab_backc_len():
    assert threewayix_symab_uniqab_backc_len(-1) == 0
    assert threewayix_symab_uniqab_backc_len(0) == 0
    assert isinstance(threewayix_symab_uniqab_backc_len(4), int)

### threewayix_symab_uniqab_uniqc_len
def test_threewayix_symab_uniqab_uniqc_len_is_concrete():
    assert_function_isconcrete(threewayix_symab_uniqab_uniqc_len)

def test_threewayix_symab_uniqab_uniqc_len():
    assert threewayix_symab_uniqab_uniqc_len(-1) == 0
    assert threewayix_symab_uniqab_uniqc_len(0) == 0
    assert isinstance(threewayix_symab_uniqab_uniqc_len(4), int)

### threewayix_asymab_anyab_anyc_len
def test_threewayix_asymab_anyab_anyc_len_is_concrete():
    assert_function_isconcrete(threewayix_asymab_anyab_anyc_len)

def test_threewayix_asymab_anyab_anyc_len():
    assert threewayix_asymab_anyab_anyc_len(-1) == 0
    assert threewayix_asymab_anyab_anyc_len(0) == 0
    assert isinstance(threewayix_asymab_anyab_anyc_len(4), int)

### threewayix_asymab_anyab_backc_len
def test_threewayix_asymab_anyab_backc_len_is_concrete():
    assert_function_isconcrete(threewayix_asymab_anyab_backc_len)

def test_threewayix_asymab_anyab_backc_len():
    assert threewayix_asymab_anyab_backc_len(-1) == 0
    assert threewayix_asymab_anyab_backc_len(0) == 0
    assert isinstance(threewayix_asymab_anyab_backc_len(4), int)

### threewayix_asymab_anyab_uniqc_len
def test_threewayix_asymab_anyab_uniqc_len_is_concrete():
    assert_function_isconcrete(threewayix_asymab_anyab_uniqc_len)

def test_threewayix_asymab_anyab_uniqc_len():
    assert threewayix_asymab_anyab_uniqc_len(-1) == 0
    assert threewayix_asymab_anyab_uniqc_len(0) == 0
    assert isinstance(threewayix_asymab_anyab_uniqc_len(4), int)

### threewayix_asymab_selfab_anyc_len
def test_threewayix_asymab_selfab_anyc_len_is_concrete():
    assert_function_isconcrete(threewayix_asymab_selfab_anyc_len)

def test_threewayix_asymab_selfab_anyc_len():
    assert threewayix_asymab_selfab_anyc_len(-1) == 0
    assert threewayix_asymab_selfab_anyc_len(0) == 0
    assert isinstance(threewayix_asymab_selfab_anyc_len(4), int)

### threewayix_asymab_selfab_backc_len
def test_threewayix_asymab_selfab_backc_len_is_concrete():
    assert_function_isconcrete(threewayix_asymab_selfab_backc_len)

def test_threewayix_asymab_selfab_backc_len():
    assert threewayix_asymab_selfab_backc_len(-1) == 0
    assert threewayix_asymab_selfab_backc_len(0) == 0
    assert isinstance(threewayix_asymab_selfab_backc_len(4), int)

### threewayix_asymab_selfab_uniqc_len
def test_threewayix_asymab_selfab_uniqc_len_is_concrete():
    assert_function_isconcrete(threewayix_asymab_selfab_uniqc_len)

def test_threewayix_asymab_selfab_uniqc_len():
    assert threewayix_asymab_selfab_uniqc_len(-1) == 0
    assert threewayix_asymab_selfab_uniqc_len(0) == 0
    assert isinstance(threewayix_asymab_selfab_uniqc_len(4), int)

### threewayix_asymab_uniqab_anyc_len
def test_threewayix_asymab_uniqab_anyc_len_is_concrete():
    assert_function_isconcrete(threewayix_asymab_uniqab_anyc_len)

def test_threewayix_asymab_uniqab_anyc_len():
    assert threewayix_asymab_uniqab_anyc_len(-1) == 0
    assert threewayix_asymab_uniqab_anyc_len(0) == 0
    assert isinstance(threewayix_asymab_uniqab_anyc_len(4), int)

### threewayix_asymab_uniqab_backc_len
def test_threewayix_asymab_uniqab_backc_len_is_concrete():
    assert_function_isconcrete(threewayix_asymab_uniqab_backc_len)

def test_threewayix_asymab_uniqab_backc_len():
    assert threewayix_asymab_uniqab_backc_len(-1) == 0
    assert threewayix_asymab_uniqab_backc_len(0) == 0
    assert isinstance(threewayix_asymab_uniqab_backc_len(4), int)

### threewayix_asymab_uniqab_uniqc_len
def test_threewayix_asymab_uniqab_uniqc_len_is_concrete():
    assert_function_isconcrete(threewayix_asymab_uniqab_uniqc_len)

def test_threewayix_asymab_uniqab_uniqc_len():
    assert threewayix_asymab_uniqab_uniqc_len(-1) == 0
    assert threewayix_asymab_uniqab_uniqc_len(0) == 0
    assert isinstance(threewayix_asymab_uniqab_uniqc_len(4), int)

### threewayix_len
def test_threewayix_len_is_concrete():
    assert_function_isconcrete(threewayix_len)
