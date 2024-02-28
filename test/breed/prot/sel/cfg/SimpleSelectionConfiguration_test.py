import pytest
from pybrops.breed.prot.sel.cfg.SimpleSelectionConfiguration import SimpleSelectionConfiguration
from pybrops.test.assert_python import assert_class_isconcrete, assert_method_isconcrete, assert_property_isconcrete, assert_class_documentation, assert_class_issemiabstract, not_raises

from .common_fixtures import *

@pytest.fixture
def selcfg(
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_pgmat,
        common_xconfig
    ):
    out = SimpleSelectionConfiguration(
        ncross = common_ncross,
        nparent = common_nparent,
        nmating = common_nmating,
        nprogeny = common_nprogeny,
        pgmat = common_pgmat,
        xconfig = common_xconfig
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_SimpleSelectionConfiguration_is_semiabstract():
    assert_class_isconcrete(SimpleSelectionConfiguration)

############################## Test class docstring ############################
def test_SimpleSelectionConfiguration_docstring():
    assert_class_documentation(SimpleSelectionConfiguration)

############################# Test class properties ############################

### __init__ ###
def test_SimpleSelectionConfiguration_init_is_concrete():
    assert_method_isconcrete(SimpleSelectionConfiguration, "__init__")

def test_init(selcfg):
    assert selcfg is not None
