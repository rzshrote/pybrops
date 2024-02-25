import pytest
from pybrops.breed.prot.sel.cfg.SimpleMateSelectionConfiguration import SimpleMateSelectionConfiguration
from pybrops.test.assert_python import assert_class_isconcrete, assert_property_isconcrete, assert_class_documentation, assert_class_issemiabstract, not_raises

from .common_fixtures import *

@pytest.fixture
def selcfg(
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_pgmat,
        common_xconfig,
        common_xconfig_xmap
    ):
    out = SimpleMateSelectionConfiguration(
        ncross = common_ncross,
        nparent = common_nparent,
        nmating = common_nmating,
        nprogeny = common_nprogeny,
        pgmat = common_pgmat,
        xconfig = common_xconfig,
        xconfig_xmap = common_xconfig_xmap
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_SimpleMateSelectionConfiguration_is_semiabstract():
    assert_class_isconcrete(SimpleMateSelectionConfiguration)

############################## Test class docstring ############################
def test_SimpleMateSelectionConfiguration_docstring():
    assert_class_documentation(SimpleMateSelectionConfiguration)

############################# Test class properties ############################

### __init__ ###
def test_SimpleMateSelectionConfiguration_init_is_concrete():
    assert_property_isconcrete(SimpleMateSelectionConfiguration, "__init__")

def test_init(selcfg):
    assert selcfg is not None
