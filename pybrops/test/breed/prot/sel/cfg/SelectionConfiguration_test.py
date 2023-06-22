from pybrops.test.assert_python import assert_abstract_class
from pybrops.breed.prot.sel.cfg.SelectionConfiguration import SelectionConfiguration
from pybrops.breed.prot.sel.cfg.SelectionConfiguration import check_is_SelectionConfiguration


################### Test class abstract/concrete properties ####################
def test_SelectionConfiguration_is_abstract():
    assert_abstract_class(SelectionConfiguration)

############################## Test class docstring ############################
