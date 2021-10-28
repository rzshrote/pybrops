from . import generic_check_is_not
from . import generic_check_dict_keys
from . import generic_check_len_eq
from . import generic_check_gteq

from . import generic_default_cond
from . import generic_cond_check_is_not
from . import generic_cond_check_dict_keys
from . import generic_cond_check_len_eq
from . import generic_cond_check_gteq
from . import generic_check_keys_in_dict_all

################################################################################
############################### check functions ################################
################################################################################
def check_is_not_None(v, vname):
    generic_check_is_not(v, vname, None, "None")

def check_len(v, vname, n):
    if len(v) != n:
        raise ValueError("the length of '{0}' is not equal to {1}".format(vname,n))

def check_len_eq(v, vname, w, wname):
    generic_check_len_eq(v, vname, w, wname)

def check_all_equal(v, vname):
    viter = iter(v)
    try:
        e0 = next(viter)
    except StopIteration:
        return # length is 1, therefore all elements are equivalent
    if any(e0 != e for e in viter):
        raise ValueError("not all elements in {0} are equal to {1}".format(vname, e0))

def check_is_positive(v, vname):
    generic_check_gteq(v, vname, 0)

##################################################
########### Dictionary check functions ###########
##################################################
def check_keys_in_dict(v, vname, *args):
    keys = v.keys()
    if any(e not in keys for e in args):
        raise ValueError("dict '{0}' must have keys: {1}".format(vname, args))

def check_keys_in_dict_all_type(v, vname, vtype):
    if any(not isinstance(e, vtype) for e in v.keys()):
        raise ValueError("not all keys in dict '{0}' are of type {1}".format(vname, str(vtype)))

def check_values_in_dict_all_type(v, vname, vtype):
    if any(not isinstance(e, vtype) for e in v.values()):
        raise ValueError("not all values in dict '{0}' are of type {1}".format(vname, str(vtype)))

def check_values_in_dict_equal_len(v, vname):
    viter = iter(v.values())
    try:
        e0 = next(viter)
    except StopIteration:
        return
    l0 = len(e0)
    if any(len(e) != l0 for e in viter):
        raise ValueError("not all values in dict '{0}' have equal length == {1}".format(vname, l0))

def check_values_in_dict_len(v, vname, l):
    viter = iter(v.values())
    try:
        e0 = next(viter)
    except StopIteration:
        raise ValueError("dict '{0}' is empty".format(vname))
    if any(len(e) != l for e in viter):
        raise ValueError("not all values in dict '{0}' have length == {1}".format(vname, l))

################################################################################
######################### conditional check functions ##########################
################################################################################
def cond_check_is_not_None(v, vname, cond = generic_default_cond):
    generic_cond_check_is_not(v, vname, None, "None", cond)

def cond_check_keys_in_dict(v, vname, *args, cond = generic_default_cond):
    generic_cond_check_dict_keys(v, vname, args, [str(e) for e in args], cond)

def cond_check_len(v, vname, n, cond = generic_default_cond):
    generic_cond_check_len(v, vname, n, cond)

def cond_check_len_eq(v, vname, w, wname, cond = generic_default_cond):
    generic_cond_check_len_eq(v, vname, w, wname, cond)

def cond_check_all_equal(v, vname, cond = generic_default_cond):
    if cond(v):
        check_all_equal(v, vname)

def cond_check_is_positive(v, vname, cond = generic_default_cond):
    generic_cond_check_gteq(v, vname, 0, cond)
