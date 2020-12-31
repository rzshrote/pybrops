from . import generic_check_is_not
from . import generic_check_dict_keys
from . import generic_check_len_eq
from . import generic_cond_check_is_not
from . import generic_cond_check_dict_keys
from . import generic_cond_check_len_eq

################################################################################
############################### check functions ################################
################################################################################
def check_is_not_None(v, vname):
    generic_check_is_not(v, vname, None, "None")

def check_keys_in_dict(v, vname, *args):
    generic_check_dict_keys(v, vname, args, [str(e) for e in args])

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

################################################################################
######################### conditional check functions ##########################
################################################################################
def cond_check_is_not_None(v, vname, cond=(lambda s: s is not None)):
    generic_cond_check_is_not(v, vname, None, "None", cond)

def cond_check_keys_in_dict(v, vname, *args, cond=(lambda s: s is not None)):
    generic_cond_check_dict_keys(v, vname, args, [str(e) for e in args], cond)

def cond_check_len_eq(v, vname, w, wname, cond=(lambda s: s is not None)):
    generic_cond_check_len_eq(v, vname, w, wname, cond)

def cond_check_all_equal(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_all_equal(v, vname)
