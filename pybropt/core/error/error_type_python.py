from . import generic_check_isinstance
from . import generic_cond_check_isinstance

################################################################################
################## basic check functions for basic data types ##################
################################################################################
def check_is_bool(v, vname):
    generic_check_isinstance(v, vname, bool)

# def check_is_bytearray(v, vname):
#     generic_check_isinstance(v, vname, bytearray)

# def check_is_bytes(v, vname):
#     generic_check_isinstance(v, vname, bytes)

# def check_is_complex(v, vname):
#     generic_check_isinstance(v, vname, complex)

def check_is_dict(v, vname):
    generic_check_isinstance(v, vname, dict)

def check_is_float(v, vname):
    generic_check_isinstance(v, vname, float)

# def check_is_frozenset(v, vname):
#     generic_check_isinstance(v, vname, frozenset)

def check_is_int(v, vname):
    generic_check_isinstance(v, vname, int)

def check_is_list(v, vname):
    generic_check_isinstance(v, vname, list)

# def check_is_memoryview(v, vname):
#     generic_check_isinstance(v, vname, memoryview)

def check_is_range(v, vname):
    generic_check_isinstance(v, vname, range)

def check_is_set(v, vname):
    generic_check_isinstance(v, vname, set)

def check_is_str(v, vname):
    generic_check_isinstance(v, vname, str)

def check_is_tuple(v, vname):
    generic_check_isinstance(v, vname, tuple)

def check_is_(v, vname):
    generic_check_isinstance(v, vname, )

################################################################################
################ compound check functions for basic data types #################
################################################################################
def check_is_str_or_iterable(v, vname):
    if not (isinstance(s, str) or hasattr(s, "__iter__")):
        raise TypeError("'{0}' must be of type str or have attribute __iter__.".format(vname))

def check_is_list_or_tuple(v, vname):
    generic_check_isinstance(v, vname, (list,tuple))


################################################################################
######################### conditional check functions ##########################
################################################################################
def cond_check_is_bool(v, vname):
    generic_cond_check_isinstance(v, vname, bool)

# def cond_check_is_bytearray(v, vname):
#     generic_cond_check_isinstance(v, vname, bytearray)

# def cond_check_is_bytes(v, vname):
#     generic_cond_check_isinstance(v, vname, bytes)

# def cond_check_is_complex(v, vname):
#     generic_cond_check_isinstance(v, vname, complex)

def cond_check_is_dict(v, vname):
    generic_cond_check_isinstance(v, vname, dict)

def cond_check_is_float(v, vname):
    generic_cond_check_isinstance(v, vname, float)

# def cond_check_is_frozenset(v, vname):
#     generic_cond_check_isinstance(v, vname, frozenset)

def cond_check_is_int(v, vname):
    generic_cond_check_isinstance(v, vname, int)

def cond_check_is_list(v, vname):
    generic_cond_check_isinstance(v, vname, list)

# def cond_check_is_memoryview(v, vname):
#     generic_cond_check_isinstance(v, vname, memoryview)

def cond_check_is_range(v, vname):
    generic_cond_check_isinstance(v, vname, range)

def cond_check_is_set(v, vname):
    generic_cond_check_isinstance(v, vname, set)

def cond_check_is_str(v, vname):
    generic_cond_check_isinstance(v, vname, str)

def cond_check_is_tuple(v, vname):
    generic_cond_check_isinstance(v, vname, tuple)


################################################################################
########## conditional compound check functions for basic data types ###########
################################################################################

def cond_check_is_str_or_iterable(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_str_or_iterable(v, vname)
