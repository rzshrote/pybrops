# import 3rd party modules
import numpy

def check_divisibility(a, aname, b, bname):
    if a % b != 0:
        raise ValueError(
            "'%s' (%d) is not divisible by '%s' (%d)." % (aname, a, bname, b)
        )

def check_keys_in_dict(d, varname, *argv):
    key_absent = [arg not in d for arg in argv]
    if any(key_absent):
        # build error string
        err_str = "'%s' does not have the required fields.\n" % varname
        for i,absent in enumerate(key_absent):
            if absent:
                err_str += "    %s is missing.\n" % argv[i]
        raise ValueError(err_str)

def check_equal_len(a, aname, b, bname):
    if len(a) != len(b):
        raise ValueError("len(%s) != len(%s)" % (aname, bname))

def check_all_equal(v, varname):
    n = v[0]
    if any(e != n for e in v):
        raise ValueError("not all elements in %s are equivalent to %s" % (varname, n))

# def check_equal_len(d):
#     for key in d.keys():
#
#     key = d.keys()
#     val = d.values()
#     pass

################################################################################

def error_readonly(varname):
    raise AttributeError("'%s' is read-only." % varname)

################################################################################
