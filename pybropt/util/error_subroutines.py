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


################################################################################

def error_readonly(varname):
    raise AttributeError("'%s' is read-only." % varname)

################################################################################
