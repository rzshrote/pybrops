import numpy
import math

def cond_str_lower(s, cond=(lambda s: isinstance(s, str))):
    return s.lower() if cond(s) else s

def srange(start, stop, step):
    yield from range(start, stop, step)
    yield stop

def matrix_is_sorted(mat):
    for i in range(mat.size-1):
        if mat[i] > mat[i+1]:
            return False
    return True

def cond_len(var, cond=(lambda var: hasattr(var, "__len__"))):
    if cond(var):
        return len(var)
    return None

def slice_to_range(s, end):
    ifnone = lambda x, y: y if x is None else x
    r = range(ifnone(s.start, 0), ifnone(s.end, end), ifnone(s.step, 1))
    return r

def slice_to_list(s, end):
    l = list(e for e in slice_to_range(s, end))
    return

# define symbol conversion coefficients
HUMAN2BYTES_DICT = {
    # for SI system
    'B'  : 10**0,     'byte' : 10**0,     'K'  : 10**3,       'kilo' : 10**3,
    'M'  : 10**6,     'mega' : 10**6,     'G'  : 10**9,       'giga' : 10**9,
    'T'  : 10**12,    'tera' : 10**12,    'P'  : 10**15,      'peta' : 10**15,
    'E'  : 10**18,    'exa'  : 10**18,    'Z'  : 10**21,      'zetta': 10**21,
    'Y'  : 10**24,    'iotta': 10**24,
    # for IEC system
    'Bi' : 2**0  /8,  'bit'  : 2**0  /8,  'Ki' : 2**10 /8,   'kibi' : 2**10 /8,
    'Mi' : 2**20 /8,  'mebi' : 2**20 /8,  'Gi' : 2**30 /8,   'gibi' : 2**30 /8,
    'Ti' : 2**40 /8,  'tebi' : 2**40 /8,  'Pi' : 2**50 /8,   'pebi' : 2**50 /8,
    'Ei' : 2**60 /8,  'exbi' : 2**60 /8,  'Zi' : 2**70 /8,   'zebi' : 2**70 /8,
    'Yi' : 2**80 /8,  'yobi' : 2**80 /8,
    # if only a number string is provided, assume bytes
    ''   : 1
}

def human2bytes(s):
    # strip any white space before or after
    s = s.strip()

    # make some variables
    sym = ''    # variable to store the symbol
    num = 0     # variable to count index and later store the parsed number

    # for each character in string 's'
    for c in s:
        # if c is not a digit nor a decimal point,
        # then we've hit the symbol
        if not c.isdigit() and c != '.':
            break
        # increment index counter
        num += 1

    # extract symbol and strip it of whitespace
    sym = s[num:].strip()

    # convert number string to floating point
    num = float(s[:num]) if num > 0 else 0

    # if the parsed symbol is not in the list of keys, raise an error
    if sym not in HUMAN2BYTES_DICT.keys():
        raise ValueError(
            "Cannot interpret string '%s':\n"\
            "    Negative numbers are not allowed.\n"\
            "    Available symbol keys:\n"\
            "        %s" % (s,HUMAN2BYTES_DICT.keys())
        )

    # return number of bytes needed rounded up
    return int(math.ceil(num * HUMAN2BYTES_DICT[sym]))
