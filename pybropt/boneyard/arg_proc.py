# The purpose of this file is to collect a common set of arguments processing
# protocols so that it is easier to write and manage genomic selection or
# genomic mating functions.

# import libraries
import numpy
import warnings
# import our libraries
from .human2bytes import human2bytes

def proc_pop_sel_sizes(pop_size, sel_size, geno):
    """
    Process selection size argument.

    If pop_size is an integer type, return it.
    Otherwise, raise an error.

    If sel_size is an integer type, return it.
    If sel_size is a floating type, multiply by pop size and return an int.
    Otherwise, raise and error.
    """
    # process 'pop_size' and make sure it's the right data type
    if not numpy.issubdtype(type(pop_size), numpy.integer):
        raise TypeError(
            "Expected 'pop_size' to be a integer type.\n"\
            "    type(pop_size) = %s" % (type(pop_size),)
        )

    # process 'sel_size' and make sure it's the right data type
    if numpy.issubdtype(type(sel_size), numpy.integer):
        pass
    elif numpy.issubdtype(type(sel_size), numpy.floating):
        sel_size = numpy.int(numpy.around(pop_size * sel_size))
    else:
        raise TypeError(
            "'sel_size' must be an integer or floating point type.\n"\
            "    type(sel_size) = %s" % (type(sel_size),)
        )

    # throw an error if the selection size is >= to the number of individuals
    if sel_size >= geno.shape[1]:
        raise ValueError(
            "'sel_size' must be less than the number of indiv:\n"\
            "    sel_size = %s\n"\
            "    geno.shape = (%s, indiv=%s, %s)\n" %
            ((sel_size,) + geno.shape)
        )

    # test if pop_size and sel_size are sane
    if pop_size < sel_size:
        raise ValueError(
            "'pop_size' must be greater than 'sel_size'.\n"\
            "    pop_size = %s\n"\
            "    sel_size = %s\n" % (pop_size, sel_size)
        )

    # return modified pop_size, sel_size as a tuple
    return pop_size, sel_size

def proc_geno_coeff(geno, coeff):
    """
    Process 'geno' and 'coeff' matrices.

    Check for right data types.
    Check for right matrix dimensions.
    Check for dimension alignment.
    Returns: Nothing
    """
    # check for 'geno' as numpy.ndarray
    if not isinstance(geno, numpy.ndarray):
        raise TypeError(
            "'geno' must be of type numpy.ndarray\n"\
            "    type(geno) = %s" % (type(geno),)
        )

    # check for 'coeff' as numpy.ndarray
    if not isinstance(coeff, numpy.ndarray):
        raise TypeError(
            "'coeff' must be of type numpy.ndarray\n"\
            "    type(coeff) = %s" % (type(coeff),)
        )

    # check for 'uint8' dtype inside 'geno'
    if geno.dtype != 'uint8':
        raise TypeError(
            "'geno' must have a dtype of 'uint8' for efficiency.\n"\
            "    geno.dtype = %s" % geno.dtype.name
        )

    # check for 'float64' dtype inside 'coeff'
    if coeff.dtype != 'float64':
        coeff = numpy.float64(coeff)
        warnings.warn(
            "'coeff' dtype is not 'float64'.\n    coeff.dtype = %s\n"\
            "Increasing precision to 'float64'." % coeff.dtype.name,
            UserWarning, stacklevel = 2
        )

    # check for correct genotype shape
    if len(geno.shape) != 3:
        raise TypeError(
            "'geno' must have 3 dimensions: (phase, indiv, locus)\n"\
            "    geno.shape = %s" % (geno.shape,)
        )

    # check for correct coeff shape
    if len(coeff.shape) != 1:
        raise TypeError(
            "'coeff' must have 1 dimension: (locus,)\n"\
            "    coeff.shape = %s" % (coeff.shape,)
        )

    # check if loci lengths are the same
    if geno.shape[2] != coeff.shape[0]:
        raise ValueError(
            "Length of 'coeff' does not match 'geno':\n"\
            "    geno.shape = (%s, %s, locus=%s)\n"\
            "    coeff.shape = (locus=%s,)\n"\
            "'locus' dimensions need to match." %
            (geno.shape + (coeff.shape[0],))
        )

    return geno, coeff

def proc_wcoeff(wcoeff, coeff, divisor = None):
    if wcoeff is None:
        wcoeff = numpy.absolute(coeff)  # take absolute value of coeff
        if divisor is None:             # if no divisor
            divisor = wcoeff.sum()      # set divisor to sum(wcoeff)
        wcoeff = wcoeff / divisor       # divide wcoeff by divisor

    return wcoeff

def proc_tfreq(tfreq):
    # check if tfreq is a numpy.ndarray, if it is increase precision if needed.
    if isinstance(tfreq, numpy.ndarray):
        if tfreq.dtype != 'float64':
            tfreq = numpy.float64(tfreq)
            # courteously warn the user
            warnings.warn(
                "'tfreq' dtype is not 'float64'.\n    tfreq.dtype = %s\n"\
                "Increasing precision to 'float64'." % tfreq.dtype.name,
                UserWarning, stacklevel = 2
            )
    # check if tfreq is a number, if it is increase precision if necessary
    elif numpy.issubdtype(type(tfreq), numpy.number):
        if not numpy.issubdtype(type(tfreq), numpy.float64):
            tfreq = numpy.float64(tfreq)
            # courteously warn the user
            warnings.warn(
                "'tfreq' is not a 'numpy.float64'.\n    type(tfreq) = %s\n"\
                "Increasing precision to 'numpy.float64'." % (type(tfreq),),
                UserWarning, stacklevel = 2
            )
    # else tfreq is not an interpretable data type.
    else:
        raise TypeError(
            "'tfreq' must be of type numpy.ndarray or numpy.number.\n"\
            "    type(tfreq) = %s" % (type(tfreq),)
        )

    return tfreq

def proc_tld(tld):
    # check if tld is a numpy.ndarray, if it is increase precision if needed.
    if isinstance(tld, numpy.ndarray):
        if tld.dtype != 'float64':
            tld = numpy.float64(tld)
            # courteously warn the user
            warnings.warn(
                "'tld' dtype is not 'float64'.\n    tld.dtype = %s\n"\
                "Increasing precision to 'float64'." % tld.dtype.name,
                UserWarning, stacklevel = 2
            )
    # check if tld is a number, if it is increase precision if necessary
    elif numpy.issubdtype(type(tld), numpy.number):
        if not numpy.issubdtype(type(tld), numpy.float64):
            tld = numpy.float64(tld)
            # courteously warn the user
            warnings.warn(
                "'tld' is not a 'numpy.float64'.\n    type(tld) = %s\n"\
                "Increasing precision to 'numpy.float64'." % (type(tld),),
                UserWarning, stacklevel = 2
            )
    # else tfreq is not an interpretable data type.
    else:
        raise TypeError(
            "'tld' must be of type numpy.ndarray or numpy.number.\n"\
            "    type(tld) = %s" % (type(tld),)
        )

    return tld

def proc_tfreq_tld(tfreq, tld):
    ################# check for correct 'tfreq' type and dtype #################
    # check if tfreq is a numpy.ndarray, if it is increase precision if needed.
    if isinstance(tfreq, numpy.ndarray):
        if tfreq.dtype != 'float64':
            tfreq = numpy.float64(tfreq)
            # courteously warn the user
            warnings.warn(
                "'tfreq' dtype is not 'float64'.\n    tfreq.dtype = %s\n"\
                "Increasing precision to 'float64'." % tfreq.dtype.name,
                UserWarning, stacklevel = 2
            )
    # check if tfreq is a number, if it is increase precision if necessary
    elif numpy.issubdtype(type(tfreq), numpy.number):
        if not numpy.issubdtype(type(tfreq), numpy.float64):
            tfreq = numpy.float64(tfreq)
            # courteously warn the user
            warnings.warn(
                "'tfreq' is not a 'numpy.float64'.\n    type(tfreq) = %s\n"\
                "Increasing precision to 'numpy.float64'." % (type(tfreq),),
                UserWarning, stacklevel = 2
            )
    # else tfreq is not an interpretable data type.
    else:
        raise TypeError(
            "'tfreq' must be of type numpy.ndarray or numpy.number.\n"\
            "    type(tfreq) = %s" % (type(tfreq),)
        )

    ################## check for correct 'tld' type and dtype ##################
    # check if tld is a numpy.ndarray, if it is increase precision if needed.
    if isinstance(tld, numpy.ndarray):
        if tld.dtype != 'float64':
            tld = numpy.float64(tld)
            # courteously warn the user
            warnings.warn(
                "'tld' dtype is not 'float64'.\n    tld.dtype = %s\n"\
                "Increasing precision to 'float64'." % tld.dtype.name,
                UserWarning, stacklevel = 2
            )
    # check if tld is a number, if it is increase precision if necessary
    elif numpy.issubdtype(type(tld), numpy.number):
        if not numpy.issubdtype(type(tld), numpy.float64):
            tld = numpy.float64(tld)
            # courteously warn the user
            warnings.warn(
                "'tld' is not a 'numpy.float64'.\n    type(tld) = %s\n"\
                "Increasing precision to 'numpy.float64'." % (type(tld),),
                UserWarning, stacklevel = 2
            )
    # else tfreq is not an interpretable data type.
    else:
        raise TypeError(
            "'tld' must be of type numpy.ndarray or numpy.number.\n"\
            "    type(tld) = %s" % (type(tld),)
        )

    return tfreq, tld

def proc_gmap_lgroup(gmap, lgroup):
    # check for 'gmap' as numpy.ndarray
    if not isinstance(gmap, numpy.ndarray):
        raise TypeError(
            "'gmap' must be of type numpy.ndarray\n"\
            "    type(gmap) = %s" % (type(gmap),)
        )
    # else 'gmap' is a numpy.ndarray; check its dtype
    elif numpy.issubdtype(gmap.dtype,numpy.number) and (gmap.dtype!='float64'):
        gmap = numpy.float64(gmap)
        warnings.warn(
            "'gmap' dtype is not 'float64'.\n    gmap.dtype = %s\n"\
            "Increasing precision to 'float64'." % gmap.dtype.name,
            UserWarning, stacklevel = 2
        )

    # check for 'lgroup' as numpy.ndarray
    if not isinstance(lgroup, numpy.ndarray):
        raise TypeError(
            "'lgroup' must be of type numpy.ndarray\n"\
            "    type(lgroup) = %s" % (type(lgroup),)
        )
    # else 'gmap' is a numpy.ndarray; check its dtype
    elif not numpy.issubdtype(lgroup.dtype, numpy.integer):
        raise TypeError(
            "'lgroup' dtype must be an integer type.\n"\
            "    lgroup.dtype = %s" % lgroup.dtype.name
        )

    # check if 'gmap' and 'lgroup' align
    if len(gmap) != lgroup.sum():
        raise ValueError(
            "The 'gmap' and 'lgroup' do not align:\n"\
            "    len(gmap)    = %s\n"\
            "    lgroup.sum() = %s\n"
            % (len(gmap), lgroup.sum())
        )

    return gmap, lgroup

def proc_dcoeff(dcoeff, n):
    # check that 'dcoeff' is a numpy.ndarray
    if not isinstance(dcoeff, numpy.ndarray):
        raise TypeError(
            "'dcoeff' must be a numpy.ndarray.\n"\
            "    type(dcoeff) = %s" % (type(dcoeff),)
        )

    # check that 'dcoeff' dtype is compatible
    if not numpy.issubdtype(dcoeff.dtype, numpy.floating):
        raise TypeError(
            "'dcoeff' must be a numpy.integer subdtype.\n"\
            "    type(dcoeff) = %s" % (type(dcoeff),)
        )

    # check that 'dcoeff' has the correct number of dimensions
    if dcoeff.ndim != 1:
        raise TypeError(
            "'dcoeff' must have 1 dimension.\n"\
            "    dcoeff.shape = %s" % (dcoeff.shape,)
        )

    # check that dcoeff has three elements
    if dcoeff.size != n:
        raise TypeError(
            "'dcoeff' must have %s elements.\n"\
            "    dcoeff.size = %s" % (n,dcoeff.size)
        )

    # dcoeff.sum() does not have to be equal to 1.0
    return dcoeff

def proc_dcoeff_t(dcoeff_t, bcycles, n):
    # check that 'dcoeff_t' is a numpy.ndarray
    if not isinstance(dcoeff_t, numpy.ndarray):
        raise TypeError(
            "'dcoeff_t' must be a numpy.ndarray.\n"\
            "    type(dcoeff_t) = %s" % (type(dcoeff_t),)
        )

    # check that 'dcoeff_t' dtype is compatible
    if not numpy.issubdtype(dcoeff_t.dtype, numpy.floating):
        raise TypeError(
            "'dcoeff_t' must be a numpy.integer subdtype.\n"\
            "    type(dcoeff_t) = %s" % (type(dcoeff_t),)
        )

    # check that 'dcoeff' has the correct number of dimensions
    if dcoeff_t.ndim != 2:
        raise TypeError(
            "'dcoeff_t' must have 2 dimensions.\n"\
            "    dcoeff_t.shape = %s" % (dcoeff_t.shape,)
        )

    # check that dcoeff has the right shape
    if dcoeff_t.shape != (bcycles, n):
        raise TypeError(
            "'dcoeff_t' must have shape %s.\n"\
            "    dcoeff_t.shape = %s" % ((bcycles,n),dcoeff_t.shape)
        )

    # dcoeff_t.sum() does not have to be equal to 1.0
    return dcoeff_t

def proc_bcycles(bcycles):
    if not numpy.issubdtype(type(bcycles), numpy.integer):
        raise TypeError(
            "Expected 'bcycles' to be a numpy.integer type.\n"\
            "    type(bcycles) = %s" % (type(bcycles),)
        )

    return bcycles

def proc_mem(mem, gmap):
    # if 'mem' is None, set to compute whole matrix
    if mem is None:
        mem = len(gmap)
    # if 'mem' is a string, parse it and calculate chunk size
    elif isinstance(mem, str):
        mem = numpy.int(numpy.floor(numpy.sqrt( # sqrt, floor, cast as int
            human2bytes(mem) /                  # parse number of required bytes
            numpy.dtype('float64').itemsize     # divide by dtype size in bytes
        )))
    # if 'mem' is not an integer type, raise an error
    elif not numpy.issubdtype(type(mem), numpy.integer):
        raise TypeError(
            "Incorrect 'mem' data type.\n"\
            "    Options: None, numpy.integer, str\n"\
            "    type(mem) = %s" % (type(mem),)
        )

    # if mem is <= 0, raise an error
    if mem <= 0:
        raise ValueError(
            "'mem' must be greater than zero:\n"\
            "    Received: %s" % mem
        )

    return mem

def proc_hc_sa_set_varg(algorithm_varg, sel_size, geno):
    # if algorithm_varg is None
    if algorithm_varg is None:
        # make a default dictionary
        algorithm_varg = {
            "k": sel_size,
            "states": numpy.arange(geno.shape[1])
        }
    # make sure we've received a dictionary
    elif not isinstance(algorithm_varg, dict):
        raise TypeError(
            "Expected 'algorithm_varg' to be a dict\n"\
            "    type(algorithm_varg) = %s" % (type(algorithm_varg),)
        )

    # if algorithm_varg does not have a 'states' key, add a default
    if "states" not in algorithm_varg.keys():
        algorithm_varg["states"] = numpy.arange(geno.shape[1])

    # if algorithm_varg does not have a 'k' key, add a default
    if "k" not in algorithm_varg.keys():
        algorithm_varg["k"] = sel_size

    return algorithm_varg

def proc_hc_sa_state_varg(algorithm_varg, sel_size, geno):
    if algorithm_varg is None:
        # make a default dictionary
        algorithm_varg = {
            "states": numpy.tile(numpy.arange(geno.shape[1]), (sel_size,)),
            "dim_sizes": numpy.repeat(geno.shape[1], sel_size)
        }
    # make sure we've received a dictionary
    elif not isinstance(algorithm_varg, dict):
        raise TypeError(
            "Expected 'algorithm_varg' to be a dict\n"\
            "    type(algorithm_varg) = %s" % (type(algorithm_varg))
        )
    # if any ["states","dim_sizes"] are not in algorithm_varg.keys
    elif any(e not in algorithm_varg.keys() for e in ["states","dim_sizes"]):
        raise ValueError(
            "Expected 'algorithm_varg' dict to have fields:\n"\
            "    states\n    dim_sizes\n"\
            "'algorithm_varg' fields present:\n    %s" % algorithm_varg.keys()
        )

    return algorithm_varg

def proc_icpso_varg(algorithm_varg):
    # define the fields we're looking for in 'algorithm_varg'
    varg_fields = ["n", "inertia_wt", "accel_coeff_pbest",  "accel_coeff_gbest",
        "scale_factor", "stpfn"]

    # make sure we've received a dictionary
    if not isinstance(algorithm_varg, dict):
        raise TypeError(
            "Expected 'algorithm_varg' to be a dict\n"\
            "    type(algorithm_varg) = %s" % (type(algorithm_varg),)
        )
    # make sure we have all of our necessary arguments.
    elif any(e not in algorithm_varg.keys() for e in varg_fields):
        raise ValueError(
            "Expected 'algorithm_varg' to have all of the following fields:\n"\
            "    Field              Received\n"\
            "    n                  %s\n"\
            "    inertia_wt         %s\n"\
            "    accel_coeff_pbest  %s\n"\
            "    accel_coeff_gbest  %s\n"\
            "    scale_factor       %s\n"\
            "    stpfn              %s\n" %
            tuple(e not in algorithm_varg.keys() for e in varg_fields)
        )

    # if the user has not specified the states && dim_sizes, default to all
    if all(e not in algorithm_varg.keys() for e in ["states", "dim_sizes"]):
        algorithm_varg["states"] = numpy.tile(
            numpy.arange(geno.shape[1]),
            (sel_size,)
        )
        algorithm_varg["dim_sizes"] = numpy.repeat(geno.shape[1], sel_size)
    # elif 'states' and 'dim_sizes' have both been provided
    elif all(e in algorithm_varg.keys() for e in ["states", "dim_sizes"]):
        pass
    else:
        raise ValueError(
            "Inconsistent 'state' and 'dim_sizes' in 'algorithm_varg'\n"\
            "    Must provide both or none"
        )


    return algorithm_varg
