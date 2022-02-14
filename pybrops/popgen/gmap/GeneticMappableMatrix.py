from pybrops.core.mat.VariantMatrix import VariantMatrix

class GeneticMappableMatrix(VariantMatrix):
    """
    Abstract class for variant matrices that can be interpolated using a
    GeneticMap.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(GeneticMappableMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################# Interpolation Methods ################
    def interp_genpos(self, gmap, **kwargs):
        """
        Interpolate genetic map postions for variants using a GeneticMap

        Parameters
        ----------
        gmap : GeneticMap
            A genetic map from which to interopolate genetic map postions for
            loci within the VariantMatrix.
        """
        raise NotImplementedError("method is abstract")

    def interp_xoprob(self, gmap, gmapfn, **kwargs):
        """
        Interpolate genetic map positions AND crossover probabilities between
        sequential markers using a GeneticMap and a GeneticMapFunction.

        Parameters
        ----------
        gmap : GeneticMap
            A genetic map from which to interopolate genetic map postions for
            loci within the VariantMatrix.
        gmapfn : GeneticMapFunction
            A genetic map function from which to interpolate crossover
            probabilities for loci within the VariantMatrix.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GeneticMappableMatrix(v):
    """
    Determine whether an object is a GeneticMappableMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a GeneticMappableMatrix object instance.
    """
    return isinstance(v, GeneticMappableMatrix)

def check_is_GeneticMappableMatrix(v, varname):
    """
    Check if object is of type GeneticMappableMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_GeneticMappableMatrix(v):
        raise TypeError("'{0}' must be a GeneticMappableMatrix".format(varname))

def cond_check_is_GeneticMappableMatrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type GeneticMappableMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        GeneticMappableMatrix.
    """
    if cond(v):
        check_is_GeneticMappableMatrix(v, varname)
