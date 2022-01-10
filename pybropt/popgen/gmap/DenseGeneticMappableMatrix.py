from pybropt.core.mat import DenseVariantMatrix
from pybropt.popgen.gmap.GeneticMappableMatrix import GeneticMappableMatrix

from pybropt.popgen.gmap.GeneticMap import check_is_GeneticMap
from pybropt.popgen.gmap.GeneticMap import check_is_GeneticMapFunction

class DenseGeneticMappableMatrix(DenseVariantMatrix,GeneticMappableMatrix):
    """docstring for DenseGeneticMappableMatrix."""

    def __init__(self, **kwargs):
        super(DenseGeneticMappableMatrix, self).__init__()

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
        # check if gmap is a GeneticMap
        check_is_GeneticMap(gmap, "gmap")

        # interpolate postions
        self.vrnt_genpos = gmap.interp_genpos(self._vrnt_chrgrp, self._vrnt_phypos)

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
        # check data types
        check_is_GeneticMap(gmap, "gmap")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")

        # check if self has been sorted and grouped
        if not self.is_grouped_vrnt():
            raise RuntimeError("must be grouped first before interpolation of crossover probabilities")

        # interpolate genetic positions
        self.vrnt_genpos = gmap.interp_genpos(self._vrnt_chrgrp, self._vrnt_phypos)

        # interpolate crossover probabilities
        self.vrnt_xoprob = gmapfn.rprob1g(gmap, self._vrnt_chrgrp, self._vrnt_genpos)



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseGeneticMappableMatrix(v):
    """
    Determine whether an object is a DenseGeneticMappableMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseGeneticMappableMatrix object instance.
    """
    return isinstance(v, DenseGeneticMappableMatrix)

def check_is_DenseGeneticMappableMatrix(v, varname):
    """
    Check if object is of type DenseGeneticMappableMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_DenseGeneticMappableMatrix(v):
        raise TypeError("'{0}' must be a DenseGeneticMappableMatrix".format(varname))

def cond_check_is_DenseGeneticMappableMatrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type DenseGeneticMappableMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        DenseGeneticMappableMatrix.
    """
    if cond(v):
        check_is_DenseGeneticMappableMatrix(v, varname)
