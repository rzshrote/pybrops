"""
Module for estimating breeding values in situations where there is no GxE and
there are replications nested within an environment.
"""

from rpy2 import robjects
from rpy2.robjects import globalenv
from rpy2.robjects import vectors
from rpy2.robjects.vectors import DataFrame as rpy2_DataFrame

from pybrops.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol
from pybrops.core.util.rpy2 import numpy_to_R_BoolVector
from pybrops.core.util.rpy2 import numpy_to_R_FloatVector
from pybrops.core.util.rpy2 import numpy_to_R_IntFactorVector
from pybrops.core.util.rpy2 import numpy_to_R_IntVector

class G_E_RnE_EstimatedBreedingValue(BreedingValueProtocol):
    """
    Class implementing protocols for estimating breeding values in situations
    where there is no GxE and there are replications nested within an
    environment. This class interfaces with R and lme4 to estimate effects.
    """

    ########################## Special Object Methods ##########################
    def __init__(self, G_efct = "fixed", E_efct = "random", RnE_efct = "random", **kwargs: dict):
        """
        Get estimated breeding values for an additive model of:
            y = G + E + R(E) + e
        Where:
            y is the response.
            G is the effect of genotype.
            E is the effect of environment.
            R(E) is the effect of replicate nested within environment.
            e is error.

        Parameters
        ----------
        G_efct : str
            A string indicating how genotype effects should be calculated.
            Option   | Description
            ---------+----------------------------------------------------------
            "fixed"  | Genotype effect is fixed.
            "random" | Genotype effect is random.
        E_efct : str
            A string indicating how environment effects should be calculated.
            Option   | Description
            ---------+----------------------------------------------------------
            "fixed"  | Environment effect is fixed.
            "random" | Environment effect is random.
        RnE_efct : str
            A string indicating how replication nested within environment
            effects should be calculated.
            Option   | Description
            ---------+----------------------------------------------------------
            "fixed"  | Replication(Environment) effect is fixed.
            "random" | Replication(Environment) effect is random.
        """
        super(G_E_RnE_EstimatedBreedingValue, self).__init__(**kwargs)
        # set data
        self.G_efct = G_efct
        self.E_efct = E_efct
        self.RnE_efct = RnE_efct

    ############################## Object Methods ##############################
    def estimate_ptdf(ptdf, gmat, **kwargs: dict):
        """
        Helper method to handle breeding value estimation when passed a
        PhenotypeDataFrame.
        """
        # get model parameters
        G_efct = self.G_efct
        E_efct = self.E_efct
        RnE_efct = self.RnE_efct

        ####################################################
        ################ Process Genotypes #################
        ####################################################

        # get taxa predictor column data
        taxa_ls = ptdf.col_data(name = "taxa")

        # make sure there aren't multiple taxa columns
        if len(taxa_ls) > 1:
            raise ValueError("Multiple 'taxa' columns detected in PhenotypeDataFrame")

        # convert taxa data to factor(int) vector
        taxa_arr = taxa_ls[0]               # get taxa vector
        taxa_uniq, taxa_int = numpy.unique( # get unique taxa
            taxa_arr,                       # input vector
            return_inverse = True           # return reverse indices
        )

        # make sure all genotypes are represented
        if not numpy.all(numpy.in1d(gmat.taxa, taxa_uniq)):
            raise ValueError("Not all genotypes in genotype matrix are in phenotype dataframe")

        # convert integer to factor(int)
        taxa_Rvec_factor = numpy_to_R_IntFactorVector(taxa_int)

        ####################################################
        ############### Process Environments ###############
        ####################################################

        # get environment predictor column data
        env_ls = ptdf.col_data(name = "env")

        # make sure there aren't multiple env columns
        if len(env_ls) > 1:
            raise ValueError("Multiple 'env' columns detected in PhenotypeDataFrame")

        # convert taxa data to factor(int) vector
        env_arr = env_ls[0]                 # get env vector
        env_uniq, env_int = numpy.unique(   # get unique env
            env_arr,                        # input vector
            return_inverse = True           # return reverse indices
        )

        # convert integer to factor(int)
        env_Rvec_factor = numpy_to_R_IntFactorVector(env_int)

        ####################################################
        ############### Process Replications ###############
        ####################################################

        # get replicate predictor column data
        rep_ls = ptdf.col_data(name = "rep")

        # make sure there aren't multiple rep columns
        if len(rep_ls) > 1:
            raise ValueError("Multiple 'rep' columns detected in PhenotypeDataFrame")

        # convert taxa data to factor(int) vector
        rep_arr = rep_ls[0]                 # get rep vector
        rep_uniq, rep_int = numpy.unique(   # get unique rep
            rep_arr,                        # input vector
            return_inverse = True           # return reverse indices
        )

        # convert integer to factor(int)
        rep_Rvec_factor = numpy_to_R_IntFactorVector(rep_int)

        ####################################################
        ################# Process Response #################
        ####################################################

        # get response columns data
        resp_ls, resp_name, resp_atype = ptdf.col_data(
            aefct = "response",
            return_name = True,
            return_atype = True
        )

        # get dictionary for mapping analysis types to conversion functions
        atype_to_R = {
            "bool" : numpy_to_R_BoolVector,
            "double" : numpy_to_R_FloatVector,
            "int" : numpy_to_R_IntVector
        }

        # make sure not all
        for atype,name in zip(resp_atype, resp_name):
            if atype not in atype_to_R.keys():
                raise ValueError("Analysis type '{0}' for '{1}' not supported. Options are: '{2}'".format(atype,name,atype_to_R.keys()))

        # extract response data and convert to R types
        # equivalent to:
        # for vec,atype,name in zip(resp_ls,resp_atype,resp_name):
        #     data_dict[name] = atype_to_R[atype](vec)
        data_dict = {n: atype_to_R[a](v) for v,a,n in zip(resp_ls,resp_atype,resp_name)}

        # add predictors
        data_dict["taxa"] = taxa_Rvec_factor
        data_dict["env"] = env_Rvec_factor
        data_dict["rep"] = rep_Rvec_factor

        # construct an R data.frame
        input_Rdf = rpy2_DataFrame(data_dict)

        # add data.frame to global environment in R
        globalenv['input_df'] = input_Rdf

        # load lme4 into R
        robjects.r('library(lme4)')

        # construct strings for defining the model
        if G_efct == "fixed":
            lme4_eqcomp_G = "taxa"
        elif G_efct == "random":
            lme4_eqcomp_G = "(1|taxa)"

        if E_efct == "fixed":
            lme4_eqcomp_E = "env"
        elif E_efct == "random":
            lme4_eqcomp_E = "(1|env)"

        if RnE_efct == "fixed":
            lme4_eqcomp_RnE = "(rep:env)"
        elif RnE_efct == "random":
            lme4_eqcomp_RnE = "(1|rep:env)"

        # run each model for each trait
        for name in resp_name:
            robjects.r("lmod = lmer({resp} ~ {G} + {E} + {RnE})".format(
                resp = name,
                G = lme4_eqcomp_G,
                E = lme4_eqcomp_E,
                RnE = lme4_eqcomp_RnE
            ))




    def estimate(self, pt_or_bv, gmat, **kwargs: dict):
        """
        Estimate breeding values.

        Parameters
        ----------
        pt_or_bv : PhenotypeDataFrame, BreedingValueMatrix
            Phenotype dataframe or breeding value matrix to use to estimate
            breeding values.
        gmat : GenotypeMatrix
            Genotype matrix used to align genotypes in estimation output.

        Returns
        -------
        bvmat : BreedingValueMatrix
            Breeding value matrix.
        """
        # check data types
        if is_PhenotypeDataFrame(pt_or_bv):
            self._estimate_ptdf(pt_or_bv, gmat, **kwargs)
        elif is_BreedingValueMatrix(pt_or_bv):
            return pt_or_bv # if we received a breeding value matrix, return it
        else:
            raise TypeError("'pt_or_bv' must be of type PhenotypeDataFrame or BreedingValueMatrix")

        raise NotImplementedError("method is abstract")
