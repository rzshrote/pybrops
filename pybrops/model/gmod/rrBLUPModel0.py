"""
Module containing RR-BLUP genomic prediction classes for a very simple 
intercept + marker effects model.
"""

import copy
from typing import Optional, Tuple, Union
from numbers import Integral, Real
import numpy
import pandas
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral, check_is_Real
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len_eq, check_ndarray_ndim
from pybrops.core.error.error_value_python import check_is_gteq
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from scipy.optimize import minimize, Bounds
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

def rrBLUP_ML0_calc_G(Z: numpy.ndarray) -> numpy.ndarray:
    """
    Calculate a genomic relationship matrix from a marker matrix.

    Parameters
    ----------
    Z : numpy.ndarray
        A genotype matrix of shape ``(nobs,nmkr)``.
    
    Returns
    -------
    G : numpy.ndarray
        A genomic relationship matrix of shape ``(nobs,nobs)``.
    """
    # copy and convert data to floating point representation
    # (nobs,nmkr) -> (nobs,nmkr)
    Z = Z.astype(float)
    # center each column around the mean marker value
    # (nobs,nmkr) - (1,nmkr) -> (nobs,nmkr)
    Z -= numpy.mean(Z, 0, keepdims = True)
    # scale each column by the standard deviation of the marker value
    # (nobs,nmkr) * (1,nmkr) -> (nobs,nmkr)
    Z *= (1.0 / numpy.std(Z, 0, keepdims = True))
    # take the outer product to get the genomic relationship matrix
    # (nobs,nmkr) @ (nmkr,nobs) -> (nobs,nobs)
    G = Z @ Z.T
    # scale the matrix by the mean diagonal value
    # (nobs,nobs) * scalar -> (nobs,nobs)
    G *= (1.0 / numpy.mean(numpy.diag(G)))
    # return genomic relationship matrix
    return G

def rrBLUP_ML0_center_y(y: numpy.ndarray) -> numpy.ndarray:
    """
    Center y values around zero.

    Parameters
    ----------
    y : numpy.ndarray
        A vector of observations of shape ``(nobs,)``.
    
    Returns
    -------
    out : numpy.ndarray
        A vector of observations of shape ``(nobs,)`` centered around zero.
    """
    # center trait data around mean
    # (nobs,) - scalar -> (nobs,)
    return y - y.mean()

def rrBLUP_ML0_calc_d_V(G: numpy.ndarray) -> Tuple[numpy.ndarray,numpy.ndarray]:
    """
    Calculate the spectral decomposition of a (symmetric) genomic relationship matrix.

    Parameters
    ----------
    G : numpy.ndarray
        A genomic relationship matrix of shape ``(nobs,nobs)``.
    
    Returns
    -------
    out : tuple
        A tuple containing ``(d,V)``.
        
        Where::
        - ``d`` is a vector of shape ``(nobs,)`` representing the diagonal of the spectral decomposition.
        - ``V`` is a matrix of shape ``(nobs,nobs)`` representing the orthonormal basis of the spectral decomposition.

        Entries are sorted from highest eigenvalue to lowest eigenvalue.
    """
    # take the spectral decomposition of the (symmetric) G matrix.
    # (nobs,nobs) -> eigenvalues = (nobs,), eigenvectors = (nobs,nobs)
    d, V = numpy.linalg.eigh(G)
    # get indices for highest to lowest eigenvalues
    # (nobs,)
    ix = numpy.argsort(d)[::-1]
    # reorder eigenvalues from highest to lowest
    # (nobs,)[(nobs,)] -> (nobs,)
    d = d[ix]
    # reorder eigenvectors from highest to lowest
    # (nobs,)[(nobs,)] -> (nobs,)
    V = V[:,ix]
    # return spectral decomposition
    return d, V

def rrBLUP_ML0_nonzero_d_V(d: numpy.ndarray, V: numpy.ndarray, tol: Real = 1e-5) -> Tuple[numpy.ndarray,numpy.ndarray]:
    """
    Extract nonzero components of eigenvalues and eigenvectors from a spectral decomposition.

    Parameters
    ----------
    d : numpy.ndarray
        A vector of shape ``(nobs,)`` representing the diagonal of the spectral decomposition.
    V : numpy.ndarray
        A matrix of shape ``(nobs,nobs)`` representing the orthonormal basis of the spectral decomposition.
    
    Returns
    -------
    out : tuple
        A tuple containing ``(d,V)``.
        
        Where::
        - ``d`` is a vector of shape ``(ncomp,)`` representing the diagonal of the spectral decomposition.
        - ``V`` is a matrix of shape ``(nobs,ncomp)`` representing the orthonormal basis of the spectral decomposition.
    """
    # for the spectral decomposition, construct a mask for eigenvalues greater than approximately zero
    # (nobs,) -> (nobs,) 
    mask = d > tol
    # extract eigenvalues that are greater than approximately zero
    # (nobs,) -> (ncomp,)
    d = d[mask]
    # extract eigenvectors for eigenvalues that are greater than approximately zero
    # (nobs,nobs) -> (nobs,ncomp)
    V = V[:,mask]
    # return values
    return d, V

def rrBLUP_ML0_calc_etasq(V: numpy.ndarray, y: numpy.ndarray) -> numpy.ndarray:
    """
    Calculate eta squared values for fast computation of likelihoods.

    Parameters
    ----------
    V : numpy.ndarray
        A matrix of shape ``(nobs,ncomp)`` containing the non-zero eigenvectors of the spectral decomposition.
    y : numpy.ndarray
        A vector of shape ``(nobs,)`` containing zero centered observations.
    
    Returns
    -------
    etasq : numpy.ndarray
        A vector of shape ``(ncomp,)`` containing squared eta values.
    """
    # calculate V.T @ y to get eta values
    # (ncomp,nobs) @ (nobs,) -> (ncomp,)
    eta = V.T @ y
    # square eta values
    # (ncomp,)^2 -> (ncomp,)
    etasq = eta**2
    # return values
    return etasq

def rrBLUP_ML0_neg2LogLik_fast(logVarComp: numpy.ndarray, etasq: numpy.ndarray, d: numpy.ndarray, n: Integral) -> Real:
    """
    -2 log-likelihood function for ML estimation.
    In optimization, this function is to be minimized.
    
    Parameters
    ----------
    logVarComp : numpy.ndarray
        A vector of shape (2,) containing parameter estimates in the log scale (log(varE),log(varU)).
        The log scale is used for more efficient search since the search space is (0,Inf).
    etasq : numpy.ndarray
        A vector of shape ``(ncomp,)`` containing squared eta values.
    d : numpy.ndarray
        A vector of shape (ncomp,) containing non-zero eigenvalues for the genomic relationship matrix.
    n : Integral
        Number of observations (nobs).
    
    Returns
    -------
    out : Real
        Scalar value proportional to the -2 log-likelihood for ML estimation.
    """
    # variance components
    varE = numpy.exp(logVarComp[0])
    varU = numpy.exp(logVarComp[1])
    # calculate the ratio between genetic variance and error variance
    lamb = varU / varE
    # calculate diagonal values for varU * G + varE * I matrix spectral decomposition
    dStar = (d * lamb + 1)
    # calculate log-determinant using sum of logs of diagonal matrix
    sumLogD = numpy.sum(numpy.log(dStar))
    # calculate -2 * log-likelihood
    n2ll = (n * numpy.log(varE) + sumLogD) + ((numpy.sum(etasq / dStar)) / varE)
    # return -2 * log-likelihood
    return n2ll

def rrBLUP_ML0_calc_ridge(varE: Real, varU: Real) -> Real:
    """
    Calculate the ridge parameter.

    Parameters
    ----------
    varE : Real
        Error variance.
    varU : Real
        Marker variance.
    
    Returns
    -------
    out : Real
        The ridge parameter.
    """
    return varE / varU

def rrBLUP_ML0_calc_ZtZplI(Z: numpy.ndarray, ridge: Real) -> numpy.ndarray:
    """
    Calculate (Z'Z + lambda * I).

    Parameters
    ----------
    Z : numpy.ndarray
        A genotype matrix of shape ``(nobs,nmkr)``.
    ridge : Real
        The ridge parameter, lambda. Must be non-negative.

    Returns
    -------
    out : numpy.ndarray
        The calculated matrix of shape ``(nmkr,nmkr)``.
    """
    # calculate marker covariance matrix: Z'Z
    # (nmkr,nobs) @ (nobs,nmkr) -> (nmkr,nmkr)
    ZtZplI = Z.T @ Z
    # extract a view of the diagonal for which to add the ridge parameter
    diagZtZplI = numpy.einsum('ii->i', ZtZplI)
    # add ridge parameter to diagonal of Z'Z
    # (nmkr,) + scalar -> (nmkr,)
    diagZtZplI += ridge
    # return values
    return ZtZplI

def rrBLUP_ML0_calc_Zty(Z: numpy.ndarray, y: numpy.ndarray) -> numpy.ndarray:
    """
    Calculate Z'y.
    
    Parameters
    ----------
    Z : numpy.ndarray
        A genotype matrix of shape ``(nobs,nmkr)``.
    y : numpy.ndarray
        A vector of shape ``(nobs,)`` containing zero centered observations.

    Returns
    -------
    out : numpy.ndarray
        A vector of shape ``(nmrk,)``.
    """
    # calculate Z'y
    # (nmkr,nobs) @ (nobs,) -> (nmkr,)
    Zty = Z.T @ y
    # return values
    return Zty

def gauss_seidel(A: numpy.ndarray, b: numpy.ndarray, atol: Real = 1e-8, maxiter: Integral = 1000) -> numpy.ndarray:
    """
    Solve the equation Ax = b using the Gauss-Seidel method.

    Parameters
    ----------
    A : numpy.ndarray
        A diagonal dominant or symmetric positive definite matrix of shape (nmkr,nmkr).
    b : numpy.ndarray
        A vector of shape (nmkr,).
    atol : Real
        Absolute tolerance. Iterate until the sum of absolute differences 
        between successive iterations is less than this value or ``maxiter``
        is reached.
        Must be non-negative.
    maxiter : Integral
        Maximum number of iterations.
    
    Returns
    -------
    x : numpy.ndarray
        Solution to the system of equations.
    """
    # get number of markers
    nmkr = len(b)
    # allocate memory for the previous x estimate
    xprev = numpy.zeros(nmkr, dtype = float)
    # allocate memory for the current x estimate
    xcurr = numpy.zeros(nmkr, dtype = float)
    # number of iterations
    niter = 0
    # absolute difference
    adiff = 2*atol
    # main loop
    while numpy.any(adiff > atol) and niter < maxiter:
        # copy current x values to previous x values without memory allocation
        xprev[:] = xcurr
        # modify current x values using the Gauss-Seidel procedure
        for i in range(nmkr):
            xcurr[i] = (b[i] - A[i,:i].dot(xcurr[:i]) - A[i,i+1:].dot(xcurr[i+1:])) / A[i,i]
        # calculate the absolute difference between iterations
        adiff = numpy.abs(xcurr - xprev)
        # increment the iteration number
        niter += 1
    # return estimates
    return xcurr

def rrBLUP_ML0(y: numpy.ndarray, Z: numpy.ndarray, varlb: Real = 1e-5, varub: Real = 1e5, gsatol: Real = 1e-8, gsmaxiter: Integral = 1000):
    """
    Ridge regression BLUP for the simple model::

    y = Zu + e

    Where::
        - ``y`` are observations.
        - ``Z`` is a design matrix for genetic markers.
        - ``u`` are marker effects which follow the distribution ``MVN(0, varU * I)``.
        - ``e`` are errors which follow the distribution ``MVN(0, varE * I)``.
    
    Uses the EMMA formulation to solve for ``varE`` and ``varU``.
    Uses the Nelder-Mead method to optimize for variance components.
    Marker effects are estimated using the Gauss-Seidel method.

    Parameters
    ----------
    y : numpy.ndarray
        A vector of observations of shape ``(nobs,)``. If not mean centered, will be centered around zero.
    Z : numpy.ndarray
        A genotype matrix of shape ``(nobs,nmkr)``.
    varlb : Real
        Lower bound permitted for variance component estimation.
        Must be non-negative.
    varub : Real
        Upper bound permitted for variance component estimation.
        Must be non-negative and greater than ``varlb``.
    gsatol : Real
        Absolute tolerance for the Gauss-Seidel method.
        Iterate until the sum of absolute differences between successive 
        iterations is less than this value or ``maxiter`` is reached.
        Must be non-negative.
    gsmaxiter : Integral
        Maximum number of iterations for the Gauss-Seidel method.
        Must be non-negative.

    Returns
    -------
    out : dict
        A dictionary of output values.
    """
    # check input types
    check_is_ndarray(y, "y")
    check_is_ndarray(Z, "Z")
    check_is_Real(varlb, "varlb")
    check_is_Real(varub, "varub")
    check_is_Real(gsatol, "gsatol")
    check_is_Integral(gsmaxiter, "gsmaxiter")

    # check input values
    check_ndarray_ndim(y, "y", 1)
    check_ndarray_ndim(Z, "Z", 2)
    check_ndarray_axis_len_eq(Z, "Z", 0, len(y))
    check_is_gteq(varlb, "varlb", 0.0)
    check_is_gteq(varub, "varub", varlb)
    check_is_gteq(gsatol, "gsatol", 0.0)
    check_is_gteq(gsmaxiter, "gsmaxiter", 0)

    # get the mean of y (the intercept)
    # (nobs,) -> scalar
    meanY = y.mean()

    # center trait data around mean
    # (nobs,)
    y = rrBLUP_ML0_center_y(y)

    # get the number of observations
    # scalar
    nobs = len(y)

    # create genomic relationship matrix
    # (nobs,nobs)
    G = rrBLUP_ML0_calc_G(Z)

    # take the spectral decomposition of the (symmetric) G matrix.
    # (nobs,nobs) -> eigenvalues = (nobs,), eigenvectors = (nobs,nobs)
    d, V = rrBLUP_ML0_calc_d_V(G)

    # remove zero components
    # eigenvalues = (nobs,), eigenvectors = (nobs,nobs) -> eigenvalues = (ncomp,), eigenvectors = (nobs,ncomp)
    d, V = rrBLUP_ML0_nonzero_d_V(d, V)

    # calculate eta squared values
    # (ncomp,)
    etasq = rrBLUP_ML0_calc_etasq(V, y)

    # calculate variance of y
    # (nobs,) -> scalar
    varY = y.var()

    # calculate initial estimates of log(varE) and log(varU); set each to half of varY
    # scalar
    logVarE0 = numpy.log(0.5 * varY)
    logVarU0 = numpy.log(0.5 * varY)

    # construct inital starting position
    # (2,)
    logVarComp0 = numpy.array([logVarE0, logVarU0])

    # construct search space boundaries
    bounds = Bounds(
        lb = numpy.repeat(numpy.log(varlb), len(logVarComp0)),
        ub = numpy.repeat(numpy.log(varub), len(logVarComp0)),
    )

    # optimize for the variance components using Nelder-Mead algorithm
    soln = minimize(
        fun = rrBLUP_ML0_neg2LogLik_fast,
        x0 = logVarComp0,
        args = (etasq, d, nobs),
        method = 'Nelder-Mead',
        bounds = bounds,
    )

    # get the solution values
    varE = numpy.exp(soln.x[0])
    varU = numpy.exp(soln.x[1])

    # calculate the ridge parameter
    ridge = rrBLUP_ML0_calc_ridge(varE, varU)

    # calculate (Z'Z + lambda * I)
    ZtZplI = rrBLUP_ML0_calc_ZtZplI(Z, ridge)

    # calculate Z'y
    Zty = rrBLUP_ML0_calc_Zty(Z, y)

    # solve for (Z'Z + lambda * I)u = Z'y using the Gauss-Seidel method
    uhat = gauss_seidel(ZtZplI, Zty, gsatol, gsmaxiter)

    # calculate heritability
    h2 = varU / (varU + varE)

    # reconstruct the log-likelihood (minus constant)
    logLik = -2 * soln.fun

    # get intercept incidence matrix
    # (nobs,1)
    X = numpy.full((nobs,1), 1.0, float)

    # get intercept
    # (1,)
    betahat = numpy.array([meanY])

    # calculate y hat
    yhat = X.dot(betahat) + Z.dot(uhat)

    # create output dictionary
    out = {
        # ridge regression elements
        "yhat": yhat,
        "X": X,
        "betahat": betahat,
        "Z": Z,
        "uhat": uhat,
        # variance estimates
        "varE": varE,
        "varU": varU,
        "h2": h2,
        "logLik": logLik,
        # optimization solution object
        "soln": soln,
    }

    return out

class rrBLUPModel0(DenseAdditiveLinearGenomicModel):
    """
    The rrBLUPModel0 class represents a simple RR-BLUP model with an intercept 
    (fixed) and marker effects (random).

    An rrBLUPModel0 is a Multivariate Multiple Linear Regression model defined as:

    .. math::
        \\mathbf{Y} = \\mathbf{XB} + \\mathbf{ZU} + \\mathbf{E}

    Where:

    - :math:`\\mathbf{Y}` is a matrix of response variables of shape ``(n,t)``.
    - :math:`\\mathbf{X}` is a matrix of ones for the intercept of shape ``(n,1)``.
    - :math:`\\mathbf{B}` is a matrix of intercept coefficients of shape ``(1,t)``.
    - :math:`\\mathbf{Z}` is a matrix of random effect predictors of shape ``(n,p)``.
    - :math:`\\mathbf{U}` is a matrix of random effect regression coefficients of shape ``(p,t)``.
    - :math:`\\mathbf{E}` is a matrix of error terms of shape ``(n,t)``.

    Block matrix modifications to :

    :math:`\\mathbf{Z}` and :math:`\\mathbf{U}` can be decomposed into block
    matrices pertaining to different sets of effects:

    .. math::
        \\mathbf{Z} = \\begin{bmatrix} \\mathbf{Z_{misc}} & \\mathbf{Z_{a}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{Z_{misc}}` is a matrix of miscellaneous random effect predictors of shape ``(n,p_misc)``
    - :math:`\\mathbf{Z_{a}}` is a matrix of additive genomic marker predictors of shape ``(n,p_a)``

    .. math::
        \\mathbf{U} = \\begin{bmatrix} \\mathbf{U_{misc}} \\\\ \\mathbf{U_{a}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{U_{misc}}` is a matrix of miscellaneous random effects of shape ``(p_misc,t)``
    - :math:`\\mathbf{U_{a}}` is a matrix of additive genomic marker effects of shape ``(p_a,t)``

    Shape definitions:

    - ``n`` is the number of individuals
    - ``p`` is the number of random effect predictors.
    - ``p_misc`` is the number of miscellaneous random effect predictors.
    - ``p_a`` is the number of additive genomic marker predictors.
    - The sum of ``p_misc`` and ``p_a`` equals ``p``.
    - ``t`` is the number of traits

    From Prototype class docstring:

    RR-BLUP model for fitting a single random effect and a single intercept fixed effect for a single trait.
    If multiple traits are provided, fit independent models for each trait.

    For a single trait, the model is::

    y = Xb + Zu + e

    Where::
        - ``y`` are observations.
        - ``X`` is a matrix of ones for the incidence of the slope.
        - ``b`` is the intercept.
        - ``Z`` is a design matrix for genetic markers.
        - ``u`` are marker effects which follow the distribution ``MVN(0, varU * I)``.
        - ``e`` are errors which follow the distribution ``MVN(0, varE * I)``.

    For a single trait, if the observations (y) are mean centered, then the model becomes::

    y = Zu + e

    Where::
        - ``y`` are observations.
        - ``Z`` is a design matrix for genetic markers.
        - ``u`` are marker effects which follow the distribution ``MVN(0, varU * I)``.
        - ``e`` are errors which follow the distribution ``MVN(0, varE * I)``.

    Multiple traits are concatenated together into the model::

    Y = XB + ZU + E
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            beta: numpy.ndarray, 
            u_misc: Union[numpy.ndarray,None], 
            u_a: Union[numpy.ndarray,None], 
            trait: Optional[numpy.ndarray] = None, 
            method: str = "ML",
            model_name: Optional[str] = None, 
            hyperparams: Optional[dict] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for rrBLUPModel0 class.

        Parameters
        ----------
        beta : numpy.ndarray
            A ``float64`` fixed effect regression coefficient matrix of shape
            ``(q,t)``.

            Where:

            - ``q`` is the number of fixed effect predictors (e.g. environments).
            - ``t`` is the number of traits.
        
        u_misc : numpy.ndarray, None
            A ``float64`` random effect regression coefficient matrix of shape
            ``(p_misc,t)`` containing miscellaneous effects.

            Where:

            - ``p_misc`` is the number of miscellaneous random effect predictors.
            - ``t`` is the number of traits.

            If ``None``, then set to an empty array of shape ``(0,t)``.
        
        u_a : numpy.ndarray, None
            A ``float64`` random effect regression coefficient matrix of shape
            ``(p_a,t)`` containing additive marker effects.

            Where:

            - ``p_a`` is the number of additive marker effect predictors.
            - ``t`` is the number of traits.

            If ``None``, then set to an empty array of shape ``(0,t)``.
        
        trait : numpy.ndarray, None
            An ``object`` array of shape ``(t,)`` containing the names of traits.

            Where:

            - ``t`` is the number of traits.
        
        method : str
            Fitting method to use. Options are ``{"ML"}``.
            
        model_name : str, None
            Name of the model.
        
        hyperparams : dict, None
            Model parameters.
        
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        # call super constructor
        super(rrBLUPModel0, self).__init__(
            beta, 
            u_misc, 
            u_a, 
            trait, 
            model_name, 
            hyperparams, 
            **kwargs
        )
        # save method metadata
        self.method = method

    #################### Model copying #####################
    def __copy__(
            self
        ) -> 'rrBLUPModel0':
        """
        Make a shallow copy of the GenomicModel.

        Returns
        -------
        out : rrBLUPModel0
            A shallow copy of the model.
        """
        out = self.__class__(
            beta        = copy.copy(self.beta),
            u_misc      = copy.copy(self.u_misc),
            u_a         = copy.copy(self.u_a),
            trait       = copy.copy(self.trait),
            method      = copy.copy(self.method),
            model_name  = copy.copy(self.model_name),
            hyperparams = copy.copy(self.hyperparams),
        )

        return out

    def __deepcopy__(
            self, 
            memo: Optional[dict] = None
        ) -> 'rrBLUPModel0':
        """
        Make a deep copy of the GenomicModel.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : rrBLUPModel0
            A deep copy of the model.
        """
        out = self.__class__(
            beta        = copy.deepcopy(self.beta, memo),
            u_misc      = copy.deepcopy(self.u_misc, memo),
            u_a         = copy.deepcopy(self.u_a, memo),
            trait       = copy.deepcopy(self.trait, memo),
            method      = copy.deepcopy(self.method, memo),
            model_name  = copy.deepcopy(self.model_name, memo),
            hyperparams = copy.deepcopy(self.hyperparams, memo),
        )

        return out

    ############################ Object Properties #############################
    @property
    def method(self) -> str:
        """Method to be used to fit the RR-BLUP model."""
        return self._method
    @method.setter
    def method(self, value: str) -> None:
        """Set the method to fit the RR-BLUP model."""
        if not isinstance(value, str):
            raise TypeError("'method' must be of type ``str``, but received type ``{0}``".format(type(value).__name__))
        options = ("ML",)
        value = value.upper()
        if value not in options:
            raise ValueError("str 'method' must be one of ``method == {0}``, but received ``method == {1}``".format(options,value))
        self._method = value

    ############################## Object Methods ##############################

    #################### Model copying #####################
    def copy(
            self
        ) -> 'rrBLUPModel0':
        """
        Make a shallow copy of the GenomicModel.

        Returns
        -------
        out : GenomicModel
            A shallow copy of the original GenomicModel
        """
        return self.__copy__()

    def deepcopy(
            self,
            memo: Optional[dict] = None
        ) -> 'rrBLUPModel0':
        """
        Make a deep copy of the GenomicModel.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : GenomicModel
            A deep copy of the original GenomicModel
        """
        return self.__deepcopy__(memo)

    ####### methods for model fitting and prediction #######
    @classmethod
    def fit_numpy(
            cls, 
            Y: numpy.ndarray, 
            X: numpy.ndarray, 
            Z: numpy.ndarray, 
            trait: Optional[numpy.ndarray] = None, 
            method: str = "ML",
            model_name: Optional[str] = None, 
            hyperparams: Optional[dict] = None, 
            **kwargs: dict
        ) -> 'rrBLUPModel0':
        """
        Fit a dense, additive linear genomic model.

        Parameters
        ----------
        Y : numpy.ndarray
            A phenotype matrix of shape (n,t).
        X : numpy.ndarray
            Not used by this model. Assumed to be ``(n,q)`` matrix of ones.
        Z : numpy.ndarray
            A genotypes matrix of shape (n,p).
        trait : numpy.ndarray
            A trait name array of shape (t,).
        method : str
            Fitting method to use. Options are ``{"ML"}``.
        model_name : str, None
            Name of the model.
        hyperparams : dict, None
            Model parameters.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : rrBLUPModel0
            An RR-BLUP model.
        """
        # type checks
        check_is_ndarray(Y, "Y")
        # check_is_ndarray(X, "X") # ignored by this model
        check_is_ndarray(Z, "Z")

        # shape checks
        check_ndarray_ndim(Y, "Y", 2)
        # check_ndarray_ndim(X, "X", 2) # ignored by this model
        check_ndarray_ndim(Z, "Z", 2)

        # convert data types to floating if needed
        if not numpy.issubdtype(Y.dtype, numpy.floating):
            Y = Y.astype(float)
        if not numpy.issubdtype(Z.dtype, numpy.floating):
            Z = Z.astype(float)

        # get number of traits
        # scalar = t
        ntrait = Y.shape[1]

        # for each trait, fit a model
        models = [rrBLUP_ML0(Y[:,i], Z) for i in range(ntrait)]

        # aggregate intercepts and stack estimates
        # (1,t)
        beta = numpy.stack([models[i]["betahat"] for i in range(ntrait)], axis = 1)

        # aggregate marker effects and stack estimates
        # (p,t)
        u_a = numpy.stack([models[i]["uhat"] for i in range(ntrait)], axis = 1)

        # create output
        out = cls(
            beta = beta,
            u_misc = None,
            u_a = u_a,
            trait = trait,
            method = method,
            model_name = model_name,
            hyperparams = hyperparams,
            **kwargs
        )

        return out

    @classmethod
    def fit(
            cls, 
            ptobj: Union[BreedingValueMatrix,numpy.ndarray], 
            cvobj: numpy.ndarray, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            trait: Optional[numpy.ndarray] = None, 
            method: str = "ML",
            model_name: Optional[str] = None, 
            hyperparams: Optional[dict] = None, 
            **kwargs: dict
        ) -> 'rrBLUPModel0':
        """
        Fit a dense, additive linear genomic model.

        Parameters
        ----------
        ptobj : BreedingValueMatrix, numpy.ndarray
            An object containing phenotype data. Must be a matrix of breeding
            values or a phenotype data frame.
        cvobj : numpy.ndarray
            An object containing covariate data.
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        trait : numpy.ndarray, None
            A trait name array of shape (t,).
        method : str
            Fitting method to use. Options are ``{"ML"}``.
            
        model_name : str, None
            Name of the model.
        
        hyperparams : dict, None
            Model parameters.
        
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : rrBLUPModel0
            An RR-BLUP model.
        """
        # process ptobj
        if isinstance(ptobj, BreedingValueMatrix):
            Y = ptobj.unscale()
        elif isinstance(ptobj, pandas.DataFrame):
            raise RuntimeError("not implmented yet")
        elif isinstance(ptobj, numpy.ndarray):
            Y = ptobj
        else:
            raise TypeError("must be BreedingValueMatrix, pandas.DataFrame, numpy.ndarray")

        # process cvobj (ignored)
        X = None

        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            Z = gtobj.mat_asformat("{0,1,2}")
        elif isinstance(gtobj, numpy.ndarray):
            Z = gtobj
        else:
            raise TypeError("must be GenotypeMatrix, numpy.ndarray")

        # fit the model
        out = cls.fit_numpy(Y, X, Z, trait, method, model_name, hyperparams, **kwargs)

        return out



################################## Utilities ###################################
def check_is_rrBLUPModel0(v: object, vname: str) -> None:
    """
    Check if object is of type rrBLUPModel0. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, rrBLUPModel0):
        raise TypeError("variable '{0}' must be a rrBLUPModel0".format(vname))
