from numbers import Integral
from numbers import Real
from typing import Callable
from typing import Optional
from typing import Union

import numpy
from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol
from pybrops.breed.prot.mate.MatingProtocol import check_is_MatingProtocol
from pybrops.breed.prot.sel.prob.RealSelectionProblem import RealSelectionProblem
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.AdditiveLinearGenomicModel import check_is_AdditiveLinearGenomicModel
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import check_is_PhasedGenotypeMatrix


class RealLookAheadGeneralizedWeightedGenomicSelectionProblem(RealSelectionProblem):
    """
    docstring for RealLookAheadGeneralizedWeightedGenomicSelectionProblem.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            fndr_pgmat: PhasedGenotypeMatrix,
            fndr_algmod: AdditiveLinearGenomicModel,
            mtprot: MatingProtocol,
            nparent: Integral,
            ncross: Integral,
            nprogeny: Integral,
            nsimul: Integral,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for RealLookAheadGeneralizedWeightedGenomicSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(RealLookAheadGeneralizedWeightedGenomicSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            **kwargs
        )
        # assignments
        self.fndr_pgmat = fndr_pgmat
        self.fndr_algmod = fndr_algmod
        self.mtprot = mtprot
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.nsimul = nsimul

    ############################ Object Properties #############################
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        return 2

    @property
    def fndr_pgmat(self) -> PhasedGenotypeMatrix:
        """Founder genotypes."""
        return self._fndr_pgmat
    @fndr_pgmat.setter
    def fndr_pgmat(self, value: PhasedGenotypeMatrix) -> None:
        """Set founder genotypes."""
        check_is_PhasedGenotypeMatrix(value, "fndr_pgmat")
        self._fndr_pgmat = value
    
    @property
    def fndr_algmod(self) -> AdditiveLinearGenomicModel:
        """Founder genomic prediction model."""
        return self._fndr_algmod
    @fndr_algmod.setter
    def fndr_algmod(self, value: AdditiveLinearGenomicModel) -> None:
        """Set founder genomic prediction model."""
        check_is_AdditiveLinearGenomicModel(value, "fndr_algmod")
        if value.ntrait > 1:
            raise ValueError("Selection method does not support >1 traits")
        self._fndr_algmod = value
    
    @property
    def mtprot(self) -> MatingProtocol:
        """mtprot."""
        return self._mtprot
    @mtprot.setter
    def mtprot(self, value: MatingProtocol) -> None:
        """Set mtprot."""
        check_is_MatingProtocol(value, "mtprot")
        self._mtprot = value
    
    @property
    def nparent(self) -> Integral:
        """nparent."""
        return self._nparent
    @nparent.setter
    def nparent(self, value: Integral) -> None:
        """Set nparent."""
        check_is_Integral(value, "nparent")
        self._nparent = value
    
    @property
    def ncross(self) -> Integral:
        """ncross."""
        return self._ncross
    @ncross.setter
    def ncross(self, value: Integral) -> None:
        """Set ncross."""
        check_is_Integral(value, "ncross")
        self._ncross = value
    
    @property
    def nprogeny(self) -> Integral:
        """nprogeny."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: Integral) -> None:
        """Set nprogeny."""
        check_is_Integral(value, "nprogeny")
        self._nprogeny = value
    
    @property
    def nsimul(self) -> Integral:
        """Number of simulations to evaluate a candidate solution."""
        return self._nsimul
    @nsimul.setter
    def nsimul(self, value: Integral) -> None:
        """Set number of simulations to evaluate a candidate solution."""
        self._nsimul = value

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
        Genomic Estimated Breeding Values (GEBV) for a population.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,)``.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A GEBV matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # gather some constants for quick access
        mtprot = self.mtprot
        nparent = self.nparent
        ncross = self.ncross
        nprogeny = self.nprogeny
        pgmat = self.fndr_pgmat
        algmod = self.fndr_algmod
        u_a = algmod.u_a
        ploidy = self.fndr_pgmat.ploidy
        # declare accumulators
        gain = 0.0
        usl = 0.0
        # run simulations
        for _ in range(self.nsimul):
            # reset pgmat to founders
            pgmat = self.fndr_pgmat.copy()
            # for each alpha value
            for alpha in x:
                # gentotype matrix (n,p)
                Z_a = pgmat.mat_asformat("{0,1,2}")
                # favorable allele frequency (p,t)
                fafreq = algmod.fafreq(pgmat)
                # prevent division by zero for fixed alleles
                fafreq[fafreq <= 0] = 1
                # calculate wGEBVs (n,t)
                wgebv = Z_a.dot(u_a * numpy.power(fafreq, -alpha))
                # take sum across trait (n,)
                wgebv = wgebv.sum(1)
                # find best parents
                sel = wgebv.argsort()[::-1][:nparent]
                # randomly shuffle indices
                numpy.random.shuffle(sel)
                # mate individuals and create progenies
                pgmat = mtprot.mate(pgmat, sel, ncross, nprogeny, None)
            # gentotype matrix (n,p)
            Z_a = pgmat.mat_asformat("{0,1,2}")
            # calculate mean GEBV and accumulate
            gain += (Z_a.dot(u_a)).sum(1).mean()
            # calculate the allele frequency (p,1)
            afreq = pgmat.afreq()[:,None]
            # get maximum attainable genotype (p,t)
            uslgeno = numpy.where(u_a > 0.0, afreq > 0.0, afreq >= 1.0)
            # calculate upper selection limit and accumulate
            usl += (float(ploidy) * u_a * uslgeno).sum()
        # take divide by number of simulations
        gain /= self.nsimul
        usl /= self.nsimul
        # create output as minimizing objectives
        out = numpy.array([-gain, -usl], dtype = float)
        return out
