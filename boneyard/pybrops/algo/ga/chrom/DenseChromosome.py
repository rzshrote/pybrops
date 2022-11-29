import copy
import numpy
from typing import Union
import numbers

from pybrops.algo.opt.ga.chrom.Chromosome import Chromosome
from pybrops.core.error import check_is_ndarray
from pybrops.core.error import check_is_bool

class DenseChromosome(Chromosome):
    """
    docstring for DenseChromosome.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, 
        chrom: numpy.ndarray, 
        feasible: Union[bool,None] = None, 
        fitness: Union[numpy.ndarray,numbers.Number,None] = None, 
        fitness_n: Union[numpy.ndarray,numbers.Number,None] = None, 
        fitness_var: Union[numpy.ndarray,numbers.Number,None] = None, 
        **kwargs):
        """
        Constructor for DenseChromosome.
        
        Parameters
        ----------
        chrom : numpy.ndarray
            Chromosome representation.
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseChromosome, self).__init__(**kwargs)
        self.chrom = chrom
        self.feasible = feasible
        self.fitness = fitness
        self.fitness_n = fitness_n
        self.fitness_var = fitness_var

    ################# Container operators ##################
    def __len__(self):
        """Length of the Chromosome"""
        return len(self._chrom)

    def __getitem__(self, key):
        """Get an element of the Chromosome"""
        return self._chrom[key]

    def __setitem__(self, key, value):
        """Set an element of the Chromosome"""
        self._chrom[key] = value

    def __delitem__(self, key):
        """Delete an element of the Chromosome"""
        del self._chrom[key]

    def __iter__(self):
        """Get an iterator for the Chromosome"""
        return iter(self._chrom)

    #################### Chromosome copying ####################
    def __copy__(self):
        """
        Make a shallow copy of the the Chromosome.

        Returns
        -------
        out : Chromosome
        """
        return self.__class__(
            chrom = copy.copy(self._chrom),
            feasible = copy.copy(self._feasible),
            fitness = copy.copy(self._fitness),
            fitness_n = copy.copy(self._fitness_n),
            fitness_var = copy.copy(self._fitness_var)
        )

    def __deepcopy__(self, memo):
        """
        Make a deep copy of the Chromosome.

        Parameters
        ----------
        memo : dict

        Returns
        -------
        out : Chromosome
        """
        return self.__class__(
            chrom = copy.deepcopy(self._chrom, memo),
            feasible = copy.deepcopy(self._feasible, memo),
            fitness = copy.deepcopy(self._fitness, memo),
            fitness_n = copy.deepcopy(self._fitness_n, memo),
            fitness_var = copy.deepcopy(self._fitness_var, memo)
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def feasible():
        doc = "The feasible property."
        def fget(self):
            """Get value for feasible."""
            return self._feasible
        def fset(self, value):
            """Set value for feasible."""
            if value is not None:
                check_is_bool(value, "feasible")
            self._feasible = value
        def fdel(self):
            """Delete value for feasible."""
            del self._feasible
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    feasible = property(**feasible())

    def fitness():
        doc = "The fitness property."
        def fget(self):
            """Get value for fitness."""
            return self._fitness
        def fset(self, value):
            """Set value for fitness."""
            if (value is not None) and (not isinstance(value, (numpy.ndarray,numbers.Number))):
                raise TypeError("property 'fitness' must be a number or numpy.ndarray")
            self._fitness = value
        def fdel(self):
            """Delete value for fitness."""
            del self._fitness
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    fitness = property(**fitness())

    def fitness_n():
        doc = "The fitness_n property."
        def fget(self):
            """Get value for fitness_n."""
            return self._fitness_n
        def fset(self, value):
            """Set value for fitness_n."""
            if (value is not None) and (not isinstance(value, (numpy.ndarray,numbers.Number))):
                raise TypeError("property 'fitness_n' must be a number or numpy.ndarray")
            self._fitness_n = value
        def fdel(self):
            """Delete value for fitness_n."""
            del self._fitness_n
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    fitness_n = property(**fitness_n())

    def fitness_var():
        doc = "The fitness_var property."
        def fget(self):
            """Get value for fitness_var."""
            return self._fitness_var
        def fset(self, value):
            """Set value for fitness_var."""
            if (value is not None) and (not isinstance(value, (numpy.ndarray,numbers.Number))):
                raise TypeError("property 'fitness_var' must be a number or numpy.ndarray")
            self._fitness_var = value
        def fdel(self):
            """Delete value for fitness_var."""
            del self._fitness_var
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    fitness_var = property(**fitness_var())

    def chrom():
        doc = "The chrom property."
        def fget(self):
            """Get value for chrom."""
            return self._chrom
        def fset(self, value):
            """Set value for chrom."""
            check_is_ndarray(value, "repr")
            self._chrom = value
        def fdel(self):
            """Delete value for chrom."""
            del self._chrom
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    chrom = property(**chrom())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    #################### Chromosome copying ####################
    def copy(self):
        """
        Make a shallow copy of the Chromosome.

        Returns
        -------
        out : Chromosome
            A shallow copy of the original Chromosome.
        """
        return copy.copy(self)

    def deepcopy(self, memo):
        """
        Make a deep copy of the Chromosome.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : Chromosome
            A deep copy of the original Chromosome.
        """
        return copy.deepcopy(self, memo)

    def decode(self):
        """
        Decode a chromosome from its native format to something which can be used by an objective/fitness function.
        
        Returns
        -------
        out : numpy.ndarray
            A decoded chromosome.
        """
        return self._chrom

