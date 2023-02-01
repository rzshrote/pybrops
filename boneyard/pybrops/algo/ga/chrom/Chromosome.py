"""
Module for defining basal Chromosome interfaces and associated error checking routines.
Chromosomes are used to represent solutions in Genetic Algorithms.
"""

class Chromosome:
    """
    docstring for Chromosome.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict):
        """
        Constructor for Chromosome class.
        """
        super(Chromosome, self).__init__(**kwargs)

    ################# Container operators ##################
    def __len__(self):
        """Length of the Chromosome"""
        raise NotImplementedError("method is abstract")

    def __getitem__(self, key):
        """Get an element of the Chromosome"""
        raise NotImplementedError("method is abstract")

    def __setitem__(self, key, value):
        """Set an element of the Chromosome"""
        raise NotImplementedError("method is abstract")

    def __delitem__(self, key):
        """Delete an element of the Chromosome"""
        raise NotImplementedError("method is abstract")

    def __iter__(self):
        """Get an iterator for the Chromosome"""
        raise NotImplementedError("method is abstract")

    #################### Chromosome copying ####################
    def __copy__(self):
        """
        Make a shallow copy of the the Chromosome.

        Returns
        -------
        out : Chromosome
        """
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def feasible():
        doc = "The feasible property."
        def fget(self):
            """Get value for feasible."""
            raise NotImplementedError("property is abstract")
        def fset(self, value):
            """Set value for feasible."""
            raise NotImplementedError("property is abstract")
        def fdel(self):
            """Delete value for feasible."""
            raise NotImplementedError("property is abstract")
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    feasible = property(**feasible())

    def fitness():
        doc = "The fitness property."
        def fget(self):
            """Get value for fitness."""
            raise NotImplementedError("property is abstract")
        def fset(self, value):
            """Set value for fitness."""
            raise NotImplementedError("property is abstract")
        def fdel(self):
            """Delete value for fitness."""
            raise NotImplementedError("property is abstract")
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    fitness = property(**fitness())

    def fitness_n():
        doc = "The fitness_n property."
        def fget(self):
            """Get value for fitness_n."""
            raise NotImplementedError("property is abstract")
        def fset(self, value):
            """Set value for fitness_n."""
            raise NotImplementedError("property is abstract")
        def fdel(self):
            """Delete value for fitness_n."""
            raise NotImplementedError("property is abstract")
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    fitness_n = property(**fitness_n())

    def fitness_var():
        doc = "The fitness_var property."
        def fget(self):
            """Get value for fitness_var."""
            raise NotImplementedError("property is abstract")
        def fset(self, value):
            """Set value for fitness_var."""
            raise NotImplementedError("property is abstract")
        def fdel(self):
            """Delete value for fitness_var."""
            raise NotImplementedError("property is abstract")
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    fitness_var = property(**fitness_var())

    def chrom():
        doc = "The chrom property."
        def fget(self):
            """Get value for chrom."""
            raise NotImplementedError("property is abstract")
        def fset(self, value):
            """Set value for chrom."""
            raise NotImplementedError("property is abstract")
        def fdel(self):
            """Delete value for chrom."""
            raise NotImplementedError("property is abstract")
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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

    def decode(self):
        """
        Decode a chromosome from its native format to something which can be used by an objective/fitness function.
        
        Returns
        -------
        out : object
            A decoded chromosome.
        """
        raise NotImplementedError("method is abstract")

