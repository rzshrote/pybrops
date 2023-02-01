class InitializationOperator:
    """
    docstring for InitializationOperator.
    """
    def __init__(self, **kwargs: dict):
        """
        Constructor for InitializationOperator.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(InitializationOperator, self).__init__(**kwargs)
    
    def initialize(self, mu, sspace):
        """
        Description
        
        Parameters
        ----------
        mu, sspace : argtype
            argdesc
        
        Returns
        -------
        out : outtype
            outdesc
        """
        raise NotImplementedError("method is abstract")