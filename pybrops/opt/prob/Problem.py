"""
Module defining optimization problems.
"""

class Problem:
    """
    docstring for Problem.
    """

    def __init__(
            self,
            **kwargs
        ) -> None:
        """
        Constructor for Problem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(Problem, self).__init__(**kwargs)
