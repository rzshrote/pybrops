import pybropt.util

class SearchSpace:
    """docstring for SearchSpace."""

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################
    def __init__(self, name = None):
        """
        name : str
            Name of the search space.
        """
        # check data types
        pybropt.util.cond_check_is_string(name, "name")

        # set private variables
        self._name = name

    ############################################################################
    ################################ Properties ################################
    ############################################################################
    def name():
        doc = "The name property."
        def fget(self):
            return self._name
        def fset(self, value):
            self._name = value
        def fdel(self):
            del self._name
        return locals()
    name = property(**name())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    def reset(self):
        raise NotImplementedError("The method 'reset' is abstract.")

    def feasible(self, pos):
        raise NotImplementedError("The method 'feasible' is abstract.")

    def feasible_vec(self, pos):
        raise NotImplementedError("The method 'feasible_vec' is abstract.")
