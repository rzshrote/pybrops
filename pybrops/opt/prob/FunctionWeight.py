"""
Module defining several helper classes for assigning weights to optimization problems.
"""

__all__ = [
    "FunctionWeight",
    "MinimizingFunctionWeight",
    "MaximizingFunctionWeight",
]

from numbers import Real

class FunctionWeight:
    """
    Class for representing weights for functions that are in their natural state.
    Not intended for general use. Please use MinimizingFunctionWeight and 
    MaximizingFunctionWeight.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            wt: Real, 
            optimization_type: str = "min"
        ) -> None:
        """
        Constructor for Weight.
        
        Parameters
        ----------
        wt : Real
            A weight value. Must be non-negative.
        optimization_type : str
            The problem optimization type for which to calculate the weight sign.
            Must be either ``"min"`` or ``"max"`` to represent a minimizing or
            maximizing optimization, respectively.
        """
        self.wt = wt
        self.optimization_type = optimization_type

    def __float__(self):
        """Convert a weight to a floating point value."""
        return float(self.wt)

    def __int__(self):
        """Convert a weight to an integer value."""
        return int(self.wt)
    
    def __str__(self) -> str:
        """Convert a weight to a string."""
        return str(self.wt)
    
    def __repr__(self) -> str:
        """Convert a weight to a string representation."""
        return str(self.wt)

    ############################ Object Properties #############################
    @property
    def wt(self) -> Real:
        """A weight value."""
        return self._wt
    @wt.setter
    def wt(self, value: Real) -> None:
        """Set weight value."""
        if not isinstance(value, Real):
            raise TypeError("'wt' must be a Real value")
        if value < 0:
            raise ValueError("'wt' must be non-negative")
        self._wt = value
    
    @property
    def optimization_type(self) -> str:
        """The problem type for which to calculate the weight sign.."""
        return self._optimization_type
    @optimization_type.setter
    def optimization_type(self, value: str) -> None:
        """Set the problem type for which to calculate the weight sign.."""
        if not isinstance(value, str):
            raise TypeError("'optimization_type' must be of type str")
        value = value.lower()
        if value != "min" and value != "max":
            raise ValueError("'optimization_type' value must be either 'min' or 'max'")
        self._optimization_type = value


class MinimizingFunctionWeight(FunctionWeight):
    """
    Class for representing weights for functions that are minimizing in their natural state.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            wt: Real, 
            optimization_type: str = "min"
        ) -> None:
        """
        Constructor for MinimizingFunctionWeight.
        
        Parameters
        ----------
        wt : Real
            A weight value. Must be non-negative.
        optimization_type : str
            The problem optimization type for which to calculate the weight sign.
            Must be either ``"min"`` or ``"max"`` to represent a minimizing or
            maximizing optimization, respectively.
            
            If the function in its original state is naturally minimizing and
            the optimization type is minimizing, then ``wt`` is represented in
            its unaltered form by this class.

            If the function in its original state is naturally minimizing and
            the optimization type is maximizing, then ``wt`` is represented in
            its negated form by this class.
        """
        super(MinimizingFunctionWeight, self).__init__(
            wt = wt,
            optimization_type = optimization_type
        )

    def __float__(self):
        """Convert a weight to a floating point value."""
        out = float(self.wt)
        return out if self.optimization_type == "min" else -out

    def __int__(self):
        """Convert a weight to an integer value."""
        out = int(self.wt)
        return out if self.optimization_type == "min" else -out

    def __str__(self) -> str:
        """Convert a weight to a string."""
        return str(self.wt if self.optimization_type == "min" else -self.wt)
    
    def __repr__(self) -> str:
        """Convert a weight to a string representation."""
        return str(self.wt if self.optimization_type == "min" else -self.wt)



class MaximizingFunctionWeight(FunctionWeight):
    """
    Class for representing weights for functions that are maximizing in their natural state.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            wt: Real, 
            optimization_type: str = "min"
        ) -> None:
        """
        Constructor for MaximizingFunctionWeight.
        
        Parameters
        ----------
        wt : Real
            A weight value. Must be non-negative.
        optimization_type : str
            The problem optimization type for which to calculate the weight sign.
            Must be either ``"min"`` or ``"max"`` to represent a minimizing or
            maximizing optimization, respectively.
            
            If the function in its original state is naturally maximizing and
            the optimization type is minimizing, then ``wt`` is represented in
            its negated form by this class.

            If the function in its original state is naturally maximizing and
            the optimization type is maximizing, then ``wt`` is represented in
            its unaltered form by this class.
        """
        super(MaximizingFunctionWeight, self).__init__(
            wt = wt,
            optimization_type = optimization_type
        )

    def __float__(self) -> float:
        """Convert a weight to a floating point value."""
        out = float(self.wt)
        return -out if self.optimization_type == "min" else out

    def __int__(self) -> int:
        """Convert a weight to an integer value."""
        out = int(self.wt)
        return -out if self.optimization_type == "min" else out

    def __str__(self) -> str:
        """Convert a weight to a string."""
        return str(-self.wt if self.optimization_type == "min" else self.wt)
    
    def __repr__(self) -> str:
        """Convert a weight to a string representation."""
        return str(-self.wt if self.optimization_type == "min" else self.wt)
