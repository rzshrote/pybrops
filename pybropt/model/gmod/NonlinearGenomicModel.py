from . import GenomicModel

class NonlinearGenomicModel(GenomicModel):
    """docstring for NonlinearGenomicModel."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for NonlinearGenomicModel class.

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(NonlinearGenomicModel, self).__init__(**kwargs)
