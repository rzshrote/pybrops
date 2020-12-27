from . import GenomicModel

class NonparametricGenomicModel(GenomicModel):
    """docstring for NonparametricGenomicModel."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, trait, mkr_name = None, model_name = "Nonparametric"):
        super(NonparametricGenomicModel, self).__init__(trait, mkr_name, model_name)
