import model

class NonparametricGenomicModel(model.GenomicModel):
    """docstring for NonparametricGenomicModel."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, trait, model_name = "Nonparametric"):
        super(NonparametricGenomicModel, self).__init__(trait, model_name)
