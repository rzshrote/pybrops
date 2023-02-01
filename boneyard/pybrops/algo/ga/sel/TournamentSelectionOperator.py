class TournamentSelectionOperator(SelectionOperator):
    """docstring for TournamentSelectionOperator."""
    def __init__(self, n, **kwargs: dict):
        super(TournamentSelectionOperator, self).__init__(
            n = n,
            **kwargs
        )
        self.n = n
    # TODO: fixme
    def select(self, t_cur, t_max, k, ppop, pscore, opop, oscore, **kwargs: dict):
        sel = []                            # declare output list
        sel_score = []                      # selection score
        ix = [k for k in range(len(score))] # create list of indices
        for i in range(k):
            tix = random.sample(ix, self.n) # sample n without replacement
            tmp = pop[ tix[0] ]
            tmp_score = score[ tix[0] ]
            for j in tix[1:]:
                if score[j] > tmp_score:
                    tmp = pop[j]
                    tmp_score = score[j]
            sel.append(tmp)
            sel_score.append(tmp_score)
        return sel, numpy.array(sel_score)
