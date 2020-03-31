class Algorithm:
    """docstring for Algorithm."""
    ############################################################################
    ############################# Class Constants ##############################
    ############################################################################

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################
    def __init__(self, name = None):
        # initialize empty lists
        self.reset()

        # set name
        self._name = name

    def __len__(self):
        """
        Get the length of the score list.
        """
        l = len(self._score)
        return l

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def algoiter():
        doc = "The algoiter property."
        def fget(self):
            return self._algoiter
        def fset(self, value):
            self._algoiter = value
        def fdel(self):
            del self._algoiter
        return locals()
    algoiter = property(**algoiter())

    def score():
        doc = "The score property."
        def fget(self):
            return self._score
        def fset(self, value):
            self._score = value
        def fdel(self):
            del self._score
        return locals()
    score = property(**score())

    def position():
        doc = "The position property."
        def fget(self):
            return self._position
        def fset(self, value):
            self._position = value
        def fdel(self):
            del self._position
        return locals()
    position = property(**position())

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
        self._algoiter = []
        self._score = []
        self._position = []

    def history_add(self, algoiter, score, position):
        """
        Add algorithm history to this object.

        Parameters
        ==========
        algoiter : int, array-like
            Iteration number.
        score : int, array-like
            Score of the position.
        position : array-like
            1D or 2D matrix of positions.
        """
        # get lengths if they have them
        algoiter_len = cond_len(algoiter)
        score_len = cond_len(score)
        #position_len = cond_len(position)   # not needed

        # convert things to arrays
        if score_len is None:
            algoiter = [algoiter]
            score = [score]
            position = [position]
        elif algoiter_len is None:
            algoiter = numpy.repeat(algoiter, score_len)

        # force conversion to numpy.array
        algoiter = numpy.array(algoiter)
        score = numpy.array(score)
        position = numpy.array(position)

        ### finally append the matrices
        self._algoiter.append(algoiter)     # append iter
        self._score.append(score)           # append score
        self._position.append(position)     # append pos
        return

    def concatenate(self):
        if len(self._iter) > 1:
            self._iter = [numpy.concatenate(self._iter, axis=0)]
        if len(self._score) > 1:
            self._score = [numpy.concatenate(self._score, axis=0)]
        if len(self._pos) > 1:
            self._pos = [numpy.concatenate(self._pos, axis=0)]

    def history_to_dict(self):
        """
        Convert history internals to a dictionary.

        Returns
        -------
        hist_dict : dict
        """
        # concatenate internal lists
        self.concatenate()

        # make a dictionary to construct the pandas.DataFrame
        hist_dict = {
            "iter" : self._algoiter[0],
            "score" : self._score[0]
        }

        # make labels for the X position headers
        xhead = [
            "x"+str(i).zfill(zfill) for i in range(self._position[0].shape[1])
        ]

        # add columns + header name to df_dict
        for i,header in enumerate(xhead):
            hist_dict[header] = self._position[0][:,i]

        return hist_dict

    def history_to_df(self, zfill = 3):
        """
        Convert internal history to pandas.DataFrame

        Parameters
        ==========
        zfill : int
            Number of zeros to fill for X position labels

        Returns
        =======
        df : pandas.DataFrame
            A pandas DataFrame of the results.
        """
        # get history as dict
        df_dict = self.history_to_dict()

        # make DataFrame
        df = pandas.DataFrame(df_dict)

        return df

    def history_to_csv(self, fname, zfill = 3, *args, **kwargs):
        """
        fname : str
            File name to write to.
        zfill : int
            Number of zeros to fill in naming variables.
        *args
            Arguments to pass to pandas.DataFrame.to_csv.
        **kwargs
            Key word arguments to pass to pandas.DataFrame.to_csv.
        """
        # make data frame
        df = self.history_to_df(zfill)

        df.to_csv(fname, *args, **kwargs)

    def gbest_ix(self, maxiter = None, cond = None, minimum = True):
        """
        Get the index of the global best position. If multiple global best
        scores are equivalent, return the first instance.

        Parameters
        ----------
        maxiter : int
            Maximum algorithm iteration to search for a global best index.
        cond : callable
            Function to apply to search for global best index.
            Example:
                cond = (lambda x: x < 9)
        minimum : boolean
            Search for a global minimum if True. Search for a global maximum if
            False.

        Returns
        =======
        ix : int
            An index of the global best score.
        """
        # concatenate everything
        self.concatenate()

        # build best function
        bestfn = numpy.argmin
        if not minimum:
            bestfn = numpy.argmax

        # build mask
        mask = None
        if maxiter is not None:
            mask = self._algoiter[0] <= maxiter
        elif callable(cond):
            mask = cond(self._algoiter[0])

        # calculate index
        ix = None
        if mask is None:
            ix = bestfn(self._score[0])
        else:
            nonzero = numpy.flatnonzero(mask)
            nz_ix = bestfn(self._score[0][nonzero])
            ix = nonzero[nz_ix]

        return ix

    def gbest(self, maxiter = None, cond = None, minimum = True):
        """
        Get global best iter, score, position as a tuple. If multiple global
        best scores are equivalent, return the first instance.

        Returns
        =======
        gbest : tuple
            A tuple of global best iter, score, pos. The tuple order is:
                (iter, score, pos)
        """
        # get best score index
        ix = self.gbest_ix(maxiter, cond, minimum)

        # construct tuple
        gbest = (
            self._algoiter[0][ix],
            self._score[0][ix],
            self._position[0][ix,:]
        )

        # return gbest
        return gbest

    def optimize(self, objfn, seed = None, nthreads = 1, verbose = False,
            *args, **kwargs):
        raise NotImplementedError("This method not implemented.")

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
