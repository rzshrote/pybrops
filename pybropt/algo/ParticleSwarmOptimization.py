# 3rd party libraries

# our libraries
from . import Algorithm

class ParticleSwarmOptimization(Algorithm):
    """docstring for ParticleSwarmOptimization."""

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################
    def __init__(self, ssize, inertia_wt, pbest_comp, gbest_comp):
        """
        Constructor for the ParticleSwarmOptimization class.

        Parameters
        ----------
        ssize : int
            Size of swarm.
        inertia_wt : float
            Velocity inertial weight.
        pbest_comp : float
            Personal best inertial component.
        gbest_comp : float
            Global best inertial component.
        """
        # call super constructor; ...should call self.reset() as well...
        super(ParticleSwarmOptimization, self).__init__()

        # type checks
        check_is_integer(ssize, "ssize")
        check_is_numeric(inertia_wt, "inertia_wt")
        check_is_numeric(pbest_comp, "pbest_comp")
        check_is_numeric(gbest_comp, "gbest_comp")

        # set private variables
        self._ssize = ssize
        self._inertia_wt = inertia_wt
        self._pbest_comp = pbest_comp
        self._gbest_comp = gbest_comp

    ############################################################################
    ################################ Properties ################################
    ############################################################################
    def velocity():
        doc = "The velocity property."
        def fget(self):
            return self._velocity
        def fset(self, value):
            self._velocity = value
        def fdel(self):
            del self._velocity
        return locals()
    velocity = property(**velocity())

    def ssize():
        doc = "The ssize property."
        def fget(self):
            return self._ssize
        def fset(self, value):
            self._ssize = value
        def fdel(self):
            del self._ssize
        return locals()
    ssize = property(**ssize())

    def inertia_wt():
        doc = "The inertia_wt property."
        def fget(self):
            return self._inertia_wt
        def fset(self, value):
            self._inertia_wt = value
        def fdel(self):
            del self._inertia_wt
        return locals()
    inertia_wt = property(**inertia_wt())

    def pbest_comp():
        doc = "The pbest_comp property."
        def fget(self):
            return self._pbest_comp
        def fset(self, value):
            self._pbest_comp = value
        def fdel(self):
            del self._pbest_comp
        return locals()
    pbest_comp = property(**pbest_comp())

    def gbest_comp():
        doc = "The gbest_comp property."
        def fget(self):
            return self._gbest_comp
        def fset(self, value):
            self._gbest_comp = value
        def fdel(self):
            del self._gbest_comp
        return locals()
    gbest_comp = property(**gbest_comp())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    def reset(self):
        # call super.reset()
        super(ParticleSwarmOptimization, self).reset()

        # set to empty list
        self._velocity = []

    def history_add(self, algoiter, score, position, velocity):
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
            velocity = [velocity]
        elif algoiter_len is None:
            algoiter = numpy.repeat(algoiter, score_len)

        # force conversion to numpy.array
        algoiter = numpy.array(algoiter)
        score = numpy.array(score)
        position = numpy.array(position)
        velocity = numpy.array(velocity)

        ### finally append the matrices
        self._algoiter.append(algoiter)     # append algoiter
        self._score.append(score)           # append score
        self._position.append(position)     # append position
        self._velocity.append(velocity)     # append velocity
        return

    def history_to_dict(self):
        """
        Convert history internals to a dictionary.

        Returns
        -------
        hist_dict : dict
        """
        hist_dict = super(ParticleSwarmOptimization, self).history_to_dict()

        # make labels for velocity headers
        vhead = [
            "v"+str(i).zfill(zfill) for i in range(self._velocity[0].shape[1])
        ]

        # add columns + header name to df_dict
        for i,header in enumerate(vhead):
            hist_dict[header] = self._velocity[0][:,i]

        return hist_dict

    def concatenate(self):
        # call super.concatenate()
        super(ParticleSwarmOptimization, self).concatenate()
        if len(self._velocity) > 1:
            self._velocity = [numpy.concatenate(self._velocity, axis=0)]

    def pbest_ix(self, ptcl = None, maxiter = None, cond = None, minimum = True):
        """
        Get personal best indices.

        Assume that internal data is divisible by self.ssize.
        Assume that iterations are grouped together.

        maxiter : int, default = None
            Assume that all iterations are in order. Calculate the personal best
            values for each particle up to, but not including the 'iter'th
            iteration.
            If None, grab pbest using all iterations.
            Example:
                iter = 0 -> grab nothing
                iter = 1 -> grab pbest from iter 0
                iter = n -> grab pbest from iter 0, 1, ..., n-1
        """
        # concatenate everything
        self.concatenate()

        # build best function
        bestfn = numpy.argmin
        if not minimum:
            bestfn = numpy.argmax

        # build particle search list
        particles = numpy.arange(self._ssize) if ptcl is None else numpy.array(ptcl)

        # build mask
        mask = None
        if maxiter is not None:
            mask = self._algoiter[0] <= maxiter
        elif callable(cond):
            mask = cond(self._algoiter[0])

        # calculate indices
        ix = None
        if mask is None:
            ix = numpy.fromiter(
                i + (bestfn(self._score[0][i::self._ssize]) * self._ssize)
                for i in particles
            )
        else:
            nonzero = numpy.flatnonzero(mask)
            subarray = self._score[0][nonzero]

            nz_ix = numpy.fromiter(
                i + (bestfn(subarray[i::self._ssize]) * self._ssize)
                for i in particles
            )
            ix = nonzero[nz_ix]

        return ix

    def pbest(self, ptcl = None, maxiter = None, cond = None, minimum = True):
        """
        Grab personal best data
        """
        # grab indices
        ix = self.pbest_ix(ptcl, maxiter, cond, minimum)

        pbest = (
            self._algoiter[0][ix],
            self._score[0][ix],
            self._position[0][ix,:],
            self._velocity[0][ix,:]
        )

        return pbest

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
            self._position[0][ix,:],
            self._velocity[0][ix,:]
        )

        # return gbest
        return gbest
