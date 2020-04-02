# 3rd party libraries
import numpy

# our libraries
from . import ParticleSwarmOptimization
import pybropt.util

class ICPSO(ParticleSwarmOptimization):
    """docstring for ICPSO."""

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################
    def __init__(self, ssize, inertia_wt, pbest_comp, gbest_comp, scale_factor, search_space):
        # call super constructor; ...should call self.reset() for us
        super(ICPSO, self).__init__(ssize, inertia_wt, pbest_comp, gbest_comp)

        # type checks
        pybropt.util.check_is_numeric(scale_factor, "scale_factor")
        pybropt.util.check_is_CategoricalSearchSpace(search_space, "search_space")

        # set private variables
        self._scale_factor = scale_factor
        self._search_space = search_space

    ############################################################################
    ################################ Properties ################################
    ############################################################################
    def scale_factor():
        doc = "The scale_factor property."
        def fget(self):
            return self._scale_factor
        def fset(self, value):
            self._scale_factor = value
        def fdel(self):
            del self._scale_factor
        return locals()
    scale_factor = property(**scale_factor())

    def sample():
        doc = "The sample property."
        def fget(self):
            return self._sample
        def fset(self, value):
            self._sample = value
        def fdel(self):
            del self._sample
        return locals()
    sample = property(**sample())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    def reset(self):
        # call super.reset()
        super(ParticleSwarmOptimization, self).reset()

        # set to empty list
        self._sample = []
        return

    def history_add(self, algoiter, score, position, velocity, sample):
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
            sample = [sample]
        elif algoiter_len is None:
            algoiter = numpy.repeat(algoiter, score_len)

        # force conversion to numpy.array
        algoiter = numpy.array(algoiter)
        score = numpy.array(score)
        position = numpy.array(position)
        velocity = numpy.array(velocity)
        sample = numpy.array(sample)

        ### finally append the matrices
        self._algoiter.append(algoiter)     # append algoiter
        self._score.append(score)           # append score
        self._position.append(position)     # append position
        self._velocity.append(velocity)     # append velocity
        self._sample.append(sample)         # append sample
        return

    def history_to_dict(self):
        """
        Convert history internals to a dictionary.

        Returns
        -------
        hist_dict : dict
        """
        hist_dict = super(ICPSO, self).history_to_dict()

        # make labels for sample headers
        shead = [
            "s"+str(i).zfill(zfill) for i in range(self._sample[0].shape[1])
        ]

        # add columns + header name to df_dict
        for i,header in enumerate(shead):
            hist_dict[header] = self._sample[0][:,i]

        return hist_dict

    def concatenate(self):
        # call super.concatenate()
        super(ParticleSwarmOptimization, self).concatenate()
        if len(self._sample) > 1:
            self._sample = [numpy.concatenate(self._sample, axis=0)]

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
            self._velocity[0][ix,:],
            self._sample[0][ix,:]
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
            self._velocity[0][ix,:],
            self._sample[0][ix,:]
        )

        # return gbest
        return gbest

    def optimize(self, objfn, stpfn, seed = None, nthreads = None,
            verbose = False, *args, **kwargs):
        """
        Perform Integer/Categorical Particle Swarm Optimization.

        Parameters
        ----------
        objfn : function
            Objective function to optimize.
        stpfn : function
            A stop function. Should accept iteration number.
            Returns False to stop the iterations.
        seed : int, default = None.
            Random number seed.
        nthreads : int, default = None
            Number of threads. If nthreads is not None, this instantiates a
            multiprocessing pool with the given number of processes.
        verbose : boolean
            Show verbose messages.
        *args
            Additional arguments for objfn.
        **kwargs
            Additional keyword arguments for objfn.

        Returns
        -------
        None
        """
        ###################################
        #### Step 0) Process arguments ####
        ###################################

        # seed RNG if a seed is provided, otherwise use previous internal state.
        pybropt.util.cond_seed_rng(seed)

        # print engine state before beginning.
        if verbose:
            print("Integer/Categorical PSO engine state:")
            print()
            print("Objective function")
            print("==================================================")
            print("Objective function =", objfn.__name__)
            print("Objective function arguments:", str(args))
            print("Objective function keyword arguments:", str(kwargs))
            print()
            print("Swarm configuration")
            print("==================================================")
            print("Number of particles =", self._ssize)
            print("Number of dimensions =", self._search_space.ndim)
            print("Possible position states =", self._search_space)
            print("Velocity inertial weight =", self._inertia_wt)
            print("Personal best acceleration coefficient =", self._pbest_comp)
            print("Global best acceleration coefficient =", self._gbest_comp)
            print("Stopping function =", stpfn.__name__)
            print()
            print("Miscellaneous")
            print("==================================================")
            print("numpy.random seed =", seed)
            print("Number of threads =", nthreads)

        ################################
        #### Step 1) Initialization ####
        ################################

        ######################################################
        # Step 1a) allocate memory, initialize select arrays #
        ######################################################

        # get dimension sizes and states as numpy arrays
        dim_sizes = self._search_space.dim_size
        states = self._search_space.state_flat

        # create arrays of dimension state start and stop indices
        dsp = dim_sizes.cumsum()    # stop indices
        dst = dsp - dim_sizes       # start indices

        #################################################
        # generate position probability array
        #################################################

        # initialize random X probability positions with n rows, len(states) columns
        pos = numpy.random.uniform(     # use uniform distribution (dtype=float64)
            low = 0.0,                  # min probability (inclusive)
            high = 1.0,                 # max probability (exclusive)
            size = (n, len(states))     # n rows by len(states) columns
        )

        #################################################
        # generate velocity vectors
        #################################################
        # initialize random X velocity with matrix shape 'pshape'
        vel = numpy.random.uniform( # use uniform distribution (dtype=float64)
            low = 0.0,              # min trajectory strength (inclusive)
            high = 1.0,             # max trajectory strength (exclusive)
            size = pos.shape        # use same shape as pos
        )

        # for each group adjust values to sum to 1
        for st,sp in zip(dst, dsp):
            pos[:,st:sp] /= pos[:,st:sp].sum(1)[:,None]
            vel[:,st:sp] /= vel[:,st:sp].sum(1)[:,None] # not needed???

        #################################################
        # generate sample array
        #################################################

        # allocate memory for an empty sample array
        smpl = numpy.empty(                 # holds samples from pos vector
            shape = (n, len(dim_sizes)),    # n rows by len(dim_sizes) columns
            dtype = states.dtype            # use same dtype as states
        )

        # take samples for each row, for each grouping
        for i in range(n):
            # for each index, start, stop (for each group)
            for j,(st,sp) in enumerate(zip(dst, dsp)):
                smpl[i,j] = numpy.random.choice(    # set sample to random choice
                    states[st:sp],                  # with these states
                    size = 1,                       # only sample 1
                    p = pos[i,st:sp]                # with these probabilities
                )

        #################################################
        # declare score matrix variable
        #################################################

        score = None

        #################################################
        # initialize personal and global best positions and score variables
        #################################################

        # personal best variables
        pbest_pos = None    # personal best positions (probability matrix)
        pbest_score = None  # personal best scores (vector; dtype depends on objfn)

        # global best variables
        gbest_pos = None    # global best position (probability vector); a matrix view
        gbest_score = None  # global best score (dtype depends on objfn)

        #################################################
        # score initial positions
        #################################################

        # score samples that have been created (threading applicable)
        if nthreads > 1:
            # https://stackoverflow.com/questions/3033952/threading-pool-similar-to-the-multiprocessing-pool
            # https://stackoverflow.com/questions/3044580/multiprocessing-vs-threading-python
            raise RuntimeError("Multithreading not supported yet...")
        else:
            # evaluate the objective function
            score = numpy.apply_along_axis(objfn, 1, smpl, *args, **kwargs)


        #################################################
        # set personal best positions and scores
        #################################################

        # copy randomly initialized position matrix to an pbest_pos matrix
        pbest_pos = pos.copy()

        # copy initial scores
        pbest_score = score.copy()

        #################################################
        # identify global best
        #################################################

        # find the index of the highest scoring particle
        gbest_ix = pbest_score.argmax()

        # get matrix of global best position
        gbest_pos = pbest_pos[gbest_ix,:].copy()

        # set the gbest_score to pbest_score[argmax]
        gbest_score = pbest_score[gbest_ix]

        #################################################
        # define iteration counter global best
        #################################################

        # define a variable to track iteration number; initialization counts as 0.
        iter = 0

        #################################################
        # append history
        #################################################
        self.history_add(iter, score, pos, vel, smpl)

        #####################################
        #### Step 2) Main iteration loop ####
        #####################################

        # while stpfn is True, we haven't been told to stop
        while stpfn(iter):
            #####################################################
            # Step 2a) update particle velocities and positions #
            #####################################################

            # increment iteration counter
            iter += 1

            # generate random acceleration coefficients for personal and global
            acoef_pbest = numpy.random.uniform(
                low = 0.0,
                high = accel_coeff_pbest,
                size = n
            )
            acoef_gbest = numpy.random.uniform(
                low = 0.0,
                high = accel_coeff_gbest,
                size = n
            )

            # calculate the new velocities
            vel = (inertia_wt * vel) + \
                    (acoef_pbest[:,None] * (pbest_pos - pos)) + \
                    (acoef_gbest[:,None] * (gbest_pos - pos))

            # map to edge
            vel[vel > 1.0] = 1.0
            vel[vel < 0.0] = 0.0

            # scale velocities
            for st,sp in zip(dst, dsp):
                vel[:,st:sp] /= vel[:,st:sp].sum(1)[:,None] # do not scale???

            # add velocity vectors to X positions
            pos += vel

            # Map position probabilities to their edge
            # Logic: if pos > 1, set to 1; if pos < 0, set to 0
            pos[pos > 1.0] = 1.0
            pos[pos < 0.0] = 0.0

            # for each group adjust pos values to sum to 1
            for st,sp in zip(dst, dsp):
                pos[:,st:sp] /= pos[:,st:sp].sum(1)[:,None]

            ################################################################
            # Step 2b) generate samples from pos and score those samples #
            ################################################################
            # take samples for each row, for each grouping
            for i in range(n):
                # for each index, start, stop (for each group)
                for j,(st,sp) in enumerate(zip(dst, dsp)):
                    smpl[i,j] = numpy.random.choice(    # set sample to random choice
                        states[st:sp],                  # with these states
                        size = 1,                       # only sample 1
                        p = pos[i,st:sp]                # with these probabilities
                    )

            # score samples (threading applicable)
            if nthreads > 1:
                raise RuntimeError("Multithreading not supported yet...")
            else:
                # evaluate the objective function for each sample
                score = numpy.apply_along_axis(objfn, 1, smpl, *args, **kwargs)

            ################################################################
            # add history
            ################################################################
            self.history_add(iter, score, pos, vel, smpl)

            ################################################################
            # determine if scores improve
            ################################################################

            improve_mask = score > pbest_score  # create an improvement mask
            for i in range(n):                  # for each score
                if improve_mask[i]:             # if score beats previous best
                    pbest_score[i] = score[i]   # update best score

                    # set new pbest_pos that is biased for values
                    for j,(st,sp) in enumerate(zip(dst,dsp)):
                        # find where sample state is not equal to potential state
                        inv_mask = states[st:sp] != smpl[i,j]

                        # calculate positive bias for sample values
                        pos_bias = ((1-scale_factor) * pos[i,st:sp] * inv_mask).sum()

                        # bias and assign new personal best position
                        pbest_pos[i,st:sp] = numpy.where(
                            inv_mask,                     # if element not sample
                            scale_factor * pos[i,st:sp],  # add neg bias
                            pos[i,st:sp] + pos_bias       # else add positive bias
                        )

            ##########################################
            # Step 2c) update the global best values #
            ##########################################

            # find the index of the highest scoring
            gbest_ix = pbest_score.argmax()

            # set the X global best position to the view of the argmax row.
            gbest_pos = pbest_pos[gbest_ix,:].copy()

            # set the gbest_score to pbest_score[argmax]
            gbest_score = pbest_score[gbest_ix]

            # if verbose:
            if verbose:
                print("Iteration", iter, "\tgBest =", gbest_score)
