# import 3rd party libraries
import numpy

# import our libraries
from . import Population
import pybropt.util

class Cross:
    """
    The Cross class.

    Supports variance calculations for 2-, 3-, 4-, dihybrid crosses.
    Simulates mating for 2-, 3-, 4-way crosses.
    """
    ############################################################################
    ############################# Class Constants ##############################
    ############################################################################
    # NOTE: these dictionary can stay here because they do not reference
    #       internal methods; they only have strings.

    # dictionary of varAfn keys
    # nesting is [varAfn][sparse]
    KEY_TO_VARAFN = {
        '2way' : {
            True :  ('varA_2way_sparse', 'varA_2way_sparse_vec'),
            False : ('varA_2way', 'varA_2way_vec')
        },
        '2wayDH' : {
            True :  ('varA_2wayDH_sparse', 'varA_2wayDH_sparse_vec'),
            False : ('varA_2wayDH', 'varA_2wayDH_vec')
        },
        '3way' : {
            True :  ('varA_3way_sparse', 'varA_3way_sparse_vec'),
            False : ('varA_3way', 'varA_3way_vec')
        },
        '3wayDH' : {
            True :  ('varA_3wayDH_sparse', 'varA_3wayDH_sparse_vec'),
            False : ('varA_3wayDH', 'varA_3wayDH_vec')
        },
        '4way' : {
            True :  ('varA_4way_sparse', 'varA_4way_sparse_vec'),
            False : ('varA_4way', 'varA_4way_vec')
        },
        '4wayDH' : {
            True :  ('varA_4wayDH_sparse', 'varA_4wayDH_sparse_vec'),
            False : ('varA_4wayDH', 'varA_4wayDH_vec')
        },
        'dihybrid' : {
            True :  ('varA_dihybrid_sparse', 'varA_dihybrid_sparse_vec'),
            False : ('varA_dihybrid', 'varA_dihybrid_vec')
        },
        'dihybridDH' : {
            True :  ('varA_dihybridDH_sparse', 'varA_dihybridDH_sparse_vec'),
            False : ('varA_dihybridDH', 'varA_dihybridDH_vec')
        },
        None : {
            None : ('varA_default', 'varA_default')
        }
    }

    # dictionary of matefn keys
    # nesting is [crossfn][matefn]
    KEY_TO_MATEFN = {
        '2way' : {
            'ctrl' :    'mate_2way_ctrl',
            'wrand' :   'mate_2way_wrand',
            'erand' :   'mate_2way_erand'
        },
        '2wayDH' : {
            'ctrl' :    'mate_2wayDH_ctrl',
            'wrand' :   'mate_2wayDH_wrand',
            'erand' :   'mate_2wayDH_erand'
        },
        '3way' : {
            'ctrl' :    'mate_3way_ctrl',
            'wrand' :   'mate_3way_wrand',
            'erand' :   'mate_3way_erand'
        },
        '3wayDH' : {
            'ctrl' :    'mate_3wayDH_ctrl',
            'wrand' :   'mate_3wayDH_wrand',
            'erand' :   'mate_3wayDH_erand'
        },
        '4way' : {
            'ctrl' :    'mate_4way_ctrl',
            'wrand' :   'mate_4way_wrand',
            'erand' :   'mate_4way_erand'
        },
        '4wayDH' : {
            'ctrl' :    'mate_4wayDH_ctrl',
            'wrand' :   'mate_4wayDH_wrand',
            'erand' :   'mate_4wayDH_erand'
        },
        'dihybrid' : {
            'ctrl' :    'mate_2way_ctrl',
            'wrand' :   'mate_2way_wrand',
            'erand' :   'mate_2way_erand'
        },
        'dihybridDH' : {
            'ctrl' :    'mate_2wayDH_ctrl',
            'wrand' :   'mate_2wayDH_wrand',
            'erand' :   'mate_2wayDH_erand'
        },
        None : {
            None : 'mate_default'
        }
    }

    # dictionary of rallocfn keys
    # nesting is [crossfn][matefn][rallocfn]
    KEY_TO_RALLOCFN = {
        '2way' : {
            'ctrl' : {
                'equal' : 'ralloc_2way_ctrl_equal'
            },
            'wrand' : {
                'equal' : 'ralloc_2way_wrand_equal'
            },
            'erand' : {
                'equal' : 'ralloc_2way_erand_equal'
            }
        },
        '2wayDH' : {
            'ctrl' : {
                'equal' : 'ralloc_2wayDH_ctrl_equal'
            },
            'wrand' : {
                'equal' : 'ralloc_2wayDH_wrand_equal'
            },
            'erand' : {
                'equal' : 'ralloc_2wayDH_erand_equal'
            }
        },
        '3way' : {
            'ctrl' : {
                'equal' : 'ralloc_3way_ctrl_equal'
            },
            'wrand' : {
                'equal' : 'ralloc_3way_wrand_equal'
            },
            'erand' : {
                'equal' : 'ralloc_3way_erand_equal'
            }
        },
        '3wayDH' : {
            'ctrl' : {
                'equal' : 'ralloc_3wayDH_ctrl_equal'
            },
            'wrand' : {
                'equal' : 'ralloc_3wayDH_wrand_equal'
            },
            'erand' : {
                'equal' : 'ralloc_3wayDH_erand_equal'
            }
        },
        '4way' : {
            'ctrl' : {
                'equal' : 'ralloc_4way_ctrl_equal'
            },
            'wrand' : {
                'equal' : 'ralloc_4way_wrand_equal'
            },
            'erand' : {
                'equal' : 'ralloc_4way_erand_equal'
            }
        },
        '4wayDH' : {
            'ctrl' : {
                'equal' : 'ralloc_4wayDH_ctrl_equal'
            },
            'wrand' : {
                'equal' : 'ralloc_4wayDH_wrand_equal'
            },
            'erand' : {
                'equal' : 'ralloc_4wayDH_erand_equal'
            }
        },
        'dihybrid' : {
            'ctrl' : {
                'equal' : 'ralloc_2way_ctrl_equal'
            },
            'wrand' : {
                'equal' : 'ralloc_2way_wrand_equal'
            },
            'erand' : {
                'equal' : 'ralloc_2way_erand_equal'
            }
        },
        'dihybridDH' : {
            'ctrl' : {
                'equal' : 'ralloc_2wayDH_ctrl_equal'
            },
            'wrand' : {
                'equal' : 'ralloc_2wayDH_wrand_equal'
            },
            'erand' : {
                'equal' : 'ralloc_2wayDH_erand_equal'
            }
        },
        None : {
            None : {
                None : 'ralloc_default'
            }
        }
    }

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    # TODO: mem accepts strings in MB, GB, etc.
    def __init__(self, population, varAfn, sparse, crossfn, matefn, rallocfn, c = 1, n = 1, s = 0, t = 0, mem = None):
        """
        population : Population
            Population from which to cross.
        varAfn : str
            Additive variance function to use in calculations of varA matrix.
            Options:
                '2way', '2wayDH',
                '3way', '3wayDH',
                '4way', '4wayDH',
                'dihybrid', 'dihybridDH'
        sparse : bool
            Boolean to indicate whether the variance matrix is to be calculated
            as a sparse matrix (for memory efficiency)
        crossfn : str
            Cross structure to use.
            Options:
                '2way', '2wayDH',
                '3way', '3wayDH',
                '4way', '4wayDH',
                'dihybrid', 'dihybridDH'
        matefn : str
            Mating format to use.
            Options:
                'ctrl'      Controlled cross
                'wrand'     Weighted random cross
                'erand'     Exact random cross
        rallocfn : str
            Resource allocation function to use.
            Options:
                'equal'     Give equal contribution to each parent.
        c : int
            Cross multiplier. Its functionality depends on the matefn chosen.
        n : int
            Number of progeny to generate per cross.
        s : int
            Number of selfing generations.
        t : int
            Number of random intermatings.

        """
        # check data types
        pybropt.util.check_is_Population(population, "population")
        pybropt.util.check_is_string(varAfn, "varAfn")
        pybropt.util.check_is_bool(sparse, "sparse")
        pybropt.util.check_is_string(crossfn, "crossfn")
        pybropt.util.check_is_string(matefn, "matefn")
        pybropt.util.check_is_string(rallocfn, "rallocfn")
        pybropt.util.check_is_integer(c, "c")
        pybropt.util.check_is_integer(n, "n")
        pybropt.util.check_is_integer_or_inf(s, "s")
        pybropt.util.check_is_integer_or_inf(t, "t")

        # set private variables
        self._population = population
        self.set_varAfn(varAfn, sparse)
        self.set_matefn(crossfn, matefn)
        self.set_rallocfn(crossfn, matefn, rallocfn)
        self._c = c
        self._n = n
        self._s = s
        self._t = t
        self._mem = len(self._population.marker_set) if mem is None else mem

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def population():
        doc = "The population property."
        def fget(self):
            return self._population
        def fset(self, value):
            self._population = value
        def fdel(self):
            del self._population
        return locals()
    population = property(**population())

    def varAfn():
        doc = "The varAfn property."
        def fget(self):
            return self._varAfn
        def fset(self, value):
            self._varAfn = value
        def fdel(self):
            del self._varAfn
        return locals()
    varAfn = property(**varAfn())

    def sparse():
        doc = "The sparse property."
        def fget(self):
            return self._sparse
        def fset(self, value):
            self._sparse = value
        def fdel(self):
            del self._sparse
        return locals()
    sparse = property(**sparse())

    def crossfn():
        doc = "The crossfn property."
        def fget(self):
            return self._crossfn
        def fset(self, value):
            self._crossfn = value
        def fdel(self):
            del self._crossfn
        return locals()
    crossfn = property(**crossfn())

    def matefn():
        doc = "The matefn property."
        def fget(self):
            return self._matefn
        def fset(self, value):
            self._matefn = value
        def fdel(self):
            del self._matefn
        return locals()
    matefn = property(**matefn())

    def rallocfn():
        doc = "The rallocfn property."
        def fget(self):
            return self._rallocfn
        def fset(self, value):
            self._rallocfn = value
        def fdel(self):
            del self._rallocfn
        return locals()
    rallocfn = property(**rallocfn())

    def c():
        doc = "The c property."
        def fget(self):
            return self._c
        def fset(self, value):
            self._c = value
        def fdel(self):
            del self._c
        return locals()
    c = property(**c())

    def n():
        doc = "The n property."
        def fget(self):
            return self._n
        def fset(self, value):
            self._n = value
        def fdel(self):
            del self._n
        return locals()
    n = property(**n())

    def s():
        doc = "The s property."
        def fget(self):
            return self._s
        def fset(self, value):
            self._s = value
        def fdel(self):
            del self._s
        return locals()
    s = property(**s())

    def t():
        doc = "The t property."
        def fget(self):
            return self._t
        def fset(self, value):
            self._t = value
        def fdel(self):
            del self._t
        return locals()
    t = property(**t())

    def mem():
        doc = "The mem property."
        def fget(self):
            return self._mem
        def fset(self, value):
            self._mem = value
        def fdel(self):
            del self._mem
        return locals()
    mem = property(**mem())

    def varA_2way_mat():
        doc = "The varA_2way_mat property."
        def fget(self):
            return self._varA_2way_mat
        def fset(self, value):
            self._varA_2way_mat = value
        def fdel(self):
            del self._varA_2way_mat
        return locals()
    varA_2way_mat = property(**varA_2way_mat())

    def varA_2wayDH_mat():
        doc = "The varA_2wayDH_mat property."
        def fget(self):
            if not hasattr(self, '_varA_2wayDH_mat'):
                self._varA_2wayDH_mat = self.calc_varA_2wayDH_mat()
            return self._varA_2wayDH_mat
        def fset(self, value):
            self._varA_2wayDH_mat = value
        def fdel(self):
            if hasattr(self, '_varA_2wayDH_mat'):
                del self._varA_2wayDH_mat
        return locals()
    varA_2wayDH_mat = property(**varA_2wayDH_mat())

    def varA_3way_mat():
        doc = "The varA_3way_mat property."
        def fget(self):
            return self._varA_3way_mat
        def fset(self, value):
            self._varA_3way_mat = value
        def fdel(self):
            del self._varA_3way_mat
        return locals()
    varA_3way_mat = property(**varA_3way_mat())

    def varA_3wayDH_mat():
        doc = "The varA_3wayDH_mat property."
        def fget(self):
            if not hasattr(self, '_varA_3wayDH_mat'):
                self._varA_3wayDH_mat = self.calc_varA_3wayDH_mat()
            return self._varA_3wayDH_mat
        def fset(self, value):
            self._varA_3wayDH_mat = value
        def fdel(self):
            if hasattr(self, '_varA_3wayDH_mat'):
                del self._varA_3wayDH_mat
        return locals()
    varA_3wayDH_mat = property(**varA_3wayDH_mat())

    def varA_4way_mat():
        doc = "The varA_4way_mat property."
        def fget(self):
            return self._varA_4way_mat
        def fset(self, value):
            self._varA_4way_mat = value
        def fdel(self):
            del self._varA_4way_mat
        return locals()
    varA_4way_mat = property(**varA_4way_mat())

    def varA_4wayDH_mat():
        doc = "The varA_4wayDH_mat property."
        def fget(self):
            if not hasattr(self, '_varA_4wayDH_mat'):
                self._varA_4wayDH_mat = self.calc_varA_4wayDH_mat()
            return self._varA_4wayDH_mat
        def fset(self, value):
            self._varA_4wayDH_mat = value
        def fdel(self):
            if hasattr(self, '_varA_4wayDH_mat'):
                del self._varA_4wayDH_mat
        return locals()
    varA_4wayDH_mat = property(**varA_4wayDH_mat())

    def varA_dihybrid_mat():
        doc = "The varA_dihybrid_mat property."
        def fget(self):
            return self._varA_dihybrid_mat
        def fset(self, value):
            self._varA_dihybrid_mat = value
        def fdel(self):
            del self._varA_dihybrid_mat
        return locals()
    varA_dihybrid_mat = property(**varA_dihybrid_mat())

    def varA_dihybridDH_mat():
        doc = "The varA_dihybridDH_mat property."
        def fget(self):
            if not hasattr(self, "_varA_dihybridDH_mat"):
                self._varA_dihybridDH_mat = self.calc_varA_dihybridDH_mat()
            return self._varA_dihybridDH_mat
        def fset(self, value):
            self._varA_dihybridDH_mat = value
        def fdel(self):
            del self._varA_dihybridDH_mat
        return locals()
    varA_dihybridDH_mat = property(**varA_dihybridDH_mat())

    def varA_2way_sparse_mat():
        doc = "The varA_2way_sparse_mat property."
        def fget(self):
            return self._varA_2way_sparse_mat
        def fset(self, value):
            self._varA_2way_sparse_mat = value
        def fdel(self):
            del self._varA_2way_sparse_mat
        return locals()
    varA_2way_sparse_mat = property(**varA_2way_sparse_mat())

    def varA_2wayDH_sparse_mat():
        doc = "The varA_2wayDH_sparse_mat property."
        def fget(self):
            # TODO: sparse matrices
            # if not hasattr(self, '_varA_2wayDH_sparse_mat'):
            #     self._varA_2wayDH_sparse_mat = None
            return self._varA_2wayDH_sparse_mat
        def fset(self, value):
            self._varA_2wayDH_sparse_mat = value
        def fdel(self):
            del self._varA_2wayDH_sparse_mat
        return locals()
    varA_2wayDH_sparse_mat = property(**varA_2wayDH_sparse_mat())

    def varA_3way_sparse_mat():
        doc = "The varA_3way_sparse_mat property."
        def fget(self):
            return self._varA_3way_sparse_mat
        def fset(self, value):
            self._varA_3way_sparse_mat = value
        def fdel(self):
            del self._varA_3way_sparse_mat
        return locals()
    varA_3way_sparse_mat = property(**varA_3way_sparse_mat())

    def varA_3wayDH_sparse_mat():
        doc = "The varA_3wayDH_sparse_mat property."
        def fget(self):
            return self._varA_3wayDH_sparse_mat
        def fset(self, value):
            self._varA_3wayDH_sparse_mat = value
        def fdel(self):
            del self._varA_3wayDH_sparse_mat
        return locals()
    varA_3wayDH_sparse_mat = property(**varA_3wayDH_sparse_mat())

    def varA_4way_sparse_mat():
        doc = "The varA_4way_sparse_mat property."
        def fget(self):
            return self._varA_4way_sparse_mat
        def fset(self, value):
            self._varA_4way_sparse_mat = value
        def fdel(self):
            del self._varA_4way_sparse_mat
        return locals()
    varA_4way_sparse_mat = property(**varA_4way_sparse_mat())

    def varA_4wayDH_sparse_mat():
        doc = "The varA_4wayDH_sparse_mat property."
        def fget(self):
            return self._varA_4wayDH_sparse_mat
        def fset(self, value):
            self._varA_4wayDH_sparse_mat = value
        def fdel(self):
            del self._varA_4wayDH_sparse_mat
        return locals()
    varA_4wayDH_sparse_mat = property(**varA_4wayDH_sparse_mat())

    def varA_dihybrid_sparse_mat():
        doc = "The varA_dihybrid_sparse_mat property."
        def fget(self):
            return self._varA_dihybrid_sparse_mat
        def fset(self, value):
            self._varA_dihybrid_sparse_mat = value
        def fdel(self):
            del self._varA_dihybrid_sparse_mat
        return locals()
    varA_dihybrid_sparse_mat = property(**varA_dihybrid_sparse_mat())

    def varA_dihybridDH_sparse_mat():
        doc = "The varA_dihybridDH_sparse_mat property."
        def fget(self):
            return self._varA_dihybridDH_sparse_mat
        def fset(self, value):
            self._varA_dihybridDH_sparse_mat = value
        def fdel(self):
            del self._varA_dihybridDH_sparse_mat
        return locals()
    varA_dihybridDH_sparse_mat = property(**varA_dihybridDH_sparse_mat())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def reset(self):
        """
        Remove all calculated matrices.
        """
        # remove matrices if they exist.
        if hasattr(self, "_varA_2wayDH_mat"):
            del self._varA_2wayDH_mat
        if hasattr(self, "_varA_3wayDH_mat"):
            del self._varA_3wayDH_mat
        if hasattr(self, "_varA_4wayDH_mat"):
            del self._varA_4wayDH_mat
        if hasattr(self, "_varA_dihybridDH_mat"):
            del self._varA_dihybridDH_mat
        if hasattr(self, "_varA_2wayDH_sparse_mat"):
            del self._varA_2wayDH_mat
        if hasattr(self, "_varA_3wayDH_sparse_mat"):
            del self._varA_3wayDH_mat
        if hasattr(self, "_varA_4wayDH_sparse_mat"):
            del self._varA_4wayDH_mat
        if hasattr(self, "_varA_dihybridDH_sparse_mat"):
            del self._varA_dihybridDH_sparse_mat

    ############################################################################
    ##################### Variance Related Object Methods ######################
    def calc_varA_2way_mat(self):
        raise NotImplementedError

    def calc_varA_2wayDH_mat(self):
        """
        Calculate a symmetrical matrix of progeny variance for each pairwise
        2-way cross. Calculations are derived from Osthushenrich et al. (2017).

        Assumes inbred individuals.
        """
        # create aliases for several variables
        geno = self._population.geno
        coeff = self._population.genomic_model.coeff
        gmap = self._population.marker_set
        mem = self._mem

        # get number of traits (t)
        ntrait = coeff.shape[1]

        # get the number of individuals (n)
        ntaxa = geno.shape[1]

        # allocate a square matrix for each pairwise variance
        varA_mat = numpy.zeros((ntrait, ntaxa, ntaxa), dtype='float64')

        # for each linkage group
        for lst, lsp in zip(gmap.chr_grp_stix, gmap.chr_grp_spix):
            # for each computational chunk
            for rst,rsp in zip(range(lst,lsp,mem),pybropt.util.srange(lst+mem,lsp,mem)):
                for cst,csp in zip(range(lst,lsp,mem),pybropt.util.srange(lst+mem,lsp,mem)):
                    # get recombination probability matrix for chunk
                    r = gmap.recomb_prob(rst, rsp, cst, csp) # (rb,cb)

                    # calculate a D1 matrix; this is specific to mating scheme
                    D1_mat = Cross.D1(r, self._s, self._t) # (rb,cb)

                    # get marker coefficients for rows and columns
                    rcoeff = coeff[rst:rsp].T # (t,rb) # rb = row block
                    ccoeff = coeff[cst:csp].T # (t,cb) # cb = column block

                    # for each mate pair (excluding selfs)
                    for female in range(1,ntaxa): # varA row index
                        for male in range(0,female): # varA col index
                            # calculate genotype differences for row, col {-1,0,1}
                            rdgeno = geno[0,female,rst:rsp] - geno[0,male,rst:rsp] # (rb,)
                            cdgeno = geno[0,female,cst:csp] - geno[0,male,cst:csp] # (cb,)

                            # calculate effect differences
                            reffect = rdgeno * rcoeff # (rb,)*(t,rb) -> (t,rb)
                            ceffect = cdgeno * ccoeff # (cb,)*(t,cb) -> (t,cb)

                            # compute dot product for each trait to get partial variance
                            # (t,rb)x(rb,cb) -> (t,cb)
                            # (t,cb)*(t,cb) -> (t,cb)
                            # (t,cb)[1] -> (t,)
                            varA_part = (reffect @ D1_mat * ceffect).sum(1)

                            # add this partial variance to the lower triangle
                            varA_mat[:,female,male] += varA_part

        # since varA matrix is symmetrical, copy lower triangle to the upper
        for female in range(1, ntaxa):
            for male in range(0, female):
                varA_mat[:,male,female] = varA_mat[:,female,male]

        return varA_mat

    def calc_varA_3way_mat(self):
        raise NotImplementedError

    def calc_varA_3wayDH_mat(self):
        """
        Calculate a symmetrical matrix of progeny variance for each pairwise
        3-way cross.

        Assumes inbred individuals.
        """
        # create aliases for several variables
        geno = self._population.geno
        coeff = self._population.genomic_model.coeff
        gmap = self._population.marker_set
        mem = self._mem

        # get number of traits (t)
        ntrait = coeff.shape[1]

        # get the number of individuals
        ntaxa = geno.shape[1]

        # allocate a cube matrix for each 3-way variance
        varA_mat = numpy.zeros((ntrait, ntaxa, ntaxa, ntaxa), dtype='float64')

        # for each linkage group
        for lst, lsp in zip(gmap.chr_grp_stix, gmap.chr_grp_spix):
            # for each computational chunk
            for rst,rsp in zip(range(lst,lsp,mem),pybropt.util.srange(lst+mem,lsp,mem)):
                for cst,csp in zip(range(lst,lsp,mem),pybropt.util.srange(lst+mem,lsp,mem)):
                    # get recombination probability matrix for chunk
                    r = gmap.recomb_prob(rst, rsp, cst, csp) # (rb,cb)

                    # calculate a D1 matrix; this is specific to mating scheme
                    D1_mat = Cross.D1(r, self._s, self._t) # (rb,cb)

                    # calculate a D2 matrix; this is specific to mating scheme
                    D2_mat = Cross.D2(r, self._s, self._t) # (rb,cb)

                    # get marker coefficients for rows and columns
                    rcoeff = coeff[rst:rsp].T # (t,rb) # rb = row block
                    ccoeff = coeff[cst:csp].T # (t,cb) # cb = column block

                    # for each 3-way cross (excluding selfs)
                    # subscript codes:
                    #   1 = recurrent parent
                    #   2 = female
                    #   3 = male
                    for recurr in range(0,ntaxa):           # varA slice index
                        for female in range(0,ntaxa):       # varA row index
                            # calculate genotype differences for row, col
                            rdgeno21 = geno[0,female,rst:rsp] - geno[0,recurr,rst:rsp] # (rb,)
                            cdgeno21 = geno[0,female,cst:csp] - geno[0,recurr,cst:csp] # (cb,)

                            # calculate effect differences
                            reffect21 = rdgeno21 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                            ceffect21 = cdgeno21 * ccoeff # (cb,)*(t,cb) -> (t,cb)

                            # compute dot product for each trait to get partial variance
                            # (t,rb)x(rb,cb) -> (t,cb)
                            # (t,cb)*(t,cb) -> (t,cb)
                            # (t,cb)[1] -> (t,)
                            varA_part21 = (reffect21 @ D1_mat * ceffect21).sum(1)

                            # only do lower triangle since symmetrical within each slice
                            for male in range(0,female):    # varA col index
                                # calculate genotype differences for row, col
                                rdgeno23 = geno[0,female,rst:rsp] - geno[0,male,rst:rsp] # (rb,)
                                cdgeno23 = geno[0,female,cst:csp] - geno[0,male,cst:csp] # (cb,)
                                rdgeno31 = geno[0,male,rst:rsp] - geno[0,recurr,rst:rsp] # (rb,)
                                cdgeno31 = geno[0,male,cst:csp] - geno[0,recurr,cst:csp] # (cb,)

                                # calculate effect differences
                                reffect23 = rdgeno23 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                                ceffect23 = cdgeno23 * ccoeff # (cb,)*(t,cb) -> (t,cb)
                                reffect31 = rdgeno31 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                                ceffect31 = cdgeno31 * ccoeff # (cb,)*(t,cb) -> (t,cb)

                                # compute dot product for each trait to get partial variance
                                # (t,rb)x(rb,cb) -> (t,cb)
                                # (t,cb)*(t,cb) -> (t,cb)
                                # (t,cb)[1] -> (t,)
                                varA_part23 = (reffect23 @ D2_mat * ceffect23).sum(1)
                                varA_part31 = (reffect31 @ D1_mat * ceffect31).sum(1)

                                # calculate varA part for this matrix chunk
                                varA_part = (2.0 * (varA_part21 + varA_part31)) + varA_part23

                                # add this partial variance to the lower triangle
                                varA_mat[:,recurr,female,male] += varA_part

        # divide entire matrix by 4 to get variance per the equation
        varA_mat /= 4.0

        # each matrix is symmetrical within a slice because exchanging female
        # and male orders is mathematically equivalent.
        # copy lower triangle to the upper since varA matrix is symmetrical within each slice
        for female in range(1, ntaxa):
            for male in range(0, female):
                varA_mat[:,:,male,female] = varA_mat[:,:,female,male]

        return varA_mat

    def calc_varA_4way_mat(self):
        raise NotImplementedError

    def calc_varA_4wayDH_mat(self):
        """
        Calculate a symmetrical matrix of progeny variance for each pairwise
        4-way cross.

        Assumes inbred individuals.
        """
        # create aliases for several variables
        geno = self._population.geno
        coeff = self._population.genomic_model.coeff
        gmap = self._population.marker_set
        mem = self._mem

        # get number of traits (t)
        ntrait = coeff.shape[1]

        # get the number of individuals
        ntaxa = geno.shape[1]

        # allocate a square matrix for each pairwise variance
        varA_mat = numpy.zeros((ntrait, ntaxa, ntaxa, ntaxa, ntaxa), dtype='float64')

        # for each linkage group
        for lst, lsp in zip(gmap.chr_grp_stix, gmap.chr_grp_spix):
            # for each computational chunk
            for rst,rsp in zip(range(lst,lsp,mem),pybropt.util.srange(lst+mem,lsp,mem)):
                for cst,csp in zip(range(lst,lsp,mem),pybropt.util.srange(lst+mem,lsp,mem)):
                    # get recombination probability matrix for chunk
                    r = gmap.recomb_prob(rst, rsp, cst, csp) # (rb,cb)

                    # calculate a D1 matrix; this is specific to mating scheme
                    D1_mat = Cross.D1(r, self._s, self._t) # (rb,cb)

                    # calculate a D2 matrix; this is specific to mating scheme
                    D2_mat = Cross.D2(r, self._s, self._t) # (rb,cb)

                    # get marker coefficients for rows and columns
                    rcoeff = coeff[rst:rsp].T # (t,rb) # rb = row block
                    ccoeff = coeff[cst:csp].T # (t,cb) # cb = column block

                    # for each 4-way cross (excluding selfs)
                    # subscript codes:
                    #   1 = female 2
                    #   2 = male 2
                    #   3 = female 1
                    #   4 = male 1
                    # TODO: make this more computationally efficient.
                    #       operations can be simplified if you mirror sections
                    #       of the 4d matrix between blocks.
                    #       how to do this is difficult to conceive
                    for female2 in range(0,ntaxa):              # varA block index (change me for efficiency?)
                        for male2 in range(0,ntaxa):            # varA slice index (change me for efficiency?)
                            # calculate genotype differences for row, col
                            rdgeno21 = geno[0,male2,rst:rsp] - geno[0,female2,rst:rsp] # (rb,)
                            cdgeno21 = geno[0,male2,cst:csp] - geno[0,female2,cst:csp] # (cb,)

                            # calculate effect differences
                            reffect21 = rdgeno21 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                            ceffect21 = cdgeno21 * ccoeff # (cb,)*(t,cb) -> (t,cb)

                            # compute dot product for each trait to get partial variance
                            # (t,rb)x(rb,cb) -> (t,cb)
                            # (t,cb)*(t,cb) -> (t,cb)
                            # (t,cb)[1] -> (t,)
                            varA_part21 = (reffect21 @ D2_mat * ceffect21).sum(1)

                            for female1 in range(0,ntaxa):      # varA row index
                                # calculate genotype differences for row, col
                                rdgeno31 = geno[0,female1,rst:rsp] - geno[0,female2,rst:rsp] # (rb,)
                                cdgeno31 = geno[0,female1,cst:csp] - geno[0,female2,cst:csp] # (cb,)
                                rdgeno32 = geno[0,female1,rst:rsp] - geno[0,male2,rst:rsp] # (rb,)
                                cdgeno32 = geno[0,female1,cst:csp] - geno[0,male2,cst:csp] # (cb,)

                                # calculate effect differences
                                reffect31 = rdgeno31 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                                ceffect31 = cdgeno31 * ccoeff # (cb,)*(t,cb) -> (t,cb)
                                reffect32 = rdgeno32 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                                ceffect32 = cdgeno32 * ccoeff # (cb,)*(t,cb) -> (t,cb)

                                # compute dot product for each trait to get partial variance
                                # (t,rb)x(rb,cb) -> (t,cb)
                                # (t,cb)*(t,cb) -> (t,cb)
                                # (t,cb)[1] -> (t,)
                                varA_part31 = (reffect31 @ D1_mat * ceffect31).sum(1)
                                varA_part32 = (reffect32 @ D1_mat * ceffect32).sum(1)

                                # only do lower triangle since symmetrical within each slice
                                for male1 in range(0,female1):  # varA col index
                                    # calculate genotype differences for row, col
                                    rdgeno41 = geno[0,male1,rst:rsp] - geno[0,female2,rst:rsp] # (rb,)
                                    cdgeno41 = geno[0,male1,cst:csp] - geno[0,female2,cst:csp] # (cb,)
                                    rdgeno42 = geno[0,male1,rst:rsp] - geno[0,male2,rst:rsp] # (rb,)
                                    cdgeno42 = geno[0,male1,cst:csp] - geno[0,male2,cst:csp] # (cb,)
                                    rdgeno43 = geno[0,male1,rst:rsp] - geno[0,female1,rst:rsp] # (rb,)
                                    cdgeno43 = geno[0,male1,cst:csp] - geno[0,female1,cst:csp] # (cb,)

                                    # calculate effect differences
                                    reffect41 = rdgeno41 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                                    ceffect41 = rdgeno41 * ccoeff # (cb,)*(t,cb) -> (t,cb)
                                    reffect42 = rdgeno42 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                                    ceffect42 = rdgeno42 * ccoeff # (cb,)*(t,cb) -> (t,cb)
                                    reffect43 = rdgeno43 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                                    ceffect43 = rdgeno43 * ccoeff # (cb,)*(t,cb) -> (t,cb)

                                    # compute dot product for each trait to get partial variance
                                    # (t,rb)x(rb,cb) -> (t,cb)
                                    # (t,cb)*(t,cb) -> (t,cb)
                                    # (t,cb)[1] -> (t,)
                                    varA_part41 = (reffect41 @ D1_mat * ceffect41).sum(1)
                                    varA_part42 = (reffect42 @ D1_mat * ceffect42).sum(1)
                                    varA_part43 = (reffect43 @ D2_mat * ceffect43).sum(1)

                                    # calculate varA part for this matrix chunk
                                    varA_part = varA_part21 + varA_part31 + varA_part32 + varA_part41 + varA_part42 + varA_part43

                                    # add this partial variance to the lower triangle
                                    varA_mat[:,female2,male2,female1,male1] += varA_part

        # divide entire matrix by 4 to get variance per the equation
        varA_mat /= 4.0

        # each matrix is symmetrical within a slice because exchanging female
        # and male orders is mathematically equivalent.
        # copy lower triangle to the upper since varA matrix is symmetrical within each slice
        for female1 in range(1, ntaxa):
            for male1 in range(0, female1):
                varA_mat[:,:,:,male1,female1] = varA_mat[:,:,:,female1,male1]

        return varA_mat

    def calc_varA_dihybrid_mat(self):
        raise NotImplementedError

    def calc_varA_dihybridDH_mat(self):
        """
        Calculate a symmetrical matrix of progeny variance for each pairwise
        dihybrid cross.

        Assumes inbred individuals.
        """
        # create aliases for several variables
        geno = self._population.geno
        coeff = self._population.genomic_model.coeff
        gmap = self._population.marker_set
        mem = self._mem

        # get number of traits (t)
        ntrait = coeff.shape[1]

        # get the number of individuals
        ntaxa = geno.shape[1]

        # allocate a square matrix for each pairwise variance
        varA_mat = numpy.zeros((ntrait, ntaxa, ntaxa), dtype='float64')

        # for each linkage group
        for lst, lsp in zip(gmap.chr_grp_stix, gmap.chr_grp_spix):
            # for each computational chunk
            for rst,rsp in zip(range(lst,lsp,mem),pybropt.util.srange(lst+mem,lsp,mem)):
                for cst,csp in zip(range(lst,lsp,mem),pybropt.util.srange(lst+mem,lsp,mem)):
                    # get recombination probability matrix for chunk
                    r = gmap.recomb_prob(rst, rsp, cst, csp) # (rb,cb)

                    # calculate a D1 matrix; this is specific to mating scheme
                    D1_mat = Cross.D1(r, self._s, self._t) # (rb,cb)

                    # calculate a D2 matrix; this is specific to mating scheme
                    D2_mat = Cross.D2(r, self._s, self._t) # (rb,cb)

                    # get marker coefficients for rows and columns
                    rcoeff = coeff[rst:rsp].T # (t,rb) # rb = row block
                    ccoeff = coeff[cst:csp].T # (t,cb) # cb = column block

                    # for each mate pair (including selfs)
                    # subscript codes:
                    #   1 = female phase 2
                    #   2 = female phase 1
                    #   3 = male phase 2
                    #   4 = male phase 1
                    for female in range(0,ntaxa):     # varA row index
                        for male in range(0,female):    # varA col index
                            # calculate genotype differences for row, col
                            rdgeno21 = geno[0,female,rst:rsp] - geno[1,female,rst:rsp] # (rb,)
                            cdgeno21 = geno[0,female,cst:csp] - geno[1,female,cst:csp] # (cb,)
                            rdgeno31 = geno[1,male,rst:rsp] - geno[1,female,rst:rsp] # (rb,)
                            cdgeno31 = geno[1,male,cst:csp] - geno[1,female,cst:csp] # (cb,)
                            rdgeno32 = geno[1,male,rst:rsp] - geno[0,female,rst:rsp] # (rb,)
                            cdgeno32 = geno[1,male,cst:csp] - geno[0,female,cst:csp] # (cb,)
                            rdgeno41 = geno[0,male,rst:rsp] - geno[1,female,rst:rsp] # (rb,)
                            cdgeno41 = geno[0,male,cst:csp] - geno[1,female,cst:csp] # (cb,)
                            rdgeno42 = geno[0,male,rst:rsp] - geno[0,female,rst:rsp] # (rb,)
                            cdgeno42 = geno[0,male,cst:csp] - geno[0,female,cst:csp] # (cb,)
                            rdgeno43 = geno[0,male,rst:rsp] - geno[1,male,rst:rsp] # (rb,)
                            cdgeno43 = geno[0,male,cst:csp] - geno[1,male,cst:csp] # (cb,)

                            # calculate effect differences
                            reffect21 = rdgeno21 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                            ceffect21 = cdgeno21 * ccoeff # (cb,)*(t,cb) -> (t,cb)
                            reffect31 = rdgeno31 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                            ceffect31 = cdgeno31 * ccoeff # (cb,)*(t,cb) -> (t,cb)
                            reffect32 = rdgeno32 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                            ceffect32 = cdgeno32 * ccoeff # (cb,)*(t,cb) -> (t,cb)
                            reffect41 = rdgeno41 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                            ceffect41 = rdgeno41 * ccoeff # (cb,)*(t,cb) -> (t,cb)
                            reffect42 = rdgeno42 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                            ceffect42 = rdgeno42 * ccoeff # (cb,)*(t,cb) -> (t,cb)
                            reffect43 = rdgeno43 * rcoeff # (rb,)*(t,rb) -> (t,rb)
                            ceffect43 = rdgeno43 * ccoeff # (cb,)*(t,cb) -> (t,cb)

                            # compute dot product for each trait to get partial variance
                            # (t,rb)x(rb,cb) -> (t,cb)
                            # (t,cb)*(t,cb) -> (t,cb)
                            # (t,cb)[1] -> (t,)
                            varA_part21 = (reffect21 @ D2_mat * ceffect21).sum(1)
                            varA_part31 = (reffect31 @ D1_mat * ceffect31).sum(1)
                            varA_part32 = (reffect32 @ D1_mat * ceffect32).sum(1)
                            varA_part41 = (reffect41 @ D1_mat * ceffect41).sum(1)
                            varA_part42 = (reffect42 @ D1_mat * ceffect42).sum(1)
                            varA_part43 = (reffect43 @ D2_mat * ceffect43).sum(1)

                            # calculate varA part for this matrix chunk
                            varA_part = varA_part21 + varA_part31 + varA_part32 + varA_part41 + varA_part42 + varA_part43

                            # add this partial variance to the lower triangle
                            varA_mat[:,female,male] += varA_part

        # divide entire matrix by 4 to get variance per the equation
        varA_mat /= 4.0

        # each matrix is symmetrical within a slice because exchanging female
        # and male orders is mathematically equivalent.
        # copy lower triangle to the upper since varA matrix is symmetrical within each slice
        for female in range(1, ntaxa):
            for male in range(0, female):
                varA_mat[:,male,female] = varA_mat[:,female,male]

        return varA_mat

    def calc_varA_2way_sparse_mat(self, sel):
        raise NotImplementedError

    def calc_varA_2wayDH_sparse_mat(self, sel):
        raise NotImplementedError("Method not implemented.")

    def calc_varA_3way_sparse_mat(self, sel):
        raise NotImplementedError

    def calc_varA_3wayDH_sparse_mat(self, sel):
        raise NotImplementedError("Method not implemented.")

    def calc_varA_4way_sparse_mat(self, sel):
        raise NotImplementedError

    def calc_varA_4wayDH_sparse_mat(self, sel):
        raise NotImplementedError("Method not implemented.")

    def calc_varA_dihybrid_sparse_mat(self, sel):
        raise NotImplementedError

    def calc_varA_dihybridDH_sparse_mat(self, sel):
        raise NotImplementedError("Method not implemented.")

    ##############################################
    def set_varAfn(self, varAfn = None, sparse = False):
        """
        Set the varA() and varA_vec() attributes. Calling these will result in
        the gathering of a default variance component.

        Parameters
        ----------
        varAfn : str, None
            A key to lookup in internal lookup tables.
            Valid keys: {'2wayDH', '3wayDH', '4wayDH', 'dihybridDH'}
        sparse : boolean, default = False
            Whether 'varAfn' is a sparse matrix calculation.
        """
        # lookup varAfn in lookup table
        attribute_varA, attribute_varA_vec = Cross.KEY_TO_VARAFN[varAfn][sparse]

        # set new function attributes
        self.varA = getattr(self, attribute_varA)
        self.varA_vec = getattr(self, attribute_varA_vec)

        # set private variables
        self._varAfn = varAfn
        self._sparse = sparse

    def varA(self, sel, *args, **kwargs):
        raise RuntimeError("varA function not set.")

    def varA_vec(self, sel, *args, **kwargs):
        raise RuntimeError("varA function not set.")

    def varA_default(self, sel, *args, **kwargs):
        """
        Default method for varA method. Raises a RuntimeError.
        """
        raise RuntimeError("varA function not set.")
    ##############################################

    ##############################################
    ########### Non-sparse matrix ops ############
    def varA_2way(self, sel):
        raise NotImplementedError()

    def varA_2wayDH(self, sel):
        """
        Retrieve additive variance components for a 2-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7

        Returns
        -------
        varA_val : numpy.ndarray
            A 1D array of variance components of shape (t,k/2)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_2wayDH_mat"):
            self._varA_2wayDH_mat = self.calc_varA_2wayDH_mat()

        # get female, male indices
        female = sel[0::2]
        male = sel[1::2]

        # get variance terms
        varA_val = self._varA_2wayDH_mat[:,female,male]

        # return variance terms
        return varA_val

    def varA_3way(self, sel):
        raise NotImplementedError()

    def varA_3wayDH(self, sel):
        """
        Retrieve additive variance components for a 3-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the recurrent parent.
                Second index is the female parent.
                Third index is the male parent.
            Example:
                [1,5,3,8,2,7]
                recurrent = 1,8
                female = 5,2
                male = 3,7

        Returns
        -------
        varA_val : numpy.ndarray
            A 1D array of variance components of shape (k/3,)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_3wayDH_mat"):
            self._varA_3wayDH_mat = self.calc_varA_3wayDH_mat()

        # get recurrent, female, male indices
        recurr = sel[0::3]
        female = sel[1::3]
        male = sel[2::3]

        # get variance terms
        varA_val = self._varA_3wayDH_mat[:,recurr,female,male]

        # return variance terms
        return varA_val

    def varA_4way(self, sel):
        raise NotImplementedError()

    def varA_4wayDH(self, sel):
        """
        Retrieve additive variance components for a 4-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals, divisible by 4.
            Indices are paired as follows:
                First index is the female parent 2.
                Second index is the male parent 2.
                Third index is the female parent 1.
                Fourth index is the male parent 1.
            Example:
                [1,5,3,8]
                female2 = 1
                male2 = 5
                female1 = 3
                male1 = 8

        Returns
        -------
        varA_val : numpy.ndarray
            A 1D array of variance components of shape (k/4,)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_4wayDH_mat"):
            self._varA_4wayDH_mat = self.calc_varA_4wayDH_mat()

        # get female2, male2, female1, male1 indices
        female2 = sel[0::4]
        male2 = sel[1::4]
        female = sel[2::4]
        male = sel[3::4]

        # get variance terms
        varA_val = self._varA_4wayDH_mat[:,female2,male2,female1,male1]

        # return variance terms
        return varA_val

    def varA_dihybrid(self, sel):
        raise NotImplementedError()

    def varA_dihybridDH(self, sel):
        """
        Retrieve additive variance components for a dihybrid DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7

        Returns
        -------
        varA_val : numpy.ndarray
            A 1D array of variance components of shape (k/2,)
        """
        if not hasattr(self, "_varA_dihybridDH_mat"):
            self._varA_dihybridDH_mat = self.calc_varA_dihybridDH_mat()

        # get female, male indices
        female = sel[0::2]
        male = sel[1::2]

        # get variance terms
        varA_val = self._varA_dihybridDH_mat[:,female,male]

        # return variance terms
        return varA_val

    def varA_2wayDH_vec(self, sel):
        """
        Retrieve additive variance components for a 2-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 2D array of indices of selected individuals of shape (j,k).
            Where:
                'j' is the number of mating configurations.
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [[1,5,3,8,2,7],
                 [2,7,3,5,4,6]]
                female = [1,3,2],[2,3,4]
                male = [5,8,7],[7,5,6]

        Returns
        -------
        varA_val : numpy.ndarray
            A 2D array of variance components of shape (j,k/2)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_2wayDH_mat"):
            self._varA_2wayDH_mat = self.calc_varA_2wayDH_mat()

        # get female, male indices
        female = sel[:,0::2]
        male = sel[:,1::2]

        # get variance terms
        varA_val = self._varA_2wayDH_mat[:,female,male]

        # return variance terms
        return varA_val

    def varA_3wayDH_vec(self, sel):
        """
        Retrieve additive variance components for a 3-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 2D array of indices of selected individuals of shape (j,k).
            Where:
                'j' is the number of mating configurations.
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the recurrent parent.
                Second index is the female parent.
                Third index is the male parent.
            Example:
                [[1,5,3,8,2,7],
                 [2,7,3,5,4,6]]
                recurrent = [1,8],[2,5]
                female = [5,2],[7,4]
                male = [3,7],[3,6]

        Returns
        -------
        varA_val : numpy.ndarray
            A 2D array of variance components of shape (j,k/3)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_3wayDH_mat"):
            self._varA_3wayDH_mat = self.calc_varA_3wayDH_mat()

        # get recurrent, female, male indices
        recurr = sel[:,0::3]
        female = sel[:,1::3]
        male = sel[:,2::3]

        # get variance terms
        varA_val = self._varA_3wayDH_mat[:,recurr,female,male]

        # return variance terms
        return varA_val

    def varA_4wayDH_vec(self, sel):
        """
        Retrieve additive variance components for a 4-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 2D array of indices of selected individuals of shape (j,k).
            Where:
                'j' is the number of mating configurations.
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the female parent 2.
                Second index is the male parent 2.
                Third index is the female parent 1.
                Fourth index is the male parent 1.
            Example:
                [[1,5,3,8]
                 [7,2,9,0]]
                female2 = [1],[7]
                male2 = [5],[2]
                female1 = [3],[9]
                male1 = [8],[0]

        Returns
        -------
        varA_val : numpy.ndarray
            A 2D array of variance components of shape (j,k/4)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_4wayDH_mat"):
            self._varA_4wayDH_mat = self.calc_varA_4wayDH_mat()

        # get female2, male2, female1, male1 indices
        female2 = sel[:,0::4]
        male2 = sel[:,1::4]
        female = sel[:,2::4]
        male = sel[:,3::4]

        # get variance terms
        varA_val = self._varA_4wayDH_mat[:,female2,male2,female1,male1]

        # return variance terms
        return varA_val

    def varA_dihybridDH_vec(self, sel):
        """
        Retrieve additive variance components for a dihybrid DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 2D array of indices of selected individuals of shape (j,k).
            Where:
                'j' is the number of mating configurations.
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [[1,5,3,8,2,7],
                 [2,7,3,5,4,6]]
                female = [1,3,2],[2,3,4]
                male = [5,8,7],[7,5,6]

        Returns
        -------
        varA_val : numpy.ndarray
            A 2D array of variance components of shape (j,k/2)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_dihybridDH_mat"):
            self._varA_dihybridDH_mat = self.calc_varA_dihybridDH_mat()

        # get female, male indices
        female = sel[:,0::2]
        male = sel[:,1::2]

        # get variance terms
        varA_val = self._varA_dihybridDH_mat[:,female,male]

        # return variance terms
        return varA_val
    ##############################################

    ##############################################
    ############# Sparse matrix ops ##############
    def varA_2wayDH_sparse(self, sel):
        raise NotImplementedError("Method not implemented.")

    def varA_3wayDH_sparse(self, sel):
        raise NotImplementedError("Method not implemented.")

    def varA_4wayDH_sparse(self, sel):
        raise NotImplementedError("Method not implemented.")

    def varA_dihybridDH_sparse(self, sel):
        raise NotImplementedError("Method not implemented.")

    def varA_2wayDH_sparse_vec(self, sel):
        raise NotImplementedError("Method not implemented.")

    def varA_3wayDH_sparse_vec(self, sel):
        raise NotImplementedError("Method not implemented.")

    def varA_4wayDH_sparse_vec(self, sel):
        raise NotImplementedError("Method not implemented.")

    def varA_dihybridDH_sparse_vec(self, sel):
        raise NotImplementedError("Method not implemented.")
    ##############################################
    ############################################################################


    ############################################################################
    ####################### Mating Related Class Methods #######################
    def meiosis(self, sel, seed = None):
        """
        Simulate meiosis. Generate gametes from a genotype matrix.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.

        Returns
        -------
        gametes : numpy.ndarray
        """
        # make gametes
        gamete = Cross.meiosis_mat(
            sel,
            self._population.geno,
            self._genetic_map,
            seed
        )

        return gamete

    def set_matefn(self, crossfn = None, matefn = None):
        # lookup matefn in lookup table
        attribute_string = Cross.KEY_TO_MATEFN[crossfn][matefn]

        # set new function attributes
        self.mate = getattr(self, attribute_string)

        # set private variables
        self._crossfn = crossfn
        self._matefn = matefn

    def mate(self, sel, seed = None, *args, **kwargs):
        raise RuntimeError("mate function not set.")

    def mate_default(self, sel, seed = None, *args, **kwargs):
        """
        Default method for mate method. Raises a RuntimeError.
        """
        raise RuntimeError("mate function not set.")

    ##############################################
    ############# Controlled mating ##############
    def mate_2way_ctrl(self, sel, c = None, n = None, seed = None):
        """
        Perform a 2-way cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7
        c : numpy.ndarray, int
            A 1D array of integers representing the number of times the
            specified cross pattern should be performed.
            shape = (k/2,)
        n : numpy.ndarray, int
            A 1D array of integers representing the number of progeny from each
            cross that should be simulated.
            Remark:
                The number of progeny is therefore c * n.
        seed : numpy.int32, None
            Random seed to initialize the random number generator. If None, do
            not alter the RNG's state.

        Returns
        -------
        population : Population
            A new population of progeny individuals.
        """
        # if all required variables are None, use ralloc to select them
        if all(v is None for v in [c, n]):
            c, n = self.ralloc(sel)

        # seed rng
        pybropt.util.cond_seed_rng(seed)

        # repeat the female, male parents 'c' times
        fsel = numpy.repeat(sel[0::2], c)
        msel = numpy.repeat(sel[1::2], c)

        # allocate empty array
        matepair = numpy.empty(len(fsel) * 2, dtype = fsel.dtype)

        # put female, male pairs in
        matepair[0::2] = fsel
        matepair[1::2] = msel

        # make genotypes
        geno = Cross.mate_mat(
            matepair,
            self._population.geno,
            self._population.marker_set,
            n
        )

        # TODO: build taxa name array

        population = Population(
            geno,
            self._population.genomic_model,
            self._population.marker_set
        )

        return population

    def mate_2wayDH_ctrl(self, sel, c = None, n = None, s = None, t = None, seed = None):
        """
        Perform a 2-way cross with double haploiding.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7
        c : numpy.ndarray, int
            A 1D array of integers representing the number of times the
            specified cross pattern should be performed.
        n : numpy.ndarray, int
            A 1D array of integers representing the number of progeny from each
            cross that should be simulated.
            Remark:
                The number of progeny is therefore c * n.
        s : int, None
            Number of times progeny from the cross should be self-fertilized.
            If 's' is None, use self.s as the value.
        t : int, None
            Not implemented yet.
            Number of times progeny from the cross should be randomly
            intermated.
            If 't' is None, use self.t as the value.
            Remark:
                's' and 't' are mutually exclusive. If 's' is nonzero, selfings
                are simulated.
        seed : numpy.int32, None
            Random seed to initialize the random number generator. If None, do
            not alter the RNG's state.

        Returns
        -------
        population : Population
            A new population of progeny individuals.
        """
        # if all required variables are None, use ralloc to select them
        if all(v is None for v in [c, n, s, t]):
            c, n, s, t = self.ralloc(sel)

        # seed rng
        pybropt.util.cond_seed_rng(seed)

        # repeat the female, male parents 'c' times
        fsel = numpy.repeat(sel[0::2], c)
        msel = numpy.repeat(sel[1::2], c)

        # allocate empty array
        matepair = numpy.empty(len(fsel) * 2, dtype = fsel.dtype)

        # put female, male pairs in
        matepair[0::2] = fsel
        matepair[1::2] = msel

        # declare 'geno' variable
        geno = None

        if t == 0:
            # make genotypes
            geno = Cross.mate_mat(
                matepair,
                self._population.geno,
                self._population.marker_set,
                1
            )

            # make doubled haploids from the hybrid genotypes
            geno = Cross.dh_mat(
                numpy.arange(geno.shape[1]), # all hybrids
                geno,
                self._population.marker_set,
                n,
                s,
            )
        else:
            raise NotImplementedError("'t' protocol not implemented yet.")

        # TODO: build taxa name array

        population = Population(
            geno,
            self._population.marker_set,
            self._population.genomic_model
        )

        return population

    def mate_3way_ctrl(self, sel, c = None, n = None, seed = None):
        """
        Perform a 3-way cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the recurrent parent.
                Second index is the female parent.
                Third index is the male parent.
            Example:
                [1,5,3,8,2,7]
                recurrent = 1,8
                female = 5,2
                male = 3,7
        c : numpy.ndarray, int
            A 1D array of integers representing the number of times the
            specified cross pattern should be performed.
        n : numpy.ndarray, int
            A 1D array of integers representing the number of progeny from each
            cross that should be simulated.
            Remark:
                The number of progeny is therefore c * n.
        seed : numpy.int32, None
            Random seed to initialize the random number generator. If None, do
            not alter the RNG's state.

        Returns
        -------
        population : Population
            A new population of progeny individuals.
        """
        # if all required variables are None, use ralloc to select them
        if all(v is None for v in [c, n]):
            c, n = self.ralloc(sel)

        # seed rng
        pybropt.util.cond_seed_rng(seed)

        # repeat the recurrent, female, male parents 'c' times
        rsel = numpy.repeat(sel[0::3], c)
        fsel = numpy.repeat(sel[1::3], c)
        msel = numpy.repeat(sel[2::3], c)

        # allocate empty array for cross 1 (female x male) specification
        cross1 = numpy.empty(len(fsel) * 2, dtype = fsel.dtype)

        # put female, male pairs in array
        cross1[0::2] = fsel
        cross1[1::2] = msel

        # make hybrid (female x male) genotypes
        geno = Cross.mate_mat(
            cross1,
            self._population.geno,
            self._population.marker_set,
            1
        )

        # concatenate hybrid and original matrices together
        geno = numpy.concatenate([geno, self._population.geno], axis=1)

        # get indices for each hybrid: there are len(fsel) hybrids
        hsel = numpy.arange(len(fsel))

        # allocate empty array for cross 2 (recurrent x hybrid) specification
        cross2 = numpy.empty(len(rsel) * 2, dtype = rsel.dtype)

        # put recurrent, hybrid pairs in array
        cross2[0::2] = (rsel + len(fsel)) # offset rsel indices by hybrid number
        cross2[1::2] = hsel

        geno = Cross.mate_mat(
            cross2,
            geno,
            self._population.marker_set,
            n
        )

        # TODO: build taxa name array

        population = Population(
            geno,
            self._population.genomic_model,
            self._population.marker_set
        )

        return population

    def mate_3wayDH_ctrl(self, sel, c = None, n = None, s = None, t = None, seed = None):
        """
        Perform a 3-way cross with double haploiding.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the recurrent parent.
                Second index is the female parent.
                Third index is the male parent.
            Example:
                [1,5,3,8,2,7]
                recurrent = 1,8
                female = 5,2
                male = 3,7
        c : numpy.ndarray, int
            A 1D array of integers representing the number of times the
            specified cross pattern should be performed.
        n : numpy.ndarray, int
            A 1D array of integers representing the number of progeny from each
            cross that should be simulated.
            Remark:
                The number of progeny is therefore c * n.
        s : int, None
            Number of times progeny from the cross should be self-fertilized.
            If 's' is None, use self.s as the value.
        t : int, None
            Number of times progeny from the cross should be randomly
            intermated.
            If 't' is None, use self.t as the value.
            Remark:
                's' and 't' are mutually exclusive. If 's' is nonzero, selfings
                are simulated.
        seed : numpy.int32, None
            Random seed to initialize the random number generator. If None, do
            not alter the RNG's state.

        Returns
        -------
        population : Population
            A new population of progeny individuals.
        """
        # if all required variables are None, use ralloc to select them
        if all(v is None for v in [c, n, s, t]):
            c, n, s, t = self.ralloc(sel)

        # seed rng
        pybropt.util.cond_seed_rng(seed)

        # repeat the recurrent, female, male parents 'c' times
        rsel = numpy.repeat(sel[0::3], c)
        fsel = numpy.repeat(sel[1::3], c)
        msel = numpy.repeat(sel[2::3], c)

        # declare 'geno' variable
        geno = None

        if t == 0:
            # allocate empty array for cross 1 (female x male) specification
            cross1 = numpy.empty(len(fsel) * 2, dtype = fsel.dtype)

            # put female, male pairs in array
            cross1[0::2] = fsel
            cross1[1::2] = msel

            # make hybrid (female x male) genotypes
            geno = Cross.mate_mat(
                cross1,
                self._population.geno,
                self._population.marker_set,
                1
            )

            # concatenate hybrid and original matrices together
            geno = numpy.concatenate([geno, self._population.geno], axis=1)

            # get indices for each hybrid: there are len(fsel) hybrids
            hsel = numpy.arange(len(fsel))

            # allocate empty array for cross 2 (recurrent x hybrid) specification
            cross2 = numpy.empty(len(rsel) * 2, dtype = rsel.dtype)

            # put recurrent, hybrid pairs in array
            cross2[0::2] = (rsel + len(fsel)) # offset rsel indices by hybrid number
            cross2[1::2] = hsel

            geno = Cross.dh_mat(
                cross2,
                geno,
                self._population.marker_set,
                n,
                s
            )
        else:
            raise NotImplementedError("'t' protocol not implemented yet.")

        # TODO: build taxa name array

        population = Population(
            geno,
            self._population.genomic_model,
            self._population.marker_set
        )

        return population

    def mate_4way_ctrl(self, sel, c = None, n = None, seed = None):
        """
        Perform a 4-way cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the female parent 2.
                Second index is the male parent 2.
                Third index is the female parent 1.
                Fourth index is the male parent 1.
            Example:
                [1,5,3,8]
                female2 = 1
                male2 = 5
                female1 = 3
                male1 = 8
        c : numpy.ndarray, int
            A 1D array of integers representing the number of times the
            specified cross pattern should be performed.
        n : numpy.ndarray, int
            A 1D array of integers representing the number of progeny from each
            cross that should be simulated.
            Remark:
                The number of progeny is therefore c * n.
        seed : numpy.int32, None
            Random seed to initialize the random number generator. If None, do
            not alter the RNG's state.

        Returns
        -------
        population : Population
            A new population of progeny individuals.
        """
        raise NotImplementedError()

    def mate_4wayDH_ctrl(self, sel, c = None, n = None, s = None, t = None, seed = None):
        """
        Perform a 4-way cross with double haploiding.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the female parent 2.
                Second index is the male parent 2.
                Third index is the female parent 1.
                Fourth index is the male parent 1.
            Example:
                [1,5,3,8]
                female2 = 1
                male2 = 5
                female1 = 3
                male1 = 8
        c : numpy.ndarray, int
            A 1D array of integers representing the number of times the
            specified cross pattern should be performed.
        n : numpy.ndarray, int
            A 1D array of integers representing the number of progeny from each
            cross that should be simulated.
            Remark:
                The number of progeny is therefore c * n.
        s : int, None
            Number of times progeny from the cross should be self-fertilized.
            If 's' is None, use self.s as the value.
        t : int, None
            Number of times progeny from the cross should be randomly
            intermated.
            If 't' is None, use self.t as the value.
            Remark:
                's' and 't' are mutually exclusive. If 's' is nonzero, selfings
                are simulated.
        seed : numpy.int32, None
            Random seed to initialize the random number generator. If None, do
            not alter the RNG's state.

        Returns
        -------
        population : Population
            A new population of progeny individuals.
        """
        raise NotImplementedError()
    ##############################################

    ##############################################
    ########### Weighted random mating ###########
    def mate_2way_wrand(self, sel, c = None, n = None, weight = None, seed = None):
        """
        Perform a 2-way cross with weighted random mating.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7
        weight : numpy.ndarray
            A 1D array of contribution weights for selected individuals of
            shape (k,). Elements in this array are floating points.
            Where:
                'k' is the number of selected individuals.
            Elements are paired as follows:
                Even indices are female weights.
                Odd indices are male weights.
            Assumptions:
                Female weights should sum to 1.
                Male weights should sum to 1.
        c : int
            An integers representing the number of random crosses that should
            be performed.
        n : numpy.ndarray, int
            A 1D array of integers representing the number of progeny from each
            cross that should be simulated.
            Remark:
                The number of progeny is therefore c * n.
        seed : numpy.int32, None
            Random seed to initialize the random number generator. If None, do
            not alter the RNG's state.

        Returns
        -------
        population : Population
            A new population of progeny individuals.
        """
        # if all required variables are None, use ralloc to select them
        if all(v is None for v in [c, n, weight]):
            c, n, weight = self.ralloc(sel)

        # seed rng if needed
        pybropt.util.cond_seed_rng(seed)

        # randomly choose indices for female, male parents
        fsel = numpy.random.choice(sel[0::2], c, True, weight[0::2])
        msel = numpy.random.choice(sel[1::2], c, True, weight[1::2])

        # allocate mate pair matrix
        matepair = numpy.empty(len(fsel) * 2, dtype = fsel.dtype)

        # copy female, male mates into matepair
        matepair[0::2] = fsel
        matepair[1::2] = msel

        # make progeny
        progeny = self.mate_2way_ctrl(matepair, 1, n)

        return progeny

    def mate_2wayDH_wrand(self, sel, c = None, n = None, weight = None, s = None, t = None, seed = None):
        # if all required variables are None, use ralloc to select them
        if all(v is None for v in [c, n, weight, s, t]):
            c, n, weight, s, t = self.ralloc(sel)

        # seed rng if needed
        pybropt.util.cond_seed_rng(seed)

        # randomly choose indices for female, male parents
        fsel = numpy.random.choice(sel[0::2], c, True, weight[0::2])
        msel = numpy.random.choice(sel[1::2], c, True, weight[1::2])

        # allocate mate pair matrix
        matepair = numpy.empty(len(fsel) * 2, dtype = fsel.dtype)

        # copy female, male mates into matepair
        matepair[0::2] = fsel
        matepair[1::2] = msel

        # make progeny
        progeny = self.mate_2wayDH_ctrl(matepair, 1, n, s, t)

        return progeny

    def mate_3way_wrand(self, sel, c = None, n = None, weight = None, seed = None):
        raise NotImplementedError()

    def mate_3wayDH_wrand(self, sel, c = None, n = None, weight = None, s = None, t = None, seed = None):
        raise NotImplementedError()

    def mate_4way_wrand(self, sel, c = None, n = None, weight = None, seed = None):
        raise NotImplementedError()

    def mate_4wayDH_wrand(self, sel, c = None, n = None, weight = None, s = None, t = None, seed = None):
        raise NotImplementedError()
    ##############################################

    ##############################################
    ############ Exact random mating #############
    def mate_2way_erand(self, sel, c = None, n = None, exact = None, seed = None):
        """
        Perform a 2-way cross with exact random mating.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7
        exact : numpy.ndarray
            A 1D array of contributions for selected individuals of shape (k,).
            Elements in this array are integers.
            Where:
                'k' is the number of selected individuals.
            Elements are paired as follows:
                Even indices are female weights.
                Odd indices are male weights.
            Assumptions:
                Female weights should sum to 1.
                Male weights should sum to 1.
        c : int
            A multiplier for 'exact' contributions.
        n : numpy.ndarray, int
            A 1D array of integers representing the number of progeny from each
            cross that should be simulated.
            Remark:
                The number of progeny is therefore c * n.
        seed : numpy.int32, None
            Random seed to initialize the random number generator. If None, do
            not alter the RNG's state.

        Returns
        -------
        population : Population
            A new population of progeny individuals.
        """
        # if all required variables are None, use ralloc to select them
        if all(v is None for v in [c, n, exact]):
            c, n, exact = self.ralloc(sel)

        # seed rng if needed
        pybropt.util.cond_seed_rng(seed)

        # replicate female, male exact number of times
        fsel = numpy.repeat(sel[0::2], c * exact[0::2])
        msel = numpy.repeat(sel[1::2], c * exact[1::2])

        # randomly shuffle parents
        numpy.random.shuffle(fsel)
        numpy.random.shuffle(msel)

        # allocate mate pair matrix
        matepair = numpy.empty(len(fsel) * 2, dtype = fsel.dtype)

        # copy female, male mates into matepair
        matepair[0::2] = fsel
        matepair[1::2] = msel

        # make progeny
        progeny = self.mate_2way_ctrl(matepair, 1, n)

        return progeny

    def mate_2wayDH_erand(self, sel, c = None, n = None, exact = None, s = None, t = None, seed = None):
        # if all required variables are None, use ralloc to select them
        if all(v is None for v in [c, n, exact, s, t]):
            c, n, exact, s, t = self.ralloc(sel)

        # seed rng if needed
        pybropt.util.cond_seed_rng(seed)

        # replicate female, male exact number of times
        fsel = numpy.repeat(sel[0::2], c * exact[0::2])
        msel = numpy.repeat(sel[1::2], c * exact[1::2])

        # randomly shuffle parents
        numpy.random.shuffle(fsel)
        numpy.random.shuffle(msel)

        # allocate mate pair matrix
        matepair = numpy.empty(len(fsel) * 2, dtype = fsel.dtype)

        # copy female, male mates into matepair
        matepair[0::2] = fsel
        matepair[1::2] = msel

        # make progeny
        progeny = self.mate_2wayDH_ctrl(matepair, 1, n, s, t)

        return progeny

    def mate_3way_erand(self, sel, c = None, n = None, exact = None, seed = None):
        raise NotImplementedError()

    def mate_3wayDH_erand(self, sel, c = None, n = None, exact = None, s = None, t = None, seed = None):
        raise NotImplementedError()

    def mate_4way_erand(self, sel, c = None, n = None, exact = None, seed = None):
        raise NotImplementedError()

    def mate_4wayDH_erand(self, sel, c = None, n = None, exact = None, s = None, t = None, seed = None):
        raise NotImplementedError()
    ##############################################
    ############################################################################


    ############################################################################
    ################# Mating Resource Allocation Class Methods #################
    def set_rallocfn(self, crossfn, matefn, rallocfn):
        # get attribute string
        attribute_string = self.KEY_TO_RALLOCFN[crossfn][matefn][rallocfn]

        # set new function attributes
        self.ralloc = getattr(self, attribute_string)

        # set private variables
        self._crossfn = crossfn
        self._matefn = matefn
        self._rallocfn = rallocfn

    def ralloc(self, sel):
        raise RuntimeError("Resource allocation function not set.")

    def ralloc_default(self, sel):
        raise RuntimeError("Resource allocation function not set.")

    ##############################################
    ############# Controlled mating ##############
    def ralloc_2way_ctrl_equal(self, sel):
        return self._c, self._n

    def ralloc_2wayDH_ctrl_equal(self, sel):
        return self._c, self._n, self._s, self._t

    def ralloc_3way_ctrl_equal(self, sel):
        return self._c, self._n

    def ralloc_3wayDH_ctrl_equal(self, sel):
        return self._c, self._n, self._s, self._t

    def ralloc_4way_ctrl_equal(self, sel):
        return self._c, self._n

    def ralloc_4wayDH_ctrl_equal(self, sel):
        return self._c, self._n, self._s, self._t
    ##############################################

    ##############################################
    ########### Weighted random mating ###########
    def ralloc_2way_wrand_equal(self, sel):
        # make a ones matrix
        weight = numpy.ones(len(sel), dtype = 'float64')

        # divide ones by number of mate sets
        weight /= (len(sel) / 2.0)

        # return tuple of resource allocations
        return self._c, self._n, weight

    def ralloc_2wayDH_wrand_equal(self, sel):
        # make a ones matrix
        weight = numpy.ones(len(sel), dtype = 'float64')

        # divide ones by number of mate sets
        weight /= (len(sel) / 2.0)

        # return tuple of resource allocations
        return self._c, self._n, weight, self._s, self._t

    def ralloc_3way_wrand_equal(self, sel):
        # make a ones matrix
        weight = numpy.ones(len(sel), dtype = 'float64')

        # divide ones by number of mate sets
        weight /= (len(sel) / 3.0)

        # return tuple of resource allocations
        return self._c, self._n, weight

    def ralloc_3wayDH_wrand_equal(self, sel):
        # make a ones matrix
        weight = numpy.ones(len(sel), dtype = 'float64')

        # divide ones by number of mate sets
        weight /= (len(sel) / 3.0)

        # return tuple of resource allocations
        return self._c, self._n, weight, self._s, self._t

    def ralloc_4way_wrand_equal(self, sel):
        # make a ones matrix
        weight = numpy.ones(len(sel), dtype = 'float64')

        # divide ones by number of mate sets
        weight /= (len(sel) / 4.0)

        # return tuple of resource allocations
        return self._c, self._n, weight

    def ralloc_4wayDH_wrand_equal(self, sel):
        # make a ones matrix
        weight = numpy.ones(len(sel), dtype = 'float64')

        # divide ones by number of mate sets
        weight /= (len(sel) / 4.0)

        # return tuple of resource allocations
        return self._c, self._n, weight, self._s, self._t
    ##############################################

    ##############################################
    ########### Weighted random mating ###########
    def ralloc_2way_erand_equal(self, sel):
        # make an exact matrix
        exact = numpy.ones(len(sel), dtype = 'int64')

        # return tuple of resource allocations
        return self._c, self._n, exact

    def ralloc_2wayDH_erand_equal(self, sel):
        # make an exact matrix
        exact = numpy.ones(len(sel), dtype = 'int64')

        # return tuple of resource allocations
        return self._c, self._n, exact, self._s, self._t

    def ralloc_3way_erand_equal(self, sel):
        # make an exact matrix
        exact = numpy.ones(len(sel), dtype = 'int64')

        # return tuple of resource allocations
        return self._c, self._n, exact

    def ralloc_3wayDH_erand_equal(self, sel):
        # make an exact matrix
        exact = numpy.ones(len(sel), dtype = 'int64')

        # return tuple of resource allocations
        return self._c, self._n, exact, self._s, self._t

    def ralloc_4way_erand_equal(self, sel):
        # make an exact matrix
        exact = numpy.ones(len(sel), dtype = 'int64')

        # return tuple of resource allocations
        return self._c, self._n, exact

    def ralloc_4wayDH_erand_equal(self, sel):
        # make an exact matrix
        exact = numpy.ones(len(sel), dtype = 'int64')

        # return tuple of resource allocations
        return self._c, self._n, exact, self._s, self._t
    ##############################################
    ############################################################################

    def copy(self):
        """
        Create copy of self.

        Returns
        -------
        cross : Cross
            A Cross object with identical properties to the originating object.
        """
        # create new Cross object (we use self.__class__ in case of inheritance)
        cross = self.__class__(
            population  = self._population,
            varAfn      = self._varAfn,
            sparse      = self._sparse,
            crossfn     = self._crossfn,
            matefn      = self._matefn,
            rallocfn    = self._rallocfn,
            c           = self._c,
            n           = self._n,
            s           = self._s,
            t           = self._t,
            mem         = self._mem
        )

        return cross

    def from_self(self, population = None, varAfn = None, sparse = None,
        crossfn = None, matefn = None, rallocfn = None, c = None, n = None,
        s = None, t = None, mem = None):
        """
        Create a new Cross object from the current Cross object with options to
        alter variables (e.g. population).

        Returns
        -------
        cross : Cross
            A Cross object.
        """
        # create new Cross object (we use self.__class__ in case of inheritance)
        cross = self.__class__(
            population  = self._population if population is None else population,
            varAfn      = self._varAfn if varAfn is None else varAfn,
            sparse      = self._sparse if sparse is None else sparse,
            crossfn     = self._crossfn if crossfn is None else crossfn,
            matefn      = self._matefn if matefn is None else matefn,
            rallocfn    = self._rallocfn if rallocfn is None else rallocfn,
            c           = self._c if c is None else c,
            n           = self._n if n is None else n,
            s           = self._s if s is None else s,
            t           = self._t if t is None else t,
            mem         = self._mem if mem is None else mem
        )

        return cross

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def rk(r, k):
        """
        Calculate the expected observed recombination rate between two loci at
        filial generation 'k' as specified by Lehermeier (2017).

        Parameters
        ----------
        r : numpy.ndarray
            Recombination probability matrix.
        k : int, inf
            Selfing filial generation number to derive gametes from.
            Example:
                k = 1    ->  Derive gametes from F1
                k = 2    ->  Derive gametes from F2
                ...
                k = inf  ->  Derive gametes from SSD

        Returns
        -------
        r_k : numpy.ndarray
            Expected observed recombination rate between two loci at filial
            generation 'k'.
        """
        # multiply matrix by two
        two_r = 2.0 * r

        # calculate first component of r_k
        r_k = two_r / (1.0 + two_r)

        # if k < inf, we do not have SSD and second term is needed
        if k < numpy.inf:
            r_k *= (1.0 - ( (0.5**k) * ((1.0 - two_r)**k) ) )

        # return result
        return r_k

    @staticmethod
    def D1(r, s, t):
        """
        Calculate a D_1 matrix. This matrix represents a linkage disequilibrium
        prototype for one recombination event prior to selfing or random
        intermating.

        Parameters
        ----------
        r : numpy.ndarray
            Recombination probability matrix.
        s : int, inf
            Selfing generation number to derive gametes from.
            Example:
                s = 0    ->  Derive gametes from F1
                s = 1    ->  Derive gametes from F2
                s = 2    ->  Derive gametes from F3
                ...
                k = inf  ->  Derive gametes from SSD
            Remark:
                's' takes priority over 't'. If s > 0, 't' calculations are not
                performed.
        t : int, inf
            Random intermating generation number to derive gametes from.
            Example:
                t = 0    ->  Derive gametes from F1
                t = 1    ->  Derive gametes from (F1 x F1)
                t = 2    ->  Derive gametes from (F1 x F1) x (F1 x F1)
                ...
                t = inf  ->  Derive gametes from unlimited random intermating
            Remark:
                's' takes priority over 't'. If s > 0, 't' calculations are not
                performed.

        Returns
        -------
        D_1 : numpy.ndarray
            A D_1 matrix.
        """
        if s == 0 and t == 0:
            return (1 - 2*r)
        elif s > 0:
            rk = Cross.rk(r, s+1)
            return (1.0 - 2.0*rk)
        elif t > 0:
            return (1.0 - 2.0*r) * ((1.0 - r)**t)
        else:
            raise ValueError("s and t must be >= 0.")

    @staticmethod
    def D2(r, s, t):
        """
        Calculate a D_2 matrix. This matrix represents a linkage disequilibrium
        prototype for two recombination events prior to selfing or random
        intermating.

        Parameters
        ----------
        r : numpy.ndarray
            Recombination probability matrix.
        s : int, inf
            Selfing generation number to derive gametes from.
            Example:
                s = 0    ->  Derive gametes from F1
                s = 1    ->  Derive gametes from F2
                s = 2    ->  Derive gametes from F3
                ...
                k = inf  ->  Derive gametes from SSD
            Remark:
                's' takes priority over 't'. If s > 0, 't' calculations are not
                performed.
        t : int, inf
            Random intermating generation number to derive gametes from.
            Example:
                t = 0    ->  Derive gametes from F1
                t = 1    ->  Derive gametes from (F1 x F1)
                t = 2    ->  Derive gametes from (F1 x F1) x (F1 x F1)
                ...
                t = inf  ->  Derive gametes from unlimited random intermating
            Remark:
                's' takes priority over 't'. If s > 0, 't' calculations are not
                performed.

        Returns
        -------
        D_2 : numpy.ndarray
            A D_2 matrix.
        """
        if s == 0 and t == 0:
            return (1.0 - 2.0*r)**2
        elif s > 0:
            rk = Cross.rk(r, s+1)
            four_r = 4.0*r
            return ( 1.0 - four_r + (four_r*rk) )
        elif t > 0:
            return ((1.0 - 2.0*r)**2) * ((1.0 - r)**t)
        else:
            raise ValueError("s and t must be >= 0.")

    @staticmethod
    def meiosis_mat(sel, geno, genetic_map, n, s = 0, seed = None):
        """
        Generate gametes from individuals. Only works with diploid data.

        Developer notes:
        The purpose of this function is to provide a "low level" (if that is
        possible in Python) function to generate gametes from individuals'
        genotype information and recombination probabilities. Thus, it is
        limited in its funtion arguments (e.g. number of gametes generated
        cannot be specified: one must alter 'sel' to increase gametes from a
        single source)

        This function is to be thought of treating individuals separately
        (i.e. knowledge about other selected individuals is not considered).
        Thus, the 't' term indicating random mating is not included because it
        requires information on other individuals.

        The 's' term indicating number of selfing generations is included
        because its calculation is a modification of recombination
        probabilities. These modified recombination probabilities can be
        applied to *single* individuals.

        Application of random mating here is undesireable because there are
        a near unlimited number of ways to apply random mating. First,
        one must consider population size. Next, how are gametes ordered in
        the final gamete array? Writing this would only serve to limit users
        of this function.

        It is the purpose of other functions to apply random mating, etc.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices for gamete generation of shape (k,).
            Where:
                'k' is the number of selected individuals.
            These indices determine from where and in what order gametes
            originate.
        geno : numpy.ndarray
            A int8 binary genotype matrix of shape (m, n, p).
            Where:
                'm' is the number of chromosome phases (only 2 will work).
                'n' is the number of individuals.
                'p' is the number of markers.
        genetic_map : GeneticMap
            A genetic map.
        n : numpy.ndarray, int
            Number of gametes to produce from a single cross. Gametes from each
            individual are placed right next to each other in the returned
            array.
        s : int, inf
            Alter recombination probability to what would be observed in
            selfing generation 's' in a single seed descent scenario.
            Example:
                s = 0    ->  Derive gametes from F1
                s = 1    ->  Derive gametes from F2
                s = 2    ->  Derive gametes from F3
                ...
                k = inf  ->  Derive gametes from SSD
        seed : numpy.int32
            A random number generator seed.

        Returns
        -------
        gamete : numpy.ndarray
            A int8 binary gamete matrix of shape (k, p).
            Where:
        """
        # if there is a RNG seed, seed the RNG
        pybropt.util.cond_seed_rng(seed)

        # make aliases for variables
        gstix = genetic_map.chr_grp_stix
        gspix = genetic_map.chr_grp_spix
        pstix = genetic_map.xo_prob_stix
        pspix = genetic_map.xo_prob_spix
        prob = genetic_map.xo_prob

        # modify recombination probabilities if needed
        if s > 0:
            prob = Cross.rk(prob, s+1)

        # replicate parent selections by the number of gametes desired.
        fsel = numpy.repeat(sel, n)

        # make an empty matrix to contain gametes
        gamete = numpy.empty(
            (len(fsel), len(genetic_map.chr_start)), # num fsel x len map
            dtype = 'int8'
        )

        # generate random numbers to determine crossover points
        rnd = numpy.random.uniform(
            0,  # minimum of 0 (inclusive)
            1,  # maximum of 1 (exclusive)
            (len(fsel), len(prob))   # num fsel x num potential crossover points
        )

        # make a matrix of random numbers indicating the starting phase
        sphase = numpy.random.binomial(
            1,      # choose phase index 0 or 1
            0.5,    # 50% prob of 0 or 1
            (len(fsel), len(genetic_map.chr_grp_len))    # num fsel x num chr
        )

        # for each gamete (i) from a female parent (s)
        for i,s in enumerate(fsel):
            # for each linkage group (j)
            for j,(gst,gsp,pst,psp) in enumerate(zip(gstix,gspix,pstix,pspix)):
                # determine where crossovers occur (add 1 for exclusive index)
                xo = numpy.flatnonzero(rnd[i,pst:psp] <= prob[pst:psp]) + 1

                # get random starting phase
                phase = sphase[i,j]

                # start point offset
                spt = 0

                # for each crossover position (end point = ept)
                for ept in xo:
                    # fill gamete genotype
                    gamete[i,gst+spt:gst+ept] = geno[phase,s,gst+spt:gst+ept]

                    # alternate phase
                    phase = 1 - phase

                    # move to the next copying segment
                    spt = ept

                # finally, copy last remaining chromosome segment.
                gamete[i,gst+spt:gsp] = geno[phase,s,gst+spt:gsp]

        return gamete

    @staticmethod
    def mate_mat(sel, geno, genetic_map, n, seed = None):
        """
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,)
            Where:
                'k' is the number of selected individuals.
            Even indices are female parents.
            Odd indices are male parents.
        n : int
            Number of offspring per cross.

        Returns
        -------
        progeny : numpy.ndarray
            A 3D genotype matrix of shape (2,k/2,p).
        """
        # seed RNG if needed
        if seed is not None:
            numpy.random.seed(seed)

        # get female and male indices
        fsel = sel[0::2]
        msel = sel[1::2]

        # generate gametes
        fgamete = Cross.meiosis_mat(fsel, geno, genetic_map, n)
        mgamete = Cross.meiosis_mat(msel, geno, genetic_map, n)

        # TODO: # OPTIMIZE:
        # generate offspring genotypes by stacking matrices to make 3d matrix
        progeny = numpy.stack([fgamete, mgamete])

        return progeny

    @staticmethod
    def dh_mat(sel, geno, genetic_map, n, s = 0, seed = None):
        """
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,)
            Where:
                'k' is the number of selected individuals.
        n : int
            Number of offspring per cross.

        Returns
        -------
        progeny : numpy.ndarray
        """
        # seed RNG if needed
        if seed is not None:
            numpy.random.seed(seed)

        # generate gametes
        fgamete = Cross.meiosis_mat(sel, geno, genetic_map, n, s)

        # generate offspring genotypes by stacking matrices to make 3d matrix
        progeny = numpy.stack([fgamete, fgamete])

        return progeny

    @staticmethod
    def pop_ld_cross(pRec, mpAB, fpAB, mfreq, ffreq):
        """
        Function from the boneyard to calculate a population x population cross.
        Calculate LD resulting from a cross between two populations.

        Parameters
        ----------
        pRec : numpy.ndarray
            Probability of recombination.
        mpAB : numpy.ndarray
            Male P(AB) matrix.
        fpAB : numpy.ndarray
            Female P(AB) matrix.
        mfreq : numpy.ndarray
            Male allele frequency matrix.
        ffreq : numpy.ndarray
            Female allele frequency matrix.

        Returns
        -------
        pop_ld : numpy.ndarray
        """
        coupling = 0.5 * (1.0 - pRec) * (mpAB + fpAB)
        repulsion = (0.5 * pRec * ((mfreq[:,None] * ffreq) + (ffreq[:,None] * mfreq)))

        pop_ld = coupling + repulsion

        return pop_ld
