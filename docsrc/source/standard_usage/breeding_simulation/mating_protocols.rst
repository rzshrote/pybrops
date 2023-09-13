Mating Protocols
################

Class Family Overview
=====================

Summary of Mating Protocol Classes
==================================

.. list-table:: Summary of classes in the ``pybrops.breed.prot.mate`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``MatingProtocol``
      - Abstract
      - Interface for all mating protocol child classes.
    * - ``SelfCross``
      - Concrete
      - Class representing self-fertilization crosses.
    * - ``TwoWayCross``
      - Concrete
      - Class representing two-way crosses.
    * - ``TwoWayDHCross``
      - Concrete
      - Class representing two-way crosses followed by doubled haploidization.
    * - ``ThreeWayCross``
      - Concrete
      - Class representing three-way crosses.
    * - ``ThreeWayDHCross``
      - Concrete
      - Class representing three-way crosses followed by doubled haploidization.
    * - ``FourWayCross``
      - Concrete
      - Class representing four-way crosses.
    * - ``FourWayDHCross``
      - Concrete
      - Class representing four-way crosses followed by doubled haploidization.

Loading Class Modules
=====================

.. code-block:: python

    # abstract interface classes
    from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol

    # concrete implementation classes
    from pybrops.breed.prot.mate.SelfCross import SelfCross
    from pybrops.breed.prot.mate.TwoWayCross import TwoWayCross
    from pybrops.breed.prot.mate.TwoWayDHCross import TwoWayDHCross
    from pybrops.breed.prot.mate.ThreeWayCross import ThreeWayCross
    from pybrops.breed.prot.mate.ThreeWayDHCross import ThreeWayDHCross
    from pybrops.breed.prot.mate.FourWayCross import FourWayCross
    from pybrops.breed.prot.mate.FourWayDHCross import FourWayDHCross

Creating Mating Protocols
=========================

.. code-block:: python

    # initialize mating operators with line and family counters at zero
    mate_self = SelfCross()
    mate_2w   = TwoWayCross()
    mate_2wdh = TwoWayDHCross()
    mate_3w   = ThreeWayCross()
    mate_3wdh = ThreeWayDHCross()
    mate_4w   = FourWayCross()
    mate_4wdh = FourWayDHCross()

    # optionally initialize counters to non-zero values
    mate_self = SelfCross(progeny_counter = 1000, family_counter = 100)
    mate_2w   = TwoWayCross(progeny_counter = 1000, family_counter = 100)
    mate_2wdh = TwoWayDHCross(progeny_counter = 1000, family_counter = 100)
    mate_3w   = ThreeWayCross(progeny_counter = 1000, family_counter = 100)
    mate_3wdh = ThreeWayDHCross(progeny_counter = 1000, family_counter = 100)
    mate_4w   = FourWayCross(progeny_counter = 1000, family_counter = 100)
    mate_4wdh = FourWayDHCross(progeny_counter = 1000, family_counter = 100)

    # make mating operator for use in rest of example
    mtprot = TwoWayDHCross()

Mating Protocol Properties
==========================

Mating protocol general properties
----------------------------------

.. list-table:: Summary of ``MatingProtocol`` general properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``nparent``
      - The number of parents required by the mating protocol.

Simulating Mating and Progeny Generation
========================================

.. code-block:: python

    #
    # Construct random genomes
    #

    # shape parameters for random genomes
    ntaxa = 100
    nvrnt = 50
    ngroup = 20
    nchrom = 10
    nphase = 2

    # create random genotypes
    mat = numpy.random.randint(0, 2, size = (nphase,ntaxa,nvrnt)).astype("int8")

    # create taxa names
    taxa = numpy.array(["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], dtype = object)

    # create taxa groups
    taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
    taxa_grp.sort()

    # create marker variant chromsome assignments
    vrnt_chrgrp = numpy.random.randint(1, nchrom+1, nvrnt)
    vrnt_chrgrp.sort()

    # create marker physical positions
    vrnt_phypos = numpy.random.choice(1000000, size = nvrnt, replace = False)
    vrnt_phypos.sort()

    # create marker genetic positions
    vrnt_genpos = numpy.random.uniform(0, 2, nvrnt)
    vrnt_genpos.sort()

    # create random crossover probabilities
    vrnt_xoprob = numpy.random.uniform(0, 0.5, nvrnt)
    step = nvrnt // nchrom
    vrnt_xoprob[0::step] = 0.5

    # create marker variant names
    vrnt_name = numpy.array(["SNP"+str(i+1).zfill(4) for i in range(nvrnt)], dtype = object)

    # create a phased genotype matrix from scratch using NumPy arrays
    pgmat = DensePhasedGenotypeMatrix(
        mat = mat,
        taxa = taxa,
        taxa_grp = taxa_grp, 
        vrnt_chrgrp = vrnt_chrgrp,
        vrnt_phypos = vrnt_phypos, 
        vrnt_name = vrnt_name, 
        vrnt_genpos = vrnt_genpos,
        vrnt_xoprob = vrnt_xoprob,
        ploidy = nphase
    )

    # create a cross configuation matrix of shape (10,2)
    xconfig = numpy.random.choice(ntaxa, (10,2), replace = False)

    # create progeny
    progeny = mtprot.mate(
        pgmat = pgmat,
        xconfig = xconfig,
        nmating = 1,
        nprogeny = 10,
        nself = 2
    )
