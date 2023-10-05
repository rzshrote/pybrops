Phenotyping Protocols
#####################

Class Family Overview
=====================

The ``PhenotypingProtocol`` family of classes is used to simulate the phenotyping of a set of individuals across one or more environments. Phenotyping typically induces errors in genotypic values, leading to increased realism in a simulation. 

Summary of Phenotyping Protocol Classes
=======================================

Phenotyping protocol classes are found in the ``pybrops.breed.prot.pt`` module in PyBrOpS. PyBrOpS provides several simple phenotyping protocols, but for more complex phenotyping scenarios, the user may implement his or her own custom phenotyping protocols using the ``PhenotypingProtocol`` interface. Phenotyping protocol classes are summarized in the table below.

.. list-table:: Summary of classes in the ``pybrops.breed.prot.pt`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``PhenotypingProtocol``
      - Abstract
      - Interface for all phenotyping protocol child classes.
    * - ``TruePhenotyping``
      - Concrete
      - Class representing phenotyping done without error.
    * - ``G_E_Phenotyping``
      - Concrete
      - Class representing phenotyping with environmental effects and **no** GxE interactions.

Phenotyping Protocol Properties
===============================

Phenotyping protocols have a couple of properties that can be grouped into two groups: genomic model properties and error variance properties. Tables summarizing these properties are in the following two subsections.

Genomic model properties
------------------------

All phenotyping protocols rely on a true genomic model to simulate phenotypes. This true genomic prediction model is accessed via the ``gpmod`` property.

.. list-table:: Summary of ``PhenotypingProtocol`` genomic model properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``gpmod``
      - The true genomic prediction model.

Error variance properties
-------------------------

All phenotyping protocols also have a basic error variance property. This error variance property quantifies the amount of pure error to be added to phenotypes. This pure error variance may be accessed via the ``var_err`` property. Other phenotyping protocols may add additional sources of error, but this is implementation dependent and no required by the ``PhenotypingProtocol`` interface.

.. list-table:: Summary of ``PhenotypingProtocol`` error variance properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``var_err``
      - The pure error variance for each trait.

Loading Class Modules
=====================

Phenotyping protocols may be imported as demonstrated in the code below.

.. code-block:: python

    # import the PhenotypingProtocol class (an abstract interface class)
    from pybrops.breed.prot.pt.PhenotypingProtocol import PhenotypingProtocol

    # import the TruePhenotyping class (a concrete implemented class)
    from pybrops.breed.prot.pt.TruePhenotyping import TruePhenotyping

    # import the G_E_Phenotyping class (a concrete implemented class)
    from pybrops.breed.prot.pt.G_E_Phenotyping import G_E_Phenotyping

Creating Phenotyping Protocols
==============================

Phenotyping protocols may be created using their constructors. Constructor definitions are implmentation dependent, but the code below demonstrates the creation of a ``G_E_Phenotyping`` object using its constructor.

.. code-block:: python

    #
    # Creating a true genomic model
    #

    # model parameters
    nfixed = 1      # number of fixed effects
    ntrait = 2      # number of traits
    nmisc = 0       # number of miscellaneous random effects
    nadditive = 50  # number of additive marker effects

    # create dummy values
    beta = numpy.random.random((nfixed,ntrait))
    u_misc = numpy.random.random((nmisc,ntrait))
    u_a = numpy.random.random((nadditive,ntrait))
    trait = numpy.array(["Trait"+str(i+1).zfill(2) for i in range(ntrait)], dtype = object)

    # create additive linear genomic model
    algmod = DenseAdditiveLinearGenomicModel(
        beta = beta,
        u_misc = u_misc,
        u_a = u_a,
        trait = trait,
        model_name = "example",
        params = None
    )

    #
    # Creating a phenotyping object
    #

    # phenotyping parameters
    nenv = 3    # number of environments
    nrep = 2    # number of replicates within each environment

    # construct phenotyping object
    ptprot = G_E_Phenotyping(
        gpmod = algmod,
        nenv = nenv,
        nrep = nrep
    )

Setting Trait Heritabilities
============================

After phenotyping protocol object construction, it is generally required to set the heritability for each trait. This may be accomplished using the ``set_h2`` or ``set_H2`` methods which set narrow-sense and broad-sense heritabilities, respectively. A founder population is needed to use these methods. ``set_h2`` and ``set_H2`` method usage is demonstrated below.

.. code-block:: python

    #
    # Construct random genomes
    #

    # shape parameters for random genomes
    ntaxa = 100
    nvrnt = nadditive
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
        ploidy = nphase
    )

    # set the heritabilities from 
    heritability = numpy.array([0.4, 0.7])

    # set the narrow sense heritability
    ptprot.set_h2(
        h2 = heritability,
        pgmat = pgmat
    )

    # set the broad sense heritability
    ptprot.set_H2(
        H2 = heritability,
        pgmat = pgmat
    )

Phenotyping Individuals
=======================

Individuals may be phenotyped using the ``phenotype`` method, which accepts a set of individual genomes and returns a phenotype ``pandas.DataFrame``. The use of this method is demonstrated below.

.. code-block:: python

    # phenotype individuals
    pheno_df = ptprot.phenotype(pgmat)
