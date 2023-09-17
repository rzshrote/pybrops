Population Genetics Data Types
##############################

The ``pybrops.popgen`` module houses modules, classes, and functions related to population genetics, as the name suggests. Population genetics data types serve as the foundations for simulations.

Genetic Maps
************

Genetic maps are used to estimate marker-wise genetic map positions in mating simulations. PyBrOpS provides functionality for the reading and writing of genetic maps and the interpolation of genetic map positions via ``GeneticMap`` classes. The ability to read a genetic map from a file allows for real genetic recombination data to be used in a breeding simulation, adding realism.

Genetic Map Functions
*********************

The ``GeneticMapFunction`` family of classes allows for the representation of genetic map functions within PyBrOpS. Genetic map functions are used to convert genetic map positions into recombination probabilities, which are useful in the simulation of genetic recombination and the estimation of expected progeny variance. Popular genetic map functions like the Haldane and Kosambi map functions are provided by PyBrOpS out-of-box.

Genotype Matrices
*****************

Perhaps the most important family of classes in PyBrOpS is the ``GenotypeMatrix`` object family. The purpose of this family of classes is to store and represent genotypic data. ``GenotypeMatrix`` classes can be used in the estimation of genomic prediction models, the estimation of breeding values, the calculation of genomic relationship matrices, and to make selection decisions. Within PyBrOpS it is possible to read genotypic data from VCF files to create ``GenotypeMatrix`` objects, allowing for real-world data to be used in breeding program simulations.

Phased Genotype Matrices
************************

The ``PhasedGenotypeMatrix`` family of classes is a derivative of the ``GenotypeMatrix`` family. The purpose of this family of classes is to store and represent genomic data for individuals. Genomic data is different from genotypic data in that the former includes allele state and chromosome phase information while the latter only includes allele state information. This family of classes is used to represent the genomes of individuals in a breeding program. Since ``PhasedGenotypeMatrix`` classes are derived from ``GenotypeMatrix`` classes, they have the same uses as a regular ``GenotypeMatrix``. Since they are phased, they can also be used to create simulated progeny. Like with ``GenotypeMatrix`` classes, it is also possible to create ``PhasedGenotypeMatrix`` objects by reading real-world data from VCF files.

Breeding Value Matrices
***********************

The ``BreedingValueMatrix`` family of classes is used to represent breeding values as its name implies. ``BreedingValueMatrix`` objects can be used in the estimation of genomic prediction models and to make selection decisions.

Coancestry Matrices
*******************

The ``CoancestryMatrix`` family of classes is used to represent coancestry relationships between individuals. ``CoancestryMatrix`` objects can be used in the estimation of genomic prediction models and to make selection decisions.
