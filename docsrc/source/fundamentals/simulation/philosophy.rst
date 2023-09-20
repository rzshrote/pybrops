PyBrOpS Simulation Philosophy
#############################

PyBrOpS adopts a script-based strategy to defining and breeding program simulations, allowing for maximum flexibility. The representation of a breeding program is accomplished using three components: Fundamental Data Types, Breeding Protocols, and Breeding Operators. These three components are used to define a breeding simulation's structure.

Fundamental Data Types
======================

PyBrOpS is equipped with a rich array of fundamental data types that play a crucial role in conducting breeding simulations. Fundamental Data Types are PyBrOpS's most basic building block and represent fundamental units of information created in a breeding program (e.g., genotypic data, genomic prediction models). Fundamental Data Types are utilized and composed by Breeding Protocols within the PyBrOpS framework. Several key data types offered by PyBrOpS are listed below:

* **Genetic Maps**: PyBrOpS provides object types to represent genetic maps, which are essential in simulating recombination and estimating marker covariances.
* **Genetic Map Functions**: These object types are designed to work in tandem with genetic maps, enabling calculations of recombination probability.
* **Genotype Matrices**: Genotype matrices are pivotal for storing genetic information on individuals in a breeding population.
* **Breeding Value Matrices**: Breeding value matrices contain estimates of an individual's genetic merit, which is vital for selecting parents with desirable traits.
* **Coancestry Matrices**: Coancestry matrices quantify the genetic relatedness between individuals in a population. These matrices are indispensable for managing genetic diversity and avoiding inbreeding in breeding programs.
* **Genomic Models**: Genomic models may be used to simulate the genetic control of one or more traits of interest. PyBrOpS supports various genomic modeling approaches, representing different modes of genetic control.
* **Progeny Variance Matrices**: Progeny variance matrices help in assessing the genetic variance and covariance among offspring. These matrices may be useful in some selection techniques which require estimates of genetic variance among progenies.
* **Optimization Problems**: PyBrOpS defines a set of optimization problem data types to facilitate the formulation and optimization of breeding program optimization problems. These problems can involve selecting optimal parents, designing mating schemes, or may be defined by the user.
* **Optimization Solutions**: After solving optimization problems, PyBrOpS records the solution(s) in this data type. Optimization solutions house data which may be used by a breeder to decide an optimal breeding decision or strategy.
* **Optimization Algorithms**: PyBrOpS offers a range of single- and multi-objective optimization algorithms designed to tackle specific optimization problems which may be faced in a breeding program.

Breeding Protocols
==================

Breeding Protocols in PyBrOpS represent small processes in a breeding program which take one or more Fundamental Data Types, perform some specialized task, and return a single Fundamental Data Type as an output.

PyBrOpS defines five Protocols: Breeding Value Estimation, Genotyping, Mating, Phenotyping, and Selection. The Breeding Value Estimation Protocol provides routines for estimating breeding values from phenotypic and genotypic data, potentially using external packages like lme4. The Genotyping Protocol simulates the genotyping of individuals. The Mating Protocol handles mating simulations for various cross configurations. The Phenotyping Protocol is responsible for simulating phenotypes resulting from field trials. Lastly, the Selection Protocol offers methods for selecting individuals within the breeding program.

Breeding Operators
==================

Operators in PyBrOpS represent processes within a breeding program. These processes take the entire state of a breeding program, transform it, and return a modified state. Operators are designed to be versatile and encompass the diverse nature of breeding programs. They can be assembled into a Universal Breeding Algorithm to simulate a breeding program effectively. On the other hand, Protocols are processes that take multiple inputs and produce a single output. They are intended to construct the processes defined by Operators.

PyBrOpS defines five Operators: Initialization, Parental Selection, Mating, Evaluation, and Survivor Selection. Initialization provides a starting point for breeding program simulations. It can involve loading genotypic and phenotypic data or defining burn-in generations. The Evaluation Operator handles tasks related to measuring the desirability of individuals in the program, including genotyping, phenotyping, genomic prediction, and breeding value estimation. Parental Selection utilizes information from the Evaluation Operator to select parents and determine cross configurations for generating offspring. The Reproduction Operator mates parents and generates progeny, mimicking meiosis and mutation. Finally, the Survivor Selection Operator uses Evaluation information to choose individuals for the next breeding cycle, merging subpopulations if needed.
