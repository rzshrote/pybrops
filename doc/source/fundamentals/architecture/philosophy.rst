PyBrOpS Architecture Philosophy
###############################

PyBrOpS uses `SOLID design principles <http://butunclebob.com/ArticleS.UncleBob.PrinciplesOfOod>`_ to guide its architecture. By using SOLID design principles, PyBrOpS is designed to be easily comprehensible, extentible, and maintainable. Many open source softwares commonly used in the plant and animal breeding community is not developed with a robust set of architecture principles, making them difficult to understand and extend. The use of SOLID architectural design principles in PyBrOpS aims to remedy this, maximizing the usefulness of the software to the user.

Briefly, SOLID design principles are summarized as:

.. list-table:: SOLID Design Principles
   :widths: 25 50
   :header-rows: 1

   * - Design Principle
     - Description
   * - Single Responsiblity Principle (SRP)
     - A software module should be responsible for only one task and have only one reason to change.
   * - Open-Closed Principle (OCP)
     - A software module should be open for exension but closed for modification.
   * - Liskov Substitutuion Principle (LSP)
     - Derived modules should be substitutable for their base modules.
   * - Interface Segregation Principle (ISP)
     - Software interfaces should have a narrow focus so as to prevent unnecessary dependencies.
   * - Dependency Inversion Principle (DIP)
     - Software modules should depend on abstractions, not concretions.


How PyBrOpS uses SOLID Design Principle
***************************************

PyBrOpS makes heavy use of abstract classes which serve as interfaces for software functionality.
