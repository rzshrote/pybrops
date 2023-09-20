PyBrOpS Software Architecture Philosophy
########################################

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


How PyBrOpS uses SOLID Design Principles
****************************************

To adhere to SOLID principles, PyBrOpS makes heavy use of `software interfaces <https://en.wikipedia.org/wiki/Interface_(computing)#In_object-oriented_languages>`_ and abstract classes using Python's `Abstract Base Class module <https://docs.python.org/3/library/abc.html>`_. Since Python uses duck typing, interfaces are manually enforced. A user seeking to create a new class implementing one of PyBrOpS's interfaces must manually ensure that the derived class adheres to the interface's function and type requirements.

Use of the Single Responsiblity Principle (SRP)
===============================================

Every family of classes in PyBrOpS has a singular, focused purpose. Commonly, these classes represent fundamental entities or data structures in breeding programs.

Use of the Open-Closed Principle (OCP)
======================================

Almost all classes in PyBrOpS implement at least one interface, and most functions and class methods depend on interface abstractions. These two design decisions allow for the adhereance to OCP. If a new feature is desired, one can simply implement a new class using one of the interfaces instead of modifying a current implementation.

Use of the Liskov Substitutuion Principle (LSP)
===============================================

All classes which implement an interface do so such that class method signatures are identical two those specified by the interface. This allows for classes within the same family to be substituted with each other to provide alternative functionality.

Use of the Interface Segregation Principle (ISP)
================================================

Interfaces within PyBrOpS have singular and focused purposes. This makes them lean and prevents classes implementing them from having unnecessary features.

Use of the Dependency Inversion Principle (DIP)
===============================================

As mentioned in the OCP section, almost all functions and class methods depend only on interface abstractions. This decouples concrete classes from each other, improving modularity.
