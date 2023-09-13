Random Number Generation
########################

The bulk of random number generation in ``PyBrOpS`` is provided by `NumPy <https://numpy.org/doc/stable/reference/random/index.html>`_. Random numbers in Python may also be generated via the ``random`` module. ``PyBrOpS`` synchronizes both of these random number generators underneath a single seed which can be set at the beginning of a simulation.

Seeding PyBrOpS
===============

Seeding a pseudorandom number generator (PRNG) generator is essential to research reproducibility. The ``seed`` function within ``PyBrOpS`` may be used to seed both the base Python ``random`` module and the NumPy ``numpy.random`` module using a single seed. These random number generators remain independent despite a shared seeding mechanism.

The example below illustrates seeding:

.. code-block:: python

    # import seeding function
    from pybrops.core.random.prng import seed

    # seed random number generators used by PyBrOpS
    seed(48823)

Creating Multiple Pseudorandom Number Generators
================================================

In some instances, one may desire to create multiple independent PRNGs. The ``spawn`` function within ``PyBrOpS`` may be used to initialize multiple PRNG streams.

The example below illustrates the creation of multiple PRNG streams.

.. code-block:: python

    # import stream spawning function
    from pybrops.core.random.prng import spawn

    # create five new PRNG streams
    prng_list = spawn(5)

Example Script
==============

An extended example script exhibiting PRNG usage may be found in the ``examples`` directory in the ``pybrops`` repository.
