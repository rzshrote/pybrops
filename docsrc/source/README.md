# Notes on package file structure for proper documentation

## `__init__.py` files
`__init__.py` files should not import objects with the same name as the module.
Instead, `__init__.py` should simply import the module.

Do this:
```
from pybrops.core.mat.Matrix import Matrix
```

Do NOT do this:
```
from . import Matrix
from .Matrix import *
```

Where `Matrix.py` contains the class `Matrix`.

# Notes on files

## index.rst
DO NOT DELETE THIS FILE! This file is used as the root file to build the entire
documentation.

## api.rst
DO NOT DELETE THIS FILE! It contains the all-important `.. autosummary::`
directive with `:recursive:` option, without which API documentation wouldn't
get extracted from docstrings by the `sphinx.ext.autosummary` engine. It is
hidden (not declared in any toctree) to remove an unnecessary intermediate
page; index.rst instead points directly to the package page.
DO NOT REMOVE THIS FILE!

# Notes on documentation generation bugs.
There was a weird bug associated with Sphinx documentation generation where
several web pages would be created multiple times. Objects imported from third-
party libraries would be documented in the website.

## Status of the bug:
* Currently squashed and not a problem, but may come back.

## Cause of the bug:
* The cause was unknown. I couldn't effectively isolate it.

## Potential causes of the bug:
* Sphinx might have searched the `PyBrOpS/dist/` directory, found the built
  package and generated documentation, erroneously.
    * Highly suspect since there were 3 duplicated documentation files generated
      (`source`, `dist`, `build`).
* Sphinx might have searched the `PyBrOpS/build/` directory, found the files
  built for packaging and generated documentation for these, erroneously.
    * Highly suspect since there were 3 duplicated documentation files generated
      (`source`, `dist`, `build`).
* `PyBrOpS/doc/source/conf.py` might have had order dependent settings.
    * Low probability of cause. Files were almost identical.
* The `.egg-info` directory might have caused the problem.
    * Low probability of cause since the packaging tutorial repo had an
      `.egg-info` directory, but still generated documentation correctly.
* The `doxygen.config` file might have caused troubles since it had the
  `.config` suffix.
    * Unknown probability. Lots of settings in this file.
* Hidden directories.
    * Low probability of cause.

Things that are definitely known NOT to cause the bug:
* The source code. There are warnings generated about malformated docstrings,
  but in isolated directories, they cause no duplication problems.
* Location of the `PyBrOpS` directory on the user file system. I thought I might
  have modified a `PATH`-type setting that was dependent on the `PyBrOpS`
  directory location, but this proved not to be an issue.
