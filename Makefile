# change Makefile shell to bash
# SHELL := /bin/bash

# pybropt:
# 	python3 -m build
#
# install:
# 	pip3 install -U dist/pybrops-1.0.0.tar.gz

# instructions for creating a virtual environment for local testing
# Line 1: if a virtual environment has not been created, create it
# Line 2: activate the virtual environment for local testing
devel-virtualenv:
	# make virtual environment if needed
	if [ ! -f env/bin/activate ]; then python3 -m venv env; fi

# instructions for installing Python libraries needed for local development
devel-dependencies: devel-virtualenv
	. env/bin/activate && python3 -m pip install numpy pandas scipy matplotlib pymoo cyvcf2 h5py pytest pytest-datadir sphinx-autodoc-typehints pydata-sphinx-theme

# instructions for installing an editable package for local testing
devel-install: devel-virtualenv devel-dependencies
	# activate virtual environment and install
	. env/bin/activate && python3 -m pip install --editable .

# build dependency tools
build-dependencies:
	python3 -m pip install --upgrade pip
	python3 -m pip install --upgrade build
	python3 -m pip install --upgrade twine

# instructions for building the package
build: build-dependencies
	python3 -m build

# instructions for distributing the package
dist: build
	python3 -m twine upload dist/*

# instructions for building the package documentation in html format
doc-html:
	cd docsrc/ && $(MAKE) html
	if [ ! -d docs ]; then mkdir docs; fi
	cp -r docsrc/build/html/* docs/
	cp -r docsrc/build/html/.* docs/

# instructions for building the package documentation in pdf format (requires LaTeX)
doc-pdf:
	cd docsrc/ && $(MAKE) latexpdf

# instructions for cleaning the virtual environment
clean-devel-virtualenv:
	rm -rf env/

clean-dist:
	rm -rf dist/
