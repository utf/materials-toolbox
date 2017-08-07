# materials-toolbox
A collection of python scripts to help with general comp chem hpc work

## Scripts

List functions here...

## Requirements

These scripts rely on a number of packages to work correctly. For symmetry
functions, `spglib`, `pymatgen`, and `seekpath` are required. For the `mp-get`
script, `ase` is an optional dependancy.

Materials-toolbox uses Pip and setuptools for installation. You *probably*
already have this; if not, your GNU/Linux package manager will be able
to oblige with a package named something like `python-setuptools`. On
Max OSX, the Python distributed with [Homebrew](http://brew.sh)
includes setuptools and Pip.

## Installation
From the directory containing this README:

    pip install --user .

will install a version of the materials-toolbox scripts in your
userspace.

## License
The materials-toolbox code is made available under the MIT License
