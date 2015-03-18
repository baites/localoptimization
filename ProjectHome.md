# LO #

_Victor E. Bazterra_

_Departament of Physics_

_University of Illinois at Chicago (2007)_


# Description #

This project is a small C++ library for local optimization. It is a general interface that implements local optimization of functions. It implements linear search algorithm with or without derivatives and Powell method and Gradient Conjugate. Moreover it has several methods for evaluating numerical derivatives.

The code is documented by using Doxygen and comes with a set of examples.

# Installation #

It is necessary to define LOROOT that should be ponting to the root directory where the library is install.

After that just type 'make' in the LOROOT directory.

Example:
~> tar xvfz lo-x.x.tar.gz
~> cd lo-x.x
lo-x.x> export LOROOT=$PWD (for bash)
        setenv LOROOT PWD (for tcsh)
lo-x.x> make

```