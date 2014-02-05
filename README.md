Non-canonical Symplectic Integrator: Basic 
==========================================

Minimalist implementation of guiding center trajectory calculation. A fourth-order Runge-Kutta and noncanonical symplectic algorithm are present for the calculation of charged particle guiding center trajectories in electric and magnetic fields. An axisymmetric tokamak field is used for demonstration. 

Installation
------------

The NCSI:Basic code uses the [BOOST::program_options](http://boost.org/) library for option specification and the [Eigen](http://eigen.tuxfamily.org) library for linear algebra tasks. Eigen is a header-only library available at http://eigen.tuxfamily.org/ . The compiler needs to know where to find the header files, so placing the Eigen package somewhere like /usr/local/include is a good idea. 

Usage
-----

After compiling and linking the executable "driver", see a summary of runtime options using 
driver --help

The options may be specified on the command line with --<option> or -<O>. Alternatively, one may specify input options using an input file. See sample_input.cfg for an example. Units are discussed in the documentation.

Upgrade to PRO Today!
---------------------
Full version of NCSI includes:
* Additional ODE systems (oscillators, magnetic field line flow)
* Additional integrators (linear multistep methods, implicit midpoint)
* Improved class hierarchy (ImplicitIntegrator, MultistepIntegrator)
* Automatic differentiation (with [ADOLC](https://projects.coin-or.org/ADOL-C))
* Field writing routines (with Python and Sympy)
* Unittests (using [gTest](http://code.google.com/p/googletest/))
