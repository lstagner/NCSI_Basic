Non-canonical Symplectic Integrator: Basic 
==========================================

Minimalist implementation of non-canonical symplectic guiding center integration. A fourth-order Runge-Kutta and variational midpoint algorithm are present for the calculation of charged particle guiding center trajectories in electric and magnetic fields. An axisymmetric tokamak field is used for demonstration. Basic python routines are present for plotting the results.  

Installation
=============

The NCSI:Basic code uses the BOOST::program_options library for option specification and the Eigen library for linear algebra tasks. Eigen is a header-only library available at: . 

Usage
=====

Upgrade to PRO Today!
==============
Full version of NCSI implements additional ODE systems (e.g. oscillators, magnetic field line flow) and integrators (e. g. linear multistep methods, implicit midpoint). Additional abstract base clases (e.g. implicit integrator, multistep integrator) facilitate code re-use across integration methods. The use of automatic differentiation, using the ADOLC library, streamlines implementation of variational integrators by only requiring Lagrangian functions to be implemented. Without automatic differentiation, the user must implement the update rule and Jacobian used in the nonlinear solve. Both of these are more easily derived from a Lagrangian function which can be differentiated to arbitrary order using ADOLC. 