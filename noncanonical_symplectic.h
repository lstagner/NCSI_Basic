/*!
*******************************************************************************
* \file noncanonical_symplectic.h
*
* \brief Header for integrator resulting from a midpoint discretization of the guiding center Lagrangian. 
*
* \date Feb 2014
* \author C. Leland Ellison
*******************************************************************************
*/

#ifndef NONCANONICAL_SYMPLECTIC_H_
#define NONCANONICAL_SYMPLECTIC_H_

#include "guiding_center.h"
#include "em_fields.h"
#include "integrator.h"
#include <iostream> 
#include <Eigen/Dense>
#include "runge-kutta4.h" // Used in InitialStep
#include "eigen_types.h"

class NoncanonicalSymplectic : public Integrator{
public:
  NoncanonicalSymplectic(const double kdt, const GuidingCenter &kGuidingCenter,
			 const double kNewtonTolerance,
			 const double kMaxIterations);
  // Method to advance (t_k, x_k) -> (t_k+1, x_k+1)
  int Step(double &t, Vector4 &x);
  // Used when running multiple initial conditions to reset state
  void Reset(){ needs_initialization_ = true; }
protected:
  // Pushes new point onto end of history. Drops oldest point.
  int StoreHistory(const Vector4 &kx);
  // Method to perform time advance while more initial conditions are needed
  int InitialStep(double &t, Vector4 &x) const;
  // Function which should evaluate to zero when algorithm is satisfied
  int UpdateRule(const double kt, const Vector4 &kx, 
		 Vector4 &error) const;
  // Method to guess a next point (t,x). Uses forward euler here.
  int NewtonGuess(double &t, Vector4 &x) const;
  // Calculates jacobian matrix used in nonlinear solve by finite difference
  int Jacobian(const double kt, const Vector4 &kx, 
	       Matrix4 &jacobian) const;
  const double kMu_;  //!< Magnetic moment of particle
  const double kNewtonTolerance_;  //!< Error threshold for nonlinear solve
  const double kMaxIterations_;  //!< Maximum allowable Newton iterations
  bool needs_initialization_;  //!< Flag for needing initial conditions
  Eigen::Matrix<double, 4, 2> x_history_;  //!< Stores previous positions
  EMFields *em_fields_;  //!< Field definitions pulled from GuidingCenter
};

#endif // NONCANONICAL_SYMPLECTIC_H_
