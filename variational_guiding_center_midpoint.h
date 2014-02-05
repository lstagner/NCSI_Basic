/*!
*******************************************************************************
* \file variational_guiding_center_midpoint.h
*
* \brief Header for integrator resulting from a midpoint discretization of the guiding center Lagrangian. 
*
* \date Feb 2014
* \author C. Leland Ellison
*******************************************************************************
*/

#ifndef VARIATIONAL_GUIDING_CENTER_MIDPOINT_H_
#define VARIATIONAL_GUIDING_CENTER_MIDPOINT_H_

#include "guiding_center.h"
#include "em_fields.h"
#include "integrator.h"
#include <stdlib.h>
#include <Eigen/Dense>
#include "runge-kutta.h" // Used in InitialStep

class VariationalGuidingCenterMidpoint : public Integrator{
public:
  VariationalGuidingCenterMidpoint(const double kdt, 
				   const GuidingCenter &kGuidingCenter,
				   const double kNewtonTolerance = 1e-12,
				   const double kMaxIterations = 15);
  // Method to advance (t_k, x_k) -> (t_k+1, x_k+1)
  int Step(double &t, Eigen::VectorXd &x);
  // Used when running multiple initial conditions to reset state
  void Reset(){ needs_initialization_ = true; }
protected:
  // Pushes new point onto end of history. Drops oldest point.
  int StoreHistory(const Eigen::VectorXd &kx);
  // Method to perform time advance while more initial conditions are needed
  int InitialStep(double &t, Eigen::VectorXd &x) const;
  // Function which should evaluate to zero when algorithm is satisfied
  int UpdateRule(const double kt, const Eigen::VectorXd &kx, 
		 Eigen::VectorXd &error) const;
  // Method to guess a next point (t,x). Uses forward euler here.
  int NewtonGuess(double &t, Eigen::VectorXd &x) const;
  // Calculates jacobian matrix used in nonlinear solve by finite difference
  int Jacobian(const double kt, const Eigen::VectorXd &kx, 
	       Eigen::MatrixXd &jacobian) const;
  const double kMu_;  //!< Magnetic moment of particle
  const double kNewtonTolerance_;  //!< Error threshold for nonlinear solve
  const double kMaxIterations_;  //!< Maximum allowable Newton iterations
  bool needs_initialization_;  //!< Flag for needing initial conditions
  Eigen::MatrixXd x_history_;  //!< Stores previous positions
  EMFields *em_fields_;
};

#endif // VARIATIONAL_GUIDING_CENTER_MIDPOINT_H_
