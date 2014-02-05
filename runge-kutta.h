/*!
*******************************************************************************
* \file runge-kutta.h
*
* \brief Header for general Runge-Kutta integrator.
*
* \date Feb 2014
* \author C. Leland Ellison
* \version 0.1
*******************************************************************************
*/
#ifndef RUNGEKUTTA_H_
#define RUNGEKUTTA_H_

#include "integrator.h"  // Abstract base class
#include <Eigen/Dense>  // Linear algebra operations 
#include "guiding_center.h"  // ODE being modeled
#include <iostream>  // Outputting messages

class RungeKutta : public Integrator {
 public:
  // General RK method constructor
  RungeKutta(const double kdt, const GuidingCenter &kGuidingCenter,  
	     const Eigen::MatrixXd a_coefficients, 
	     const Eigen::VectorXd b_coefficients, 
	     const Eigen::VectorXd c_coefficients);
  // Convenience constructor which sets coefficients for order 2, 4
  RungeKutta(const double kdt, const GuidingCenter &kGuidingCenter, 
	     const int kOrder);
  ~RungeKutta();
  // Map from (t_k, x_k) -> (t_{k+1}, x_{k+1})
  int Step(double &t, Eigen::VectorXd &x);  
  //// Accessors
  //! Access a_ coefficients matrix
  Eigen::MatrixXd a_coefficients() const { return a_; }
  //! Access b_ coefficients vector
  Eigen::VectorXd b_coefficients() const { return b_; }
  //! Access c_ coefficients vector
  Eigen::VectorXd c_coefficients() const { return c_; }
 private:
  Eigen::MatrixXd a_;  //!< Specifies for a_ij coefficients of the RK method
  Eigen::VectorXd b_;  //!< Specifies b_i coefficients of the RK method
  Eigen::VectorXd c_;  //!< Specifies c_i coefficients of the RK method
  Eigen::MatrixXd k_;  //!< Temporary space for k_i evaluations
  Eigen::VectorXd xtemp_; //!< Temporary space used at internal stages
  Eigen::VectorXd ftemp_; //!< Temporary space for function evaluations
};

#endif  // RUNGEKUTTA_H_
