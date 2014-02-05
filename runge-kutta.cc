/*!
*******************************************************************************
* \file runge-kutta.cc
*
* \brief Implementation of general Runge-Kutta. 
*
* \date October 2013
* \author C. Leland Ellison
* \version 0.1
*
*******************************************************************************
*/

#include "runge-kutta.h"

 /*!
  * 
  * \brief Constructor for general runge-kutta methods. Specify the coefficients in eigen vectors/matrices. 
  *
  * @param[in] kdt Numerical step size
  * @param[in] kGuidingCenter ODE system being modeled
  * @param[in] a_coeffs Matrix with a_{ij} coefficients in Butcher tableau
  * @param[in] b_coeffs Vector with b_i coefficients in Butcher tableau
  * @param[in] c_coeffs Vector with c_i coefficients in Butcher tableau
  * 
  */
RungeKutta::RungeKutta(const double kdt, const GuidingCenter &kGuidingCenter, 
		       const Eigen::MatrixXd a_coefficients, 
		       const Eigen::VectorXd b_coefficients, 
		       const Eigen::VectorXd c_coefficients)
  :   Integrator(kdt, kGuidingCenter), a_(a_coefficients), b_(b_coefficients),
      c_(c_coefficients) { 
  k_.resize(kDimen_, a_.cols());
  xtemp_.resize(kDimen_);
  ftemp_.resize(kDimen_);
}

 /*!
  * 
  * \brief Constructor which implements common methods for convenience. Specify an order with an integer, and this constructor will set the corresponding coefficients.
  *
  * @param[in] kdt Numerical step size
  * @param[in] model ODE system being solved
  * 
  */
RungeKutta::RungeKutta(const double kdt, const GuidingCenter &kGuidingCenter, 
		       const int kOrder) : Integrator(kdt, kGuidingCenter) {
  //// Set size of matrices
  a_ = Eigen::MatrixXd::Zero(kOrder,kOrder);
  b_ = Eigen::MatrixXd::Zero(kOrder, 1);
  c_ = Eigen::MatrixXd::Zero(kOrder, 1);

  k_.resize(kDimen_, a_.cols());
  xtemp_.resize(kDimen_);
  ftemp_.resize(kDimen_);

  //// Setting of the various coefficients
  // RK2
  if (kOrder == 2){
    a_(1,0) = 0.5;
    b_(1) = 1.0;
    c_(1) = 0.5;
  }
  // Rk4
  else if (kOrder == 4){
    a_(1,0) = 0.5;
    a_(2,1) = 0.5;
    a_(3,2) = 1.0;
    b_(0) = 1./6.;
    b_(1) = 1./3.;
    b_(2) = 1./3.;
    b_(3) = 1./6.;
    c_(1) = 0.5;
    c_(2) = 0.5;
    c_(3) = 1.0;
  }
  else {
    std::cout << "Unknown/Unimplemented RK order" << std::endl;
  }
}


/*!
 * \brief Explicit advance x(t) forward in time by step size kdt. RK:
 *        x_{n+1} = x_n + sum_{i=1}^s b_i k_i
 *            k_i = h f(t_n + c_i h, y_n + sum_{j=1}^{i-1} a_ij k_j 
 * 
 * @param[in, out] t Simulation time. Advanced by kdt_
 * @param[in, out] x Position. At in: x(t=t_k) At out: x(t=t_{k+1})
 */
int RungeKutta::Step(double &t, Eigen::VectorXd &x) {
  
  // Assign k_i values
  for (int i=0; i<k_.cols(); ++i){
    // Get temporary position
    xtemp_ = x;
    for (int j=0; j<i; ++j){
      xtemp_ += a_(i,j)*k_.col(j);
    }
    // Evaluate the rhs at x_temp to determine k
    guiding_center_.VectorField(t + c_(i)*kdt_, xtemp_, ftemp_);
    k_.col(i) = kdt_*ftemp_;
  }

  // Update x
  x = x + k_*b_;
  // Update t
  t = t + kdt_;

  return 0;
}
