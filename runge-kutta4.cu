/*!
*******************************************************************************
* \file runge-kutta.cc
*
* \brief Implementation of Runge-Kutta4 - adapted from general version and for CUDA. 
*
* \date July 2014
* \author C. Leland Ellison
* \version 0.1
*
*******************************************************************************
*/

#include "runge-kutta4.h"

 /*!
  * 
  * \brief Constructor which implements common methods for convenience. Specify an order with an integer, and this constructor will set the corresponding coefficients.
  *
  * @param[in] kdt Numerical step size
  * @param[in] kGuidingCenter ODE system being solved
  * @param[in] kOrder Order of Runge-Kutta method. 2 and 4 are implemented.
  * 
  */
RungeKutta4::RungeKutta4(const double kdt, 
			 const GuidingCenter &kGuidingCenter) : 
  Integrator(kdt, kGuidingCenter), kOrder_(4) {
 
  //// Setting of the various coefficients
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

RungeKutta4::~RungeKutta4(){
}

/*!
 * \brief Explicit advance x(t) forward in time by step size kdt. RK:
 *        x_{n+1} = x_n + sum_{i=1}^s b_i k_i
 *            k_i = h f(t_n + c_i h, y_n + sum_{j=1}^{i-1} a_ij k_j 
 * 
 * @param[in, out] t Simulation time. Advanced by kdt_
 * @param[in, out] x Position. At in: x(t=t_k) At out: x(t=t_{k+1})
 */
int RungeKutta4::Step(double &t, Vector4 &x) {

  // Assign k_i values
  for (int i=0; i<k_.cols(); ++i){
    // Get temporary position
    xtemp_ = x;
    for (int j=0; j<i; ++j){
      xtemp_ += a_(i,j)*k_.col(j);
    }
    // Evaluate the rhs at x_temp to determine k
    kGuidingCenter_.VectorField(t + c_(i)*kdt_, xtemp_, ftemp_);
    k_.col(i) = kdt_*ftemp_;
  }

  // Update x
  x = x + k_*b_;
  // Update t
  t = t + kdt_;

  /////////////////// DOUBLE * COMPATIBLE VERSION ////////////////////////
  // //// RECALL: the i-th row, j-th column of matrix m is m[j*n_rows + i]
  // // Assign k_i values
  // for (int i=0; i<kOrder_; ++i){
  //   // Get temporary position
  //   for(int j=0; j<kDimen_; ++j){
  //     xtemp_[j] = x[j];
  //   }
  //   for (int j=0; j<i; ++j){
  //     for (int k=0; k<kDimen_; ++k){
  // 	xtemp_[k] += a_(i,j)*k_[j*kDimen_ + k]; // Right??
  //     }
  //   }
  //   // Evaluate the rhs at x_temp to determine k
  //   model_.VectorField(t + c_(i)*kdt_, xtemp_, ftemp_);
  //   for (int j=0; j<kDimen_; ++j){
  //     k_[i*kDimen_ + j] = kdt_*ftemp_[j];
  //   }
  // }

  // // Update x
  // for(int j=0; j<kOrder_; ++j){
  //   for(int i=0; i<kDimen_; ++i){
  //     // x^i = k^i_j*b^j
  //     x[i] += k_[j*kDimen_ + i]*b_[j]; // Unintuitive and likely slow
  //   }
  // }
  // // Update t
  // t = t + kdt_;

  // Finally, check theta coordinate for exceeding 2pi rad
  // Assumes coordinates are cylindrical! 
  double pi=3.141592653589793;
  if(std::abs(x[1]) > 2*pi){
    int sign;
    if(x[1]>0){
      sign=1;
    }
    else{
      sign=-1;
    }
    x[1] -= sign*2*pi;
  }

  return 0;
}
