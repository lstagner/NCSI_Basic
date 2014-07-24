/*!
*******************************************************************************
* \file runge-kutta.h
*
* \brief Header for general Runge-Kutta integrator. Modified for CUDA.
*
* \date July 2014
* \author C. Leland Ellison
* \version 0.1
*
*******************************************************************************
*/
#ifndef RUNGEKUTTA4_H_
#define RUNGEKUTTA4_H_

#include "integrator.h"
#include <Eigen/Dense>
#include "guiding_center.h"
#include "cuda_member.h"
#include "eigen_types.h"

class RungeKutta4 : public Integrator {
 public:
  CUDA_MEMBER RungeKutta4(const double kdt, const Model &model);
  CUDA_MEMBER ~RungeKutta4();
  // Map from (t_k, x_k) -> (t_{k+1}, x_{k+1})
  CUDA_MEMBER int Step(double &t, Vector4 &x);  
  //// Accessors
  //! Return name of integrator: runge-kutta4
  HOST_MEMBER std::string IntegratorName() const { return "runge-kutta4";  }
  //! Access a_ coefficients matrix
  HOST_MEMBER Matrix4 a_coefficients() const { return a_; }
  //! Access b_ coefficients vector
  HOST_MEMBER Vector4 b_coefficients() const { return b_; }
  //! Access c_ coefficients vector
  HOST_MEMBER Vector4 c_coefficients() const { return c_; }
 private:
  Matrix4 a_;  //!< Specifies for a_ij coefficients of the RK method
  Vector4 b_;  //!< Specifies b_i coefficients of the RK method
  Vector4 c_;  //!< Specifies c_i coefficients of the RK method
  Matrix4 k_;  //!< Temporary space for k_i evaluations
  Vector4 xtemp_; //!< Temporary space used at internal stages
  Vector4 ftemp_; //!< Temporary space for function evaluations
  const int kOrder_;
};

#endif  // RUNGEKUTTA4_H_
