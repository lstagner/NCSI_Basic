/*!
*******************************************************************************
* \file integrator.h
*
* \brief Header for Integrator abstract base class. Step is pure virtual. Declares members kdt_, model_ and kDimen_.
*
* \date Feb 2014
* \author C. Leland Ellison
*******************************************************************************
*/

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "guiding_center.h" // ODE system being modeled
#include <Eigen/Dense> // Provides linear algebra capabilities

class Integrator {
 public:
  /*!
   * \brief Constructor - Save model, stepsize, and dimension
   * @param[in] kdt Numerical Step Size
   * @param[in] kGuidingCenter ODE Model instance
  */
  Integrator(const double kdt, const GuidingCenter &kGuidingCenter):
      kdt_(kdt), kGuidingCenter_(kGuidingCenter), 
	kDimen_(kGuidingCenter.kDimen()) {}
  virtual ~Integrator() {}
  //! Time-advance map from (t_k, x_k) -> (t_{k+1}, x_{k+1}).  
  virtual int Step(double &t, Eigen::VectorXd &x) = 0;
  //! Method used in multistep integrators
  virtual void Reset(){}; 
 protected:
  const double kdt_;  //!< Numerical Step Size
  const GuidingCenter &kGuidingCenter_;  //!< ODE model
  const int kDimen_;  //!< Dimension of the ODE system
};

#endif  // INTEGRATOR_H_
