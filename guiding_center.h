/*!
*******************************************************************************
* \file guiding_center.h
*
* \brief Header for GuidingCenter system describing the motion of charged magnetized particles in magnetic fields.
*
* Guiding center equations of motion
*   dot{x}_i = (B^dag_i u - (b x E^dag)_i) / B^dag_\par 
*   dot{u} = (B^dag dot E^dag)/B^dag_par
*
* \date Feb 2014
* \author C. Leland Ellison
*********************************************************************************/

#ifndef GUIDING_CENTER_H_
#define GUIDING_CENTER_H_

#include <Eigen/Dense>
#include "em_fields.h"  // Defines electric/magnetic fields

// Guiding Center System
class GuidingCenter {
 public:
  GuidingCenter(EMFields *em_fields, const double kMu);
  // Evaluates vector field of ODE. f in \dot{x} = f(x)
  int VectorField(const double kt, const Eigen::VectorXd &kx, 
		  Eigen::VectorXd &fx) const;

  // Accessors
  //! Return pointer to em_fields instance
  EMFields *em_fields() const { return em_fields_; }
  int kDimen() const {return kDimen_; }  //!< Return dimension of ODE system
  double kMu() const {return kMu_;}  //!< Return value of mu
  
 private:
  static const int kDimen_ = 4;  //!< Dimension of ODE system
  EMFields *em_fields_;  //!< Pointer to class defining electromagnetic fields
  const double kMu_;  //!< Magnetic moment
};

#endif  // GUIDING_CENTER_H_
