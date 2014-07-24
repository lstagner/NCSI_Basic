/*!
*******************************************************************************
* \file em_fields.h
*
* \brief Header for EMFields abstract base class which defines the interfaces for electromagnetic field quantities used in calculating, for instance,  guiding center trajectories. 
*
* \date Feb 2014
* \author C. Leland Ellison
*********************************************************************************/

#ifndef EM_FIELDS_H_
#define EM_FIELDS_H_

#include <Eigen/Dense>
#include "eigen_types.h"
#include "cuda_member.h"

class EMFields{
 public:
  //! Fetch vector potential A 
  CUDA_MEMBER virtual void VectorPotentialA(const double kt, 
					    const Vector3 &kx,  
					    Vector3 &a) const = 0; 
  //! Fetch gradient matrix of the vector potential A  
  CUDA_MEMBER virtual void GradA(const double kt, const Vector3 &kx,
				 Matrix3 &grad_a) const = 0;
  //! Fetch unit vector in direction of magnetic field
  CUDA_MEMBER virtual void BHat(const double kt, const Vector3 &kx, 
				Vector3 &b_hat) const = 0;
  //! Fetch gradient matrix of magnetic field unit vector
  CUDA_MEMBER virtual void GradBHat(const double kt, const Vector3 &kx, 
				    Matrix3 &grad_b_hat) const = 0;
  //! Fetch gradient of scalar potential
  CUDA_MEMBER virtual void GradPhi(const double kt, const Vector3 &kx, 
				   Vector3 &grad_phi) const = 0;
  //! Fetch gradient of magnetic field magnitude
  CUDA_MEMBER virtual void GradModB(const double kt, const Vector3 &kx,
				    Vector3 &grad_mod_b) const = 0;
  // Accessors for parameters which may be set in derived classes
  CUDA_MEMBER virtual double kR0() const { return 0.0; }
  CUDA_MEMBER virtual double kB0() const { return 0.0; }
};


#endif // EM_FIELDS_H_
