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

class EMFields{
 public:
  //! Fetch vector potential A 
  virtual void VectorPotentialA(const double kt, const Eigen::VectorXd &kx,  
    				Eigen::Vector3d &a) const = 0; 
  //! Fetch gradient matrix of the vector potential A  
  virtual void GradA(const double kt, const Eigen::VectorXd &kx, 
		     Eigen::MatrixXd &grad_a) const = 0;
  //! Fetch unit vector in direction of magnetic field
  virtual void BHat(const double kt, const Eigen::VectorXd &kx, 
		    Eigen::Vector3d &b_hat) const = 0;
  //! Fetch gradient matrix of magnetic field unit vector
  virtual void GradBHat(const double kt, const Eigen::VectorXd &kx, 
			Eigen::MatrixXd &grad_b_hat) const = 0;
  //! Fetch gradient of scalar potential
  virtual void GradPhi(const double kt, const Eigen::VectorXd &kx, 
		       Eigen::Vector3d &grad_phi) const = 0;
  // Fetch gradient of magnetic field magnitude
  virtual void GradModB(const double kt, const Eigen::VectorXd &kx,
			Eigen::Vector3d &grad_mod_b) const = 0;
};


#endif // EM_FIELDS_H_
