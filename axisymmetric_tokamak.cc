/*!
*******************************************************************************
* \file axisymmetric_tokamak.cc
*
* \brief EMFields derived class which implements the axisymmetric fields in Qin_2009.
*
* \date Feb 2014
* \author C. Leland Ellison
*******************************************************************************
*/

#include "axisymmetric_tokamak.h"

/*!
 * @brief Construct AxisymmetricTokamak which derives from EMFields. Set constants b0, r0.
 * 
 * @param[in] kB0 Magnetic field strength in [Tesla] (see documentation)
 * @param[in] kR0 Major radius of tokamak in [cm]
 */
AxisymmetricTokamak::AxisymmetricTokamak(const double kB0, const double kR0): 
  kB0_(kB0), kR0_(kR0){}

/*! 
 * @brief Evaluates Vector potential A
 *                                                                   
 * @param[in] kUnusedt Current time (fields time independent)              
 * @param[in] kx Position in cartesian coordinates                   
 * @param[out] a Vector potential A 
 */ 
void AxisymmetricTokamak::VectorPotentialA(const double kUnusedt, 
					   const Eigen::VectorXd &kx, 
					   Eigen::Vector3d &a) const{ 
 
  a(0) = kB0_*kR0_*kx[2]/(2*kx[0]);
 
  a(1) = sqrt(2)*kB0_*(pow(kx[2], 2) + pow(-kR0_ + kx[0], 2))/(4*kx[0]);
 
  a(2) = -kB0_*kR0_*log(kx[0]/kR0_)/2;
 
} 
 

/*! 
 * @brief Evaluates Matrix of derivatives of vector potential A
 *                                                                   
 * @param[in] kUnusedt Current time (fields time independent)              
 * @param[in] kx Position in cartesian coordinates                   
 * @param[out] grad_a Matrix of derivatives of vector potential A 
 */ 
void AxisymmetricTokamak::GradA(const double kUnusedt, 
				const Eigen::VectorXd &kx, 
				Eigen::MatrixXd &grad_a) const{ 
 
  grad_a(0,0) = -kB0_*kR0_*kx[2]/(2*pow(kx[0], 2));
 
  grad_a(0,1) = 0;
 
  grad_a(0,2) = kB0_*kR0_/(2*kx[0]);
 
  grad_a(1,0) = sqrt(2)*kB0_*(-2*kR0_ + 2*kx[0])/(4*kx[0]) - sqrt(2)*kB0_*(pow(kx[2], 2) + pow(-kR0_ + kx[0], 2))/(4*pow(kx[0], 2));
 
  grad_a(1,1) = 0;
 
  grad_a(1,2) = sqrt(2)*kB0_*kx[2]/(2*kx[0]);
 
  grad_a(2,0) = -kB0_*kR0_/(2*kx[0]);
 
  grad_a(2,1) = 0;
 
  grad_a(2,2) = 0;
 
} 
 

/*! 
 * @brief Evaluates Unit vector in the direction of the magnetic field
 *                                                                   
 * @param[in] kUnusedt Current time (fields time independent)              
 * @param[in] kx Position in cartesian coordinates                   
 * @param[out] b_hat Unit vector in the direction of the magnetic field 
 */ 
void AxisymmetricTokamak::BHat(const double kUnusedt, 
			       const Eigen::VectorXd &kx, 
			       Eigen::Vector3d &b_hat) const{ 
 
  b_hat(0) = -sqrt(2)*kB0_*kx[2]/(2*kx[0]*sqrt(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2))));
 
  b_hat(1) = kB0_*kR0_/(kx[0]*sqrt(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2))));
 
  b_hat(2) = sqrt(2)*kB0_*(-2*kR0_ + 2*kx[0])/(4*kx[0]*sqrt(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2))));
 
} 
 

/*! 
 * @brief Evaluates Gradient matrix of magnetic field unit vector
 *                                                                   
 * @param[in] kUnusedt Current time (fields time independent)              
 * @param[in] kx Position in cartesian coordinates                   
 * @param[out] grad_b_hat Gradient matrix of magnetic field unit vector 
 */ 
void AxisymmetricTokamak::GradBHat(const double kUnusedt, 
				   const Eigen::VectorXd &kx, 
				   Eigen::MatrixXd &grad_b_hat) const{ 
 
  grad_b_hat(0,0) = -sqrt(2)*kB0_*kx[2]*(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 3) - pow(kB0_, 2)*(-8*kR0_ + 8*kx[0])/(16*pow(kx[0], 2)) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 3)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 3)))/(2*kx[0]*pow(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2)), 3.0/2.0)) + sqrt(2)*kB0_*kx[2]/(2*pow(kx[0], 2)*sqrt(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2))));
 
  grad_b_hat(0,1) = 0;
 
  grad_b_hat(0,2) = sqrt(2)*pow(kB0_, 3)*pow(kx[2], 2)/(4*pow(kx[0], 3)*pow(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2)), 3.0/2.0)) - sqrt(2)*kB0_/(2*kx[0]*sqrt(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2))));
 
  grad_b_hat(1,0) = kB0_*kR0_*(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 3) - pow(kB0_, 2)*(-8*kR0_ + 8*kx[0])/(16*pow(kx[0], 2)) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 3)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 3)))/(kx[0]*pow(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2)), 3.0/2.0)) - kB0_*kR0_/(pow(kx[0], 2)*sqrt(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2))));
 
  grad_b_hat(1,1) = 0;
 
  grad_b_hat(1,2) = -pow(kB0_, 3)*kR0_*kx[2]/(2*pow(kx[0], 3)*pow(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2)), 3.0/2.0));
 
  grad_b_hat(2,0) = sqrt(2)*kB0_*(-2*kR0_ + 2*kx[0])*(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 3) - pow(kB0_, 2)*(-8*kR0_ + 8*kx[0])/(16*pow(kx[0], 2)) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 3)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 3)))/(4*kx[0]*pow(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2)), 3.0/2.0)) + sqrt(2)*kB0_/(2*kx[0]*sqrt(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2)))) - sqrt(2)*kB0_*(-2*kR0_ + 2*kx[0])/(4*pow(kx[0], 2)*sqrt(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2))));
 
  grad_b_hat(2,1) = 0;
 
  grad_b_hat(2,2) = -sqrt(2)*pow(kB0_, 3)*kx[2]*(-2*kR0_ + 2*kx[0])/(8*pow(kx[0], 3)*pow(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2)), 3.0/2.0));
 
} 
 

/*! 
 * @brief Evaluates Gradient of scalar potential phi
 *                                                                   
 * @param[in] kUnusedt Current time (fields time independent)              
 * @param[in] kUnusedx Position in cartesian coordinates                   
 * @param[out] grad_phi Gradient of scalar potential phi 
 */ 
void AxisymmetricTokamak::GradPhi(const double kUnusedt, 
				  const Eigen::VectorXd &kUnusedx, 
				  Eigen::Vector3d &grad_phi) const{ 
 
  grad_phi(0) = 0;
 
  grad_phi(1) = 0;
 
  grad_phi(2) = 0;
 
} 
 

/*! 
 * @brief Evaluates Gradient of magnetic field magnitude
 *                                                                   
 * @param[in] kUnusedt Current time (fields time independent)              
 * @param[in] kx Position in cartesian coordinates                   
 * @param[out] grad_mod_b Gradient of magnetic field magnitude 
 */ 
void AxisymmetricTokamak::GradModB(const double kUnusedt, 
				   const Eigen::VectorXd &kx, 
				   Eigen::Vector3d &grad_mod_b) const{ 
 
  grad_mod_b(0) = (-pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 3) + pow(kB0_, 2)*(-8*kR0_ + 8*kx[0])/(16*pow(kx[0], 2)) - pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 3)) - pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 3)))/sqrt(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2)));
 
  grad_mod_b(1) = 0;
 
  grad_mod_b(2) = pow(kB0_, 2)*kx[2]/(2*pow(kx[0], 2)*sqrt(pow(kB0_, 2)*pow(kR0_, 2)/pow(kx[0], 2) + pow(kB0_, 2)*pow(kx[2], 2)/(2*pow(kx[0], 2)) + pow(kB0_, 2)*pow(-2*kR0_ + 2*kx[0], 2)/(8*pow(kx[0], 2))));
 
} 

