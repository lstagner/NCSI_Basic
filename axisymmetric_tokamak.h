/*!
*******************************************************************************
* \file axisymmetric_tokamak.h
*
* \brief EMFields derived class which implements the axisymmetric fields in Qin_2009 in cylindrical coordinates.
*
* \date Feb 2014
* \author C. Leland Ellison
*******************************************************************************
*/

#ifndef AXISYMMETRIC_TOKAMAK_H_
#define AXISYMMETRIC_TOKAMAK_H_

#include "em_fields.h"
#include <Eigen/Dense>

class AxisymmetricTokamak : public EMFields{
 public:
  AxisymmetricTokamak(const double kB0 = 1.0, const double kR0 = 1.0);
  void VectorPotentialA(const double kt, const Eigen::VectorXd &kx,
			Eigen::Vector3d &a) const;
  void GradA(const double kt, const Eigen::VectorXd &kx, 
	     Eigen::MatrixXd &grad_a) const;
  void BHat(const double kt, const Eigen::VectorXd &kx, 
	    Eigen::Vector3d &b_hat) const;
  void GradBHat(const double kt, const Eigen::VectorXd &kx, 
		Eigen::MatrixXd &grad_b_hat) const;
  void GradPhi(const double kt, const Eigen::VectorXd &kx, 
	       Eigen::Vector3d &grad_phi) const;
  void GradModB(const double kt, const Eigen::VectorXd &kx, 
		Eigen::Vector3d &grad_mod_b) const;
 private:
  const double kB0_;
  const double kR0_;
};
#endif // AXISYMMETRIC_TOKAMAK_H_
