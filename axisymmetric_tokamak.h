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
#include "eigen_types.h"

class AxisymmetricTokamak : public EMFields{
 public:
  AxisymmetricTokamak(const double kB0, const double kR0);
  void VectorPotentialA(const double kt, const Vector3 &kx, Vector3 &a) const;
  void GradA(const double kt, const Vector3 &kx, Matrix3 &grad_a) const;
  void BHat(const double kt, const Vector3 &kx, Vector3 &b_hat) const;
  void GradBHat(const double kt, const Vector3 &kx, Matrix3 &grad_b_hat) const;
  void GradPhi(const double kt, const Vector3 &kx, Vector3 &grad_phi) const;
  void GradModB(const double kt, const Vector3 &kx, Vector3 &grad_mod_b) const;
  double kR0() const { return kR0_;}
  double kB0() const { return kB0_;}
 private:
  const double kB0_; //!< Magnetic field on-axis amplitude in [Tesla]
  const double kR0_; //!< Major radius of tokamak in [cm]
};
#endif // AXISYMMETRIC_TOKAMAK_H_
