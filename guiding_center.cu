/*! 
****************************************************************************
* @file guiding_center.cc 
* 
* @brief Implementation of GuidingCenter system. 
* 
* @date Feb 2014 
* @author C. Leland Ellison 
****************************************************************************/
 

#include "guiding_center.h" 
 
/*! 
 * @brief Constructor which sets magnetic moment and initializes fields 
 * 
 * @param[in] em_fields
 * @param[in] kMu Magnetic moment in [cm^2/10^-8 s] 
 */
GuidingCenter::GuidingCenter(EMFields *em_fields, const double kMu) 
  : kDimen_(4), em_fields_(em_fields), kMu_(kMu){}


/*!  
 * @brief Evaluation of guiding center equations of motion.  
 *  
 * Letting x be particle position and u the velocity:  
 *   dot{x}_i = (B^dag_i u - (b x E^dag)_i) / B^dag_\par  
 *   dot{u} = (B^dag dot E^dag)/B^dag_par  
 *  
 * @param[in] kt Time  
 * @param[in] kx Position of guiding center particle [x u]  
 * @param[out] fx Right hand side of dot{x} = f(x)  
 * @return zero if success  
 */
int GuidingCenter::VectorField(const double kt, const Vector4 &kx, 
			       Vector4 &fx) const {

  // Declare intermediate convenience variables
  Vector3 b_hat;  // b_hat = magnetic field unit vector
  Matrix3 da;  // Gradient matrix of vector potential
  Matrix3 db_hat;  // Gradient matrix of Bfield unit vector
  Vector3 grad_phi;  // Gradient of scalar potential phi
  Vector3 grad_mod_b;  // Gradient of magnitute of B field
  Vector3 b_dag;  // B^dag = curl(A + ub)
  Vector3 e_dag;  // E^dag = -grad(phi + mu mod_B)
  double b_dag_par; // B^dag_parallel =  B^dag dot b
  Vector3 b_cross_e_dag;  // variable for b x E^dag

  // Retrieve fields at current position, time from em_fields_
  //   kx.head(3) picks off the first three elements of the vector kx = [x u]
  // Fill in b_hat
  em_fields_->BHat(kt, kx.head(3), b_hat); 
  // Fill in da
  em_fields_->GradA(kt, kx.head(3), da);
  // Fill in dbhat
  em_fields_->GradBHat(kt, kx.head(3), db_hat);
  // Fill in grad_phi
  em_fields_->GradPhi(kt, kx.head(3), grad_phi);
  // Fill in grad_mod_b
  em_fields_->GradModB(kt, kx.head(3), grad_mod_b);
  
  // Update intermediate variables
  // Fill in b_dag = curl(A + ub)
  b_dag[0] = da(2,1) + kx[3]*db_hat(2,1) - da(1,2) - kx[3]*db_hat(1,2);
  b_dag[1] = da(0,2) + kx[3]*db_hat(0,2) - da(2,0) - kx[3]*db_hat(2,0);
  b_dag[2] = da(1,0) + kx[3]*db_hat(1,0) - da(0,1) - kx[3]*db_hat(0,1);
  // Fill in e_dag = -grad(phi + mu mod_B)
  e_dag = -1.0*(grad_phi + kMu_*grad_mod_b);
  // Remaining intermediate variables
  b_dag_par = b_dag.dot(b_hat);
  b_cross_e_dag = b_hat.cross(e_dag);

  // Evaluate ODE vector field
  for (int i=0; i<3; ++i){
    fx[i] = (b_dag[i]*kx[3] - b_cross_e_dag[i])/b_dag_par;
  }
  fx[3] = b_dag.dot(e_dag)/b_dag_par;

  return 0;
} 
 

