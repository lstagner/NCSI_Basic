/*!
*******************************************************************************
* \file noncanonical_symplectic.cc
*
* \brief Implementation of integrator resulting from a midpoint discretization of the guiding center Lagrangian. Implements a step which interfaces with the base class integrators and defines the update rule. 
*
* \date Jan 2014
* \author C. Leland Ellison
*******************************************************************************
*/

#include "noncanonical_symplectic.h"

/*!
 * Constructor - Initializes members and base class Integrator. Sets size of
 x_history_ to have a number of rows equal to the dimension of the ODE system and number of columns equal to the number of steps in the multistep method. In this case, two.
 *
 * @param[in] kdt Numerical Step Size
 * @param[in] kGuidingCenter ODE being modeled
 * @param[in] kNewtonTolerance Error threshold for nonlinear solve
 * @param[in] kMaxIterations Maximum number of nonlinear solve iterations
 */
NoncanonicalSymplectic::NoncanonicalSymplectic(const double kdt, 
					  const GuidingCenter &kGuidingCenter, 
					  const double kNewtonTolerance,
					  const double kMaxIterations) 
  : Integrator(kdt, kGuidingCenter), kMu_(kGuidingCenter.kMu()), 
    kNewtonTolerance_(kNewtonTolerance), kMaxIterations_(kMaxIterations),
    needs_initialization_(true), x_history_(kGuidingCenter.kDimen(), 2){
  x_history_.setZero(); // Just to be safe.
  em_fields_ = kGuidingCenter.em_fields();
}


/*!
 * \brief Advance t and x forward in time.
 *
 * If this is the first step, some additional initial conditions need to be generated. At present, this is performed using RK4. Otherwise, solve the implicit map defined by the discrete Euler-Lagrange equations.
 * 
 * @param[in] t Physical time at this step. For time-dependent systems.
 * @param[in,out] x Vector of current coordinates. Will be updated to x(t+h). 
 * @return Integer code which is 0 if successful. 
 */
int NoncanonicalSymplectic::Step(double &t, Eigen::VectorXd &x){
  if(needs_initialization_){

    StoreHistory(x); // x_history_ --> [0, x0]

    needs_initialization_ = false; // No longer need init

    return InitialStep(t, x); // RK4 advance
  }
  else{
    StoreHistory(x); // x_history_ --> [x_k-1, x_k]
   
    // Declare local variables used in nonlinear solve
    Eigen::VectorXd error(kDimen_);
    int n_iterations = 0;
    Eigen::MatrixXd jacobian(kDimen_, kDimen_);

    NewtonGuess(t, x); // Get estimate of new position. Advances (t,x).
    UpdateRule(t, x, error); // Check how far t,x are from solution
    
    // Newton-Rhapson nonlinear solve
    while((error.norm() > kNewtonTolerance_) && 
	  (n_iterations < kMaxIterations_)){
      n_iterations++; // Increment how many times we've done this
      Jacobian(t, x, jacobian);  // Fetch new Jacobian
      x -= jacobian.inverse()*error; // Adjust x
      UpdateRule(t, x, error); // Check new error
    }
    if (n_iterations == kMaxIterations_){
      std::cout << "Newton solve did not converge!!" << std::endl;
    }
  } // end else 

  // Finally, check theta coordinate for exceeding 2pi rad
  // Assumes coordinates are cylindrical! 
  double pi=3.141592653589793;
  double r0 = em_fields_->kR0();
  if(abs(x[1]) > 2*pi*r0){
    int sign;
    if(x[1]>0){
      sign=1;
    }
    else{
      sign=-1;
    }
    x[1] -= sign*2*pi*r0;
    x_history_(1,1) -= sign*2*pi*r0;
  }

  return 0;
}


/*! 
 * @brief StoreHistory updates member x_history_  by knocking off the oldest data and putting on the newest. 
 * 
 * @param[in] kx New position to store
 * @returns 0 upon success
 */
int NoncanonicalSymplectic::StoreHistory(const Eigen::VectorXd &kx){
  // Shift everything back one
  x_history_.col(0) = x_history_.col(1);
  
  // Next, write the new data at the end
  x_history_.col(1) = kx;

  return 0;
}


/*! InitialStep advances t and x using RK4. To be used while needs_initialization_ is true.
 * 
 * @param[in,out] t Simulation time. Advanced by kdt_. 
 * @param[in,out] x Position. At in: x(t=t_k) At out: x(t=t_k+1)
 * @returns 0 upon success
 */
int NoncanonicalSymplectic::InitialStep(double &t, 
						  Eigen::VectorXd &x) const{
  RungeKutta rk4(kdt_, kGuidingCenter_, 4);
  return rk4.Step(t,x);
}


/*!
 * \brief Function which should evaluate to zero(vector) when the algorithm is satisfied.
 *
 * @param[in] kt Simulation time at proposed new position
 * @param[in] kx New position to test
 * @param[out] error Error vector which is zero for satisfied update rule.
 * @return Zero for success
 */
int NoncanonicalSymplectic::UpdateRule(const double kt,
				       const Eigen::VectorXd &kx,
				       Eigen::VectorXd &error) const
{
  //// Define intermediate variables
  // Conventions: kp1 --> ``k plus 1"
  //              km1 --> ``k minus 1"
  //              kphalf --> ``k plus 1/2"
  //              kmhalf --> ``k minus 1/2"
  double t_kphalf;
  double t_kmhalf;
  Eigen::Vector3d x_kphalf;
  Eigen::Vector3d x_kmhalf;
  Eigen::Vector3d delta_x_kphalf;
  Eigen::Vector3d delta_x_kmhalf;
  Eigen::Vector3d a_kphalf;
  Eigen::Vector3d a_kmhalf;
  Eigen::Vector3d b_hat_kphalf;
  Eigen::Vector3d b_hat_kmhalf;
  Eigen::MatrixXd grad_a_kphalf(3,3);
  Eigen::MatrixXd grad_a_kmhalf(3,3);
  Eigen::MatrixXd grad_b_hat_kphalf(3,3);
  Eigen::MatrixXd grad_b_hat_kmhalf(3,3);
  Eigen::Vector3d grad_mod_b_kphalf;
  Eigen::Vector3d grad_mod_b_kmhalf;
  Eigen::Vector3d grad_phi_kphalf;
  Eigen::Vector3d grad_phi_kmhalf;
  Eigen::MatrixXd grad_a_dag_kphalf(3,3);
  Eigen::MatrixXd grad_a_dag_kmhalf(3,3);
  Eigen::Vector3d a_dag_kphalf;
  Eigen::Vector3d a_dag_kmhalf;
  
  //// Update intermediate variables
  // t coming in is t_kp1. Define t_kphalf, tkmhalf
  t_kphalf = kt - 0.5*kdt_;
  t_kmhalf = kt - 1.5*kdt_;

  // Calculate x from history and kx (which is x_kp1)
  // Recall x_history_ contains x_km1 in column 0 and x_k in column 1
  // block<a,b>(c,d) selects a rows and b cols of a matrix starting from (c,d)
  // head<a> picks first a elements of a vector
  x_kphalf = 0.5*(x_history_.block<3,1>(0,1) + kx.head<3>());
  x_kmhalf = 0.5*(x_history_.block<3,1>(0,0) + x_history_.block<3,1>(0,1));
  delta_x_kphalf = kx.head<3>() - x_history_.block<3,1>(0,1);
  delta_x_kmhalf = x_history_.block<3,1>(0,1) - x_history_.block<3,1>(0,0);

  // Fields
  em_fields_->VectorPotentialA(t_kphalf, x_kphalf, a_kphalf);
  em_fields_->VectorPotentialA(t_kmhalf, x_kmhalf, a_kmhalf);
  em_fields_->BHat(t_kphalf, x_kphalf, b_hat_kphalf);
  em_fields_->BHat(t_kmhalf, x_kmhalf, b_hat_kmhalf);
  em_fields_->GradA(t_kphalf, x_kphalf, grad_a_kphalf);
  em_fields_->GradA(t_kmhalf, x_kmhalf, grad_a_kmhalf);
  em_fields_->GradBHat(t_kphalf, x_kphalf, grad_b_hat_kphalf);
  em_fields_->GradBHat(t_kmhalf, x_kmhalf, grad_b_hat_kmhalf);
  em_fields_->GradModB(t_kphalf, x_kphalf, grad_mod_b_kphalf);
  em_fields_->GradModB(t_kmhalf, x_kmhalf, grad_mod_b_kmhalf);
  em_fields_->GradPhi(t_kphalf, x_kphalf, grad_phi_kphalf);
  em_fields_->GradPhi(t_kmhalf, x_kmhalf, grad_phi_kmhalf);

  // Convenience variables: A^\dagger = [A + u b]
  a_dag_kphalf = a_kphalf + kx[3]*b_hat_kphalf;
  a_dag_kmhalf = a_kmhalf + x_history_(3,1)*b_hat_kmhalf;
  grad_a_dag_kphalf = grad_a_kphalf + kx[3]*grad_b_hat_kphalf;
  grad_a_dag_kmhalf = grad_a_kmhalf + x_history_(3,1)*grad_b_hat_kmhalf;

  //// Evaluate variational integrator update rule
  // j = 0, 1, 2
  // 1/2[A_{i,j}(k+1/2) + u_{k+1/2} b_{i,j}(k+1/2)]*[x_{k+1}^i - x_k^i] +
  // 1/2[A_{i,j}(k-1/2) + u_{k-1/2} b_{i,j}(k-1/2)]*[x_{k}^i - x_{k-1}^i] - 
  // [A_j(k+1/2) + u_{k+1/2} b_j(k+1/2) - A_j(k-1/2) - u_{k-1/2} b_j(k-1/2)] -
  // h/2[muB_{,j}(k+1/2) + muB_{,j}(k-1/2) + phi_{,j}(k+1/2) + phi_{,j}(k-1/2)]
  // NOTE: Transpose is to sum over proper index
  error.head(3) = 0.5*grad_a_dag_kphalf.transpose()*delta_x_kphalf + 
    0.5*grad_a_dag_kmhalf.transpose()*delta_x_kmhalf - 
    (a_dag_kphalf - a_dag_kmhalf) - kdt_*0.5*(kMu_*(grad_mod_b_kphalf +
    grad_mod_b_kmhalf) + grad_phi_kphalf + grad_phi_kmhalf);
  
  // j = 3
  // b_i(k+1/2)*(x_{k+1}^i - x_k^i) - h*u_{k+1/2} = 0
  error[3] = b_hat_kphalf.dot(delta_x_kphalf) - kdt_*kx[3];

  return 0;
}


/*!
 * \brief Guess for initializing NewtonSolver. Uses ForwardEuler.
 *
 * @param[in,out] t Simulation time. Advanced by step size kdt
 * @param[in,out] x Position to advance to initial guess for solver
 * @return Zero for success
 */
int NoncanonicalSymplectic::NewtonGuess(double &t, 
						  Eigen::VectorXd &x) const {
  // Stores output of VectorField
  Eigen::VectorXd fx(kDimen_);

  // Fetch f(x) in \dot{x} = f(x)
  kGuidingCenter_.VectorField(t,x,fx); 
  t += kdt_;
  x = x + kdt_*fx; 
  return 0;
}


/*!
 * \brief Calculates Jacobian matrix used in nonlinear solve. This base-class version uses a centered finite difference on the update rule.
 * 
 * @param[in] kt Simulation time at newest point.
 * @param[in] kx Position at newest point.
 * @param[out] jacobian derivative of the update rule w.r.t the new position
 * @return Integer code which is 0 if successful. 
 */
int NoncanonicalSymplectic::Jacobian(const double kt, 
					       const Eigen::VectorXd &kx, 
					       Eigen::MatrixXd &jacobian) 
  const{

  double delta = 1e-5; // Step size for finite difference

  Eigen::VectorXd left_evaluation(kDimen_);
  Eigen::VectorXd right_evaluation(kDimen_);
  Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(kDimen_, kDimen_);

  // jacobian needs to be kDimen by kDimen!
  for (int i=0; i<kDimen_; ++i){
    UpdateRule(kt, kx + 0.5*delta*eye.col(i), right_evaluation);
    UpdateRule(kt, kx - 0.5*delta*eye.col(i), left_evaluation);
    jacobian.col(i) = (right_evaluation - left_evaluation)/delta;
  }

  return 0;
}
