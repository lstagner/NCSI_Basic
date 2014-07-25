/*!
*******************************************************************************
* \file driver.cc
*
* \brief Implements driver for integrating guiding center trajectories
*
* \date Feb 2014
* \author C. Leland Ellison
*******************************************************************************
*/

#include <stdlib.h>
#include <iostream>
#include <Eigen/Dense>
#include <ctime> // For timing runtime
#include "input_parser.h"
#include "guiding_center.h"
#include "em_fields.h"
#include "axisymmetric_tokamak.h"
#include "integrator.h"
#include "runge-kutta4.h"
#include "noncanonical_symplectic.h"
#include "eigen_types.h"

#define GC_DIM 4 // Dimension of guiding center system

/*!
 * \brief Prints a line to standard out giving [time x[0] x[1] ....]
 *
 * @param[in] t Time
 * @param[in] x Position vector
 * @param[in] n_digits Number of digits to display in output
 */
void PrintState(double t, const Vector4 &x, int n_digits){
  // Formatting options: n_digits precision, all vector elements on one line
  static Eigen::IOFormat OneLineNDeep(n_digits, 0, "     ", "     ",
				      "", "", "", "");
  std::cout.setf( std::ios::fixed, std:: ios::floatfield );
  std::cout << t << "     " << x.format(OneLineNDeep) << std::endl;
}

/*!
 * \brief Body of the driver. Use program options to specify ode, integrator, dt, and n_steps.
 *
 */
int main(int argc, char *argv[]) {
  
  //// Read and extract input
  InputParser input_parser;
  int read_result;
  read_result = input_parser.ReadInput(argc, argv);
  if (read_result){return read_result; } // Quits on, e.g. --help

  // Initialize runtime parameters determined by input
  // double dt, b0, r0, mu, newton_tolerance;
  // int n_steps, save_nth, print_precision, max_iterations;
  // std::vector<double> initial_conditions; 
  // bool time_flag;
  // std::string integrator_name;

  // Retrieve the values from the input parser
  const double kdt = input_parser.GetValue<double>("dt");
  const int kNSteps = input_parser.GetValue<int>("n_steps");
  const int kSaveNth = input_parser.GetValue<int>("save_nth");
  std::vector<double> initial_conditions = 
    input_parser.GetValue<std::vector<double> >("initial_conditions");
  const bool kTimeFlag = input_parser.GetValue<bool>("time");
  const int kPrintPrecision = input_parser.GetValue<int>("precision");
  const double kB0 = input_parser.GetValue<double>("b0");
  const double kR0 = input_parser.GetValue<double>("r0");
  const std::string kIntegratorName = 
    input_parser.GetValue<std::string>("integrator");
  const double kMu = input_parser.GetValue<double>("mu");
  const double kSolveTolerance = input_parser.GetValue<double>("tol");
  const int kMaxIterations = input_parser.GetValue<int>("max_iter");

  //// Initialize model and integrator
  AxisymmetricTokamak em_fields(kB0, kR0);
  GuidingCenter guiding_center((EMFields *) &em_fields, kMu);
  Integrator *integrator;
  if (kIntegratorName.compare("rk4")==0){
    integrator = new RungeKutta4(kdt, guiding_center);
  }
  else if (kIntegratorName.compare("ncsi")==0){
    integrator = new NoncanonicalSymplectic(kdt, guiding_center, 
					    kSolveTolerance, kMaxIterations);
  }
  else{
    std::cout << "Unrecognized integrator. Try rk4 or ncsi" 
	      << std::endl;
    return 1;
  }

  //// Initial conditions
  // Make a vector matching the size of the ode system
  Vector4 x; 
  int n_initial_conditions;
  
  // If no initial conditions specified or they are of the wrong dimension
  if((!initial_conditions.size()) || 
     (initial_conditions.size() % (GC_DIM) )){
    // Use the default initial conditions: Vector of ones
    n_initial_conditions = 1;
    for (int i = 0; i < x.size(); ++i) {
      initial_conditions.push_back(1.0);
    }
  }
  else{
   // Otherwise, we have more than one initial condition to simulate
   n_initial_conditions = initial_conditions.size()/GC_DIM;
  }

  //// Record time?
  std::clock_t run_time;
  if(kTimeFlag){
    run_time = std::clock();
  }

  //// Time advance each initial condition
  for (int j = 0; j < n_initial_conditions; ++j){
    integrator->Reset(); // Resets temporary variables in integrators
    double t = 0;
    // Set initial condition
    for (int i=0; i<GC_DIM; ++i){
      x[i] = initial_conditions[j*GC_DIM + i];
    }
    PrintState(t, x, kPrintPrecision); // Print initial position
    // Run standard stepping
    for (int i = 1; i <= kNSteps; ++i) {
      integrator->Step(t, x);
      if(!(i%kSaveNth)){
	PrintState( t, x, kPrintPrecision);
      }
    }
  }

  //// Print run time?
  if(kTimeFlag){
    run_time = std::clock() - run_time;
    std::cout << "Run time: " << (double)run_time/CLOCKS_PER_SEC << std::endl;
  }

  //// Clean up
  delete integrator;
  return 0;
}
