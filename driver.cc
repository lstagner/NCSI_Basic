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


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <Eigen/Dense>
#include <ctime> // For timing runtime
#include "input_parser.h"
#include "guiding_center.h"
#include "em_fields.h"
#include "axisymmetric_tokamak.h"
#include "integrator.h"
#include "runge-kutta.h"
#include "variational_guiding_center_midpoint.h"

/*!
 * \brief Prints a line to standard out giving [time x[0] x[1] ....]
 *
 * @param[in] t Time
 * @param[in] x Position vector
 */
void PrintState(double t, const Eigen::VectorXd &x, int n_digits){
  // Formatting options: n_digits precision, all vector elements on one line
  static Eigen::IOFormat OneLineNDeep(n_digits, 0, "     ", "     ",
				      "", "", "", "");
  std::cout.setf( std::ios::fixed, std:: ios::floatfield ); // static?
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
  double dt, b0, r0, mu;
  int n_steps, save_nth, print_precision;
  std::vector<double> initial_conditions; 
  bool time_flag;
  std::string integrator_name;

  // Retrieve the values from the input parser
  input_parser.GetValue("dt", dt);
  input_parser.GetValue("n_steps", n_steps);
  input_parser.GetValue("save_nth", save_nth);
  input_parser.GetValue("initial_conditions", initial_conditions);
  input_parser.GetValue("time", time_flag);
  input_parser.GetValue("precision", print_precision);
  input_parser.GetValue("b0", b0);
  input_parser.GetValue("r0", r0);
  input_parser.GetValue("integrator", integrator_name);
  input_parser.GetValue("mu", mu);

  //// Initialize model and integrator
  EMFields *em_fields = new AxisymmetricTokamak(b0, r0);
  GuidingCenter *guiding_center = new GuidingCenter(em_fields, mu);
  Integrator *integrator;
  if (integrator_name.compare("runge-kutta4")==0){
    integrator = new RungeKutta(dt, *guiding_center, 4);
  }
  else if (integrator_name.compare("variational")==0){
    integrator = new VariationalGuidingCenterMidpoint(dt, *guiding_center);
  }
  else{
    std::cout << "Unrecognized integrator. Try runge-kutta4 or variational" 
	      << std::endl;
    return 1;
  }

  //// Initial conditions
  // Make a vector matching the size of the ode system
  Eigen::VectorXd x(guiding_center->kDimen()); 
  int n_initial_conditions;
  
  // If no initial conditions specified or they are of the wrong dimension
  if((!initial_conditions.size()) || 
     (initial_conditions.size() % (guiding_center->kDimen()) )){
    // Use the default initial conditions: Vector of ones
    n_initial_conditions = 1;
    for (int i = 0; i < x.size(); ++i) {
      initial_conditions.push_back(1.0);
    }
  }
  else{
   // Otherwise, we have more than one initial condition to simulate
   n_initial_conditions = initial_conditions.size()/guiding_center->kDimen();
  }

  //// Record time?
  std::clock_t run_time;
  if(time_flag){
    run_time = std::clock();
  }

  //// Time advance each initial condition
  for (int j = 0; j < n_initial_conditions; ++j){
    integrator->Reset(); // Resets temporary variables in integrators
    double t = 0;
    // Set initial condition
    for (int i=0; i<guiding_center->kDimen(); ++i){
      x[i] = initial_conditions[j*guiding_center->kDimen() + i];
    }
    PrintState(t, x, print_precision); // Print initial position
    // Run standard stepping
    for (int i = 1; i <= n_steps; ++i) {
      integrator->Step(t, x);
      if(!(i%save_nth)){
	PrintState( t, x, print_precision);
      }
    }
  }

  //// Print run time?
  if(time_flag){
    run_time = std::clock() - run_time;
    std::cout << "Run time: " << (double)run_time/CLOCKS_PER_SEC << std::endl;
  }

  //// Clean up
  delete integrator;
  delete guiding_center;
  delete em_fields;
  return 0;
}
