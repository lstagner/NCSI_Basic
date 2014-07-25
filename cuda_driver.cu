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
#include <ctime> // For timing run time
#include "input_parser.h"
#include "guiding_center.h"
#include "em_fields.h"
#include "axisymmetric_tokamak.h"
#include "integrator.h"
#include "runge-kutta4.h"
#include "noncanonical_symplectic.h"
#include "eigen_types.h" // Typedef for Vector4, etc.
#include "cuda_error.h" // HANDLE_ERROR macro


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


// Kernel for time advance
// template <class I>
__global__ void step_positions(Vector4 *x, double t, const double kdt,  
			       const int kNSteps, const int kNParticles, 
			       const double kB0, const double kR0,
			       const double kMu){
  // Thread identification
  int idx=blockIdx.x*blockDim.x + threadIdx.x;

  // Integrator initialization
  AxisymmetricTokamak em_fields(kB0, kR0);
  GuidingCenter model((EMFields *) &em_fields, kMu);
  RungeKutta4 integrator(kdt, model);

  // Time advance
  for(int i=0; i<kNSteps; ++i){
    if (idx < kNParticles){
      integrator.Step(t, x[idx]);
    }
    __syncthreads(); // Likely not necessary, but doesn't slow down
  }
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

  //// Retrieve the values from the input parser
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
  const int kNParticles = input_parser.GetValue<int>("n_particles"); 
  const int kBlockSize = input_parser.GetValue<int>("block_size");  

  //// Check for valid parameters
  // Integrator
  if (kIntegratorName.compare("rk4") && kIntegratorName.compare("ncsi") ){
    std::cout << "Unrecognized integrator. Try rk4 or ncsi" << std::endl;
    return 1;
  }
  // Initial conditions
  if (initial_conditions.size() % (GC_DIM) ){
    std:: cout << "Wrong number of initial conditions. Use 4." << std::endl;
    return 1;
  }

  //// Set up CUDA parameters
  //       Should  check for cuda-capable device, here
  cudaSetDevice(1); // Some people hop on first device while running use
  dim3 dimBlock(kBlockSize);
  dim3 dimGrid(ceil(kNParticles/(float)kBlockSize));
  
  //// Data set up
  // Host side
  Vector4 x_host[kNParticles];
  for(int i=0; i<kNParticles; ++i){
    for(int j=0; j<GC_DIM; ++j){
      x_host[i][j] = initial_conditions[j];
    }
  }
  double t = 0.0;
  // Device side
  Vector4 *x_device;
  HANDLE_ERROR( cudaMalloc((void **)&x_device, 
			   kNParticles*sizeof(Vector4)) );
  HANDLE_ERROR( cudaMemcpy(x_device, x_host, kNParticles*sizeof(Vector4), 
			  cudaMemcpyHostToDevice) );

  //// Record time?
  std::clock_t run_time;
  if(kTimeFlag){
    run_time = std::clock();
  }

  // Print initial state
  PrintState(t, x_host[0], kPrintPrecision);

  //// Time advance
  for(int i=0; i<kNSteps/kSaveNth; ++i){
    // step_positions<RungeKutta4><<<dimGrid, dimBlock>>>(x_device, t, kdt,
    // 						      kSaveNth, kNParticles,
    // 						      kB0, kR0, kMu);
    step_positions<<<dimGrid, dimBlock>>>(x_device, t, kdt, kSaveNth, 
					  kNParticles, kB0, kR0, kMu);

    HANDLE_ERROR( cudaGetLastError() ); // Check for kernel errors
    t += kSaveNth*kdt; // Advance host time
    // Advance host positions
    HANDLE_ERROR( cudaMemcpy(x_host, x_device, kNParticles*sizeof(Vector4),
			     cudaMemcpyDeviceToHost) );
    // Print state
    PrintState(t, x_host[0], kPrintPrecision);
  }

  //// Print run time?
  if(kTimeFlag){
    run_time = std::clock() - run_time;
    std::cout << "Run time: " << (double)run_time/CLOCKS_PER_SEC 
	      << std::endl;
  }

  //// Clean up
  HANDLE_ERROR( cudaFree(x_device) );
  return 0;
}
