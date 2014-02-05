/*!
*******************************************************************************
* \file input_parser.cc
*
* \brief Implementation of InputParser class. Sets runtime options.
*
* \date Jan 2014
* \author C. Leland Ellison
* \version 0.1
*******************************************************************************
*/

#include "input_parser.h"

/*!
 * Read command line arguments and notify variables map of runtime options.
 *
 * @param[in] argc Number of input arguments
 * @param[in] argv Input arguments
 * @returns Zero for success
 */
int InputParser::ReadInput(int argc, char** argv){
 
  //// Initialize different options lists ////
  // Manditory Options - Failure if not specified
  po::options_description manditory_options("Manditory Options");
  manditory_options.add_options()
    ("integrator,I", po::value<std::string>()->required(), 
     "Integration method")
    ("dt,h", po::value<double>()->required(), "Numerical step size dt")
    ("n_steps,N", po::value<int>()->required(), "Number of time steps")
  ;

  // Descriptive Options - Helps the user 
  po::options_description descriptive_options("Descriptive Options");
  descriptive_options.add_options()
    ("help", "Produce help message")
  ;

  // Declare hidden options - Not manditory, but helpful on command line
  po::options_description hidden_options("Hidden Options");
  hidden_options.add_options()
    ("input_file", po::value<std::string>(), "Input file defining run options")
    ("save_nth", po::value<int>()->default_value(1),
     "Save every save_nth step")
    ("initial_conditions,x", po::value<std::vector<double> >()
     ->multitoken()->default_value(std::vector<double>(),""),
     "Initial conditions")
  ;

  // Declare output options
  po::options_description output_options("Output Options");
  output_options.add_options()
    ("time,T", po::value<bool>()->default_value(false), 
     "Print total runtime")
    ("precision,P", po::value<int>()->default_value(8),
     "Number of printed digits")
  ;

  // Declare model options
  po::options_description model_options("Model Options");
  model_options.add_options()
    ("mu", po::value<double>()->default_value(0.21), "mu parameter")
    ("b0", po::value<double>()->default_value(1.0), "magnetic field strength")
    ("r0", po::value<double>()->default_value(100.0), "radius scale")
  ;

  // Declare integrator options
  po::options_description integrator_options("Integrator Options");
  integrator_options.add_options()
    ("tol", po::value<double>()->default_value(1e-12), 
     "Nonlinear solve tolerance")
    ("max_iter", po::value<int>()->default_value(15),
     "Maximum nonlinear solve iterations")
  ;

  // Declare all options
  po::options_description all_options("All Options");
  all_options.add(descriptive_options).add(manditory_options).
    add(hidden_options).add(output_options).add(model_options).
    add(integrator_options);

  // Positional option for input file
  po::positional_options_description positional_options;
  positional_options.add("input_file", -1);
  //// END OPTIONS SPECIFICATION

  //// Set up/notify variables map
  po::store(po::command_line_parser(argc, argv).
	    options(all_options).positional(positional_options).run(),
	    variables_map_);

  if (variables_map_.count("help")){
    std::cout << descriptive_options << std::endl
	      << manditory_options << std::endl
              << model_options << std::endl
	      << integrator_options << std::endl
	      << output_options << std::endl
	      << "USAGE:" << std::endl
   	      << "driver -I <integrator> -h <double> -N <int>"
   	      << std::endl << "-OR-" << std::endl 
   	      << "driver <input_file>" << std::endl << std::endl; 
    return 1;
  }  

  // Read in options from input file, if necessary
  std::string input_file;
  // If there was an input file, read the options in there.
  if (variables_map_.count("input_file")) {
    input_file = variables_map_["input_file"].as<std::string>();
    std::ifstream ifs(input_file.c_str());
    if (!ifs){
      std::cout << "Cannot open input file: " << input_file << std::endl;
      return 1;
    }
    else{
      po::store(po::parse_config_file(ifs, all_options), variables_map_);
    }
  }
  po::notify(variables_map_);
  //// End variables map setup
  
  return 0; // Success
}
