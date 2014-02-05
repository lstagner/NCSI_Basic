/*!
*******************************************************************************
* \file input_parser.h
*
* \brief Interface for input parser class. Recieves command line input from driver and outputs essential information for running the driver.
*
* \date Jan 2014
* \author C. Leland Ellison
* \version 0.1
*
*******************************************************************************
*/

#ifndef INPUT_PARSER_H_
#define INPUT_PARSER_H_

#include <stdlib.h>
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

class InputParser{
 public:
  int ReadInput(int argc, char** argv);
  //// Accessors
  //! Retrieve the value at key and output it to value_out 
  template <typename T>
  int GetValue(const char* key, T &value_out) const{
    value_out = variables_map_[key].as<T>();
    return 0;
  }

 private:
  po::variables_map variables_map_; //!< Program Options variables map
};

#endif //INPUT_PARSER_H_
