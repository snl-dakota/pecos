#ifndef MATLAB_ENGINE_HPP
#define MATLAB_ENGINE_HPP

#ifdef ENABLE_MATLAB

#include "LinearAlgebra.hpp"

#include "engine.h"
#include "matrix.h"

class MatlabEngine
{
public:
  Engine *ep_;

  char *buffer_;

  std::vector<mxArray*> arrays_;

  MatlabEngine();

  ~MatlabEngine();

  /**
   * \brief Start the MATLAB engine locally by executing the string
   * "matlab"
   *
   * To start the session on a remote host, use the name of
   * the host as the string rather than \0
   *
   * For more complicated cases, use any string with whitespace,
   * and that string will be executed literally to start MATLAB
   */
  void start();

  /**
   * \brief Stop the MATLAB engine
   */
  void stop();

  /**
   * \brief Pass a vector to the MATLAB engine
   */
  void send_vector_double( RealVector &vector, 
			   const std::string &var_name );

  void send_matrix_double( RealMatrix &matrix, 
			   const std::string &var_name );

  void send_string( const std::string &string,
		    const std::string &var_name );
  
  void receive_array_double( RealVector &result, 
			     const std::string &var_name,
			     bool copy );

  void receive_matrix_double( RealMatrix &result, 
			      const std::string &var_name,
			      bool copy );

  void execute_command( const std::string &command );

  /**
   * \brief Tell c++ to capture the matlab output

   * Tells any subsequent calls to engEvalString to save the first size 
   * characters of output in the character buffer. Must be called before
   * print_matlab_buffer()
   */
  void capture_matlab_buffer( int size );

  void print_matlab_buffer();

  void delete_arrays();
};

#endif // ENABLE_MATLAB

#endif // MATLAB_ENGINE_HPP
