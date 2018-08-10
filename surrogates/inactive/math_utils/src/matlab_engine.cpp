#include "matlab_engine.hpp"

#ifdef ENABLE_MATLAB

MatlabEngine::MatlabEngine(){buffer_ = NULL;};

MatlabEngine::~MatlabEngine(){buffer_ = NULL;};

void MatlabEngine::start()
{
  if ( !( ep_ = engOpen( "\0" ) ) ) {
    throw( std::runtime_error( "Can't start MATLAB engine" ) ) ;
  }
};

void MatlabEngine::stop()
{
  engClose( ep_ );
};

void MatlabEngine::send_vector_double( RealVector &vector, 
				       const std::string& var_name )
{
  /* 
   * Create a variable for our data
   */
  int n = vector.length();
  mxArray *array = mxCreateDoubleMatrix( n, 1, mxREAL );
  memcpy( (void *)mxGetPr( array ), (void *)vector.values(), n * sizeof(Real) );
  engPutVariable( ep_, var_name.c_str(), array );
  
  // Warning: must always delete matrix variables when we are done with 
  // them 
  //mxDestroyArray( array );
  arrays_.push_back( array );
};

void MatlabEngine::send_matrix_double( RealMatrix &matrix, 
				       const std::string& var_name )
{
  /* 
   * Create a variable for our data
   */
  int m = matrix.numRows(), n = matrix.numCols();
  mxArray *array = mxCreateDoubleMatrix( m, n, mxREAL );
  memcpy( (void *)mxGetPr( array ), (void *)matrix.values(), 
	  m * n * sizeof( Real ) );
  //printf("matrix size: %d,%d\n",m,n);
  engPutVariable( ep_, var_name.c_str(), array );
  
  // Warning: must always delete matrix variables when we are done with 
  // them 
  //mxDestroyArray( array );
  arrays_.push_back( array );
};

void MatlabEngine::send_string( const std::string &string,
				const std::string& var_name )
{
   mxArray *str = mxCreateString( string.c_str() );
   engPutVariable( ep_, var_name.c_str(), str );
   // Warning: must always delete matrix variables when we are done with 
   // them 
   // mxDestroyArray( str );
   arrays_.push_back( str );
};

void MatlabEngine::receive_array_double( RealVector &vector, 
					 const std::string& var_name,
					 bool copy )
{
  mxArray *result = NULL;
  if ( ( result = engGetVariable( ep_, var_name.c_str() ) ) == NULL )
    throw( std::runtime_error( "Variable " + var_name + " is not defined." ) ) ;
  else {
    if ( !mxIsNumeric( result ) )
      throw( std::runtime_error( "Variable "+var_name+" is not numeric." ));
    if ( mxGetNumberOfDimensions( result ) != 2 )
      // matlab always returns array with dimension 2 or greater
      throw( std::runtime_error( "Variable "+var_name+" is not a 2d vector." ));
    const mwSize *size_array = mxGetDimensions( result );
    //printf("receive array %s size: %d,%d\n",var_name.c_str(), (int)size_array[0],(int)size_array[1]);
    
    if ( size_array[0] != 1 && size_array[1] != 1)
      throw( std::runtime_error( "Variable "+var_name+
				 " must be a 1d vector." ) );
    if ( mxGetClassID( result ) != mxDOUBLE_CLASS )
      throw( std::runtime_error( "Variable "+var_name+
				 " is not a 1d vector of doubles." ) );
    int n = mxGetNumberOfElements( result );
    if ( copy )
      {
	// Memory will persist for life of vector
	vector.sizeUninitialized( n );
	memcpy( (void*)vector.values(), (void*)mxGetData( result ), 
		n * sizeof( Real ) );
      }
    else
      {
	// Memory will be deleted when matlab deletes var_name.
	RealVector temp( Teuchos::View, (double*)mxGetData( result ), n );
	vector = temp;
      }    
  }
  // Warning must always delete matrix variables when we are done with them
  // mxDestroyArray( result );
};

void MatlabEngine::receive_matrix_double( RealMatrix &matrix, 
					  const std::string& var_name,
					  bool copy )
{
  mxArray *result = NULL;
  if ( ( result = engGetVariable( ep_, var_name.c_str() ) ) == NULL )
    throw( std::runtime_error( "Variable " + var_name + " is not defined." ) ) ;
  else {
    if ( !mxIsNumeric( result ) )
      throw( std::runtime_error( "Variable "+var_name+" is not numeric." ));
    if ( mxGetNumberOfDimensions( result ) != 2 )
      // matlab always returns array with dimension 2 or greater
      throw( std::runtime_error( "Variable "+var_name+" is not a 2d vector." ));
    const mwSize *size_array = mxGetDimensions( result );
    if ( mxGetClassID( result ) != mxDOUBLE_CLASS )
      throw( std::runtime_error( "Variable "+var_name+
				 " is not a 1d vector of doubles." ) );
    
    int m = size_array[0], n = size_array[1];
    if ( copy )
      {
	// Memory will persist for life of vector
	matrix.shapeUninitialized( m, n );
	memcpy( (void*)matrix.values(), (void*)mxGetData( result ), 
		m * n * sizeof( Real ) );
      }
    else
      {
	// Memory will be deleted when matlab deletes var_name.
	RealMatrix temp( Teuchos::View, (double*)mxGetData( result ), m, m, n );
	matrix = temp;
      }    
  }
  // Warning must always delete matrix variables when we are done with them
  // mxDestroyArray( result );
};

  
void MatlabEngine::execute_command( const std::string& command )
{
  engEvalString( ep_, command.c_str() );
};

void MatlabEngine::capture_matlab_buffer( int size )
{
  buffer_ = new char[size+1];
  buffer_[size] = '\0';
  engOutputBuffer( ep_, buffer_, size );
};

void MatlabEngine::print_matlab_buffer()
{
  std::printf( "%s\n", buffer_ );
}

void MatlabEngine::delete_arrays()
{
  for ( int i = 0; i < (int)arrays_.size(); i++ )
    {
      mxDestroyArray( arrays_[i] );
    }
  arrays_.clear();
}
 


#endif // ENABLE_MATLAB
