#include "ParallelObject.hpp"

#ifdef ENABLE_LIBHEAT_MPI

namespace Pecos {

void send( RealVector &M, int dest, MPI_Comm mpi_comm )
{
  int size( M.length() );

  // Send array size
  MPI_Send( &size, 1, MPI_INT, dest, TAG_SENDING_MATRIX, MPI_COMM_WORLD );

  // Send array
  Real *data = M.values();
  MPI_Send( data,               // message buffer
	    size,               // number of items
	    MPI_DOUBLE,         // items are of type Real
	    dest,               // destination process rank
	    TAG_SENDING_MATRIX, // user defined message tag
	    mpi_comm );         // communicator
};

void receive( RealVector &M, int source, MPI_Comm mpi_comm )
{
  int size;
  MPI_Status status;

  // Receive array size
  MPI_Recv( &size, 1, MPI_INT, source, TAG_SENDING_MATRIX, mpi_comm, 
	    &status );

  // Allocate memory for message buffer
  resize( M, size );
  Real *data = M.values();

  // Receive array
  MPI_Recv( data,               // message buffer
	    size,               // number of items
	    MPI_DOUBLE,         // items are of type Real
	    source,             // receive from source
	    TAG_SENDING_MATRIX, // user defined message tag
	    mpi_comm,           // communicator
	    &status );          // info about the received message
};

void send( RealMatrix &M, int dest, MPI_Comm mpi_comm )
{
  int *size;
  size = new int [2];
  size[0] = M.numRows(); size[1] = M.numCols();
  int len = size[0]*size[1];

  Real *data = M.values();
  // Send array size
  MPI_Send( size, 2, MPI_INT, dest, TAG_SENDING_MATRIX, mpi_comm );

  // Send array
  MPI_Send( data,                // message buffer
	    len,                 // number of items
	    MPI_DOUBLE,          // items are of type Real
	    dest,                // destination process rank
	    TAG_SENDING_MATRIX,  // user defined message tag
	    mpi_comm );          // communicator

  delete [] size;
};

/**
 * Recieve an entire array from another process
 */
void receive( RealMatrix &M, int source, MPI_Comm mpi_comm )
{
  // Receive array size
  MPI_Status status;
  int *size = new int [2];
  MPI_Recv( size, 2, MPI_INT, source, TAG_SENDING_MATRIX, mpi_comm, 
	    &status );

  // Allocate memory for message buffer
  reshape( M, size[0], size[1] );
  Real *data = M.values();

  // Receive array
  MPI_Recv( data,               // message buffer
	    size[0]*size[1],    // number of items
	    MPI_DOUBLE,         // items are of type Real
	    source,             // receive from source
	    TAG_SENDING_MATRIX, // user defined message tag
	    mpi_comm,           // communicator
	    &status );          // info about the received message

  delete [] size;
};

void send( const std::string &s, int destination, MPI_Comm mpi_comm )
{
  int size = (int)s.size();
  MPI_Send( &size, 1, MPI_INT, destination, TAG_SENDING_SCALAR, 
	    mpi_comm );
  MPI_Send( (void*)s.c_str(), size, MPI_CHAR, destination, TAG_SENDING_STRING, 
	    mpi_comm );
};

void receive( std::string &s, int source, MPI_Comm mpi_comm )
{
  int size;
  MPI_Status status;
  MPI_Recv( &size, 1, MPI_INT, source, TAG_SENDING_SCALAR, 
	    mpi_comm, &status );
  char *buffer = new char [size];
  MPI_Recv( buffer, size, MPI_CHAR, source, TAG_SENDING_STRING, 
	    mpi_comm, &status );
  s.assign( buffer, size );
  delete [] buffer;
};

void send( RealMatrixList &list, int destination, MPI_Comm mpi_comm )
{
  // Send list length
  int len( list.size() );
  MPI_Send( &len, 1, MPI_INT, destination, TAG_SENDING_MATRIX_LIST, mpi_comm );

  for ( int i = 0; i < len; i++ )
    {
      send( list[i], destination, mpi_comm );
    }
};

void receive( RealMatrixList &list, int source, MPI_Comm mpi_comm )
{
  // Send list length
  int len;
  MPI_Status status;
  MPI_Recv( &len, 1, MPI_INT, source, TAG_SENDING_MATRIX_LIST, 
	    mpi_comm, &status );

  list.resize( len );
  for ( int i = 0; i < len; i++ )
    {
      receive( list[i], source, mpi_comm );
    }
};

} // namespace Pecos

#endif // ENABLE_LIBHEAT_MPI
