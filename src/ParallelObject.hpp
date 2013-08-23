#ifndef RUNTIME_ENVIRONMENT_HPP
#define RUNTIME_ENVIRONMENT_HPP

namespace Pecos {

#ifndef ENABLE_LIBHEAT_MPI

typedef int MPI_Comm;
#define MPI_COMM_WORLD   0
#define MPI_COMM_NULL    0

#else

#include <mpi.h>

const int TAG_SENDING_SCALAR = 0, TAG_SENDING_ARRAY = 1, TAG_SENDING_STRING = 2, 
  TAG_SENDING_MATRIX = 3, TAG_SENDING_MATRIX_LIST = 3;

void send( RealVector& M, int dest, MPI_Comm mpi_comm );

void receive( RealVector& M, int source, MPI_Comm mpi_comm );

void send( RealMatrix &M, int dest, MPI_Comm mpi_comm );

/**
 * Recieve an entire array from another process
 */
void receive( RealMatrix &M, int source, MPI_Comm mpi_comm );

void send( const std::string &s, int destination, MPI_Comm mpi_comm );

void receive( std::string &s, int source, MPI_Comm mpi_comm );

void send( RealMatrixArray &list, int destination, MPI_Comm mpi_comm );

void receive( RealMatrixArray &list, int source, MPI_Comm mpi_comm );

#endif // ENABLE_LIBHEAT_MPI

class ParallelObject
{
 
protected:
  /// The unqiue identifier of the processor
  int processorId_;

  /// The unqiue identifier of the master processor
  int masterProcessorId_;

  /// The number of processors in MPICommunicator_
  int numProcessors_;

  /// MPI Communicator on which ParallelObject is running
  MPI_Comm MPICommunicator_;

  /// Defines the verbosity level of i/o
  int verbosity_;

public:
  /// Default constructor
  ParallelObject() :
    processorId_( 0 ), masterProcessorId_( 0 ), numProcessors_( 1 ),
    MPICommunicator_( MPI_COMM_NULL )
  {};

 ~ParallelObject()
  {};

  /// Set the MPI communicator
  void mpi_communicator( MPI_Comm mpi_comm )
  {
    MPICommunicator_ =  mpi_comm;
#ifdef ENABLE_LIBHEAT_MPI
    // Get total number of processes
    MPI_Comm_size( MPICommunicator_, &numProcessors_ ); 
    
    // Get process number for this process
    MPI_Comm_rank( MPICommunicator_, &processorId_ );
#endif
  };

  int num_processors()
  {
    return numProcessors_;
  };

  /// Set the verbosity level of i/o
  void verbosity( int verbosity_in )
  {
    verbosity_ = verbosity_in;
  }
  
  /// Determine if the current processor is the master
  bool is_master()
  {
    return ( processorId_ == masterProcessorId_ );
  };

  /// Return the id of the master processor
  bool master_processor_id()
  {
    return masterProcessorId_;
  };

  /// Return the unqiue identifier of the current processor. 
  int processor_id() 
  {
    return processorId_;
  };

};

} // namespace Pecos 

#endif //RUNTIME_ENVIRONMENT_HPP
