#ifndef SERIALIZATION_HPP
#define SERIALIZATION_HPP

#ifdef ENABLE_SERIALIZATION

#include "LinearAlgebra.hpp"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/shared_ptr.hpp>

namespace boost { 
  namespace serialization {
    template<class Archive, typename O, typename T >
    void save( Archive & ar, const Teuchos::SerialDenseMatrix<O,T>& matrix, 
	       const unsigned int version )
    {
      int M = matrix.numRows(), N = matrix.numCols();
      ar & M; ar & N;
      for ( int n = 0 ; n < N; n++ )
	{
	  for ( int m = 0 ; m < M; m++ )
	    ar & matrix(m,n);
	}
    } 

    template<class Archive, typename O, typename T >
    void load( Archive & ar, Teuchos::SerialDenseMatrix<O,T>& matrix, 
	       const unsigned int version )
    {
      int M = 0, N = 0;
      ar >> M; ar >> N;
      matrix.shapeUninitialized( M, N );
      for ( int n = 0 ; n < N; n++ )
	{
	  for ( int m = 0 ; m < M; m++ )
	    ar >> matrix(m,n);
	}
    } 

    template<class Archive, typename O, typename T>
    inline void serialize( Archive & ar, 
			   Teuchos::SerialDenseMatrix<O,T> & matrix,
			   const unsigned int version)
    {
      split_free( ar, matrix, version ); 
    }

    template<class Archive, typename O, typename T>
    inline void serialize( Archive & ar, 
			   Teuchos::SerialDenseVector<O,T> & vector,
			   const unsigned int version)
    {
      ar & boost::serialization::base_object< Teuchos::SerialDenseMatrix<O,T> >( vector );
    }

  }
}

template<typename O,typename T>
void save_matrix( Teuchos::SerialDenseMatrix<O,T>& matrix, 
		  const std::string& filename )

{
  //std::ofstream ofs(filename.c_str());
  //boost::archive::text_oarchive oa(ofs);
  std::ofstream ofs(filename.c_str(), std::ios::binary );
  boost::archive::binary_oarchive oa(ofs,boost::archive::no_header);
  oa << matrix;
}

// can only have one of either boost::archive::text_oarchive or boost::archive::text_iarchive created at a time. By creating two functions load and save matrix
// the archives become out of scope and are deleted soon as the funcion is exited
template<typename O,typename T>
void load_matrix( Teuchos::SerialDenseMatrix<O,T>& result, 
		  const std::string& filename )
{
  //std::ifstream ifs(filename.c_str());
  //boost::archive::text_iarchive ia(ifs);
  std::ifstream ifs(filename.c_str(), std::ios::binary );
  boost::archive::binary_iarchive ia(ifs, boost::archive::no_header);
  // restore the matrix from the archive
  ia >> result;
}

template<typename T>
void save_grid( T& grid, const std::string& filename )

{
  std::ofstream ofs(filename.c_str(), std::ios::binary );
  boost::archive::binary_oarchive oa(ofs,boost::archive::no_header);
  // store the grid to the archive
  oa << const_cast<T&>(grid);
}

template<typename T>
void load_grid( T& grid, const std::string& filename )
{
  std::ifstream ifs(filename.c_str(), std::ios::binary );
  boost::archive::binary_iarchive ia(ifs, boost::archive::no_header);
  // restore the grid from the archive
  ia >> grid;
}

template<typename T>
void save_pce( T& pce, const std::string& filename, int archive_type )

{
  std::ios_base::openmode mode;
  if ( archive_type == 0 )
    {
      std::ofstream ofs( filename.c_str(), std::ios::binary );
      boost::archive::binary_oarchive oa( ofs, boost::archive::no_header );
      // store the pce to the archive
      oa << const_cast<T&>(pce);
    }
  else
    {
      std::ofstream ofs( filename.c_str() );
      boost::archive::text_oarchive oa( ofs, boost::archive::no_header );
      // store the pce to the archive
      oa << const_cast<T&>(pce);
    }
}

template<typename T>
void load_pce( T& pce, const std::string& filename, int archive_type )
{
 if ( archive_type == 0 )
    {
      std::ifstream ifs( filename.c_str(), std::ios::binary );
      boost::archive::binary_iarchive ia( ifs, boost::archive::no_header );
      // restore the pce from the archive
      ia >> pce;
    }
 else
   {
     std::ifstream ifs( filename.c_str() );
     boost::archive::text_iarchive ia( ifs, boost::archive::no_header );
     // restore the pce from the archive
     ia >> pce;
   }
}

#endif // ENABLE_SERIALIZATION

#endif // SERIALIZATION_HPP
