/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Description:  Functions for solving undertermined linear systems using 
// compressed sensing

#ifndef COMPRESSED_SENSING_HPP
#define COMPRESSED_SENSING_HPP

#include "pecos_data_types.hpp"

namespace Pecos {

class CompressedSensing
{
public:

  /**
   * @breif Default constructor. This class contains no member variables
   * so nothing is done here.
   */
  CompressedSensing(){};

  /**
   *@brief Deconstructor.  This class contains no member variables
   * so nothing is done here.
   */
  ~CompressedSensing(){};

  /**
   * @brief Write a real matrix to file.
   */
  void write_matrix_to_file( RealMatrix& M, std::string filename );

  /**
   * @brief Perform a deep copy of data in an array to a matrix
   */
  void copy_data_to_matrix( Real *A_matrix, int m, int n, RealMatrix& A );

  /**
   * @brief Perform a deep copy of data in a matrix to an array
   */
  void copy_data_from_matrix( RealMatrix& A, Real *A_matrix  );

  /**
   * @brief Find the argument max of a vector
   */
  int ArgMaxMagnitude(RealVector& v);
  

  /**
   * @brief Compute the dot product of two vectors
   */
  Real DotProduct ( RealVector &v1, RealVector& v2 );

  /**
   * @brief Compute a matrix vector product. That is
   * y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
   */
  void Product(RealMatrix& A, RealVector& x, RealVector& y,
	       Teuchos::ETransp trans, Real alpha, Real beta);

  /**
   * @brief Solve AX=B using the singular value decomposition.
   * The matrix A can be square or rectangular.
   * under or over determined systems are both permissible.
   */
  int SVDSolve( RealMatrix& A, RealMatrix& B,
		RealMatrix& X );  

  /**
   * @brief Solve AX = B  using the Cholesky factorization,
   * where A is a symmetric positive definite matrix
   */
  int CholeskySolve( RealMatrix& A, RealMatrix& B,
		     RealMatrix& X );

  /**
   * @brief Compute the duality gap of the primal-dual interior point 
   * optimization algorithm. This is called by the function BasisPursuit().
   */
  double GetSurrogateDualityGapAndSlackness( RealMatrix& primal_residual,
					     RealMatrix& fu1, 
					     RealMatrix& fu2, 
					     RealMatrix& lamu1, 
					     RealMatrix& lamu2, 
					     RealMatrix& AtV,
					     double mu,
					     double tol,
					     RealVector& sdg,
					     RealVector& tau,
					     RealVector& slackness_norm,
					     IntVector& complete);

  /**
   * @brief Compute a sparse solution to AX=B using the Basis Pursuit (BP) 
   * algorithm
   */
  void BasisPursuit( RealMatrix& A, RealMatrix& B, RealMatrix& X );

  /* @brief Use the modified gram-schmidt algorithm to add a new column to a
   * matrix A and ensure all columns remain orthogonal.
   * See A Wavelet Tour of Signal Processing By Stephane Mallat p 428
   */
  void UpdateOrthogonalMatrix(RealMatrix& Q, int n, RealVector& new_column);

  /**
   * @brief Compute a solution to AX=B using the OrthogonalMatchingPursuit (OMP) 
   * algorithm
   */
  void OrthogonalMatchingPursuit( RealMatrix& A, RealMatrix& B, RealMatrix& X );

};// class CompressedSensing

} // namespace Pecos

#endif
