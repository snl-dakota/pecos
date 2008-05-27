/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- File:        Xform tester
//- Description: Fires assert if any Transformation tests fail
//- Checked by:
//- Version:
//- $Id$

// WJB: wait for MikesLATEST vs. HACK something?? #include "Transformation.H"
#include "Teuchos_SerialSymDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

#include "fftw3.h"
#include <cassert>
#include <iostream>
#include <string>

//using namespace Pecos;


int main(int, char** argv)
{
//std::cout << "Using " << FFTW_CONCAT(fftw_, version) << std::endl;

  assert(   !strcmp(FFTW_CONCAT(fftw_, version), "fftw-3.1.2")
         && "Incorrect version of FFTW" );

  assert(0 == 1 && "helloPecosWorld");  // WJB ToDo: enable Teuchos testing

  // Creating a double-precision matrix can be done in several ways:
  // Create an empty matrix with no dimension
  Teuchos::SerialSymDenseMatrix<int,double> Empty_Matrix;
  // Create an empty 4x4 matrix
  Teuchos::SerialSymDenseMatrix<int,double> My_Matrix( 4 );
  // Basic copy of My_Matrix
  Teuchos::SerialSymDenseMatrix<int,double> My_Copy1( My_Matrix ),
    // (Deep) Copy of principle 3x3 submatrix of My_Matrix
    My_Copy2( Teuchos::Copy, My_Matrix, 3 ),
    // (Shallow) Copy of 3x3 submatrix of My_Matrix
    My_Copy3( Teuchos::View, My_Matrix, 3, 1 );

  // The matrix dimensions and strided storage information can be obtained:
  int rows, cols, stride;
  rows = My_Copy3.numRows();  // number of rows
  cols = My_Copy3.numCols();  // number of columns
  stride = My_Copy3.stride(); // storage stride

  // Matrices can change dimension:
  Empty_Matrix.shape( 3 );      // size non-dimensional matrices
  My_Matrix.reshape( 3 );       // resize matrices and save values

  // Filling matrices with numbers can be done in several ways:
  My_Matrix.random();             // random numbers
  My_Copy1.putScalar( 1.0 );      // every entry is 1.0
  My_Copy2(1,1) = 10.0;           // individual element access
  Empty_Matrix = My_Matrix;       // copy My_Matrix to Empty_Matrix

  // Basic matrix arithmetic can be performed:
  Teuchos::SerialDenseMatrix<int,double> My_Prod( 4, 3 ), My_GenMatrix( 4, 3 );
  My_GenMatrix.putScalar(1.0);
  // Matrix multiplication ( My_Prod = 1.0*My_GenMatrix*My_Matrix )
  My_Prod.multiply( Teuchos::RIGHT_SIDE, 1.0, My_Matrix, My_GenMatrix, 0.0 );
  My_Copy2 += My_Matrix;   // Matrix addition
  My_Copy2 *= 0.5;         // Matrix scaling

  // Matrices can be compared:
  // Check if the matrices are equal in dimension and values
  if (Empty_Matrix == My_Matrix) {
    std::cout<< "The matrices are the same!" <<std::endl;
  }
  // Check if the matrices are different in dimension or values
  if (My_Copy2 != My_Matrix) {
    std::cout<< "The matrices are different!" <<std::endl;
  }

  return 0;
}

