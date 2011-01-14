/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        PiecewiseInterpPolynomial
//- Description:  Implementation code for PiecewiseInterpPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "PiecewiseInterpPolynomial.hpp"

namespace Pecos {


//void PiecewiseInterpPolynomial::precompute_data()
//{
//}


/** Compute value of Piecewise polynomial for interpolation point i. */
const Real& PiecewiseInterpPolynomial::
get_value(const Real& x, unsigned short i)
{
  basisPolyValue = 1.;
  //for (size_t j=0; j<numInterpPts; j++)
  //  ;
  return basisPolyValue;
}


/** Compute derivative with respect to x of Piecewise polynomial for
    interpolation point i. */
const Real& PiecewiseInterpPolynomial::
get_gradient(const Real& x, unsigned short i)
{ 
  basisPolyGradient = 1.;
  //for (size_t j=0; j<numInterpPts; j++)
  //  ;
  return basisPolyGradient;
}

} // namespace Pecos
