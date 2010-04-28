/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HermiteOrthogPolynomial
//- Description:  Implementation code for HermiteOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "HermiteOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.H"
#endif


namespace Pecos {


const Real& HermiteOrthogPolynomial::
get_value(const Real& x, unsigned short order)
{
  switch (order) {
  case 0:
    basisPolyValue = 1.;
    break;
  case 1:
    basisPolyValue = x;
    break;
  case 2:
    basisPolyValue = x*x - 1.;
    break;
  case 3:
    basisPolyValue = x*(x*x - 3.);
    break;
  case 4: {
    Real x2 = x*x;
    basisPolyValue = x2*x2 - 6.*x2 + 3.;
    break;
  }
  case 5: {
    Real x2 = x*x;
    basisPolyValue = x*(x2*x2 - 10.*x2 + 15.);
    break;
  }
  case 6: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue = x2*x4 - 15.*x4 + 45.*x2 - 15.;
    break;
  }
  case 7: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue = x*(x2*x4 - 21.*x4 + 105.*x2 - 105.);
    break;
  }
  case 8: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue = x4*x4 - 28.*x4*x2 + 210.*x4 - 420.*x2 + 105.;
    break;
  }
  case 9: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue = x*(x4*x4 - 36.*x4*x2 + 378.*x4 - 1260.*x2 + 945.);
    break;
  }
  case 10: {
    Real x2 = x*x, x4 = x2*x2, x8 = x4*x4;
    basisPolyValue = x2*x8 - 45.*x8 + 630.*x2*x4 - 3150.*x4 + 4725.*x2 - 945.;
    break;
  }
  default:
    // Support higher order polynomials using the 3 point recursion formula:
    Real x2 = x*x, x4 = x2*x2, x8 = x4*x4,
      He_n = x2*x8 - 45.*x8 + 630.*x2*x4 - 3150.*x4 + 4725.*x2 - 945., // He_10
      He_nminus1 = x*(x8 - 36.*x4*x2 + 378.*x4 - 1260.*x2 + 945.);     // He_9
    for (size_t i=10; i<order; i++) {
      basisPolyValue = x*He_n - i*He_nminus1; // He_nplus1
      if (i != order-1) {
	He_nminus1 = He_n;
	He_n       = basisPolyValue;
      }
    }
    break;
  }

  return basisPolyValue;
}


const Real& HermiteOrthogPolynomial::
get_gradient(const Real& x, unsigned short order)
{ 
  basisPolyGradient = (order) ? order*get_value(x, order-1): 0;
  return basisPolyGradient;
}


const Real& HermiteOrthogPolynomial::norm_squared(unsigned short order)
{ 
  orthogPolyNormSq = factorial(order);
  return orthogPolyNormSq;
}


const RealArray& HermiteOrthogPolynomial::gauss_points(unsigned short order)
{
  // pull this outside block below since order=0 is initial gauss pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "HermiteOrthogPolynomial::gauss_points()." << std::endl;
    abort_handler(-1);
  }

  if (gaussPoints.size() != order) { // if not already computed
    gaussPoints.resize(order);
#ifdef HAVE_SPARSE_GRID
    if (order <= 20) { // retrieve full precision tabulated values
      webbur::hermite_lookup_points(order, &gaussPoints[0]);
      for (size_t i=0; i<order; i++)
	gaussPoints[i] *= ptFactor; // scale H_n roots by sr2 to get He_n roots
    }
    else { // calculates points/weights together
      if (gaussWeights.size() != order)
	gaussWeights.resize(order);
      webbur::hermite_compute(order, &gaussPoints[0], &gaussWeights[0]);
      for (size_t i=0; i<order; i++) {
	gaussPoints[i]  *= ptFactor; // scale H_n roots by sr2 to get He_n roots
	gaussWeights[i] *= wtFactor; // polynomial weight fn -> PDF
      }
    }
#else
    switch (order) {
    case 1: // zeros of He_1(x) for one Gauss-Hermite point:
      gaussPoints[0] = 0.0;  break;
    case 2: // zeros of He_2(x) for two Gauss-Hermite points:
      gaussPoints[0] = -1.0;
      gaussPoints[1] =  1.0; break;
    case 3: { // zeros of He_3(x) for three Gauss-Hermite points:
      Real sr3 = std::sqrt(3.);
      gaussPoints[0] = -sr3;
      gaussPoints[1] =  0.0;
      gaussPoints[2] =  sr3; break;
    }
    case 4: { // zeros of He_4(x) for four Gauss-Hermite points:
      Real sr3 = std::sqrt(3.), sr6 = std::sqrt(6.),
	sr3psr6 = std::sqrt(3.+sr6), sr3msr6 = std::sqrt(3.-sr6);
      gaussPoints[0] = -sr3psr6;
      gaussPoints[1] = -sr3msr6;
      gaussPoints[2] =  sr3msr6;
      gaussPoints[3] =  sr3psr6; break;
    }
    case 5: { // zeros of He_5(x) for five Gauss-Hermite points:
      Real sr10 = std::sqrt(10.), sr5psr10 = std::sqrt(5.+sr10),
	sr5msr10 = std::sqrt(5.-sr10);
      gaussPoints[0] = -sr5psr10;
      gaussPoints[1] = -sr5msr10;
      gaussPoints[2] =  0.0;
      gaussPoints[3] =  sr5msr10;
      gaussPoints[4] =  sr5psr10; break;
    }
    // tabulated values from Abramowitz & Stegun have limited precision
    case 6:
      gaussPoints[0] = -2.350604973674492 * ptFactor;
      gaussPoints[1] = -1.335849074013697 * ptFactor;
      gaussPoints[2] = -0.436077411927617 * ptFactor;
      gaussPoints[3] = -gaussPoints[2];
      gaussPoints[4] = -gaussPoints[1];
      gaussPoints[5] = -gaussPoints[0]; break;
    case 7:
      gaussPoints[0] = -2.651961356835233 * ptFactor;
      gaussPoints[1] = -1.673551628767471 * ptFactor;
      gaussPoints[2] = -0.816287882858965 * ptFactor;
      gaussPoints[3] =  0.0;
      gaussPoints[4] = -gaussPoints[2];
      gaussPoints[5] = -gaussPoints[1];
      gaussPoints[6] = -gaussPoints[0]; break;
    case 8:
      gaussPoints[0] = -2.930637420257244 * ptFactor;
      gaussPoints[1] = -1.981656756695843 * ptFactor;
      gaussPoints[2] = -1.157193712446780 * ptFactor;
      gaussPoints[3] = -0.381186990207322 * ptFactor;
      gaussPoints[4] = -gaussPoints[3];
      gaussPoints[5] = -gaussPoints[2];
      gaussPoints[6] = -gaussPoints[1];
      gaussPoints[7] = -gaussPoints[0]; break;
    case 9:
      gaussPoints[0] = -3.190993201781528 * ptFactor;
      gaussPoints[1] = -2.266580584531843 * ptFactor;
      gaussPoints[2] = -1.468553289216668 * ptFactor;
      gaussPoints[3] = -0.723551018752838 * ptFactor;
      gaussPoints[4] =  0.0;
      gaussPoints[5] = -gaussPoints[3];
      gaussPoints[6] = -gaussPoints[2];
      gaussPoints[7] = -gaussPoints[1];
      gaussPoints[8] = -gaussPoints[0]; break;
    case 10:
      gaussPoints[0] = -3.436159118837738 * ptFactor;
      gaussPoints[1] = -2.532731674232790 * ptFactor;
      gaussPoints[2] = -1.756683649299882 * ptFactor;
      gaussPoints[3] = -1.036610829789514 * ptFactor;
      gaussPoints[4] = -0.342901327223705 * ptFactor;
      gaussPoints[5] = -gaussPoints[4];
      gaussPoints[6] = -gaussPoints[3];
      gaussPoints[7] = -gaussPoints[2];
      gaussPoints[8] = -gaussPoints[1];
      gaussPoints[9] = -gaussPoints[0]; break;
    default:
      PCerr << "Error: overflow in maximum quadrature order limit (10) in "
	    << "HermiteOrthogPolynomial::gauss_points().  Configure with "
	    << "VPISparseGrid to extend range." << std::endl;
      abort_handler(-1); break;
    }
#endif
  }

  return gaussPoints;
}


const RealArray& HermiteOrthogPolynomial::gauss_weights(unsigned short order)
{
  // pull this outside block below since order=0 is initial gauss pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "HermiteOrthogPolynomial::gauss_weights()." << std::endl;
    abort_handler(-1);
  }

  // The sums of the weights = 1, which is the integral of the density function
  // 1/sqrt(2*PI) exp(-x^2/2) over the support range of [-infinity,+infinity]
  // (the std normal CDF for +infinity).

  if (gaussWeights.size() != order) { // if not already computed
    gaussWeights.resize(order);
#ifdef HAVE_SPARSE_GRID
    if (order <= 20) { // tabulated values from sandia_rules have full precision
      webbur::hermite_lookup_weights(order, &gaussWeights[0]);
      for (size_t i=0; i<order; i++)
	gaussWeights[i] *= wtFactor; // polynomial weight fn -> PDF
    }
    else { // sandia_rules calculates points/weights together
      if (gaussPoints.size() != order)
	gaussPoints.resize(order);
      webbur::hermite_compute(order, &gaussPoints[0], &gaussWeights[0]);
      for (size_t i=0; i<order; i++) {
	gaussPoints[i]  *= ptFactor; // scale H_n roots by sr2 to get He_n roots
	gaussWeights[i] *= wtFactor; // polynomial weight fn -> PDF
      }
    }
#else
    switch (order) {
    case 1: // weights for one Gauss-Hermite point:
      gaussWeights[0] = 1.0; break;
    case 2: // weights for two Gauss-Hermite points:
      gaussWeights[0] = gaussWeights[1] = 0.5; break;
    case 3: // weights for three Gauss-Hermite points:
      gaussWeights[0] = gaussWeights[2] = 1./6.;
      gaussWeights[1] = 2./3.; break;
    case 4: { // weights for four Gauss-Hermite points:
      Real sr6 = std::sqrt(6.);
      gaussWeights[0] = gaussWeights[3] = 1./4./(3.+sr6);
      gaussWeights[1] = gaussWeights[2] = 1./4./(3.-sr6); break;
    }
    case 5: { // weights for five Gauss-Hermite points:
      Real w2sr10 = 2.*std::sqrt(10.);
      gaussWeights[0] = gaussWeights[4] = 3./20./(7.+w2sr10);
      gaussWeights[1] = gaussWeights[3] = 3./20./(7.-w2sr10);
      gaussWeights[2] = 8./15.; break;
    }
    // tabulated values from Abramowitz & Stegun have limited precision
    case 6:
      gaussWeights[0] = gaussWeights[5] = 4.530009905509e-3 * wtFactor;
      gaussWeights[1] = gaussWeights[4] = 0.1570673203229 * wtFactor;
      gaussWeights[2] = gaussWeights[3] = 0.7246295952244 * wtFactor; break;
    case 7:
      gaussWeights[0] = gaussWeights[6] = 9.717812450995e-4 * wtFactor;
      gaussWeights[1] = gaussWeights[5] = 5.451558281913e-2 * wtFactor;
      gaussWeights[2] = gaussWeights[4] = 0.4256072526101 * wtFactor;
      gaussWeights[3] = 0.8102646175568 * wtFactor; break;
    case 8:
      gaussWeights[0] = gaussWeights[7] = 1.996040722114e-4 * wtFactor;
      gaussWeights[1] = gaussWeights[6] = 1.707798300741e-2 * wtFactor;
      gaussWeights[2] = gaussWeights[5] = 0.2078023258149 * wtFactor;
      gaussWeights[3] = gaussWeights[4] = 0.6611470125582 * wtFactor; break;
    case 9:
      gaussWeights[0] = gaussWeights[8] = 3.960697726326e-5 * wtFactor;
      gaussWeights[1] = gaussWeights[7] = 4.943624275537e-3 * wtFactor;
      gaussWeights[2] = gaussWeights[6] = 8.847452739438e-2 * wtFactor;
      gaussWeights[3] = gaussWeights[5] = 0.4326515590026 * wtFactor;
      gaussWeights[4] = 0.7202352156061 * wtFactor; break;
    case 10:
      gaussWeights[0] = gaussWeights[9] = 7.640432855233e-6 * wtFactor;
      gaussWeights[1] = gaussWeights[8] = 1.343645746781e-3 * wtFactor;
      gaussWeights[2] = gaussWeights[7] = 3.387439445548e-2 * wtFactor;
      gaussWeights[3] = gaussWeights[6] = 0.2401386110823 * wtFactor;
      gaussWeights[4] = gaussWeights[5] = 0.6108626337353 * wtFactor; break;
    default:
      // define Gauss wts from Gauss pts using
      // -(A_{n+1} gamma_n)/(A_n Phi_n'(x_i) Phi_{n+1}(x_i)),
      // which for He(x), is n!/(n He_{n-1}(x_i))^2.
      const RealArray& gauss_pts = gauss_points(order);
      for (size_t i=0; i<order; i++)
	gaussWeights[i]
	  = factorial(order)/std::pow(order*get_value(gauss_pts[i], order-1),2);
      break;
    }
#endif
  }

  return gaussWeights;
}

} // namespace Pecos
