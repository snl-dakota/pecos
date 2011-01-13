/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LegendreOrthogPolynomial
//- Description:  Implementation code for LegendreOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "LegendreOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.H"
#endif

//#define DEBUG


namespace Pecos {

const Real& LegendreOrthogPolynomial::
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
    basisPolyValue = (3.*x*x - 1.)/2.;
    break;
  case 3:
    basisPolyValue = x*(5.*x*x - 3.)/2.;
    break;
  case 4: {
    Real x2 = x*x;
    basisPolyValue = (35.*x2*x2 - 30.*x2 + 3.)/8.;
    break;
  }
  case 5: {
    Real x2 = x*x;
    basisPolyValue = x*(63.*x2*x2 - 70.*x2 + 15.)/8.;
    break;
  }
  case 6: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue = (231.*x2*x4 - 315.*x4 + 105.*x2 - 5.)/16.;
    break;
  }
  case 7: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue = x*(429.*x2*x4 - 693.*x4 + 315.*x2 - 35.)/16.;
    break;
  }
  case 8: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue
      = (6435.*x4*x4 - 12012.*x2*x4 + 6930.*x4 - 1260.*x2 + 35.)/128.;
    break;
  }
  case 9: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue
      = x*(12155.*x4*x4 - 25740.*x2*x4 + 18018.*x4 - 4620.*x2 + 315.)/128.;
    break;
  }
  case 10: {
    Real x2 = x*x, x4 = x2*x2, x8 = x4*x4;
    basisPolyValue  = (46189.*x2*x8 - 109395.*x8 + 90090.*x2*x4 - 30030.*x4 +
			3465.*x2 - 63.)/256.;
    break;
  }
  default:
    // Support higher order polynomials using the 3 point recursion formula:
    Real x2 = x*x, x4 = x2*x2, x8 = x4*x4,
      P_n = (46189.*x2*x8 - 109395.*x8 + 90090.*x2*x4 - 30030.*x4 + 3465.*x2 -
	     63.)/256., // P_10
      P_nminus1 = x*(12155.*x8 - 25740.*x2*x4 + 18018.*x4 - 4620.*x2 + 315.)
                / 128.; // P_9
    for (size_t i=10; i<order; i++) {
      basisPolyValue = ( (2.*i+1.)*x*P_n - i*P_nminus1 ) / (i+1.); // P_nplus1
      if (i != order-1) {
	P_nminus1 = P_n;
	P_n       = basisPolyValue;
      }
    }
    break;
  }

  return basisPolyValue;
}


const Real& LegendreOrthogPolynomial::
get_gradient(const Real& x, unsigned short order)
{
#ifdef DEBUG
  // See Abramowitz & Stegun, Section 22.8, p.783
  //basisPolyGradient = (order) ?
  //  order*(x*get_value(x, order) - get_value(x, order-1))/(x*x - 1.) : 0.;
  if (order) { // be careful with reference to changing basisPolyValue
    Real P_n = get_value(x, order), P_nminus1 = get_value(x, order-1);
    basisPolyGradient = order*(x*P_n - P_nminus1)/(x*x - 1.);
  }
  else
    basisPolyGradient = 0.;
  PCout << "Legendre gradient approach 1: " << basisPolyGradient << '\n';
#endif // DEBUG

  // The previous approach, while very compact, produces 0/0 = NaN at x = +/-1.
  // It is therefore only used for cross-validation purposes.  To avoid NaN
  // issues at bounds, differentiate the 3 pt value recursion to get a 3 point
  // gradient recursion
  switch (order) {
  case 0:
    basisPolyGradient = 0.;
    break;
  case 1:
    basisPolyGradient = 1;
    break;
  case 2:
    basisPolyGradient = 3.*x;
    break;
  case 3:
    basisPolyGradient = (15.*x*x - 3.)/2.;
    break;
  case 4:
    basisPolyGradient = x*(35.*x*x - 15.)/2.;
    break;
  case 5: {
    Real x2 = x*x;
    basisPolyGradient = (315.*x2*x2 - 210.*x2 + 15.)/8.;
    break;
  }
  case 6: {
    Real x2 = x*x;
    basisPolyGradient = x*(693.*x2*x2 - 630.*x2 + 105.)/8.;
    break;
  }
  default:
    // Support higher order polynomials using a 3 point recursion formula:
    Real x2 = x*x, x4 = x2*x2,
      dPdx_n       = x*(693.*x4 - 630.*x2 + 105.)/8., // P'_6
      dPdx_nminus1 = (315.*x4 - 210.*x2 + 15.)/8.;    // P'_5
    for (size_t i=6; i<order; i++) {
      basisPolyGradient // dPdx_nplus1
	= ( (2.*i+1.)*(x*dPdx_n + get_value(x,i)) - i*dPdx_nminus1 ) / (i+1.);
      if (i != order-1) {
	dPdx_nminus1 = dPdx_n;
	dPdx_n       = basisPolyGradient;
      }
    }
    break;
  }
#ifdef DEBUG
  PCout << "Legendre gradient approach 2: " << basisPolyGradient << '\n';
#endif // DEBUG

  return basisPolyGradient;
}


const Real& LegendreOrthogPolynomial::norm_squared(unsigned short order)
{
  // Abramowitz & Stegun: w(x) = 1
  //orthogPolyNormSq = 2./(2.*order + 1.);

  // sampling density f(x) = 1/(U-L) = 1/2 for [L,U] = [-1,1]
  orthogPolyNormSq = 1./(2.*order + 1.);

  return orthogPolyNormSq;
}


const RealArray& LegendreOrthogPolynomial::gauss_points(unsigned short order)
{
  // pull this outside block below since order=0 is initial gauss pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "LegendreOrthogPolynomial::gauss_points()." << std::endl;
    abort_handler(-1);
  }

  bool mode_err = false;
  if (gaussPoints.size() != order) { // if not already computed
    gaussPoints.resize(order);
#ifdef HAVE_SPARSE_GRID
    if (gaussMode == GAUSS_PATTERSON)
      webbur::patterson_lookup_points(order, &gaussPoints[0]);
    else if (gaussMode == GAUSS_LEGENDRE) {
      if (order <= 33) // retrieve full precision tabulated values
	webbur::legendre_lookup_points(order, &gaussPoints[0]);
      else { // sandia_rules calculates points/weights together
	if (gaussWeights.size() != order)
	  gaussWeights.resize(order);
	webbur::legendre_compute(order, &gaussPoints[0], &gaussWeights[0]);
	for (size_t i=0; i<order; i++)
	  gaussWeights[i] *= wtFactor; // polynomial weight fn -> PDF
      }
    }
    else
      mode_err = true;
#else
    if (gaussMode == GAUSS_PATTERSON) {
      PCerr << "Error: VPISparseGrid required for Gauss-Patterson points in "
	    << "LegendreOrthogPolynomial::gauss_points()." << std::endl;
      abort_handler(-1);
    }
    else if (gaussMode == GAUSS_LEGENDRE) {
      switch (order) {
      case 1: // zeros of P_1(x) for one Gauss-Legendre point:
	gaussPoints[0] = 0.0;
	break;
      case 2: { // zeros of P_2(x) for two Gauss-Legendre points:
	Real z1osr3 = 1./std::sqrt(3.);
	gaussPoints[0] = -z1osr3;
	gaussPoints[1] =  z1osr3;
	break;
      }
      case 3: { // zeros of P_3(x) for three Gauss-Legendre points:
	Real sr3o5 = std::sqrt(3./5.);
	gaussPoints[0] = -sr3o5;
	gaussPoints[1] =  0.0;
	gaussPoints[2] =  sr3o5;
	break;
      }
      case 4: { // zeros of P_4(x) for four Gauss-Legendre points:
	Real sr30 = std::sqrt(30.), sr525p70sr30 = std::sqrt(525.+70.*sr30)/35.,
	  sr525m70sr30 = std::sqrt(525.-70.*sr30)/35.;
	gaussPoints[0] = -sr525p70sr30;
	gaussPoints[1] = -sr525m70sr30;
	gaussPoints[2] =  sr525m70sr30;
	gaussPoints[3] =  sr525p70sr30;
	break;
      }
      case 5: { // zeros of P_5(x) for five Gauss-Legendre points:
	Real sr70 = std::sqrt(70.), sr245p14sr70 = std::sqrt(245.+14.*sr70)/21.,
	  sr245m14sr70 = std::sqrt(245.-14.*sr70)/21.;
	gaussPoints[0] = -sr245p14sr70;
	gaussPoints[1] = -sr245m14sr70;
	gaussPoints[2] =  0.0;
	gaussPoints[3] =  sr245m14sr70;
	gaussPoints[4] =  sr245p14sr70;
	break;
      }
      // tabulated values from Abramowitz & Stegun have limited precision
      case 6:
	gaussPoints[0] = -0.932469514203152;
	gaussPoints[1] = -0.661209386466265;
	gaussPoints[2] = -0.238619186083197;
	gaussPoints[3] = -gaussPoints[2];
	gaussPoints[4] = -gaussPoints[1];
	gaussPoints[5] = -gaussPoints[0]; break;
      case 7:
	gaussPoints[0] = -0.949107912342759;
	gaussPoints[1] = -0.741531185599394;
	gaussPoints[2] = -0.405845151377397;
	gaussPoints[3] =  0.0;
	gaussPoints[4] = -gaussPoints[2];
	gaussPoints[5] = -gaussPoints[1];
	gaussPoints[6] = -gaussPoints[0]; break;
      case 8:
	gaussPoints[0] = -0.960289856497536;
	gaussPoints[1] = -0.796666477413627;
	gaussPoints[2] = -0.525532409916329;
	gaussPoints[3] = -0.183434642495650;
	gaussPoints[4] = -gaussPoints[3];
	gaussPoints[5] = -gaussPoints[2];
	gaussPoints[6] = -gaussPoints[1];
	gaussPoints[7] = -gaussPoints[0]; break;
      case 9:
	gaussPoints[0] = -0.968160239507626;
	gaussPoints[1] = -0.836031107326636;
	gaussPoints[2] = -0.613371432700590;
	gaussPoints[3] = -0.324253423403809;
	gaussPoints[4] =  0.0;
	gaussPoints[5] = -gaussPoints[3];
	gaussPoints[6] = -gaussPoints[2];
	gaussPoints[7] = -gaussPoints[1];
	gaussPoints[8] = -gaussPoints[0]; break;
      case 10:
	gaussPoints[0] = -0.973906528517172;
	gaussPoints[1] = -0.865063366688985;
	gaussPoints[2] = -0.679409568299024;
	gaussPoints[3] = -0.433395394129247;
	gaussPoints[4] = -0.148874338981631;
	gaussPoints[5] = -gaussPoints[4];
	gaussPoints[6] = -gaussPoints[3];
	gaussPoints[7] = -gaussPoints[2];
	gaussPoints[8] = -gaussPoints[1];
	gaussPoints[9] = -gaussPoints[0]; break;
      default:
	PCerr << "Error: overflow in maximum quadrature order limit (10) in "
	      << "LegendreOrthogPolynomial::gauss_points().  Configure with "
	      << "VPISparseGrid to extend range." << std::endl;
	abort_handler(-1); break;
      }
    }
    else
      mode_err = true;
#endif
  }

  if (mode_err) {
    PCerr << "Error: unsupported Gauss point mode in "
	  << "LegendreOrthogPolynomial::gauss_points()." << std::endl;
    abort_handler(-1);
  }

  return gaussPoints;
}


const RealArray& LegendreOrthogPolynomial::gauss_weights(unsigned short order)
{
  // pull this outside block below since order=0 is initial gauss pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "LegendreOrthogPolynomial::gauss_weights()." << std::endl;
    abort_handler(-1);
  }

  // The sums of the weights = 1, which is the integral of the density
  // function 1/2 over the support range of [-1,+1].  These differ from
  // Abramowitz & Stegun by a constant factor of 1/2.

  bool mode_err = false;
  if (gaussWeights.size() != order) { // if not already computed
    gaussWeights.resize(order);
#ifdef HAVE_SPARSE_GRID
    if (gaussMode == GAUSS_PATTERSON)
      webbur::patterson_lookup_weights(order, &gaussWeights[0]);
    else if (gaussMode == GAUSS_LEGENDRE) {
      if (order <= 33) // retrieve full precision tabulated values
	webbur::legendre_lookup_weights(order, &gaussWeights[0]);
      else { // sandia_rules calculates points/weights together
	if (gaussPoints.size() != order)
	  gaussPoints.resize(order);
	webbur::legendre_compute(order, &gaussPoints[0], &gaussWeights[0]);
      }
    }
    else
      mode_err = true;
    for (size_t i=0; i<order; i++)
      gaussWeights[i] *= wtFactor;
#else
    if (gaussMode == GAUSS_PATTERSON) {
      PCerr << "Error: VPISparseGrid required for Gauss-Patterson weights in "
	    << "LegendreOrthogPolynomial::gauss_weights()." << std::endl;
      abort_handler(-1);
    }
    else if (gaussMode == GAUSS_LEGENDRE) {
      switch (order) {
      case 1: // weights for one Gauss-Legendre point:
	gaussWeights[0] = 1.0; break;
      case 2: // weights for two Gauss-Legendre points:
	gaussWeights[0] = gaussWeights[1] = 0.5; break;
      case 3: // weights for three Gauss-Legendre points:
	gaussWeights[0] = gaussWeights[2] = 5./18.;
	gaussWeights[1] = 4./9.; break;
      case 4: { // weights for four Gauss-Legendre points:
	Real sr30 = std::sqrt(30.);
	gaussWeights[0] = gaussWeights[3] = (18.-sr30)/72.;
	gaussWeights[1] = gaussWeights[2] = (18.+sr30)/72.; break;
      }
      case 5: { // weights for five Gauss-Legendre points:
	Real sr70 = std::sqrt(70.);
	gaussWeights[0] = gaussWeights[4] = (322.-13.*sr70)/1800.;
	gaussWeights[1] = gaussWeights[3] = (322.+13.*sr70)/1800.;
	gaussWeights[2] = 64./225.; break;
      }
      // tabulated values from Abramowitz & Stegun have limited precision
      case 6:
	gaussWeights[0] = gaussWeights[5] = 0.171324492379170 * wtFactor;
	gaussWeights[1] = gaussWeights[4] = 0.360761573048139 * wtFactor;
	gaussWeights[2] = gaussWeights[3] = 0.467913934572691 * wtFactor; break;
      case 7:
	gaussWeights[0] = gaussWeights[6] = 0.129484966168870 * wtFactor;
	gaussWeights[1] = gaussWeights[5] = 0.279705391489277 * wtFactor;
	gaussWeights[2] = gaussWeights[4] = 0.381830050505119 * wtFactor;
	gaussWeights[3] = 0.417959183673469 * wtFactor;	break;
      case 8:
	gaussWeights[0] = gaussWeights[7] = 0.101228536290376 * wtFactor;
	gaussWeights[1] = gaussWeights[6] = 0.222381034453374 * wtFactor;
	gaussWeights[2] = gaussWeights[5] = 0.313706645877887 * wtFactor;
	gaussWeights[3] = gaussWeights[4] = 0.362683783378362 * wtFactor; break;
      case 9:
	gaussWeights[0] = gaussWeights[8] = 0.081274388361574 * wtFactor;
	gaussWeights[1] = gaussWeights[7] = 0.180648160694857 * wtFactor;
	gaussWeights[2] = gaussWeights[6] = 0.260610696402935 * wtFactor;
	gaussWeights[3] = gaussWeights[5] = 0.312347077040003 * wtFactor;
	gaussWeights[4] = 0.330239355001260 * wtFactor;	break;
      case 10:
	gaussWeights[0] = gaussWeights[9] = 0.066671344308688 * wtFactor;
	gaussWeights[1] = gaussWeights[8] = 0.149451349150581 * wtFactor;
	gaussWeights[2] = gaussWeights[7] = 0.219086362515982 * wtFactor;
	gaussWeights[3] = gaussWeights[6] = 0.269266719309996 * wtFactor;
	gaussWeights[4] = gaussWeights[5] = 0.295524224714753 * wtFactor; break;
      default:
	// define Gauss wts from Gauss pts using
	// -(A_{n+1} gamma_n)/(A_n Phi_n'(x_i) Phi_{n+1}(x_i)),
	// which for P(x) with w(x) = 1/2 is (1-x_i^2)/(n P_{n-1}(x_i))^2.
	const RealArray& gauss_pts = gauss_points(order);
	for (size_t i=0; i<order; i++) {
	  const Real& x_i = gauss_pts[i];
	  gaussWeights[i] = (1.-x_i*x_i)
	                  / std::pow(order*get_value(x_i,order-1),2);
	}
	break;
      }
    }
    else
      mode_err = true;
#endif
  }

  if (mode_err) {
    PCerr << "Error: unsupported Gauss weight mode in "
	  << "LegendreOrthogPolynomial::gauss_weights()." << std::endl;
    abort_handler(-1);
  }

  return gaussWeights;
}

} // namespace Pecos
