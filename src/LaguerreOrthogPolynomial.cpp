/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LaguerreOrthogPolynomial
//- Description:  Implementation code for LaguerreOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "LaguerreOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.H"
#endif

//#define DEBUG


namespace Pecos {

const Real& LaguerreOrthogPolynomial::
get_value(const Real& x, unsigned short order)
{
  switch (order) {
  case 0:
    basisPolyValue = 1.;
    break;
  case 1:
    basisPolyValue = -x + 1.;
    break;
  case 2:
    basisPolyValue = (x*x - 4.*x + 2.)/2.;
    break;
  case 3: {
    Real x2 = x*x;
    basisPolyValue = (-x*x2 + 9.*x2 - 18.*x + 6.)/6.;
    break;
  }
  case 4: {
    Real x2 = x*x;
    basisPolyValue = (x2*x2 - 16.*x*x2 + 72.*x2 - 96.*x + 24.)/24.;
    break;
  }
  case 5: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue
      = (-x*x4 + 25.*x4 - 200.*x*x2 + 600.*x2 - 600.*x + 120.) / 120.;
    break;
  }
  case 6: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue = (x4*x2 - 36.*x*x4 + 450.*x4 - 2400.*x*x2 + 5400.*x2 -
		      4320.*x + 720.)/720.;
    break;
  }
  case 7: {
    Real x2 = x*x, x4 = x2*x2, x6 = x4*x2;
    basisPolyValue = (-x*x6 + 49.*x6 - 882.*x*x4 + 7350.*x4 - 29400.*x*x2 +
		      52920.*x2 - 35280.*x + 5040.)/5040.;
    break;
  }
  case 8: {
    Real x2 = x*x, x4 = x2*x2, x6 = x4*x2;
    basisPolyValue = (x2*x6 - 64.*x*x6 + 1568.*x6 - 18816.*x*x4 + 117600.*x4 -
		      376320.*x*x2 + 564480.*x2 - 322560.*x + 40320.)/40320.;
    break;
  }
  case 9: {
    Real x2 = x*x, x4 = x2*x2, x6 = x4*x2, x8 = x4*x4;
    basisPolyValue = (-x*x8 + 81.*x8 - 2592.*x*x6 + 42336.*x6 - 381024.*x*x4 +
		      1905120.*x4 - 5080320.*x*x2 + 6531840.*x2 - 3265920.*x +
		      362880.)/362880.;
    break;
  }
  case 10: {
    Real x2 = x*x, x4 = x2*x2, x6 = x4*x2, x8 = x4*x4;
    basisPolyValue = (x2*x8 - 100.*x*x8 + 4050.*x8 - 86400.*x*x6 +
		      1058400.*x6 - 7620480.*x*x4 + 31752000.*x4 -
		      72576000.*x*x2 + 81648000.*x2 - 36288000.*x +
		      3628800.)/3628800.;
    break;
  }
  default:
    // Support higher order polynomials using the 3 point recursion formula
    Real x2 = x*x, x4 = x2*x2, x6 = x4*x2, x8 = x4*x4,
      L_n = (x2*x8 - 100.*x*x8 + 4050.*x8 - 86400.*x*x6 + 1058400.*x6 -
	     7620480.*x*x4 + 31752000.*x4 - 72576000.*x*x2 + 81648000.*x2 -
	     36288000.*x + 3628800.)/3628800., // L_10
      L_nminus1 = (-x*x8 + 81.*x8 - 2592.*x*x6 + 42336.*x6 - 381024.*x*x4 +
		   1905120.*x4 - 5080320.*x*x2 + 6531840.*x2 - 3265920.*x +
		   362880.)/362880.;           // L_9
    for (size_t i=10; i<order; i++) {
      basisPolyValue = ( (2.*i+1.-x)*L_n - i*L_nminus1 ) / (i+1.); // L_nplus1
      if (i != order-1) {
	L_nminus1 = L_n;
	L_n       = basisPolyValue;
      }
    }
    break;
  }

  return basisPolyValue;
}


const Real& LaguerreOrthogPolynomial::
get_gradient(const Real& x, unsigned short order)
{ 
#ifdef DEBUG
  // See Abramowitz & Stegun, Section 22.8, p.783
  //basisPolyGradient = (order) ?
  //  order*(get_value(x, order) - get_value(x, order-1))/x : 0.;
  if (order) { // be careful with reference to changing basisPolyValue
    Real L_n = get_value(x, order), L_nminus1 = get_value(x, order-1);
    basisPolyGradient = order*(L_n - L_nminus1)/x;
  }
  else
    basisPolyGradient = 0.;
  PCout << "Laguerre gradient approach 1: " << basisPolyGradient << '\n';
#endif // DEBUG

  // The previous approach, while very compact, produces 0/0 = NaN at x = 0.
  // To avoid NaN issue at lower bound, differentiate the 3 pt value recursion
  // to get a 3 point gradient recursion
  switch (order) {

  case 0:
    basisPolyGradient = 0.;
    break;
  case 1:
    basisPolyGradient = -1.;
    break;
  case 2:
    basisPolyGradient = x - 2.;
    break;
  case 3:
    basisPolyGradient = (-x*x + 6.*x - 6.)/2.;
    break;
  case 4: {
    Real x2 = x*x;
    basisPolyGradient = (x2*x - 12.*x2 + 36.*x - 24.)/6.;
    break;
  }
  case 5: {
    Real x2 = x*x;
    basisPolyGradient = (-x2*x2 + 20.*x2*x - 120.*x2 + 240.*x - 120.) / 24.;
    break;
  }
  case 6: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyGradient
      = (x4*x - 30.*x4 + 300.*x2*x - 1200.*x2 + 1800.*x - 720.) / 120.;
    break;
  }
  default:
    // Support higher order polynomials using the 3 point recursion formula
    Real x2 = x*x, x4 = x2*x2,
      dLdx_n
        = (x4*x - 30.*x4 + 300.*x2*x - 1200.*x2 + 1800.*x - 720.)/120., // L'_6
      dLdx_nminus1 = (-x2*x2 + 20.*x2*x - 120.*x2 + 240.*x - 120.)/24.; // L'_5
    for (size_t i=6; i<order; i++) {
      basisPolyGradient // dLdx_nplus1
	= ( (2.*i+1.-x)*dLdx_n - get_value(x,i) - i*dLdx_nminus1 ) / (i+1.);
      if (i != order-1) {
	dLdx_nminus1 = dLdx_n;
	dLdx_n       = basisPolyGradient;
      }
    }
    break;
  }
#ifdef DEBUG
  PCout << "Laguerre gradient approach 2: " << basisPolyGradient << '\n';
#endif // DEBUG

  return basisPolyGradient;
}


const Real& LaguerreOrthogPolynomial::norm_squared(unsigned short order)
{ 
  orthogPolyNormSq = 1.;
  return orthogPolyNormSq;
}


const RealArray& LaguerreOrthogPolynomial::gauss_points(unsigned short order)
{
  // pull this outside block below since order=0 is initial gauss pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "LaguerreOrthogPolynomial::gauss_points()." << std::endl;
    abort_handler(-1);
  }

  if (gaussPoints.size() != order) { // if not already computed
    gaussPoints.resize(order);
#ifdef HAVE_SPARSE_GRID
    if (order <= 20) // retrieve full precision tabulated values
      webbur::laguerre_lookup_points(order, &gaussPoints[0]);
    else { // calculates points/weights together
      if (gaussWeights.size() != order)
	gaussWeights.resize(order);
      webbur::laguerre_compute(order, &gaussPoints[0], &gaussWeights[0]);
    }
#else
    switch (order) {
    case 1: // zeros of L_1(x) for one Gauss-Laguerre point:
      gaussPoints[0] =  1.0; break;
    case 2: { // zeros of L_2(x) for two Gauss-Laguerre points:
      Real sr2 = std::sqrt(2.);
      gaussPoints[0] =  2. - sr2;
      gaussPoints[1] =  2. + sr2; break;
    }
    // Only ~12 digits of precision in Abramowitz & Stegun tabulated values
    case 3:
      gaussPoints[0] =  0.415774556783;
      gaussPoints[1] =  2.294280360279;
      gaussPoints[2] =  6.289945082937; break;
    case 4:
      gaussPoints[0] =  0.322547689619;
      gaussPoints[1] =  1.745761101158;
      gaussPoints[2] =  4.536620296921;
      gaussPoints[3] =  9.395070912301; break;
    case 5:
      gaussPoints[0] =  0.263560319718;
      gaussPoints[1] =  1.413403059107;
      gaussPoints[2] =  3.596425771041;
      gaussPoints[3] =  7.085810005859;
      gaussPoints[4] = 12.640800844276; break;
    case 6:
      gaussPoints[0] =  0.222846604179;
      gaussPoints[1] =  1.188932101673;
      gaussPoints[2] =  2.992736326059;
      gaussPoints[3] =  5.775143569105;
      gaussPoints[4] =  9.837467418383;
      gaussPoints[5] = 15.982873980602; break;
    case 7:
      gaussPoints[0] =  0.193043676560;
      gaussPoints[1] =  1.026664895339;
      gaussPoints[2] =  2.567876744951;
      gaussPoints[3] =  4.900353084526;
      gaussPoints[4] =  8.182153444563;
      gaussPoints[5] = 12.734180291798;
      gaussPoints[6] = 19.395727862263; break;
    case 8:
      gaussPoints[0] =  0.170279632305;
      gaussPoints[1] =  0.903701776799;
      gaussPoints[2] =  2.251086629866;
      gaussPoints[3] =  4.266700170288;
      gaussPoints[4] =  7.045905402393;
      gaussPoints[5] = 10.758516010181;
      gaussPoints[6] = 15.740678641278;
      gaussPoints[7] = 22.863131736889; break;
    case 9:
      gaussPoints[0] =  0.152322227732;
      gaussPoints[1] =  0.807220022742;
      gaussPoints[2] =  2.005135155619;
      gaussPoints[3] =  3.783473973331;
      gaussPoints[4] =  6.204956777877;
      gaussPoints[5] =  9.372985251688;
      gaussPoints[6] = 13.466236911092;
      gaussPoints[7] = 18.833597788992;
      gaussPoints[8] = 26.374071890927; break;
    case 10:
      gaussPoints[0] =  0.137793470540;
      gaussPoints[1] =  0.729454549503;
      gaussPoints[2] =  1.808342901740;
      gaussPoints[3] =  3.401433697855;
      gaussPoints[4] =  5.552496140064;
      gaussPoints[5] =  8.330152746764;
      gaussPoints[6] = 11.843785837900;
      gaussPoints[7] = 16.279257831378;
      gaussPoints[8] = 21.996585811981;
      gaussPoints[9] = 29.920697012274; break;
    default:
      PCerr << "Error: overflow in maximum quadrature order limit (10) in "
	    << "LaguerreOrthogPolynomial::gauss_points().  Configure with "
	    << "VPISparseGrid to extend range." << std::endl;
      abort_handler(-1); break;
    }
#endif
  }

  return gaussPoints;
}


const RealArray& LaguerreOrthogPolynomial::gauss_weights(unsigned short order)
{
  // pull this outside block below since order=0 is initial gauss pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "LaguerreOrthogPolynomial::gauss_weights()." << std::endl;
    abort_handler(-1);
  }

  // The sums of the weights = 1, which is the integral of the density function
  // exp(-x) over the support range of [0,+infinity].

  if (gaussWeights.size() != order) { // if not already computed
    gaussWeights.resize(order);
#ifdef HAVE_SPARSE_GRID
    if (order <= 20) // tabulated values from sandia_rules have full precision
      webbur::laguerre_lookup_weights(order, &gaussWeights[0]);
    else { // sandia_rules calculates points/weights together
      if (gaussPoints.size() != order)
	gaussPoints.resize(order);
      webbur::laguerre_compute(order, &gaussPoints[0], &gaussWeights[0]);
    }
#else
    switch (order) {
    case 1: // weights for one Gauss-Laguerre point:
      gaussWeights[0] = 1.0; break;
    case 2: { // weights for two Gauss-Laguerre points:
      Real sr2 = std::sqrt(2.);
      gaussWeights[0] = (2. + sr2)/4.;
      gaussWeights[1] = (2. - sr2)/4.; break;
    }
    // Only ~12 digits of precision in Abramowitz & Stegun tabulated values
    case 3:
      gaussWeights[0] = 0.711093009929;
      gaussWeights[1] = 0.278517733569;
      gaussWeights[2] = 0.0103892565016; break;
    case 4:
      gaussWeights[0] = 0.603154104342;
      gaussWeights[1] = 0.357418692438;
      gaussWeights[2] = 0.0388879085150;
      gaussWeights[3] = 0.000539294705561; break;
    case 5:
      gaussWeights[0] = 0.521755610583;
      gaussWeights[1] = 0.398666811083;
      gaussWeights[2] = 0.0759424496817;
      gaussWeights[3] = 0.00361175867992;
      gaussWeights[4] = 2.33699723858e-5; break;
    case 6:
      gaussWeights[0] = 0.458964673950;
      gaussWeights[1] = 0.417000830772;
      gaussWeights[2] = 0.113373382074;
      gaussWeights[3] = 0.0103991974531;
      gaussWeights[4] = 0.000261017202815;
      gaussWeights[5] = 8.98547906430e-7; break;
    case 7:
      gaussWeights[0] = 0.409318951701;
      gaussWeights[1] = 0.421831277862;
      gaussWeights[2] = 0.147126348658;
      gaussWeights[3] = 0.0206335144687;
      gaussWeights[4] = 0.00107401014328;
      gaussWeights[5] = 1.58654643486e-5;
      gaussWeights[6] = 3.17031547900e-8; break;
    case 8:
      gaussWeights[0] = 0.369188589342;
      gaussWeights[1] = 0.418786780814;
      gaussWeights[2] = 0.175794986637;
      gaussWeights[3] = 0.0333434922612;
      gaussWeights[4] = 0.00279453623523;
      gaussWeights[5] = 9.07650877336e-5;
      gaussWeights[6] = 8.48574671627e-7;
      gaussWeights[7] = 1.04800117487e-9; break;
    case 9:
      gaussWeights[0] = 0.336126421798;
      gaussWeights[1] = 0.411213980424;
      gaussWeights[2] = 0.199287525371;
      gaussWeights[3] = 0.0474605627657;
      gaussWeights[4] = 0.00559962661079;
      gaussWeights[5] = 0.000305249767093;
      gaussWeights[6] = 6.59212302608e-6;
      gaussWeights[7] = 4.11076933035e-8;
      gaussWeights[8] = 3.29087403035e-11; break;
    case 10:
      gaussWeights[0] = 0.308441115765;
      gaussWeights[1] = 0.401119929155;
      gaussWeights[2] = 0.218068287612;
      gaussWeights[3] = 0.0620874560987;
      gaussWeights[4] = 0.00950151697518;
      gaussWeights[5] = 0.000753008388588;
      gaussWeights[6] = 2.82592334960e-5;
      gaussWeights[7] = 4.24931398496e-7;
      gaussWeights[8] = 1.83956482398e-9;
      gaussWeights[9] = 9.91182721961e-13; break;
    default:
      // define Gauss wts from Gauss pts using
      // -(A_{n+1} gamma_n)/(A_n Phi_n'(x_i) Phi_{n+1}(x_i)),
      // which for L(x), is x_i/(n L_{n-1}(x_i))^2.
      const RealArray& gauss_pts = gauss_points(order);
      for (size_t i=0; i<order; i++) {
	const Real& x_i = gauss_pts[i];
	gaussWeights[i] = x_i/std::pow(order*get_value(x_i, order-1), 2);
      }
      break;
    }
#endif
  }

  return gaussWeights;
}

} // namespace Pecos
