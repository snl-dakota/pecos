/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HermiteInterpPolynomial
//- Description:  Implementation code for HermiteInterpPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "HermiteInterpPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.H"
#endif


namespace Pecos {

/** Pre-compute denominator products that are only a function of the
    interpolation points. */
void HermiteInterpPolynomial::precompute_data()
{
  // set up divided difference tables for type1/type2 matching conditions

  int i, num_pts = interpPts.size(), n2 = 2*num_pts, n2m1 = n2 - 1;
  RealArray values(num_pts, 0.), derivs(num_pts, 0.);

  // type 1
  xT1ValDiffTab.resize(num_pts);  yT1ValDiffTab.resize(num_pts);
  xT1GradDiffTab.resize(num_pts); yT1GradDiffTab.resize(num_pts);
  for (i=0; i<num_pts; ++i) {
    RealArray& xt1v_tab_i = xT1ValDiffTab[i];  xt1v_tab_i.resize(n2);
    RealArray& yt1v_tab_i = yT1ValDiffTab[i];  yt1v_tab_i.resize(n2);
    RealArray& xt1g_tab_i = xT1GradDiffTab[i]; xt1g_tab_i.resize(n2m1);
    RealArray& yt1g_tab_i = yT1GradDiffTab[i]; yt1g_tab_i.resize(n2m1);
    values[i] = 1.; // i-th type1 matching condition
    webbur::hermite_interpolant(num_pts, &interpPts[0], &values[0],
				&derivs[0], &xt1v_tab_i[0], &yt1v_tab_i[0],
				&xt1g_tab_i[0], &yt1g_tab_i[0]);
    values[i] = 0.; // restore
  }

  // type 2
  xT2ValDiffTab.resize(num_pts);  yT2ValDiffTab.resize(num_pts);
  xT2GradDiffTab.resize(num_pts); yT2GradDiffTab.resize(num_pts);
  for (i=0; i<num_pts; ++i) {
    RealArray& xt2v_tab_i = xT2ValDiffTab[i];  xt2v_tab_i.resize(n2);
    RealArray& yt2v_tab_i = yT2ValDiffTab[i];  yt2v_tab_i.resize(n2);
    RealArray& xt2g_tab_i = xT2GradDiffTab[i]; xt2g_tab_i.resize(n2m1);
    RealArray& yt2g_tab_i = yT2GradDiffTab[i]; yt2g_tab_i.resize(n2m1);
    derivs[i] = 1.; // i-th type1 matching condition
    webbur::hermite_interpolant(num_pts, &interpPts[0], &values[0],
				&derivs[0], &xt2v_tab_i[0], &yt2v_tab_i[0],
				&xt2g_tab_i[0], &yt2g_tab_i[0]);
    derivs[i] = 0.; // restore
  }
}


/** Compute value of the Hermite type 1 polynomial corresponding to
    interpolation point i. */
const Real& HermiteInterpPolynomial::
type1_value(const Real& x, unsigned short i)
{
  int num_pts = interpPts.size(); Real x_copy = x;
  webbur::hermite_interpolant_value(num_pts, &xT1ValDiffTab[i][0],
    &yT1ValDiffTab[i][0], &xT1GradDiffTab[i][0], &yT1GradDiffTab[i][0],
    1, &x_copy, &basisPolyValue, &basisPolyGradient);
  return basisPolyValue;
}


/** Compute value of the Hermite type 2 polynomial corresponding to
    interpolation point i. */
const Real& HermiteInterpPolynomial::
type2_value(const Real& x, unsigned short i)
{
  int num_pts = interpPts.size(); Real x_copy = x;
  webbur::hermite_interpolant_value(num_pts, &xT2ValDiffTab[i][0],
    &yT2ValDiffTab[i][0], &xT2GradDiffTab[i][0], &yT2GradDiffTab[i][0],
    1, &x_copy, &basisPolyValue, &basisPolyGradient);
  return basisPolyValue;
}


/** Compute derivative with respect to x of the Hermite type 1
    polynomial corresponding to interpolation point i. */
const Real& HermiteInterpPolynomial::
type1_gradient(const Real& x, unsigned short i)
{ 
  int num_pts = interpPts.size(); Real x_copy = x;
  webbur::hermite_interpolant_value(num_pts, &xT1ValDiffTab[i][0],
    &yT1ValDiffTab[i][0], &xT1GradDiffTab[i][0], &yT1GradDiffTab[i][0],
    1, &x_copy, &basisPolyValue, &basisPolyGradient);
  return basisPolyGradient;
}


/** Compute derivative with respect to x of the Hermite type 2
    polynomial corresponding to interpolation point i. */
const Real& HermiteInterpPolynomial::
type2_gradient(const Real& x, unsigned short i)
{ 
  int num_pts = interpPts.size(); Real x_copy = x;
  webbur::hermite_interpolant_value(num_pts, &xT2ValDiffTab[i][0],
    &yT2ValDiffTab[i][0], &xT2GradDiffTab[i][0], &yT2GradDiffTab[i][0],
    1, &x_copy, &basisPolyValue, &basisPolyGradient);
  return basisPolyGradient;
}


const RealArray& HermiteInterpPolynomial::
type1_collocation_weights(unsigned short order)
{
  // pull this outside block below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum order (1) in HermiteInterpPolynomial"
	  << "::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  // TO DO: integration bounds! (semi-infinite and infinite domains)
  // TO DO: interpPts matching order ???
  int i, num_pts = interpPts.size();
  RealArray wts(2*num_pts);
  webbur::hermite_interpolant_rule(num_pts, -1., 1., &interpPts[0], &wts[0]);
  for (i=0; i<num_pts; ++i)
    type1InterpWts[i] = wts[2*i];
  return type1InterpWts;
}


const RealArray& HermiteInterpPolynomial::
type2_collocation_weights(unsigned short order)
{
  // pull this outside block below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum order (1) in HermiteInterpPolynomial"
	  << "::type2_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  // TO DO: integration bounds! (semi-infinite and infinite domains)
  // TO DO: interpPts matching order ???
  int i, num_pts = interpPts.size();
  RealArray wts(2*num_pts);
  webbur::hermite_interpolant_rule(num_pts, -1., 1., &interpPts[0], &wts[0]);
  for (i=0; i<num_pts; ++i)
    type2InterpWts[i] = wts[2*i+1];
  return type2InterpWts;
}

} // namespace Pecos
