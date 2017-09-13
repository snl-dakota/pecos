/* univariate_polynomials.i */
%module(directors=1,package="PyDakota",autodoc=1) univariate_polynomials
%{
#include "numpy_include.hpp"

#include "pecos_data_types.hpp"

#include "CharlierOrthogPolynomial.hpp"
#include "ChebyshevOrthogPolynomial.hpp"
#include "GenLaguerreOrthogPolynomial.hpp"
#include "HahnOrthogPolynomial.hpp"
#include "HermiteInterpPolynomial.hpp"
#include "HermiteOrthogPolynomial.hpp"
#include "InterpolationPolynomial.hpp"
#include "JacobiOrthogPolynomial.hpp"
#include "LagrangeInterpPolynomial.hpp"
#include "LaguerreOrthogPolynomial.hpp"
#include "LegendreOrthogPolynomial.hpp"
#include "MeixnerOrthogPolynomial.hpp"
#include "NumericGenOrthogPolynomial.hpp"

  using namespace Pecos;
%}
%include "fundamentals.i"

%include "CharlierOrthogPolynomial.hpp"
%include "ChebyshevOrthogPolynomial.hpp"
%include "GenLaguerreOrthogPolynomial.hpp"
%include "HahnOrthogPolynomial.hpp"
%include "HermiteInterpPolynomial.hpp"
%include "HermiteOrthogPolynomial.hpp"
%include "InterpolationPolynomial.hpp"
%include "JacobiOrthogPolynomial.hpp"
%include "LagrangeInterpPolynomial.hpp"
%include "LaguerreOrthogPolynomial.hpp"
%include "LegendreOrthogPolynomial.hpp"
%include "MeixnerOrthogPolynomial.hpp"
%include "NumericGenOrthogPolynomial.hpp"
