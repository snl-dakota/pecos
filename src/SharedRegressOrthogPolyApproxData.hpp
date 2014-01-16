/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedRegressPolyApproxData
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        John Jakeman

#ifndef SHARED_REGRESS_ORTHOG_POLY_APPROX_DATA_HPP
#define SHARED_REGRESS_ORTHOG_POLY_APPROX_DATA_HPP

#include "SharedOrthogPolyApproxData.hpp"
#include "CompressedSensing.hpp"

namespace Pecos {


/// Derived approximation class for multivariate orthogonal polynomial
/// approximation with coefficient estimation via regression.

/** The SharedRegressOrthogPolyApproxData class provides a global
    approximation based on multivariate orthogonal polynomials, where
    the coefficients are computed using regression approaches such as
    least squares (L2) or compressed sensing (L1).  It is used
    primarily for polynomial chaos expansion aproaches to UQ. */

class SharedRegressOrthogPolyApproxData: public SharedOrthogPolyApproxData
{
  //
  //- Heading: Friends
  //

  friend class RegressOrthogPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// lightweight constructor
  SharedRegressOrthogPolyApproxData(short basis_type,
				    const UShortArray& approx_order,
				    size_t num_vars);
  /// full constructor
  SharedRegressOrthogPolyApproxData(short basis_type,
				    const UShortArray& approx_order,
				    size_t num_vars,
				    const ExpansionConfigOptions& ec_options,
				    const BasisConfigOptions& bc_options);
  /// destructor
  ~SharedRegressOrthogPolyApproxData();

  //
  //- Heading: Member functions
  //

  /// set crossValidation flag
  void cross_validation(bool flag);
  /// set the noise tolerance(s) for compressed sensing approaches
  void noise_tolerance(const RealVector& noise_tol);
  /// set the L2 penalty parameter for LASSO (elastic net variant)
  void l2_penalty(Real l2_pen);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_data();
  //void increment_data();
  //void decrement_data();
  //void restore_data();
  //void finalize_data();

private:

  //
  //- Heading: Member functions
  //

  /// pack polynomial contributions to Psi matrix for regression
  void pack_polynomial_data(const RealVector& c_vars, const UShortArray& mi,
			    bool add_val,  double* pack_val,  size_t& pv_cntr,
			    bool add_grad, double* pack_grad, size_t& pg_cntr);
  /// pack response contributions to RHS for regression
  void pack_response_data(const SurrogateDataResp& sdr,
			  bool add_val,  double* pack_val,  size_t& pv_cntr,
			  bool add_grad, double* pack_grad, size_t& pg_cntr);

  //
  //- Heading: Data
  //

  /// Wrapper class that is used to solve regression problems
  CompressedSensingTool CSTool;

  /// flag for use of automatic cross-validation for parameter
  /// selection in regression approaches
  bool crossValidation;
  /// noise tolerance(s) for compressed sensing algorithms; vector form
  /// used in cross-validation
  RealVector noiseTols;
  /// the L2 penalty parameter for LASSO (elastic net variant)
  Real l2Penalty;

  /// lower matrix factor in factorization
  RealMatrix lowerFactor;
  /// upper matrix factor in factorization
  RealMatrix upperFactor;
  /// pivoting history of block-LU factorization
  RealMatrix pivotHistory;
  /// pivoting vector
  IntVector pivotVect;
};


inline SharedRegressOrthogPolyApproxData::
SharedRegressOrthogPolyApproxData(short basis_type,
				  const UShortArray& approx_order,
				  size_t num_vars):
  SharedOrthogPolyApproxData(basis_type, approx_order, num_vars), l2Penalty(0.)
{ }


inline SharedRegressOrthogPolyApproxData::
SharedRegressOrthogPolyApproxData(short basis_type,
				  const UShortArray& approx_order,
				  size_t num_vars,
				  const ExpansionConfigOptions& ec_options,
				  const BasisConfigOptions& bc_options):
  SharedOrthogPolyApproxData(basis_type, approx_order, num_vars,
			     ec_options, bc_options), l2Penalty(0.)
{ }


inline SharedRegressOrthogPolyApproxData::~SharedRegressOrthogPolyApproxData()
{ }


inline void SharedRegressOrthogPolyApproxData::cross_validation(bool flag)
{ crossValidation = flag; }


inline void SharedRegressOrthogPolyApproxData::
noise_tolerance(const RealVector& noise_tol)
{ noiseTols = noise_tol; }


inline void SharedRegressOrthogPolyApproxData::l2_penalty(Real l2_pen)
{ l2Penalty = l2_pen; }

} // namespace Pecos

#endif
