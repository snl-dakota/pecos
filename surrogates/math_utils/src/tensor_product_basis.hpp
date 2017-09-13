/**
 * \file TensorProductBasis.hpp
 * \author John D. Jakeman
 * \date 21 March 2013
 */

#ifndef TENSOR_PRODUCT_BASIS_HPP
#define TENSOR_PRODUCT_BASIS_HPP

#include "OrthogonalPolynomial1D.hpp"
#include "PolynomialIndex.hpp"

/**
 * \class TensorProductBasis
 * \brief Evaluation methods for tensor products of 1D bases
 */
class TensorProductBasis
{
private:
  std::vector< OrthogPoly_ptr > bases1D_;
  int numDims_;
  std::vector< IntVector > dimensionList_;

  void compute_max_degree_for_each_dimension( 
	       const std::vector<PolyIndex_ptr>& indices,
	       IntVector &max_degree_1d ) const;

  void compute_1d_basis_values_for_each_dimension( 
               const RealMatrix &pts,
	       const IntVector &max_degree_1d,
	       std::vector< std::vector<RealVector> > &basis_values_1d ) const;

  void compute_1d_basis_derivatives_for_single_dimension( 
               const RealMatrix &pts,
	       int derivative_dim,
	       int max_degree,
	       std::vector<RealVector> &basis_derivatives_1d ) const;

  void derivative_set( const RealMatrix &pts,
	       const std::vector<boost::shared_ptr<PolynomialIndex> >& indices,
	       int derivative_dim,
	       const IntVector &max_degree_1d,
	       const std::vector<std::vector<RealVector> > &basis_values_1d,
	       RealMatrix &result ) const;

public:
  
  TensorProductBasis();
  
  ~TensorProductBasis();

  void clear();

  void set_num_dims( int num_dims );

  OrthogPoly_ptr get_basis_1d( int dim );

  /**
   * \brief Initialise a set of one-dimensional bases.
   */
  void set_bases( const std::vector< OrthogPoly_ptr >& bases, 
		  const std::vector< IntVector >& dimension_list, 
		  int set_size );

    /** 
   * \brief Compute the value of the d-dimensional basis of 
   * degree index = [i1,...,id] at a single point
   */
  Real value( const RealVector &pt, PolyIndex_ptr index ) const;

  /** 
   * \brief Compute the value of the d-dimensional basis of 
   * degree index = [i1,...,id] at a set of points
   */
  void value_set( const RealMatrix &pts,
		  boost::shared_ptr<PolynomialIndex> index,
		  RealVector &result ) const;

  void value_set( const RealMatrix &pts,
		  std::vector< boost::shared_ptr< PolynomialIndex> >& indices,
		  RealMatrix &result ) const;

  /** 
   * \brief Compute the derivative of the d-dimensional basis of 
   * degree index = [i1,...,id] at a single point
   */
  Real derivative( const RealVector &pt, PolyIndex_ptr index,
		   int derivative_dim ) const;

  /** 
   * \brief Compute the derivative of the d-dimensional basis of 
   * degree index = [i1,...,id] at a set of points
   */
  void derivative_set( const RealMatrix &pts,
		       const boost::shared_ptr<PolynomialIndex> index,
		       int derivative_dim,
		       RealVector &result ) const;

  void derivative_set_deprecated( const RealMatrix &pts,
		     std::vector<boost::shared_ptr<PolynomialIndex> >& indices,
		     int derivative_dim,
		     RealMatrix &result ) const;

  void derivative_set( const RealMatrix &pts,
	       const std::vector<boost::shared_ptr<PolynomialIndex> >& indices,
	       int derivative_dim,
	       RealMatrix &result ) const;

  inline OrthogPoly_ptr operator[]( int dim ) const
  {
    return bases1D_[dim];
  };

  void gradient_set( const RealMatrix &pts,
	     const std::vector< boost::shared_ptr< PolynomialIndex> >& indices,
	     RealMatrix &result ) const;
};

typedef boost::shared_ptr<TensorProductBasis> TensorProductBasis_ptr;

#endif // TENSOR_PRODUCT_BASIS_HPP
