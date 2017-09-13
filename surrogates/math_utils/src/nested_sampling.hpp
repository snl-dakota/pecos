#ifndef SAMPLER_HPP
#define SAMPLER_HPP

#include "LinearAlgebra.hpp"
#include "PolynomialChaosExpansion.hpp"

class Sampler
{
 
public:
  
  /// Default onstrcutor
  Sampler(){};
  
  /// Destructor
  ~Sampler(){};

  /// Generate a set of samples.
  void generate_samples( int num_dims, int num_samples, 
			 RealMatrix &result ) const;

  /// Enrinch a set of M samples with N new samples
  virtual void enrich_samples( int num_dims, const RealMatrix &initial_samples, 
			       int num_new_samples,
			       RealMatrix &result ) const = 0;


  /// Get the best indices of a set of candidate samples to enrinch a 
  /// set of M samples with N new samples from a set of K candidate
  /// samples, K>N. This function assumes the candidate samples are unique.
  virtual void get_enriched_sample_indices( int num_dims, 
					    const RealMatrix &initial_samples,
					    int num_new_samples,
					    const RealMatrix & candidate_samples,
					    IntVector &result ) const = 0;

  /// Enrinch a set of M samples with N new samples from a set of K candidate
  /// samples, K>N
  void enrich_samples( int num_dims, const RealMatrix &initial_samples,
		       int num_new_samples,
		       const RealMatrix & candidate_samples,
		       RealMatrix &result ) const;

  /// Remove any candidate samples that are already in an initial set of samples
  static void get_unique_samples( const RealMatrix &initial_samples,
				  int num_new_samples,
				  const RealMatrix &candidate_samples,
				  RealMatrix &result );

};

class LejaSampler : public Sampler
{
private:
  PolynomialChaosExpansion pce_;

  int numCandidateSamples_;
  int seed_;
  bool precondition_;

public:
  /// Default constructor
  LejaSampler() : numCandidateSamples_(10000), seed_(0),
		  precondition_(true){};

  /// Destructor
  ~LejaSampler(){};

  /** \brief  Use pivoted LU factorization to get the indices of a set of 
   * candidate samples to enrinch a set of M samples with N new samples 
   * from a set of K candidate samples, K>N. 
   *
   * The selected samples with be a greedy optimization of the determinant 
   * of the polynomial vandermonde. The resulting set of points is often
   * referred to as a Leja sequence
   */
  void get_enriched_sample_indices( int num_dims, 
				    const RealMatrix &initial_samples,
				    int num_new_samples,
				    const RealMatrix & candidate_samples,
				    IntVector &result ) const;

  /** \brief Draw samples from the equilibrium measure. For all bounded variables
   * the equilibrium measure is the arcsine measure, i.e. "Chebyshev" density
   */
  static void get_candidate_samples( int num_dims, int num_samples, int seed,
				     RealMatrix &result );

  /// Enrinch a set of M samples with N new samples
  void enrich_samples( int num_dims, const RealMatrix &initial_samples,
		       int num_new_samples,
		       RealMatrix &result ) const;

  /// Set the seed of the random number generator
  void set_seed( int seed );

  /// Set the number of canidadates used to select new points
  void set_num_candidate_samples( int num_candidate_samples );

  /// Set the basis for which we want to maximize the determinant
  void set_pce( PolynomialChaosExpansion & pce );

  /// Specify whether to use preconditioning when selecting points
  void set_precondition( bool precondition );

  /**
   * \brief Use the Christoffel function to precondition the polynomial
   * basis matrix in place.
   */
  static void apply_preconditioning( RealMatrix &basisMatrix );

  /**
   * Evaluate the basis functions at a set of samples.
   */
  void build_basis_matrix( const RealMatrix &samples, 
			   RealMatrix &result ) const;

};


#endif //SAMPLER_HPP
