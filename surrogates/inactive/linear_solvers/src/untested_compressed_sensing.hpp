#ifndef UNTESTED_COMPRESSED_SENSING_HPP
#define UNTESTED_COMPRESSED_SENSING_HPP

#include "IndexHashSet.hpp"
#include "linear_solvers.hpp"
#include "PolynomialChaosExpansion.hpp"
#include "least_interpolation.hpp"

void get_forward_neighbours( PolyIndex_ptr index, 
		     int num_dims,
		     bool expansion_type,
		     std::set<int> &non_zero_indices,
		     IndexHashSet<PolynomialIndex>::hashSet &indices,
 		     IndexHashSet<PolynomialIndex>::hashSet &neighbours );

void get_ancestors( PolyIndex_ptr index, int num_dims, 
		    std::set<int> &non_zero_indices,
		    IndexHashSet<PolynomialIndex>::hashSet &indices,
		    IndexHashSet<PolynomialIndex>::hashSet &ancestors );

void get_children( PolyIndex_ptr index, int num_dims, int depth, 
		   int expansion_type,
		   std::set<int> &non_zero_indices,
		   IndexHashSet<PolynomialIndex>::hashSet &indices,
		   IndexHashSet<PolynomialIndex>::hashSet &children );

void expand_basis( PolyIndex_ptr index,
		   std::vector<PolyIndex_ptr> &indices, 
		   IndexHashSet<PolynomialIndex>::hashSet &indices_set,
		   int num_dims, int depth,
		   RealMatrix &build_points,
		   RealMatrix &A, RealVector &b,
		   TensorProductBasis_ptr basis,
		   bool normalise,
		   RealVector &column_norms,
		   RealMatrix &Atb,
		   std::set<int> &non_zero_indices, 
		   int expansion_type,
		   IndexHashSet<PolynomialIndex>::hashSet &ancestors );

#ifdef QR_UPDATED_TO_USE_NEW_QR_FACTORIZTION_UPDATE
void tree_orthogonal_matching_pursuit( RealMatrix &A,
				       RealVector &b,
				       PolynomialChaosExpansion *pce,
				       RealMatrix &solutions,
				       RealMatrix &solution_metrics,
				       Real epsilon, 
				       Real ratio,
				       int max_num_non_zero_entries,
				       int verbosity,
				       bool add_ancestors );

class TOMPSolver : public LinearSolver
{
private:
  PolynomialChaosExpansion *pce_;
  std::vector<PolyIndex_ptr> initialIndices_;
  Real ratio_;
  bool addAncestors_;

public:
  TOMPSolver() : pce_( NULL ), ratio_( 0.5 ),
		 addAncestors_( true ) {};

  ~TOMPSolver(){};

  void set_initial_indices( std::vector<PolyIndex_ptr> &initial_indices )
  {
    initialIndices_ = initial_indices;
  }

  void set_pce( PolynomialChaosExpansion *pce )
  {
    pce_ = pce;
  }

  void set_ratio( Real ratio )
  {
    ratio_ = ratio;
  }

  void set_add_ancestors( bool add_ancestors )
  {
    addAncestors_ = add_ancestors;
  }

  /**
   * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
   */
  virtual void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
		      RealMatrix &result_1 )
  {
    if ( !pce_ ) 
      throw( std::runtime_error( "TOMPSolver::solve() Must set pce" ) );
    
    if ( B.numCols() != 1 )
      throw( std::runtime_error(" TOMPSolver::solve() B must be a vector") );

    RealMatrix A_copy( A );
    RealVector column_norms;
    if ( normaliseInputs_ )
      normalise_columns( A_copy, column_norms );

    RealVector b( Teuchos::View, B[0], B.numRows() );

    tree_orthogonal_matching_pursuit( A_copy, b, pce_, 
				      result_0, result_1, 
				      residualTols_[0], 
				      ratio_,
				      maxIters_, 
				      verbosity_,
				      addAncestors_ );

    if ( normaliseInputs_ )
      adjust_coefficients( column_norms, result_0 );
  };

  void build_matrix( const RealMatrix &build_points, RealMatrix &result_0 )
  {
    // return A using the all the indices computed by adaptive omp
    pce_->build_vandermonde( build_points, result_0 );
  }
};
#endif


class OMPCrossValidationSolver : public LinearSolver
{
  // todo: Currently this does only one at a time cross validation and 
  // trainingIndices and validationIndices are never used.
private:
  std::vector< IntVector > trainingIndices_;
  std::vector< IntVector > validationIndices_;
  std::vector< std::vector < RealVector > > cvErrors_;

  RealVector cvScores_;
  int thinningRatio_;

public:
  OMPCrossValidationSolver(){ thinningRatio_ = 1; };

  ~OMPCrossValidationSolver(){};

  OMPCrossValidationSolver( int thinning_ratio )
  { thinningRatio_ = thinning_ratio; };

  void set_cross_validation_indices( std::vector<IntVector> &training_indices,
				     std::vector<IntVector> &validation_indices )
  {
    trainingIndices_ = training_indices;
    validationIndices_ = validation_indices;
  }

  void get_cross_validation_errors( std::vector< RealMatrix > &result )
  {
    int num_folds = validationIndices_.size();
    result.resize( num_folds );
    for ( int i = 0; i < num_folds; i++ )
      {
	int num_valid_pts = validationIndices_[i].length(),
	  num_sols = cvErrors_.size();
	result[i].shapeUninitialized( num_valid_pts, num_sols );
	for ( int j = 0; j < num_valid_pts; j++ )
	  {
	    for ( int k = 0; k < num_sols; k++ )	      
	      result[i](j,k) = cvErrors_[k][i][j];
	  } 
      }   
  }

  /**
   * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
   */
  void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
	      RealMatrix &result_1 )
  {
    if ( B.numCols() != 1 )
      throw( std::runtime_error(" OMPCrossValidationSolver::solve() B must be a vector") );

    RealVector b( Teuchos::View, B[0], B.numRows() );
    orthogonal_matching_pursuit_cholesky( A, b, result_0, result_1, 
					  trainingIndices_, validationIndices_,
					  cvErrors_, cvScores_,
					  residualTols_[0], maxIters_, 
					  verbosity_ );
  };

  void get_cross_validation_scores( RealVector &result )
  {
    result = cvScores_;
  }
};

class BlatmanSolver : public LinearSolver
{
private:

  int solver_;
  RealVector cvScores_;
  
public:
  BlatmanSolver() : solver_( LEAST_ANGLE_REGRESSION ){};

  ~BlatmanSolver(){clear();};

  void clear()
  {
    LinearSolver::clear();
    solver_ = LEAST_ANGLE_REGRESSION;
  };

  void set_sub_solver( int solver_id )
  {
    if ( ( solver_id != LEAST_ANGLE_REGRESSION ) )
      {
	std::stringstream msg;
	msg << "set_sub_solver() solver id must be: " << LEAST_ANGLE_REGRESSION
	    << "\n";
	throw( std::runtime_error( msg.str() ) );
      }
    solver_ = solver_id;
  };

  /**
   * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
   */
  void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
	      RealMatrix &result_1 )
  {
    if ( B.numCols() != 1 )
      throw( std::runtime_error(" LARSSolver::solve() B must be a vector") );

    RealVector b( Teuchos::View, B[0], B.numRows() );
    RealMatrix A_copy( A );
    RealMatrix sols, metrics;
    
    RealVector column_norms;
    if ( normaliseInputs_ )
      normalise_columns( A_copy, column_norms );

    least_angle_regression( A_copy, b, sols, metrics, 
			    residualTols_[0], solver_, 0., 
			    maxIters_, maxNNZ_, verbosity_, 
			    normaliseChoice_ );

    if ( normaliseInputs_ )
      adjust_coefficients( column_norms, sols );

    int M = A.numRows(), N = A.numCols();
    IntVector ordering( std::min( metrics.numCols(), M-1 ) );
    for ( int i = 0; i < ordering.length(); i++ )
      ordering[i] = metrics(1,i);

    RealMatrix best_sol;
    RealVector cvScores_;
    int best_sol_index = loo_step_lsq_cross_validation( A, b, ordering, sols, 
							cvScores_, verbosity_,
							true );

    result_0.shapeUninitialized( N, 1 );
    for ( int i = 0; i < N; i++ )
      result_0(i,0) = sols(i,best_sol_index);

    result_1.shapeUninitialized( 2, 1 );
    result_1(0,0) = cvScores_[best_sol_index];
    result_1(1,0) = best_sol_index;
  }

  void get_cross_valdiations_scores( RealVector &result )
  {
    result = cvScores_;
  }
};


#endif //UNTESTED_COMPRESSED_SENSING_HPP
