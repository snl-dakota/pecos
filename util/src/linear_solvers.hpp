#ifndef LINEAR_SOLVERS_HPP
#define LINEAR_SOLVERS_HPP

#include "LinearSystemSolver.hpp"
#include "OMPSolver.hpp"
#include "LARSolver.hpp"
#include "LSQSolver.hpp"
#include "EqConstrainedLSQSolver.hpp"
#include "CrossValidatedSolver.hpp"

namespace Surrogates {
  
/**\
 * \brief Initialize a linear system solver from a list of options and return 
 * a shared pointer to the base class.
 *
 * \param[in] opts a list of options
 *
 * opts (required parameters)
 * -------------------------
 * 
 * "regression_type" : RegressionType
 *     The regression-type. Accepted values are SVD_LEAST_SQ_REGRESSION, 
 *     QR_LEAST_SQ_REGRESSION, LU_LEAST_SQ_REGRESSION, ORTHOG_MATCH_PURSUIT
 *     LEAST_ANGLE_REGRESSION, LASSO_REGRESSION, EQ_CONS_LEAST_SQ_REGRESSION
 *
 * opts (optional parameters)
 * -------------------------
 * 
 * "use-cross-validation" : boolean default=false
 *     If true return a CrossValidatedSolver which wraps a standard linear 
 *     solver of "regression_type"
 */
boost::shared_ptr<LinearSystemSolver> regression_solver_factory(OptionsList &opts);

boost::shared_ptr<OMPSolver> cast_linear_system_solver_to_ompsolver(boost::shared_ptr<LinearSystemSolver> &solver);

boost::shared_ptr<LARSolver> cast_linear_system_solver_to_larssolver(boost::shared_ptr<LinearSystemSolver> &solver);

boost::shared_ptr<LSQSolver> cast_linear_system_solver_to_lsqsolver(boost::shared_ptr<LinearSystemSolver> &solver);

boost::shared_ptr<EqConstrainedLSQSolver> cast_linear_system_solver_to_eqconstrainedlsqsolver(boost::shared_ptr<LinearSystemSolver> &solver);
  
} // namespace Surrogates

#endif // include guard
