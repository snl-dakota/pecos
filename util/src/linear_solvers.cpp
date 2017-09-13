#include "linear_solvers.hpp"

//#include "linear_algebra.hpp"
//#include "OptionsList.hpp"

namespace Surrogates {

boost::shared_ptr<LinearSystemSolver> regression_solver_factory(OptionsList &opts){
  std::string name = "regression_type";
  RegressionType regression_type =
    get_enum_enforce_existance<RegressionType>(opts, name);
  
  bool use_cross_validation = opts.get("use-cross-validation", false);
  if (use_cross_validation){
    boost::shared_ptr<CrossValidatedSolver> cv_solver(new CrossValidatedSolver);
    cv_solver->set_linear_system_solver(regression_type);
    return cv_solver;
  }
  
  switch (regression_type){
  case ORTHOG_MATCH_PURSUIT : {
    boost::shared_ptr<OMPSolver> omp_solver(new OMPSolver);
      return omp_solver;
  }
  case LEAST_ANGLE_REGRESSION : {
    boost::shared_ptr<LARSolver> lars_solver(new LARSolver);
    lars_solver->set_sub_solver(LEAST_ANGLE_REGRESSION);
    return lars_solver;
  }
  case LASSO_REGRESSION : {
    boost::shared_ptr<LARSolver> lars_solver(new LARSolver);
    lars_solver->set_sub_solver(LASSO_REGRESSION);
    return lars_solver;
  }
  case EQ_CONS_LEAST_SQ_REGRESSION : {
    boost::shared_ptr<EqConstrainedLSQSolver>
      eqlsq_solver(new EqConstrainedLSQSolver);
    return eqlsq_solver;
  }
  case SVD_LEAST_SQ_REGRESSION: case LU_LEAST_SQ_REGRESSION:
  case QR_LEAST_SQ_REGRESSION:
  {
    //\todo add set_lsq_solver to lsqsolver class so we can switch
    // between svd, qr and lu factorization methods
    boost::shared_ptr<LSQSolver> lsq_solver(new LSQSolver);
    return lsq_solver;
  }
  default: {
    throw(std::runtime_error("Incorrect \"regression-type\""));
  }
  }
}

boost::shared_ptr<OMPSolver> cast_linear_system_solver_to_ompsolver(boost::shared_ptr<LinearSystemSolver> &solver){
  boost::shared_ptr<OMPSolver> solver_cast =
    boost::dynamic_pointer_cast<OMPSolver>(solver);
  if (!solver_cast)
    throw(std::runtime_error("Could not cast to OMPSolver shared_ptr"));
  return solver_cast;
}

boost::shared_ptr<LARSolver> cast_linear_system_solver_to_larssolver(boost::shared_ptr<LinearSystemSolver> &solver){
  boost::shared_ptr<LARSolver> solver_cast =
    boost::dynamic_pointer_cast<LARSolver>(solver);
  if (!solver_cast)
    throw(std::runtime_error("Could not cast to LARSolver shared_ptr"));
  return solver_cast;
}

boost::shared_ptr<LSQSolver> cast_linear_system_solver_to_lsqsolver(boost::shared_ptr<LinearSystemSolver> &solver){
  boost::shared_ptr<LSQSolver> solver_cast =
    boost::dynamic_pointer_cast<LSQSolver>(solver);
  if (!solver_cast)
    throw(std::runtime_error("Could not cast to LSQSolver shared_ptr"));
  return solver_cast;
}

boost::shared_ptr<EqConstrainedLSQSolver> cast_linear_system_solver_to_eqconstrainedlsqsolver(boost::shared_ptr<LinearSystemSolver> &solver){
  boost::shared_ptr<EqConstrainedLSQSolver> solver_cast =
    boost::dynamic_pointer_cast<EqConstrainedLSQSolver>(solver);
  if (!solver_cast)
    throw(std::runtime_error("Could not cast to EqConstrainedLSQSolver shared_ptr"));
  return solver_cast;
}


} // namespace surrogates
