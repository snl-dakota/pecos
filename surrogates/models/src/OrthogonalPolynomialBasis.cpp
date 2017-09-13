
namespace Surrogates {

OrthogonalPolynomialBasis::OrthogonalPolynomialBasis() nestedRules_(true){};
OrthogonalPolynomialBasis::~OrthogonalPolynomialBasis(){};

OrthogonalPolynomialBasis::
void get_basis_type_and_colloc_rule(short type, bool nested_rules,
				    short &basis_type, short &colloc_rule){
  switch (type) {
  case STD_NORMAL:
    basis_type  = HERMITE_ORTHOG;
    colloc_rule = (nested_rules) ? GENZ_KEISTER : GAUSS_HERMITE;
    break;
  case STD_UNIFORM:
    basis_type = LEGENDRE_ORTHOG;
    colloc_rule = (nested_rules) ? GAUSS_PATTERSON : GAUSS_LEGENDRE;
    break;
  case STD_EXPONENTIAL:
    basis_type = LAGUERRE_ORTHOG; colloc_rule = GAUSS_LAGUERRE; break;
  case STD_BETA:
    basis_type = JACOBI_ORTHOG; colloc_rule = GAUSS_JACOBI; break;
  case STD_GAMMA:
    basis_type = GEN_LAGUERRE_ORTHOG; colloc_rule = GEN_GAUSS_LAGUERRE; break;
  case POISSON:
    basis_type = CHARLIER_DISCRETE; colloc_rule = GAUSS_CHARLIER; break;
  case BINOMIAL:
    basis_type = KRAWTCHOUK_DISCRETE; colloc_rule = GAUSS_KRAWTCHOUK; break;
  case NEGATIVE_BINOMIAL:
    basis_type = MEIXNER_DISCRETE; colloc_rule = GAUSS_MEIXNER; break;
  default:
    basis_type = NUM_GEN_ORTHOG; colloc_rule = GOLUB_WELSCH; break;
  }
}

OrthogonalPolynomialBasis::
void get_basis_types_and_colloc_rules(const boost::shared_ptr<Surrogates::AleatoryVariables> vars, bool nested_rules, ShortArray& colloc_rules, ShortArray& basis_types){
  size_t num_vars = vars.num_vars();
  basis_types.resize(num_vars);
  colloc_rules.resize(num_vars);
  for (size_t i=0; i<num_vars; ++i){
    get_basis_type_and_colloc_rule(vars.type(i), nested_rules, basis_types[i],
				   colloc_rules[i]);
  }
}

OrthogonalPolynomialBasis::
void initialize_polynomial_basis(const boost::shared_ptr<Surrogates::AleatoryVariables> vars, std::vector<BasisPolynomial> &polynomial_basis){

  boost::shared_ptr<Surrogates::AleatoryVariables> aleatory_vars = Teuchos::rcp_dynamic_cast<Surrogates::AleatoryVariables>(vars,true);
  if (aleatory_vars.is_null())
    throw(std::runtime_error("vars is not an object of type AleatoryVariables"));

  ShortArray basis_types, colloc_rules;
  get_basis_types_and_colloc_rules(vars, nestedRules_, basis_types,
				   colloc_rules);

  size_t num_vars = varTransform_.num_vars();
  polynomial_basis.resize(num_vars);
  for (size_t i=0; i<num_vars; ++i)
    polynomial_basis_[i] = BasisPolynomial(basis_types[i], colloc_rules[i]);
}

} // namespace Surrogates
