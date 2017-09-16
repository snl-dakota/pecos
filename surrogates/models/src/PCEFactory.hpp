#ifndef PCE_FACTORY_HPP
#define PCE_FACTORY_HPP


namespace surrogates {

/**
   \class PCEFactory
   \brief The object used to build Polynomial Chaos Expansions (PCE).
*/
class PCEFactory {
  public:

  /** \brief Build a PCE using the specifications in opts
   *
   * \param[in] opts Options used to build the surrogate
   *
   */  
  boost::shared_ptr<Approximation> build(const OptionsList opts){
    std::string construction_type =
      opts.get<std::string>("construction method");
    switch (construction_type){
      case "regression" : {
        // boost::shared_ptr<PolyApproximation> PCE
  //       break;
  //     }
  //     case "spectral": { 
  // 	throw(std::runtime_error("spectral pce not yet implemented"));
  //       break;
  //     }
  //     default : {
  //       throw(std::runtime_error("Incorrect approximation type"));
  //     }
  //   }
  // }

  // /// Constructor
  // PCEFactory(){};
  
  /// Destructor
  virtual ~PCEFactory(){};
}

}; // namespace surrogates

#endif //PCE_FACTORY_HPP

  
