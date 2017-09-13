#ifndef SURROGATE_FACTORY_HPP
#define SURROGATE_FACTORY_HPP


namespace Surrogates {
  
enum ApproxType {PCE,GP};

/**
   \class SurrogateFactory
   \brief The object used to build surrogates.
*/
class SurrogateFactory {
  public:

  /** \brief Build a surrogate using the specifications in opts
   *
   * \param[in] opts Options used to build the surrogate
   *
   */  
  build(const OptionsList opts){
    ApproxType approx_type = opts.get<ApproxType>("approximation type");
    switch (approx_type){
      case PCE : {
        PCEFactory(opts, approx);
        break;
      }
      case GP : { 
        GPFactory(opts, approx);
        break;
      }
      default : {
        throw(std::runtime_error("Incorrect approximation type"));
      }
    }
  }

  /// Constructor
  SurrogateFactory(){};
  
  /// Destructor
  virtual ~SurrogateFactory(){};
}

}; // namespace Surrogates

#endif //SURROGATE_FACTORY_HPP

  