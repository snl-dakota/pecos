#ifndef VARIABLES_HPP
#define VARIABLES_HPP

#include <OptionsList.hpp>
namespace Surrogates {

  /**
     \class Variables
     \brief The object representing the meta data of a vector of variables.
  */
  class Variables {
  public:
    Variables();

    virtual ~Variables();

    /**\brief Set options specific to the model
     */
    virtual void set_options(const OptionsList &opts);

    /**\brief Get the model specific options
     */
    virtual void get_options(OptionsList &opts);

    /**\brief Set the number of variables
     */
    
    void set_num_vars(int num_vars);

    /**\brief Get the number of variables
     */
    int num_vars();

  protected:
    /// The number of variables
    int numVars_;

  }; // class Variables


}; // namespace Surrogates

#endif // VARIABLES_HPP
