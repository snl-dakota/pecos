/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef SURROGATE_DATA_HPP
#define SURROGATE_DATA_HPP

#include "pecos_data_types.hpp"
#include <boost/math/special_functions/fpclassify.hpp> //for boostmath::isfinite


namespace Pecos {

/// The representation of a SurrogateDataVars instance.  This representation,
/// or body, may be shared by multiple SurrogateDataVars handle instances.

/** The SurrogateDataVars/SurrogateDataVarsRep pairs utilize a
    handle-body idiom (Coplien, Advanced C++). */

class SurrogateDataVarsRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class SurrogateDataVars;

private:

  //
  //- Heading: Private member functions
  //

  /// constructor
  SurrogateDataVarsRep(const RealVector& c_vars, const IntVector& di_vars,
		       const RealVector& dr_vars, short mode);
  /// alternate constructor (data sizing only)
  SurrogateDataVarsRep(size_t num_c_vars, size_t num_di_vars,
		       size_t num_dr_vars);
  /// destructor
  ~SurrogateDataVarsRep();

  //
  //- Heading: Private data members
  //

  RealVector continuousVars;   ///< continuous variables
  IntVector  discreteIntVars;  ///< discrete integer variables
  RealVector discreteRealVars; ///< discrete real variables

  int referenceCount;        ///< number of handle objects sharing sdvRep
};


inline SurrogateDataVarsRep::
SurrogateDataVarsRep(const RealVector& c_vars, const IntVector& di_vars,
		     const RealVector& dr_vars, short mode): referenceCount(1)
{
  // Note: provided a way to query DataAccess mode for c_vars, could make
  // greater use of operator= for {DEEP,SHALLOW}_COPY modes
  if (mode == DEEP_COPY) {         // enforce deep vector copy
    if (c_vars.length())  copy_data( c_vars, continuousVars);
    if (di_vars.length()) copy_data(di_vars, discreteIntVars);
    if (dr_vars.length()) copy_data(dr_vars, discreteRealVars);
  }
  else if (mode == SHALLOW_COPY) { // enforce shallow vector copy
    if (c_vars.length())
      continuousVars
	= RealVector(Teuchos::View,  c_vars.values(),  c_vars.length());
    if (di_vars.length())
      discreteIntVars
	= IntVector(Teuchos::View,  di_vars.values(), di_vars.length());
    if (dr_vars.length())
      discreteRealVars
	= RealVector(Teuchos::View, dr_vars.values(), dr_vars.length());
  }
  else {                           // default: assume existing Copy/View state
    if (c_vars.length())    continuousVars = c_vars;
    if (di_vars.length())  discreteIntVars = di_vars;
    if (dr_vars.length()) discreteRealVars = dr_vars;
  }
}


inline SurrogateDataVarsRep::
SurrogateDataVarsRep(size_t num_c_vars, size_t num_di_vars, size_t num_dr_vars):
  referenceCount(1)
{
  continuousVars.sizeUninitialized(num_c_vars);
  discreteIntVars.sizeUninitialized(num_di_vars);
  discreteRealVars.sizeUninitialized(num_dr_vars);
}


inline SurrogateDataVarsRep::~SurrogateDataVarsRep()
{ }


/// Container class encapsulating basic parameter data for defining a
/// "truth" data point.

/** A set of these input data points is contained in SurrogateData and
    provides the data to build the approximation.  A handle-body idiom
    is used to avoid excessive data copying overhead. */

class SurrogateDataVars
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  /// default constructor
  SurrogateDataVars();
  /// standard constructor
  SurrogateDataVars(const RealVector& c_vars, const IntVector& di_vars,
		    const RealVector& dr_vars, short mode = DEFAULT_COPY);
  /// alternate constructor (data sizing only)
  SurrogateDataVars(size_t num_c_vars, size_t num_di_vars, size_t num_dr_vars);
  /// copy constructor
  SurrogateDataVars(const SurrogateDataVars& sdv);
  /// destructor
  ~SurrogateDataVars();

  /// assignment operator
  SurrogateDataVars& operator=(const SurrogateDataVars& sdv);
  // equality operator
  //bool operator==(const SurrogateDataVars& sdv) const;

  //
  //- Heading: member functions
  //

  /// return deep copy of SurrogateDataVars instance
  SurrogateDataVars copy() const;

  /// set i^{th} entry within continuousVars
  void continuous_variable(Real c_var, size_t i);
  /// set continuousVars
  void continuous_variables(const RealVector& c_vars,
			    short mode = DEFAULT_COPY);
  /// get continuousVars
  const RealVector& continuous_variables() const;
  /// get view of continuousVars for updating in place
  RealVector continuous_variables_view();

  /// set i^{th} entry within discreteIntVars
  void discrete_int_variable(int di_var, size_t i);
  /// set discreteIntVars
  void discrete_int_variables(const IntVector& di_vars,
			      short mode = DEFAULT_COPY);
  /// get discreteIntVars
  const IntVector& discrete_int_variables() const;
  /// get view of discreteIntVars for updating in place
  IntVector discrete_int_variables_view();

  /// set i^{th} entry within discreteRealVars
  void discrete_real_variable(Real dr_var, size_t i);
  /// set discreteRealVars
  void discrete_real_variables(const RealVector& dr_vars,
			       short mode = DEFAULT_COPY);
  /// get discreteRealVars
  const RealVector& discrete_real_variables() const;
  /// get view of discreteRealVars for updating in place
  RealVector discrete_real_variables_view();

  /// function to check sdvRep (does this handle contain a body)
  bool is_null() const;

private:

  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  SurrogateDataVarsRep* sdvRep;
};


inline SurrogateDataVars::SurrogateDataVars(): sdvRep(NULL)
{ }


inline SurrogateDataVars::
SurrogateDataVars(const RealVector& c_vars, const IntVector& di_vars,
		  const RealVector& dr_vars, short mode):
  sdvRep(new SurrogateDataVarsRep(c_vars, di_vars, dr_vars, mode))
{ }


inline SurrogateDataVars::
SurrogateDataVars(size_t num_c_vars, size_t num_di_vars, size_t num_dr_vars):
  sdvRep(new SurrogateDataVarsRep(num_c_vars, num_di_vars, num_dr_vars))
{ }


inline SurrogateDataVars::SurrogateDataVars(const SurrogateDataVars& sdv)
{
  // Increment new (no old to decrement)
  sdvRep = sdv.sdvRep;
  if (sdvRep) // Check for an assignment of NULL
    ++sdvRep->referenceCount;
}


inline SurrogateDataVars::~SurrogateDataVars()
{
  if (sdvRep) { // Check for NULL
    --sdvRep->referenceCount; // decrement
    if (sdvRep->referenceCount == 0)
      delete sdvRep;
  }
}


inline SurrogateDataVars& SurrogateDataVars::
operator=(const SurrogateDataVars& sdv)
{
  if (sdvRep != sdv.sdvRep) { // prevent re-assignment of same rep
    // Decrement old
    if (sdvRep) // Check for NULL
      if ( --sdvRep->referenceCount == 0 ) 
	delete sdvRep;
    // Increment new
    sdvRep = sdv.sdvRep;
    if (sdvRep) // Check for an assignment of NULL
      ++sdvRep->referenceCount;
  }
  // else if assigning same rep, then leave referenceCount as is

  return *this;
}


//inline bool SurrogateDataVars::operator==(const SurrogateDataVars& sdv) const
//{
//  return (sdvRep->continuousVars   == sdv.sdvRep->continuousVars  &&
//          sdvRep->discreteIntVars  == sdv.sdvRep->discreteIntVars &&
//          sdvRep->discreteRealVars == sdv.sdvRep->discreteRealVars) ?
//    true : false;
//}


/// deep copy of SurrogateDataVars instance
inline SurrogateDataVars SurrogateDataVars::copy() const
{
  SurrogateDataVars sdv(sdvRep->continuousVars,   sdvRep->discreteIntVars,
			sdvRep->discreteRealVars, DEEP_COPY);
  return sdv;
}


inline void SurrogateDataVars::continuous_variable(Real c_var, size_t i)
{ sdvRep->continuousVars[i] = c_var; }


inline void SurrogateDataVars::
continuous_variables(const RealVector& c_vars, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(c_vars, sdvRep->continuousVars);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    sdvRep->continuousVars
      = RealVector(Teuchos::View, c_vars.values(), c_vars.length());
  else                           // default: assume existing Copy/View state
    sdvRep->continuousVars = c_vars;
}


inline const RealVector& SurrogateDataVars::continuous_variables() const
{ return sdvRep->continuousVars; }


inline RealVector SurrogateDataVars::continuous_variables_view()
{
  return RealVector(Teuchos::View, sdvRep->continuousVars.values(),
		    sdvRep->continuousVars.length());
}


inline void SurrogateDataVars::discrete_int_variable(int di_var, size_t i)
{ sdvRep->discreteIntVars[i] = di_var; }


inline void SurrogateDataVars::
discrete_int_variables(const IntVector& di_vars, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(di_vars, sdvRep->discreteIntVars);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    sdvRep->discreteIntVars
      = IntVector(Teuchos::View, di_vars.values(), di_vars.length());
  else                           // default: assume existing Copy/View state
    sdvRep->discreteIntVars = di_vars;
}


inline const IntVector& SurrogateDataVars::discrete_int_variables() const
{ return sdvRep->discreteIntVars; }


inline IntVector SurrogateDataVars::discrete_int_variables_view()
{
  return IntVector(Teuchos::View, sdvRep->discreteIntVars.values(),
		   sdvRep->discreteIntVars.length());
}


inline void SurrogateDataVars::discrete_real_variable(Real dr_var, size_t i)
{ sdvRep->discreteRealVars[i] = dr_var; }


inline void SurrogateDataVars::
discrete_real_variables(const RealVector& dr_vars, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(dr_vars, sdvRep->discreteRealVars);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    sdvRep->discreteRealVars
      = RealVector(Teuchos::View, dr_vars.values(), dr_vars.length());
  else                           // default: assume existing Copy/View state
    sdvRep->discreteRealVars = dr_vars;
}


inline const RealVector& SurrogateDataVars::discrete_real_variables() const
{ return sdvRep->discreteRealVars; }


inline RealVector SurrogateDataVars::discrete_real_variables_view()
{
  return RealVector(Teuchos::View, sdvRep->discreteRealVars.values(),
		    sdvRep->discreteRealVars.length());
}


inline bool SurrogateDataVars::is_null() const
{ return (sdvRep) ? false : true; }


/// The representation of a surrogate data response.  This representation,
/// or body, may be shared by multiple SurrogateDataResp handle instances.

/** The SurrogateDataResp/SurrogateDataRespRep pairs utilize a
    handle-body idiom (Coplien, Advanced C++). */

class SurrogateDataRespRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class SurrogateDataResp;

private:

  //
  //- Heading: Private member functions
  //

  /// constructor
  SurrogateDataRespRep(Real fn_val, const RealVector& fn_grad,
		       const RealSymMatrix& fn_hess, short bits, short mode);
  /// alternate constructor (data sizing only)
  SurrogateDataRespRep(short data_order, size_t num_vars);
  /// destructor
  ~SurrogateDataRespRep();

  //
  //- Heading: Private data members
  //

  short           activeBits; ///< active data bits: 1 (fn), 2 (grad), 4 (hess)
  Real            responseFn; ///< truth response function value
  RealVector    responseGrad; ///< truth response function gradient
  RealSymMatrix responseHess; ///< truth response function Hessian
  int         referenceCount; ///< number of handle objects sharing sdrRep

};


inline SurrogateDataRespRep::
SurrogateDataRespRep(Real fn_val, const RealVector& fn_grad,
		     const RealSymMatrix& fn_hess, short bits, short mode):
  responseFn(fn_val), // always deep copy for scalars
  activeBits(bits), referenceCount(1)
{
  // Note: provided a way to query incoming grad/hess DataAccess modes,
  // could make greater use of operator= for {DEEP,SHALLOW}_COPY modes
  if (mode == DEEP_COPY) {          // enforce vector/matrix deep copy
    if (activeBits & 2) copy_data(fn_grad, responseGrad);
    if (activeBits & 4) copy_data(fn_hess, responseHess);
  }
  else if (mode == SHALLOW_COPY) {  // enforce vector/matrix shallow copy
    if (activeBits & 2)
      responseGrad = RealVector(Teuchos::View, fn_grad.values(),
				fn_grad.length());
    if (activeBits & 4)
      responseHess = RealSymMatrix(Teuchos::View, fn_hess, fn_hess.numRows());
  }
  else {                            // default: assume existing Copy/View state
    if (activeBits & 2) responseGrad = fn_grad;
    if (activeBits & 4) responseHess = fn_hess;
  }
}


inline SurrogateDataRespRep::
SurrogateDataRespRep(short data_order, size_t num_vars):
  activeBits(data_order), referenceCount(1)
{
  if (data_order & 2)
    responseGrad.sizeUninitialized(num_vars);
  if (data_order & 4)
    responseHess.shapeUninitialized(num_vars);
}


inline SurrogateDataRespRep::~SurrogateDataRespRep()
{ }


/// Container class encapsulating basic parameter and response data
/// for defining a "truth" data point.

/** A list of these data points is contained in Approximation instances
    (e.g., Dakota::Approximation::currentPoints) and provides the data
    to build the approximation.  A handle-body idiom is used to avoid
    excessive data copying overhead. */

class SurrogateDataResp
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  /// default constructor
  SurrogateDataResp();
  /// standard constructor
  SurrogateDataResp(Real fn_val, const RealVector& fn_grad,
		    const RealSymMatrix& fn_hess, short bits,
		    short mode = DEFAULT_COPY);
  /// alternate constructor (data sizing only)
  SurrogateDataResp(short data_order, size_t num_vars);
  /// copy constructor
  SurrogateDataResp(const SurrogateDataResp& sdr);
  /// destructor
  ~SurrogateDataResp();

  /// assignment operator
  SurrogateDataResp& operator=(const SurrogateDataResp& sdr);
  // equality operator
  //bool operator==(const SurrogateDataResp& sdr) const;

  //
  //- Heading: member functions
  //

  /// return deep copy of SurrogateDataResp instance
  SurrogateDataResp copy() const;

  /// set activeBits
  void active_bits(short val);
  /// get activeBits
  short active_bits() const;

  /// set responseFn
  void response_function(Real fn);
  /// get responseFn
  Real response_function() const;
  /// get "view" of responseFn for updating in place
  Real& response_function_view();

  /// set i^{th} entry within responseGrad
  void response_gradient(Real grad_i, size_t i);
  /// set responseGrad
  void response_gradient(const RealVector& grad, short mode = DEFAULT_COPY);
  /// get responseGrad
  const RealVector& response_gradient() const;
  /// get view of responseGrad for updating in place
  RealVector response_gradient_view();

  /// set i-j^{th} entry within responseHess
  void response_hessian(Real hess_ij, size_t i, size_t j);
  /// set responseHess
  void response_hessian(const RealSymMatrix& hess, short mode = DEFAULT_COPY);
  /// get responseHess
  const RealSymMatrix& response_hessian() const;
  /// get view of responseHess for updating in place
  RealSymMatrix response_hessian_view();

  /// function to check sdrRep (does this handle contain a body)
  bool is_null() const;
  /// output response function, gradient, and Hessian data
  void write(std::ostream& s) const;

private:

  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  SurrogateDataRespRep* sdrRep;
};


inline SurrogateDataResp::SurrogateDataResp(): sdrRep(NULL)
{ }


inline SurrogateDataResp::
SurrogateDataResp(Real fn_val, const RealVector& fn_grad,
		  const RealSymMatrix& fn_hess, short bits, short mode):
  sdrRep(new SurrogateDataRespRep(fn_val, fn_grad, fn_hess, bits, mode))
{ }


inline SurrogateDataResp::
SurrogateDataResp(short data_order, size_t num_vars):
  sdrRep(new SurrogateDataRespRep(data_order, num_vars))
{ }


inline SurrogateDataResp::SurrogateDataResp(const SurrogateDataResp& sdr)
{
  // Increment new (no old to decrement)
  sdrRep = sdr.sdrRep;
  if (sdrRep) // Check for an assignment of NULL
    ++sdrRep->referenceCount;
}


inline SurrogateDataResp::~SurrogateDataResp()
{
  if (sdrRep) { // Check for NULL
    --sdrRep->referenceCount; // decrement
    if (sdrRep->referenceCount == 0)
      delete sdrRep;
  }
}


inline SurrogateDataResp& SurrogateDataResp::
operator=(const SurrogateDataResp& sdr)
{
  if (sdrRep != sdr.sdrRep) { // prevent re-assignment of same rep
    // Decrement old
    if (sdrRep) // Check for NULL
      if ( --sdrRep->referenceCount == 0 ) 
	delete sdrRep;
    // Increment new
    sdrRep = sdr.sdrRep;
    if (sdrRep) // Check for an assignment of NULL
      ++sdrRep->referenceCount;
  }
  // else if assigning same rep, then leave referenceCount as is

  return *this;
}


//inline bool SurrogateDataResp::operator==(const SurrogateDataResp& sdr) const
//{
//  return ( sdrRep->responseFn   == sdr.sdrRep->responseFn   &&
//	     sdrRep->responseGrad == sdr.sdrRep->responseGrad &&
//	     sdrRep->responseHess == sdr.sdrRep->responseHess ) ? true : false;
//}


/// deep copy of SurrogateDataResp instance
inline SurrogateDataResp SurrogateDataResp::copy() const
{
  SurrogateDataResp sdr(sdrRep->responseFn,   sdrRep->responseGrad,
			sdrRep->responseHess, sdrRep->activeBits, DEEP_COPY);
  return sdr;
}


inline void SurrogateDataResp::active_bits(short val)
{ sdrRep->activeBits = val; }


inline short SurrogateDataResp::active_bits() const
{ return sdrRep->activeBits; }


inline void SurrogateDataResp::response_function(Real fn)
{ sdrRep->responseFn = fn; }


inline Real SurrogateDataResp::response_function() const
{ return sdrRep->responseFn; }


inline Real& SurrogateDataResp::response_function_view()
{ return sdrRep->responseFn; }


inline void SurrogateDataResp::response_gradient(Real grad_i, size_t i)
{ sdrRep->responseGrad[i] = grad_i; }


inline void SurrogateDataResp::
response_gradient(const RealVector& grad, short mode)
{
  if (mode == DEEP_COPY)          // enforce vector deep copy
    copy_data(grad, sdrRep->responseGrad);
  else if (mode == SHALLOW_COPY)  // enforce vector shallow copy
    sdrRep->responseGrad
      = RealVector(Teuchos::View, grad.values(), grad.length());
  else                            // default: assume existing Copy/View state
    sdrRep->responseGrad = grad;
}


inline const RealVector& SurrogateDataResp::response_gradient() const
{ return sdrRep->responseGrad; }


inline RealVector SurrogateDataResp::response_gradient_view()
{
  return RealVector(Teuchos::View, sdrRep->responseGrad.values(),
		    sdrRep->responseGrad.length());
}


inline void SurrogateDataResp::
response_hessian(Real hess_ij, size_t i, size_t j)
{ sdrRep->responseHess(i,j) = hess_ij; }


inline void SurrogateDataResp::
response_hessian(const RealSymMatrix& hess, short mode)
{
  if (mode == DEEP_COPY)          // enforce matrix deep copy
    copy_data(hess, sdrRep->responseHess);
  else if (mode == SHALLOW_COPY)  // enforce matrix shallow copy
    sdrRep->responseHess = RealSymMatrix(Teuchos::View, hess, hess.numRows());
  else                            // default: assume existing Copy/View state
    sdrRep->responseHess = hess;
}


inline const RealSymMatrix& SurrogateDataResp::response_hessian() const
{ return sdrRep->responseHess; }


inline RealSymMatrix SurrogateDataResp::response_hessian_view()
{
  return RealSymMatrix(Teuchos::View, sdrRep->responseHess,
		       sdrRep->responseHess.numRows());
}


inline bool SurrogateDataResp::is_null() const
{ return (sdrRep) ? false : true; }


inline void SurrogateDataResp::write(std::ostream& s) const
{
  if (sdrRep->activeBits & 1)
    s << "function value    =  " << std::setw(WRITE_PRECISION+7)
      << sdrRep->responseFn << '\n';
  if (sdrRep->activeBits & 2) {
    s << "function gradient =\n";
    write_data_trans(s, sdrRep->responseGrad, true, true, true);
  }
  if (sdrRep->activeBits & 4) {
    s << "function Hessian  =\n";
    write_data(s, sdrRep->responseHess, true, true, true);
  }
}


/// std::ostream insertion operator for SurrogateDataResp
inline std::ostream& operator<<(std::ostream& s, const SurrogateDataResp& sdr)
{ sdr.write(s); return s; }


/// Representation of management class for surrogate data defined from
/// input variable data and output response data.

/** Response sets are generally unique for each approximated response
    function, but variable sets are often (but not always) replicated.
    Thus the design accommodates the case where inputs and/or outputs
    involve either shallow or deep copies. */

class SurrogateDataRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class SurrogateData;

public:

private:

  //
  //- Heading: Constructors and destructor
  //

  SurrogateDataRep();  ///< constructor
  ~SurrogateDataRep(); ///< destructor

  //
  //- Heading: Private data members
  //

  /// database of reference variable data sets, with lookup by model/level index
  std::map<UShortArray, SDVArray> varsData;
  /// database of reference response data sets, with lookup by model/level index
  std::map<UShortArray, SDRArray> respData;

  /// sets of popped variables data sets, with lookup by model/level index.
  /// Each popped set is an SDVArray extracted from varsData.
  std::map<UShortArray, SDVArrayDeque> poppedVarsData;
  /// sets of popped response data sets, with lookup by model/level index.
  /// Each popped set is an SDRArray extracted from respData.
  std::map<UShortArray, SDRArrayDeque> poppedRespData;
  /// a stack managing the number of points previously appended that
  /// can be removed by calls to pop()
  std::map<UShortArray, SizetArray> popCountStack;

  /// database key indicating the currently active {SDV,SDR}Arrays.
  /// the key is a multi-index managing multiple modeling dimensions
  /// such as model form, doscretization level, etc.
  UShortArray activeKey;
  /// index of anchor point within {vars,resp}Data, _NPOS if none; for now,
  /// we restrict anchor to reference data to simplify bookkeeping (assume
  /// anchor does not migrate within pushed/popped data)
  std::map<UShortArray, size_t> anchorIndex;

  /// map from failed respData indices to failed data bits; defined
  /// in sample_checks() and used for fault tolerance
  std::map<UShortArray, SizetShortMap> failedRespData;

  /*
  /// a set of variables samples used to build the approximation.
  /// These sample points are fit approximately (e.g., using least
  /// squares regression); exact matching is not enforced.
  SDVArray varsData;
  /// a set of response samples used to build the approximation.  These
  /// sample points are fit approximately (e.g., using least squares
  /// regression); exact matching is not enforced.
  SDRArray respData;

  /// sets of variables samples that have been popped off varsData but
  /// which are available for future restoration.  The granularity of
  /// this 2D array corresponds to multiple trial sets, each
  /// contributing an SDVArray.
  SDVArrayDeque poppedVarsTrials;
  /// sets of response samples that have been popped off respData but
  /// which are available for future restoration.  The granularity of
  /// this 2D array corresponds to multiple trial sets, each
  /// contributing a SDRArray.
  SDRArrayDeque poppedRespTrials;
  /// a stack managing the number of points previously appended that
  /// can be removed by calls to pop()
  SizetArray popCountStack;

  /// a set of variables samples that have been stored for future
  /// restoration.  The granularity of this 2D array corresponds to a
  /// wholesale caching of the current varsData state, where each of
  /// multiple cachings contributes an SDVArray.
  SDVArrayDeque storedVarsData;
  /// a set of response samples that have been stored for future
  /// restoration.  The granularity of this 2D array corresponds to a
  /// wholesale caching of the current respData state, where each of
  /// multiple cachings contributes an SDRArray.
  SDRArrayDeque storedRespData;

  /// map from failed respData indices to failed data bits; defined
  /// in sample_checks() and used for fault tolerance
  SizetShortMap failedRespData;
  */

  /// number of handle objects sharing sdRep
  int referenceCount;
};


inline SurrogateDataRep::SurrogateDataRep(): referenceCount(1)
{ }


inline SurrogateDataRep::~SurrogateDataRep()
{ }


/// Management class for surrogate data defined from input variable
/// data and output response data.

/** Response sets are generally unique for each approximated response
    function, but variable sets are often (but not always) replicated.
    Thus the design accommodates the case where inputs and/or outputs
    involve either shallow or deep copies. */

class SurrogateData
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  SurrogateData();                        ///< default constructor
  SurrogateData(const SurrogateData& sd); ///< copy constructor
  ~SurrogateData();                       ///< destructor

  /// assignment operator
  SurrogateData& operator=(const SurrogateData& sdv);

  //
  //- Heading: Member functions
  //

  /// deep copy of SurrogateData instance with options for shallow copies
  /// of the SurrogateData{Vars,Resp} objects
  SurrogateData copy(short sdv_mode = DEEP_COPY,
		     short sdr_mode = DEEP_COPY) const;

  /// set anchor{Vars,Resp}
  void anchor_point(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr);
  /// set {vars,resp}Data
  void data_points(const SDVArray& sdv_array, const SDRArray& sdr_array);

  /// set anchorVars
  void anchor_variables(const SurrogateDataVars& sdv);
  /// get anchorVars
  const SurrogateDataVars& anchor_variables() const;
  /// set anchorResp
  void anchor_response(const SurrogateDataResp& sdr);
  /// get anchorResp
  const SurrogateDataResp& anchor_response() const;

  /// set varsData
  void variables_data(const SDVArray& sdv_array);
  /// get varsData
  const SDVArray& variables_data() const;
  /// set respData
  void response_data(const SDRArray& sdr_array);
  /// get respData
  const SDRArray& response_data() const;

  /*
  /// get anchorVars.continuous_variables()
  const RealVector& anchor_continuous_variables() const;
  /// get anchorVars.discrete_int_variables()
  const IntVector& anchor_discrete_int_variables() const;
  /// get anchorVars.discrete_real_variables()
  const RealVector& anchor_discrete_real_variables() const;
  /// get anchorResp.active_bits()
  short anchor_active_bits() const;
  /// get anchorResp.response_function()
  Real anchor_function() const;
  /// get anchorResp.response_gradient()
  const RealVector& anchor_gradient() const;
  /// get anchorResp.response_hessian()
  const RealSymMatrix& anchor_hessian() const;

  /// get varsData[i].continuous_variables()
  const RealVector& continuous_variables(size_t i) const;
  /// get varsData[i].discrete_int_variables()
  const IntVector& discrete_int_variables(size_t i) const;
  /// get varsData[i].discrete_real_variables()
  const RealVector& discrete_real_variables(size_t i) const;
  /// get respData[i].active_bits()
  short response_active_bits(size_t i) const;

  /// set respData[i].response_function()
  void response_function(Real fn_val, size_t i);
  /// get respData[i].response_function()
  Real response_function(size_t i) const;
  /// set respData[i].response_gradient()
  void response_gradient(const RealVector& fn_grad, size_t i);
  /// get respData[i].response_gradient()
  const RealVector& response_gradient(size_t i) const;
  /// set respData[i].response_hessian()
  void response_hessian(const RealSymMatrix& fn_hess, size_t i);
  /// get respData[i].response_hessian()
  const RealSymMatrix& response_hessian(size_t i) const;
  */

  /// push sdv onto end of varsData
  void push_back(const SurrogateDataVars& sdv);
  /// push sdr onto end of respData
  void push_back(const SurrogateDataResp& sdr);
  /// push {sdv,sdr} onto ends of {vars,resp}Data
  void push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr);

  /// remove num_pop_pts entries from ends of {vars,resp}Data
  void pop(bool save_data = true);
  /// return a previously popped data set (identified by index) to the
  /// ends of {vars,resp}Data
  void push(size_t index, bool erase_popped = true);

  /// append count to popCountStack
  void pop_count(size_t count);
  /// return popCountStack.back()
  size_t pop_count() const;

  /*
  /// move all entries from {vars,resp}Data to stored{Vars,Resp}Data
  /// (default is push_back)
  void store(size_t index = _NPOS);
  /// return an entry set from stored{Vars,Resp}Data to {vars,resp}Data
  void restore(size_t index = _NPOS);
  /// remove an entry from stored{Vars,Resp}Data (default is pop_back)
  void remove_stored(size_t index = _NPOS);
  /// swap stored and active data sets (variables and response)
  void swap(size_t index);
  */

  /// query presence of anchor{Vars,Resp}
  bool anchor() const;
  /// assign anchorIndex[activeKey] to incoming index
  void anchor_index(size_t index)
  /// return size of {vars,resp}Data arrays (neglecting anchor point)
  size_t points() const;

  /// return total number of available data components
  size_t response_size() const;
  /// return number of failed data components
  size_t failed_response_size() const;
  /// return net number of active data components (total minus failed)
  size_t active_response_size() const;

  // return number of 1D arrays within stored{Vars,Resp}Data 2D arrays
  //size_t stored_sets() const;
  /// return number of 1D arrays within popped{Vars,Resp}Trials 2D arrays
  size_t popped_sets() const;
  /// return number of derivative variables as indicated by size of
  /// gradient/Hessian arrays
  size_t num_derivative_variables() const;

  /// convenience function used by data_checks() for anchorResp and respData
  void response_check(const SurrogateDataResp& sdr, short& failed_data);
  /// screen data sets for samples with Inf/Nan that should be excluded;
  /// defines failedRespData
  void data_checks();
  /// return failedRespData corresponding to active anchorIndex
  short failed_anchor_data() const;
  /// return active failedRespData
  const SizetShortMap& failed_response_data() const;

  /// clear {vars,resp}Data
  void clear_data();
  /// clear popped{Vars,Resp}Trials
  void clear_popped();
  // clear stored{Vars,Resp}Data
  //void clear_stored();

  /// return sdRep
  SurrogateDataRep* data_rep() const;

private:

  //
  //- Heading: Member functions
  //

  /// define or retrieve anchorIndex[activeKey]
  size_t assign_anchor_index();
  /// retrieve anchorIndex[activeKey]
  size_t SurrogateData::retrieve_anchor_index(bool hard_fail);
  
  // set failedAnchorData
  //void failed_anchor_data(short fail_anchor);
  /// get failedRespData
  const std::map<UShortArray, SizetShortMap>& failed_response_data_map() const;
  /// set failedRespData
  void failed_response_data_map(
    const std::map<UShortArray, SizetShortMap>&	fail_resp);

  /*
  /// get storedVarsData
  const SDVArrayDeque& stored_variables_data() const;
  /// set storedVarsData
  void stored_variables_data(const SDVArrayDeque& stored_vars);
  /// get storedRespData
  const SDRArrayDeque& stored_response_data() const;
  /// set storedRespData
  void stored_response_data(const SDRArrayDeque& stored_resp);
  */

  /// get poppedVarsTrials
  const std::map<UShortArray, SDVArrayDeque>& popped_variables_trials() const;
  /// set poppedVarsTrials
  void popped_variables_trials(
    const std::map<UShortArray, SDVArrayDeque>& popped_vars);
  /// get poppedRespTrials
  const std::map<UShortArray, SDRArrayDeque>& popped_response_trials() const;
  /// set poppedRespTrials
  void popped_response_trials(
    const std::map<UShortArray, SDRArrayDeque>& popped_resp);

  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  SurrogateDataRep* sdRep;
};


inline SurrogateData::SurrogateData(): sdRep(new SurrogateDataRep())
{ }


inline SurrogateData::SurrogateData(const SurrogateData& sd)
{
  // Increment new (no old to decrement)
  sdRep = sd.sdRep;
  if (sdRep) // Check for an assignment of NULL
    ++sdRep->referenceCount;
}


inline SurrogateData::~SurrogateData()
{
  if (sdRep) { // Check for NULL
    --sdRep->referenceCount; // decrement
    if (sdRep->referenceCount == 0)
      delete sdRep;
  }
}


inline SurrogateData& SurrogateData::operator=(const SurrogateData& sd)
{
  if (sdRep != sd.sdRep) { // prevent re-assignment of same rep
    // Decrement old
    if (sdRep) // Check for NULL
      if ( --sdRep->referenceCount == 0 ) 
	delete sdRep;
    // Increment new
    sdRep = sd.sdRep;
    if (sdRep) // Check for an assignment of NULL
      ++sdRep->referenceCount;
  }
  // else if assigning same rep, then leave referenceCount as is

  return *this;
}


inline void SurrogateData::
data_points(const SDVArray& sdv_array, const SDRArray& sdr_array)
{
  sdRep->varsData[activeKey] = sdv_array;
  sdRep->respData[activeKey] = sdr_array;
}


inline void SurrogateData::anchor_index(size_t index)
{ if (index != _NPOS) anchorIndex[activeKey] = index; }


inline size_t SurrogateData::assign_anchor_index()
{
  std::map<UShortArray, size_t>::iterator anchor_it
    = anchorIndex.find(activeKey);
  if (anchor_it == anchorIndex.end()) { // no anchor defined
    size_t index = sdv_array.size();
    anchorIndex[activeKey] = index;
    return index;
  }
  else
    return anchor_it->second;
}


inline size_t SurrogateData::retrieve_anchor_index(bool hard_fail)
{
  std::map<UShortArray, size_t>::iterator anchor_it
    = anchorIndex.find(activeKey);
  if (anchor_it == anchorIndex.end() || anchor_it->second == _NPOS) {
    if (hard_fail) {
      PCerr << "Error: lookup failure in SurrogateData::retrieve_anchor()."
	    << std::endl;
      abort_handler(-1);
    }
    return _NPOS;
  }
  else
    return anchor_it->second;
}


inline void SurrogateData::assign_anchor_variables(size_t index)
{
  SDVArray& sdv_array = sdRep->varsData[activeKey];
  if (index == sdv_array.size() || index == _NPOS) sdv_array.push_back(sdv);
  else sdv_array[index] = sdv;
}


inline void SurrogateData::assign_anchor_response(size_t index)
{
  SDRArray& sdr_array = sdRep->respData[activeKey];
  if (index == sdr_array.size() || index == _NPOS) sdr_array.push_back(sdr);
  else sdr_array[index] = sdr;
}


inline void SurrogateData::
anchor_point(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr)
{
  size_t index = assign_anchor_index();
  assign_anchor_variables(index);
  assign_anchor_response(index);
}


inline void SurrogateData::anchor_variables(const SurrogateDataVars& sdv)
{
  size_t index = assign_anchor_index();
  assign_anchor_variables(index);
}


inline const SurrogateDataVars& SurrogateData::anchor_variables() const
{
  size_t index = retrieve_anchor_index(true); // abort on index error
  return sdRep->varsData[activeKey][index];
}


inline void SurrogateData::anchor_response(const SurrogateDataResp& sdr)
{
  size_t index = assign_anchor_index();
  assign_anchor_response(index);
}


inline const SurrogateDataResp& SurrogateData::anchor_response() const
{
  size_t index = retrieve_anchor_index(true); // abort on index error
  return sdRep->respData[activeKey][index];
}


inline void SurrogateData::variables_data(const SDVArray& sdv_array)
{ sdRep->varsData[activeKey] = sdv_array; }


inline const SDVArray& SurrogateData::variables_data() const
{ return sdRep->varsData[activeKey]; }


inline void SurrogateData::response_data(const SDRArray& sdr_array)
{ sdRep->respData[activeKey] = sdr_array; }


inline const SDRArray& SurrogateData::response_data() const
{ return sdRep->respData[activeKey]; }


/*
inline const RealVector& SurrogateData::anchor_continuous_variables() const
{ return anchor_variables().continuous_variables(); }


inline const IntVector& SurrogateData::anchor_discrete_int_variables() const
{ return anchor_variables().discrete_int_variables(); }


inline const RealVector& SurrogateData::anchor_discrete_real_variables() const
{ return anchor_variables().discrete_real_variables(); }


inline const RealVector& SurrogateData::continuous_variables(size_t i) const
{ return sdRep->varsData[activeKey][i].continuous_variables(); }


inline const IntVector& SurrogateData::discrete_int_variables(size_t i) const
{ return sdRep->varsData[activeKey][i].discrete_int_variables(); }


inline const RealVector& SurrogateData::discrete_real_variables(size_t i) const
{ return sdRep->varsData[activeKey][i].discrete_real_variables(); }


inline short SurrogateData::anchor_active_bits() const
{ return anchor_response().active_bits(); }


inline short SurrogateData::response_active_bits(size_t i) const
{ return sdRep->respData[activeKey][i].active_bits(); }


inline Real SurrogateData::anchor_function() const
{ return anchor_response().response_function(); }


inline const RealVector& SurrogateData::anchor_gradient() const
{ return anchor_response().response_gradient(); }


inline const RealSymMatrix& SurrogateData::anchor_hessian() const
{ return anchor_response().response_hessian(); }


inline void SurrogateData::response_function(Real fn_val, size_t i)
{ sdRep->respData[activeKey][i].response_function(fn_val); }


inline Real SurrogateData::response_function(size_t i) const
{ return sdRep->respData[activeKey][i].response_function(); }


inline void SurrogateData::
response_gradient(const RealVector& fn_grad, size_t i)
{ sdRep->respData[activeKey][i].response_gradient(fn_grad); }


inline const RealVector& SurrogateData::response_gradient(size_t i) const
{ return sdRep->respData[activeKey][i].response_gradient(); }


inline void SurrogateData::
response_hessian(const RealSymMatrix& fn_hess, size_t i)
{ sdRep->respData[activeKey][i].response_hessian(fn_hess); }


inline const RealSymMatrix& SurrogateData::response_hessian(size_t i) const
{ return sdRep->respData[activeKey][i].response_hessian(); }
*/


inline void SurrogateData::push_back(const SurrogateDataVars& sdv)
{ sdRep->varsData[activeKey].push_back(sdv); }


inline void SurrogateData::push_back(const SurrogateDataResp& sdr)
{ sdRep->respData[activeKey].push_back(sdr); }


inline void SurrogateData::
push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr)
{
  sdRep->varsData[activeKey].push_back(sdv);
  sdRep->respData[activeKey].push_back(sdr);
}


inline void SurrogateData::pop(bool save_data)
{
  if (sdRep->popCountStack.empty()) {
    PCerr << "\nError: empty count stack in SurrogateData::pop()." << std::endl;
    abort_handler(-1);
  }
  SizetArray& pop_count_stack = sdRep->popCountStack[activeKey];
  size_t num_pop_pts = pop_count_stack.back();
  if (num_pop_pts) {
    size_t data_size = points();
    if (data_size < num_pop_pts) {
      PCerr << "Error: pop count (" << num_pop_pts << ") exceeds data size ("
	    << data_size << ") in SurrogateData::pop(size_t)." << std::endl;
      abort_handler(-1);
    }
    SDVArray& sdv_array_ref = sdRep->varsData[activeKey];
    SDRArray& sdr_array_ref = sdRep->respData[activeKey];
    SDVArrayDeque& popped_sdv_arrays = sdRep->poppedVarsData[activeKey];
    SDRArrayDeque& popped_sdr_arrays = sdRep->poppedRespData[activeKey];
    if (save_data) {
      // append empty arrays and then update them in place
      popped_sdv_arrays.push_back(SDVArray());
      popped_sdr_arrays.push_back(SDRArray());
      SDVArray& last_popped_sdv_array = popped_sdv_arrays.back();
      SDRArray& last_popped_sdr_array = popped_sdr_arrays.back();
      SDVArray::iterator v_end = sdv_array_ref.end();
      SDRArray::iterator r_end = sdr_array_ref.end();
      last_sdv_array.insert(last_popped_sdv_array.begin(),
			    v_end - num_pop_pts, v_end);
      last_sdr_array.insert(last_popped_sdr_array.begin(),
			    r_end - num_pop_pts, r_end);
    }
    size_t new_size = data_size - num_pop_pts;
    sdv_array_ref.resize(new_size); sdr_array_ref.resize(new_size);

    // TO DO: prune failedRespData[activeKey] or leave in map ?
    data_checks(); // from scratch for now...
  }
  pop_count_stack.pop_back();
}


inline void SurrogateData::push(size_t index, bool erase_popped)
{
  SDVArray& sdv_array_ref = sdRep->varsData[activeKey];
  SDRArray& sdr_array_ref = sdRep->respData[activeKey];
  SDVArrayDeque& popped_sdv_arrays = sdRep->poppedVarsData[activeKey];
  SDRArrayDeque& popped_sdr_arrays = sdRep->poppedRespData[activeKey];

  SDVArrayDeque::iterator vit = popped_sdv_arrays.begin();
  SDRArrayDeque::iterator rit = popped_sdr_arrays.begin();
  std::advance(vit, index); std::advance(rit, index);
  size_t num_pts = std::min(vit->size(), rit->size());

  sdv_array_ref.insert(sdv_array_ref.end(), vit->begin(), vit->end());
  sdr_array_ref.insert(sdr_array_ref.end(), rit->begin(), rit->end());

  // TO DO: update failedRespData[activeKey] ?
  data_checks(); // from scratch for now...

  if (erase_popped)
    { popped_sdv_arrays.erase(vit); popped_sdr_arrays.erase(rit); }

  sdRep->popCountStack[activeKey].push_back(num_pts);
}


inline void SurrogateData::pop_count(size_t count)
{ sdRep->popCountStack[activeKey].push_back(count); }


inline size_t SurrogateData::pop_count() const
{
  std::map<UShortArray, SizetArray>::iterator pop_it
    = sdRep->popCountStack.find(activeKey);
  return (pop_it == sdRep->popCountStack.end() || pop_it->second.empty()) ?
    _NPOS : pop_it->second.back();
}


/*
inline void SurrogateData::store(size_t index)
{
  size_t stored_len = std::min(sdRep->storedVarsData.size(),
			       sdRep->storedRespData.size());
  if (index == _NPOS || index == stored_len) { // append
    sdRep->storedVarsData.push_back(sdRep->varsData); // shallow copies
    sdRep->storedRespData.push_back(sdRep->respData); // shallow copies
  }
  else if (index < stored_len) { // replace
    sdRep->storedVarsData[index] = sdRep->varsData; // shallow copies
    sdRep->storedRespData[index] = sdRep->respData; // shallow copies
  }
  else {
    PCerr << "Error: bad index (" << index
	  << ") passed in SurrogateData::store()" << std::endl;
    abort_handler(-1);
  }
  clear_data();
}


inline void SurrogateData::restore(size_t index)
{
  size_t stored_len = std::min(sdRep->storedVarsData.size(),
			       sdRep->storedRespData.size());
  if (index == _NPOS) {
    sdRep->varsData = sdRep->storedVarsData.back(); // shallow copies
    sdRep->respData = sdRep->storedRespData.back(); // shallow copies
  }
  else if (index < stored_len) { // replace
    sdRep->varsData = sdRep->storedVarsData[index]; // shallow copies
    sdRep->respData = sdRep->storedRespData[index]; // shallow copies
  }
  else {
    PCerr << "Error: bad index (" << index
	  << ") passed in SurrogateData::store()" << std::endl;
    abort_handler(-1);
  }
}


inline void SurrogateData::remove_stored(size_t index)
{
  size_t stored_len = std::min(sdRep->storedVarsData.size(),
			       sdRep->storedRespData.size());
  if (index == _NPOS || index == stored_len - 1) {
    sdRep->storedVarsData.pop_back(); // shallow copies
    sdRep->storedRespData.pop_back(); // shallow copies
  }
  else if (index < stored_len) { // replace
    SDVArrayDeque::iterator vit = sdRep->storedVarsData.begin();
    std::advance(vit, index);  sdRep->storedVarsData.erase(vit);
    SDRArrayDeque::iterator rit = sdRep->storedRespData.begin();
    std::advance(rit, index);  sdRep->storedRespData.erase(rit);
  }
  else {
    PCerr << "Error: bad index (" << index
	  << ") passed in SurrogateData::remove_stored()" << std::endl;
    abort_handler(-1);
  }
  clear_data();
}


inline void SurrogateData::swap(size_t index)
{
  if (index == _NPOS)
    return;
  else if (index >= sdRep->storedVarsData.size()) {
    PCerr << "Error: index out of range in SurrogateData::swap()" << std::endl;
    abort_handler(-1);
  }

  // swap stored and active using shallow copies
  SDVArray tmp_vars_data = sdRep->varsData;
  SDRArray tmp_resp_data = sdRep->respData;

  sdRep->varsData = sdRep->storedVarsData[index];
  sdRep->respData = sdRep->storedRespData[index];

  sdRep->storedVarsData[index] = tmp_vars_data;
  sdRep->storedRespData[index] = tmp_resp_data;
}
*/


inline bool SurrogateData::anchor() const
{
  std::map<UShortArray, size_t>::iterator anchor_it
    = anchorIndex.find(activeKey);
  return (anchor_it == anchorIndex.end() || anchor_it->second == _NPOS) ?
    false : true;
}


inline size_t SurrogateData::points() const
{
  return std::min(sdRep->varsData[activeKey].size(),
		  sdRep->respData[activeKey].size());
}


inline size_t SurrogateData::response_size() const
{
  const SDRArray& sdr_array = sdRep->respData[activeKey];
  size_t i, data_size = 0, num_resp = sdr_array.size(), nh;
  short active_bits;
  for (i=0; i<num_resp; ++i) {
    const SurrogateDataResp& sdr = sdr_array[i];
    active_bits = sdr.active_bits(i);
    if (active_bits & 1) ++data_size;
    if (active_bits & 2) data_size += sdr.response_gradient().length();
    if (active_bits & 4) {
      nh = sdr.response_hessian().numRows();
      if (nh) data_size += nh * (nh + 1) / 2;
    }
  }
  return data_size;
}


inline size_t SurrogateData::failed_response_size() const
{
  const SDRArray& sdr_array = sdRep->respData[activeKey];
  size_t failed_size = 0, nh; short fail_bits;
  const SizetShortMap& failed_resp = failed_response_data();
  for (StShMCIter cit=failed_resp.begin(); cit!=failed_resp.end(); ++cit) {
    fail_bits = cit->second;
    const SurrogateDataResp& sdr = sdr_array[cit->first];
    if (fail_bits & 1) ++failed_size;
    if (fail_bits & 2) failed_size += sdr.response_gradient().length();
    if (fail_bits & 4) {
      nh = sdr.response_hessian().numRows();
      if (nh) failed_size += nh * (nh + 1) / 2;
    }
  }
  return failed_size;
}


inline size_t SurrogateData::active_response_size() const
{ return response_size() - failed_response_size(); }


inline size_t SurrogateData::popped_sets() const
{
  return std::min(sdRep->poppedVarsTrials.size(),
		  sdRep->poppedRespTrials.size());
}


inline size_t SurrogateData::num_derivative_variables() const
{
  return //(sdRep->anchorResp.is_null()) ?
    sdRep->respData[activeKey][0].response_gradient().length();// :
    //sdRep->anchorResp.response_gradient().length();
}


inline void SurrogateData::
response_check(const SurrogateDataResp& sdr, short& failed_data)
{
  // We take a conservative approach of rejecting all data of derivative
  // order greater than or equal to a detected failure:

  short resp_bits = sdr.active_bits();
  failed_data = 0;
  if (resp_bits & 1) {
    if (!boost::math::isfinite<Real>(sdr.response_function()))
      failed_data = resp_bits;       // all data for this & higher deriv orders
  }
  if ( (resp_bits & 2) && !failed_data ) {
    const RealVector& grad = sdr.response_gradient();
    size_t j, num_deriv_vars = grad.length();
    for (j=0; j<num_deriv_vars; ++j)
      if (!boost::math::isfinite<Real>(grad[j]))
	{ failed_data = (resp_bits & 6); break; } // this & higher deriv orders
  }
  if ( (resp_bits & 4) && !failed_data ) {
    const RealSymMatrix& hess = sdr.response_hessian();
    size_t j, k, num_deriv_vars = hess.numRows();
    for (j=0; j<num_deriv_vars; ++j)
      for (k=0; k<=j; ++k)
	if (!boost::math::isfinite<Real>(hess(j,k)))
	  { failed_data = 4; break; }             // this & higher deriv orders
  }
}


inline void SurrogateData::data_checks()
{
  SizetShortMap& failed_resp_data = sdRep->failedRespData[activeKey];
  failed_resp_data.clear();
  const SDRArray& resp_data = sdRep->respData[activeKey];
  size_t i, num_resp = resp_data.size(); short failed_data;
  for (i=0; i<num_resp; ++i) {
    response_check(resp_data[i], failed_data);
    if (failed_data)
      failed_resp_data[i] = failed_data;
  }

#ifdef DEBUG
  if (!failed_resp_data.empty()) {
    PCout << "failedRespData:\n";
    for (SizetShortMap::iterator it=failed_resp_data.begin();
	 it!=failed_resp_data.end(); ++it)
      PCout << "index: " << std::setw(6) << it->first
	    << " data: " << it->second << '\n';
  }
#endif // DEBUG
}


inline short SurrogateData::failed_anchor_data() const
{
  SizetShortMap& failed_resp_data = sdRep->failedRespData[activeKey];
  size_t index = retrieve_anchor_index(false); // no hard failure
  SizetShortMap::iterator it = failed_resp_data.find(index);
  return (it == failed_resp_data.end()) ? 0 : it->second;    
}


//inline void SurrogateData::failed_anchor_data(short fail_anchor)
//{ sdRep->failedAnchorData = fail_anchor; }


inline const SizetShortMap& SurrogateData::failed_response_data() const
{ return sdRep->failedRespData[activeKey]; }


inline void SurrogateData::
failed_response_data_map(const std::map<UShortArray, SizetShortMap>& fail_resp)
{ sdRep->failedRespData = fail_resp; }


inline const std::map<UShortArray, SizetShortMap>& SurrogateData::
failed_response_data_map() const
{ return sdRep->failedRespData; }


/*
inline size_t SurrogateData::stored_sets() const
{ return std::min(sdRep->storedVarsData.size(), sdRep->storedRespData.size()); }


inline const SDVArrayDeque& SurrogateData::stored_variables_data() const
{ return sdRep->storedVarsData; }


inline void SurrogateData::
stored_variables_data(const SDVArrayDeque& stored_vars)
{ sdRep->storedVarsData = stored_vars; }


inline const SDRArrayDeque& SurrogateData::stored_response_data() const
{ return sdRep->storedRespData; }


inline void SurrogateData::
stored_response_data(const SDRArrayDeque& stored_resp)
{ sdRep->storedRespData = stored_resp; }
*/


inline const std::map<UShortArray, SDVArrayDeque>& SurrogateData::
popped_variables_trials() const
{ return sdRep->poppedVarsTrials; }


inline void SurrogateData::
popped_variables_trials(const std::map<UShortArray, SDVArrayDeque>& popped_vars)
{ sdRep->poppedVarsTrials = popped_vars; }


inline const std::map<UShortArray, SDRArrayDeque>& SurrogateData::
popped_response_trials() const
{ return sdRep->poppedRespTrials; }


inline void SurrogateData::
popped_response_trials(const std::map<UShortArray, SDRArrayDeque>& popped_resp)
{ sdRep->poppedRespTrials = popped_resp; }


inline SurrogateData SurrogateData::copy(short sdv_mode, short sdr_mode) const
{
  SurrogateData sd;
  //bool anchor_pt = anchor();

  if (sdv_mode == DEEP_COPY) {
    //if (anchor_pt) sd.anchor_variables(sdRep->anchorVars.copy());

    std::map<UShortArray, SDVArray>&    vars_map = sdRep->varsData;
    std::map<UShortArray, SDVArray> new_vars_map;
    std::map<UShortArray, SDVArray>::iterator v_iter;
    for (v_iter=vars_map.begin(); v_iter!=vars_map.end(); ++v_iter) {
      const SDVArray& sdv_array = v_iter->second;
      size_t i, num_pts = sdv_array.size();
      SDVArray new_sdv_array(num_pts);
      for (i=0; i<num_pts; ++i)
	new_sdv_array[i] = sdv_array[i].copy();
      new_vars_map[v_iter->first] = sdv_array;
    }
    sd.variables_data_map(new_vars_map);

    /*
    size_t j, num_sdva = sdRep->storedVarsData.size();
    SDVArrayDeque stored_sdv(num_sdva);
    for (i=0; i<num_sdva; ++i) {
      const SDVArray& rep_stored_i = sdRep->storedVarsData[i];
      num_pts = rep_stored_i.size();
      stored_sdv[i].resize(num_pts);
      for (j=0; j<num_pts; ++j)
	stored_sdv[i][j] = rep_stored_i[j].copy();
    }
    sd.stored_variables_data(stored_sdv);
    */

    std::map<UShortArray, SDVArrayDeque>& popped_vars_map
      = sdRep->poppedVarsTrials;
    std::map<UShortArray, SDVArrayDeque> new_popped_vars_map;
    std::map<UShortArray, SDVArrayDeque>::iterator v2_iter;
    for (v2_iter  = popped_vars_map.begin();
	 v2_iter != popped_vars_map.end(); ++v2_iter) {
      const SDVArrayDeque& sdv_2d_array = v2_iter->second;
      num_sdva = sdv_2d_array.size();
      SDVArrayDeque new_popped_sdv(num_sdva);
      for (i=0; i<num_sdva; ++i) {
	const SDVArray& sdva_i = sdv_2d_array[i];
	num_pts = sdva_i.size();
	new_popped_sdv[i].resize(num_pts);
	for (j=0; j<num_pts; ++j)
	  new_popped_sdv[i][j] = sdva_i[j].copy();
      }
      new_popped_vars_map[v2_iter->first] = new_sdv_2d_array;
    }
    sd.popped_variables_trials(new_popped_vars_map);
  }
  else { // shallow SDV copies based on operator=
    //if (anchor_pt) sd.anchor_variables(sdRep->anchorVars);
    sd.variables_data(sdRep->varsData);
    //sd.stored_variables_data(sdRep->storedVarsData);
    sd.popped_variables_trials(sdRep->poppedVarsTrials);
  }

  if (sdr_mode == DEEP_COPY) {
    //if (anchor_pt) sd.anchor_response(sdRep->anchorResp.copy());

    std::map<UShortArray, SDRArray>&    resp_map = sdRep->respData;
    std::map<UShortArray, SDRArray> new_resp_map;
    std::map<UShortArray, SDRArray>::iterator r_iter;
    for (r_iter = resp_map.begin(); r_iter != resp_map.end(); ++r_iter) {
      const SDRArray& sdr_array = r_iter->second;
      size_t i, num_pts = sdr_array.size();
      SDRArray new_sdr_array(num_pts);
      for (i=0; i<num_pts; ++i)
	new_sdr_array[i] = sdr_array[i].copy();
      new_resp_map[r_iter->first] = sdr_array;
    }
    sd.response_data_map(new_resp_map);

    /*
    size_t j, num_sdra = sdRep->storedRespData.size();
    SDRArrayDeque stored_sdr(num_sdra);
    for (i=0; i<num_sdra; ++i) {
      const SDRArray& rep_stored_i = sdRep->storedRespData[i];
      num_pts = rep_stored_i.size();
      stored_sdr[i].resize(num_pts);
      for (j=0; j<num_pts; ++j)
	stored_sdr[i][j] = rep_stored_i[j].copy();
    }
    sd.stored_response_data(stored_sdr);
    */

    std::map<UShortArray, SDRArrayDeque>& popped_resp_map
      = sdRep->poppedRespTrials;
    std::map<UShortArray, SDRArrayDeque> new_popped_resp_map;
    std::map<UShortArray, SDRArrayDeque>::iterator r2_iter;
    for (r2_iter  = popped_resp_map.begin();
	 r2_iter != popped_resp_map.end(); ++r2_iter) {
      const SDRArrayDeque& sdr_2d_array = r2_iter->second;
      num_sdra = sdr_2d_array.size();
      SDRArrayDeque new_popped_sdr(num_sdra);
      for (i=0; i<num_sdra; ++i) {
	const SDRArray& sdra_i = sdr_2d_array[i];
	num_pts = sdra_i.size();
	new_popped_sdr[i].resize(num_pts);
	for (j=0; j<num_pts; ++j)
	  new_popped_sdr[i][j] = sdra_i[j].copy();
      }
      new_popped_resp_map[r2_iter->first] = new_sdr_2d_array;
    }
    sd.popped_response_trials(new_popped_resp_map);
  }
  else { // shallow SDR copies based on operator=
    sd.response_data(sdRep->respData);
    //sd.stored_response_data(sdRep->storedRespData);
    sd.popped_response_trials(sdRep->poppedRespTrials);
  }

  sd.failed_response_data_map(sd.failed_response_data_map());

  return sd;
}


/*
inline void SurrogateData::clear_stored()
{ sdRep->storedVarsData.clear(); sdRep->storedRespData.clear(); }
*/


inline void SurrogateData::clear_data()
{ sdRep->varsData.clear(); sdRep->respData.clear(); }


inline void SurrogateData::clear_popped()
{
  sdRep->poppedVarsTrials.clear(); sdRep->poppedRespTrials.clear();
  sdRep->popCountStack.clear();
}


inline SurrogateDataRep* SurrogateData::data_rep() const
{ return sdRep; }

} // namespace Pecos

#endif
