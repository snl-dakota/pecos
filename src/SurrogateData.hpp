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
  //- Heading: Member functions
  //

  /// update {varsData,respData}Iter from activeKey
  void update_active_iterators();

  //
  //- Heading: Private data members
  //

  /// database of reference variable data sets, with lookup by model/level index
  std::map<UShortArray, SDVArray> varsData;
  /// iterator to active entry within varsData
  std::map<UShortArray, SDVArray>::iterator varsDataIter;

  /// database of reference response data sets, with lookup by model/level index
  std::map<UShortArray, SDRArray> respData;
  /// iterator to active entry within respData
  std::map<UShortArray, SDRArray>::iterator respDataIter;

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


inline void SurrogateDataRep::update_active_iterators()
{
  varsDataIter = varsData.find(activeKey);
  if (varsDataIter == varsData.end()) {
    std::pair<UShortArray, SDVArray> sdv_pair(activeKey, SDVArray());
    varsDataIter = varsData.insert(sdv_pair).first;
  }
  respDataIter = respData.find(activeKey);
  if (respDataIter == respData.end()) {
    std::pair<UShortArray, SDRArray> sdr_pair(activeKey, SDRArray());
    respDataIter = respData.insert(sdr_pair).first;
  }
}


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
  SurrogateData(bool handle);             ///< handle constructor
  SurrogateData(const SurrogateData& sd); ///< copy constructor
  ~SurrogateData();                       ///< destructor

  /// assignment operator
  SurrogateData& operator=(const SurrogateData& sdv);

  //
  //- Heading: Member functions
  //

  /// copy incoming sd instance with options for shallow or deep copies
  /// of the SurrogateData{Vars,Resp} objects
  void copy(const SurrogateData& sd, short sdv_mode = DEEP_COPY,
	    short sdr_mode = DEEP_COPY) const;
  /// returns a new SurrogateData instance with options for shallow or
  /// deep copies of the SurrogateData{Vars,Resp} objects
  SurrogateData copy(short sdv_mode = DEEP_COPY,
		     short sdr_mode = DEEP_COPY) const;
  /// copy active data for incoming sd instance with options for
  /// shallow or deep copies of the SurrogateData{Vars,Resp} objects
  void copy_active(const SurrogateData& sd, short sdv_mode,
		   short sdr_mode) const;
  
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

  /// set varsData[activeKey]
  void variables_data(const SDVArray& sdv_array);
  /// get varsData[activeKey]
  const SDVArray& variables_data() const;
  /// get varsData[activeKey]
  SDVArray& variables_data();

  /// set respData[activeKey]
  void response_data(const SDRArray& sdr_array);
  /// get respData[activeKey]
  const SDRArray& response_data() const;
  /// get respData[activeKey]
  SDRArray& response_data();

  /// get varsData
  const std::map<UShortArray, SDVArray>& variables_data_map() const;
  /// set varsData
  void variables_data_map(const std::map<UShortArray, SDVArray>& vars_map);
  /// get respData
  const std::map<UShortArray, SDRArray>& response_data_map() const;
  /// set respData
  void response_data_map(const std::map<UShortArray, SDRArray>& resp_map);

  /// push sdv onto end of varsData
  void push_back(const SurrogateDataVars& sdv);
  /// push sdr onto end of respData
  void push_back(const SurrogateDataResp& sdr);
  /// push {sdv,sdr} onto ends of {vars,resp}Data
  void push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr);

  /// remove the first entry from {vars,resp}Data, managing anchorIndex
  /// Note: inefficient for std::vector's, but needed in rare cases.
  void pop_front();
  /// remove the last entry from {vars,resp}Data, managing anchorIndex
  void pop_back();

  /// remove num_pop_pts entries from ends of {vars,resp}Data
  void pop(bool save_data = true);
  /// return a previously popped data set (identified by index) to the
  /// ends of {vars,resp}Data
  void push(size_t index, bool erase_popped = true);

  /// append count to popCountStack
  void pop_count(size_t count);
  /// return popCountStack.back()
  size_t pop_count() const;

  /// query presence of anchor{Vars,Resp}
  bool anchor() const;
  /// assign anchorIndex[activeKey] to incoming index
  void anchor_index(size_t index);
  /// return anchorIndex[activeKey], if defined
  size_t anchor_index() const;
  /// erase anchorIndex[activeKey]
  void clear_anchor_index();

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
  /// return number of 1D arrays within popped{Vars,Resp}Data 2D arrays
  size_t popped_sets() const;

  /// return number of gradient variables from size of gradient arrays
  size_t num_gradient_variables() const;
  /// return number of Hessian variables from size of Hessian arrays
  size_t num_hessian_variables() const;
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
  /// set active failedRespData
  void failed_response_data(const SizetShortMap& fail_data);
  /// return active failedRespData
  const SizetShortMap& failed_response_data() const;

  /// assign activeKey and update active iterators
  void active_key(const UShortArray& key);
  /// return activeKey
  const UShortArray& active_key() const;
  /// reset initial state by clearing all model keys (empties all maps)
  void clear_keys();

  /// clear active key within {vars,resp}Data
  void clear_active();
  /// clear all inactive data within {vars,resp}Data
  void clear_inactive();
  /// clear active key within popped{Vars,Resp}Data
  void clear_active_popped();
  /// clear popped{Vars,Resp}Data
  void clear_popped();

  /// return sdRep
  SurrogateDataRep* data_rep() const;

  /// function to check sdRep (does this handle contain a body)
  bool is_null() const;

private:

  //
  //- Heading: Member functions
  //

  /// define or retrieve anchorIndex[activeKey]
  size_t assign_anchor_index();
  /// retrieve anchorIndex[activeKey]
  size_t retrieve_anchor_index(bool hard_fail = false) const;

  /// assign sdv within varsData[activeKey] at indicated index
  void assign_variables(const SurrogateDataVars& sdv, size_t index);
  /// assign sdr within respData[activeKey] at indicated index
  void assign_response(const SurrogateDataResp& sdr, size_t index);

  // set failedAnchorData
  //void failed_anchor_data(short fail_anchor);
  /// get failedRespData
  const std::map<UShortArray, SizetShortMap>& failed_response_data_map() const;
  /// set failedRespData
  void failed_response_data_map(
    const std::map<UShortArray, SizetShortMap>&	fail_resp);

  /// get active poppedVarsData
  const SDVArrayDeque& popped_variables() const;
  /// set active poppedVarsData
  void popped_variables(const SDVArrayDeque& popped_vars);
  /// get active poppedRespData
  const SDRArrayDeque& popped_response() const;
  /// set active poppedRespData
  void popped_response(const SDRArrayDeque& popped_resp);

  /// get poppedVarsData
  const std::map<UShortArray, SDVArrayDeque>& popped_variables_map() const;
  /// set poppedVarsData
  void popped_variables_map(
    const std::map<UShortArray, SDVArrayDeque>& popped_vars);
  /// get poppedRespData
  const std::map<UShortArray, SDRArrayDeque>& popped_response_map() const;
  /// set poppedRespData
  void popped_response_map(
    const std::map<UShortArray, SDRArrayDeque>& popped_resp);

  //
  //- Heading: Private data members
  //

 
  /// pointer to the body (handle-body idiom)
  SurrogateDataRep* sdRep;
};


inline SurrogateData::SurrogateData(): sdRep(NULL)
{ }


inline SurrogateData::SurrogateData(bool handle):
  sdRep(new SurrogateDataRep())
{ sdRep->update_active_iterators(); } // default activeKey is empty array


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


inline void SurrogateData::active_key(const UShortArray& key)
{
  if (sdRep->activeKey != key) {
    sdRep->activeKey = key;
    sdRep->update_active_iterators();
  }
}


inline const UShortArray& SurrogateData::active_key() const
{ return sdRep->activeKey; }


inline void SurrogateData::
data_points(const SDVArray& sdv_array, const SDRArray& sdr_array)
{
  sdRep->varsDataIter->second = sdv_array;
  sdRep->respDataIter->second = sdr_array;
}


inline void SurrogateData::anchor_index(size_t index)
{ if (index != _NPOS) sdRep->anchorIndex[sdRep->activeKey] = index; }


inline size_t SurrogateData::anchor_index() const
{ return retrieve_anchor_index(false); }


inline void SurrogateData::clear_anchor_index()
{ sdRep->anchorIndex.erase(sdRep->activeKey); }


inline size_t SurrogateData::assign_anchor_index()
{
  std::map<UShortArray, size_t>& anchor_index = sdRep->anchorIndex;
  const UShortArray& key = sdRep->activeKey;
  std::map<UShortArray, size_t>::iterator anchor_it = anchor_index.find(key);
  size_t index;
  if (anchor_it == anchor_index.end()) { // no anchor defined
    index = sdRep->varsDataIter->second.size();
    anchor_index[key] = index;
  }
  else {
    index = anchor_it->second;
    if (index == _NPOS)
      anchor_it->second = index = sdRep->varsDataIter->second.size();
  }
  return index;
}


inline size_t SurrogateData::retrieve_anchor_index(bool hard_fail) const
{
  std::map<UShortArray, size_t>& anchor_index = sdRep->anchorIndex;
  std::map<UShortArray, size_t>::iterator anchor_it
    = anchor_index.find(sdRep->activeKey);
  if (anchor_it == anchor_index.end() || anchor_it->second == _NPOS) {
    if (hard_fail) {
      PCerr << "Error: lookup failure in SurrogateData::retrieve_anchor_index"
	    << "()." << std::endl;
      abort_handler(-1);
    }
    return _NPOS;
  }
  else
    return anchor_it->second;
}


inline void SurrogateData::
assign_variables(const SurrogateDataVars& sdv, size_t index)
{
  SDVArray& sdv_array = sdRep->varsDataIter->second;
  if (index == sdv_array.size() || index == _NPOS) sdv_array.push_back(sdv);
  else sdv_array[index] = sdv;
}


inline void SurrogateData::
assign_response(const SurrogateDataResp& sdr, size_t index)
{
  SDRArray& sdr_array = sdRep->respDataIter->second;
  if (index == sdr_array.size() || index == _NPOS) sdr_array.push_back(sdr);
  else sdr_array[index] = sdr;
}


inline void SurrogateData::
anchor_point(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr)
{
  size_t index = assign_anchor_index();
  assign_variables(sdv, index);
  assign_response(sdr, index);
}


inline void SurrogateData::anchor_variables(const SurrogateDataVars& sdv)
{
  size_t index = assign_anchor_index();
  assign_variables(sdv, index);
}


inline const SurrogateDataVars& SurrogateData::anchor_variables() const
{
  size_t index = retrieve_anchor_index(true); // abort on index error
  return sdRep->varsDataIter->second[index];
}


inline void SurrogateData::anchor_response(const SurrogateDataResp& sdr)
{
  size_t index = assign_anchor_index();
  assign_response(sdr, index);
}


inline const SurrogateDataResp& SurrogateData::anchor_response() const
{
  size_t index = retrieve_anchor_index(true); // abort on index error
  return sdRep->respDataIter->second[index];
}


inline void SurrogateData::variables_data(const SDVArray& sdv_array)
{ sdRep->varsDataIter->second = sdv_array; }


inline const SDVArray& SurrogateData::variables_data() const
{ return sdRep->varsDataIter->second; }


inline SDVArray& SurrogateData::variables_data()
{ return sdRep->varsDataIter->second; }


inline void SurrogateData::response_data(const SDRArray& sdr_array)
{ sdRep->respDataIter->second = sdr_array; }


inline const SDRArray& SurrogateData::response_data() const
{ return sdRep->respDataIter->second; }


inline SDRArray& SurrogateData::response_data()
{ return sdRep->respDataIter->second; }


inline const std::map<UShortArray, SDVArray>& SurrogateData::
variables_data_map() const
{ return sdRep->varsData; }


inline void SurrogateData::
variables_data_map(const std::map<UShortArray, SDVArray>& vars_map)
{ sdRep->varsData = vars_map; }


inline const std::map<UShortArray, SDRArray>& SurrogateData::
response_data_map() const
{ return sdRep->respData; }


inline void SurrogateData::
response_data_map(const std::map<UShortArray, SDRArray>& resp_map)
{ sdRep->respData = resp_map; }


inline void SurrogateData::push_back(const SurrogateDataVars& sdv)
{ sdRep->varsDataIter->second.push_back(sdv); }


inline void SurrogateData::push_back(const SurrogateDataResp& sdr)
{ sdRep->respDataIter->second.push_back(sdr); }


inline void SurrogateData::
push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr)
{
  sdRep->varsDataIter->second.push_back(sdv);
  sdRep->respDataIter->second.push_back(sdr);
}


inline void SurrogateData::pop_back()
{
  SDVArray& sdv_array = sdRep->varsDataIter->second;
  SDRArray& sdr_array = sdRep->respDataIter->second;
  if (!sdv_array.empty()) sdv_array.pop_back();
  if (!sdr_array.empty()) sdr_array.pop_back();
  if (retrieve_anchor_index() >= points()) // popped point was anchor
    clear_anchor_index();
}


inline void SurrogateData::pop_front()
{
  if (points()) {
    SDVArray& sdv_array = sdRep->varsDataIter->second;
    SDRArray& sdr_array = sdRep->respDataIter->second;
    if (!sdv_array.empty()) sdv_array.erase(sdv_array.begin());
    if (!sdr_array.empty()) sdr_array.erase(sdr_array.begin());
    size_t index = retrieve_anchor_index();
    if (index == 0) // popped point was anchor
      clear_anchor_index();
    else if (index != _NPOS) // anchor (still) exists, decrement its index
      anchor_index(index - 1);
  }
}


inline void SurrogateData::pop(bool save_data)
{
  SDVArray& sdv_array_ref = sdRep->varsDataIter->second;
  SDRArray& sdr_array_ref = sdRep->respDataIter->second;
  size_t num_ref_pts = std::min(sdv_array_ref.size(), sdr_array_ref.size());

  // harden logic for case of an empty SurrogateData (e.g.,
  // distinct discrepancy at level 0)
  if (sdRep->popCountStack.empty()) {
    if (num_ref_pts == 0)
      return; // assume inactive SurrogateData -> ignore pop request
    else {
      PCerr << "\nError: empty count stack in SurrogateData::pop()."<<std::endl;
      abort_handler(-1);
    }
  }

  const UShortArray& key = sdRep->activeKey;
  SizetArray& pop_count_stack = sdRep->popCountStack[key];
  size_t num_pop_pts = pop_count_stack.back();
  if (num_pop_pts) {
    if (num_ref_pts < num_pop_pts) {
      PCerr << "Error: pop count (" << num_pop_pts << ") exceeds data size ("
	    << num_ref_pts << ") in SurrogateData::pop(size_t)." << std::endl;
      abort_handler(-1);
    }
    SDVArrayDeque& popped_sdv_arrays = sdRep->poppedVarsData[key];
    SDRArrayDeque& popped_sdr_arrays = sdRep->poppedRespData[key];
    if (save_data) {
      // append empty arrays and then update them in place
      popped_sdv_arrays.push_back(SDVArray());
      popped_sdr_arrays.push_back(SDRArray());
      SDVArray& last_popped_sdv_array = popped_sdv_arrays.back();
      SDRArray& last_popped_sdr_array = popped_sdr_arrays.back();
      SDVArray::iterator v_end = sdv_array_ref.end();
      SDRArray::iterator r_end = sdr_array_ref.end();
      last_popped_sdv_array.insert(last_popped_sdv_array.begin(),
				   v_end - num_pop_pts, v_end);
      last_popped_sdr_array.insert(last_popped_sdr_array.begin(),
				   r_end - num_pop_pts, r_end);
    }
    size_t new_size = num_ref_pts - num_pop_pts;
    sdv_array_ref.resize(new_size); sdr_array_ref.resize(new_size);

    // TO DO: prune failedRespData[key] or leave in map ?
    data_checks(); // from scratch for now...
  }
  pop_count_stack.pop_back();
}


inline void SurrogateData::push(size_t index, bool erase_popped)
{
  const UShortArray& key = sdRep->activeKey;
  SDVArrayDeque& popped_sdv_arrays = sdRep->poppedVarsData[key];
  SDRArrayDeque& popped_sdr_arrays = sdRep->poppedRespData[key];
  size_t num_popped
    = std::min(popped_sdv_arrays.size(), popped_sdr_arrays.size());

  // harden logic for case of an empty SurrogateData (e.g.,
  // distinct discrepancy at level 0)
  if (num_popped > index) {
    SDVArrayDeque::iterator vit = popped_sdv_arrays.begin();
    SDRArrayDeque::iterator rit = popped_sdr_arrays.begin();
    std::advance(vit, index); std::advance(rit, index);
    size_t num_pts = std::min(vit->size(), rit->size());

    SDVArray& sdv_array_ref = sdRep->varsDataIter->second;
    SDRArray& sdr_array_ref = sdRep->respDataIter->second;
    sdv_array_ref.insert(sdv_array_ref.end(), vit->begin(), vit->end());
    sdr_array_ref.insert(sdr_array_ref.end(), rit->begin(), rit->end());

    // TO DO: update failedRespData[activeKey] ?
    data_checks(); // from scratch for now...

    if (erase_popped)
      { popped_sdv_arrays.erase(vit); popped_sdr_arrays.erase(rit); }

    sdRep->popCountStack[key].push_back(num_pts);
  }
  else if (num_popped) { // not empty
    PCerr << "Error: index out of range for active popped arrays in "
	  << "SurrogateData::push()." << std::endl;
    abort_handler(-1);
  }
  // else ignore push request for empty popped (SurrogateData assumed inactive)
}


inline void SurrogateData::pop_count(size_t count)
{ sdRep->popCountStack[sdRep->activeKey].push_back(count); }


inline size_t SurrogateData::pop_count() const
{
  std::map<UShortArray, SizetArray>::iterator pop_it
    = sdRep->popCountStack.find(sdRep->activeKey);
  return (pop_it == sdRep->popCountStack.end() || pop_it->second.empty()) ?
    _NPOS : pop_it->second.back();
}


inline bool SurrogateData::anchor() const
{
  std::map<UShortArray, size_t>::iterator anchor_it
    = sdRep->anchorIndex.find(sdRep->activeKey);
  return (anchor_it == sdRep->anchorIndex.end() || anchor_it->second == _NPOS)
    ? false : true;
}


inline size_t SurrogateData::points() const
{
  return std::min(sdRep->varsDataIter->second.size(),
		  sdRep->respDataIter->second.size());
}


inline size_t SurrogateData::response_size() const
{
  const SDRArray& sdr_array = sdRep->respDataIter->second;
  size_t i, data_size = 0, num_resp = sdr_array.size(), nh;
  short active_bits;
  for (i=0; i<num_resp; ++i) {
    const SurrogateDataResp& sdr = sdr_array[i];
    active_bits = sdr.active_bits();
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
  const SDRArray& sdr_array = sdRep->respDataIter->second;
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
  const UShortArray& key = sdRep->activeKey;
  return std::min(sdRep->poppedVarsData[key].size(),
		  sdRep->poppedRespData[key].size());
}


inline size_t SurrogateData::num_gradient_variables() const
{
  const SDRArray& sdr_array = sdRep->respDataIter->second;
  return (sdr_array.empty()) ? 0 : sdr_array[0].response_gradient().length();
}


inline size_t SurrogateData::num_hessian_variables() const
{
  const SDRArray& sdr_array = sdRep->respDataIter->second;
  return (sdr_array.empty()) ? 0 : sdr_array[0].response_hessian().numRows();
}


inline size_t SurrogateData::num_derivative_variables() const
{
  size_t num_grad_vars = num_gradient_variables();
  if (num_grad_vars) return num_grad_vars;           // precedence
  else               return num_hessian_variables(); // fall-back
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
  SizetShortMap failed_resp;
  const SDRArray& resp_data = sdRep->respDataIter->second;
  size_t i, num_resp = resp_data.size(); short failed_data;
  for (i=0; i<num_resp; ++i) {
    response_check(resp_data[i], failed_data);
    if (failed_data)
      failed_resp[i] = failed_data;
  }

  if (!failed_resp.empty()) {
    failed_response_data(failed_resp);
#ifdef DEBUG
    PCout << "failedRespData:\n";
    for (SizetShortMap::iterator it=failed_resp_map.begin();
	 it!=failed_resp_map.end(); ++it)
      PCout << "index: " << std::setw(6) << it->first
	    << " data: " << it->second << '\n';
#endif // DEBUG
  }
}


inline short SurrogateData::failed_anchor_data() const
{
  std::map<UShortArray, SizetShortMap>::const_iterator cit1
    = sdRep->failedRespData.find(sdRep->activeKey);
  if (cit1 == sdRep->failedRespData.end()) return 0;
  else {
    const SizetShortMap& failed_resp_data = cit1->second;
    SizetShortMap::const_iterator cit2 = failed_resp_data.find(anchor_index());
    return (cit2 == failed_resp_data.end()) ? 0 : cit2->second;
  }
}


//inline void SurrogateData::failed_anchor_data(short fail_anchor)
//{ sdRep->failedAnchorData = fail_anchor; }


inline void SurrogateData::failed_response_data(const SizetShortMap& fail_data)
{ sdRep->failedRespData[sdRep->activeKey] = fail_data; }


inline const SizetShortMap& SurrogateData::failed_response_data() const
{ return sdRep->failedRespData[sdRep->activeKey]; }


inline void SurrogateData::
failed_response_data_map(const std::map<UShortArray, SizetShortMap>& fail_resp)
{ sdRep->failedRespData = fail_resp; }


inline const std::map<UShortArray, SizetShortMap>& SurrogateData::
failed_response_data_map() const
{ return sdRep->failedRespData; }


inline const SDVArrayDeque& SurrogateData::popped_variables() const
{
  std::map<UShortArray, SDVArrayDeque>::const_iterator cit
    = sdRep->poppedVarsData.find(sdRep->activeKey);
  if (cit == sdRep->poppedVarsData.end()) {
    PCerr << "Error: active key lookup failure in SurrogateData::"
	  << "popped_variables()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline void SurrogateData::popped_variables(const SDVArrayDeque& popped_vars)
{ sdRep->poppedVarsData[sdRep->activeKey] = popped_vars; }


inline const SDRArrayDeque& SurrogateData::popped_response() const
{
  std::map<UShortArray, SDRArrayDeque>::const_iterator cit
    = sdRep->poppedRespData.find(sdRep->activeKey);
  if (cit == sdRep->poppedRespData.end()) {
    PCerr << "Error: active key lookup failure in SurrogateData::"
	  << "popped_response()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline void SurrogateData::
popped_response(const SDRArrayDeque& popped_resp)
{ sdRep->poppedRespData[sdRep->activeKey] = popped_resp; }


inline const std::map<UShortArray, SDVArrayDeque>& SurrogateData::
popped_variables_map() const
{ return sdRep->poppedVarsData; }


inline void SurrogateData::
popped_variables_map(const std::map<UShortArray, SDVArrayDeque>& popped_vars)
{ sdRep->poppedVarsData = popped_vars; }


inline const std::map<UShortArray, SDRArrayDeque>& SurrogateData::
popped_response_map() const
{ return sdRep->poppedRespData; }


inline void SurrogateData::
popped_response_map(const std::map<UShortArray, SDRArrayDeque>& popped_resp)
{ sdRep->poppedRespData = popped_resp; }


inline void SurrogateData::
copy(const SurrogateData& sd, short sdv_mode, short sdr_mode) const
{
  if (sdv_mode == DEEP_COPY) {
    size_t i, j, num_pts, num_sdva;
    const std::map<UShortArray, SDVArray>& vars_map = sd.variables_data_map();
    std::map<UShortArray, SDVArray>::const_iterator v_cit;
    std::map<UShortArray, SDVArray>& new_vars_map = sdRep->varsData;
    new_vars_map.clear();
    for (v_cit = vars_map.begin(); v_cit != vars_map.end(); ++v_cit) {
      const SDVArray& sdv_array = v_cit->second;
      num_pts = sdv_array.size();
      SDVArray new_sdv_array(num_pts);
      for (i=0; i<num_pts; ++i)
	new_sdv_array[i] = sdv_array[i].copy();
      new_vars_map[v_cit->first] = sdv_array;
    }

    const std::map<UShortArray, SDVArrayDeque>& popped_vars_map
      = sd.popped_variables_map();
    std::map<UShortArray, SDVArrayDeque>::const_iterator v2_cit;
    std::map<UShortArray, SDVArrayDeque>& new_popped_vars_map
      = sdRep->poppedVarsData;
    new_popped_vars_map.clear();
    for (v2_cit = popped_vars_map.begin(); v2_cit != popped_vars_map.end();
	 ++v2_cit) {
      const SDVArrayDeque& sdv_2d_array = v2_cit->second;
      num_sdva = sdv_2d_array.size();
      SDVArrayDeque new_popped_sdv_2d(num_sdva);
      for (i=0; i<num_sdva; ++i) {
	const SDVArray& sdva_i = sdv_2d_array[i];
	num_pts = sdva_i.size();
	new_popped_sdv_2d[i].resize(num_pts);
	for (j=0; j<num_pts; ++j)
	  new_popped_sdv_2d[i][j] = sdva_i[j].copy();
      }
      new_popped_vars_map[v2_cit->first] = new_popped_sdv_2d;
    }
  }
  else { // shallow SDV copies based on operator=
    sdRep->varsData       = sd.variables_data_map();
    sdRep->poppedVarsData = sd.popped_variables_map();
  }

  if (sdr_mode == DEEP_COPY) {
    size_t i, j, num_pts, num_sdra;
    const std::map<UShortArray, SDRArray>& resp_map = sd.response_data_map();
    std::map<UShortArray, SDRArray>::const_iterator r_cit;
    std::map<UShortArray, SDRArray>& new_resp_map = sdRep->respData;
    new_resp_map.clear();
    for (r_cit = resp_map.begin(); r_cit != resp_map.end(); ++r_cit) {
      const SDRArray& sdr_array = r_cit->second;
      num_pts = sdr_array.size();
      SDRArray new_sdr_array(num_pts);
      for (i=0; i<num_pts; ++i)
	new_sdr_array[i] = sdr_array[i].copy();
      new_resp_map[r_cit->first] = new_sdr_array;
    }

    const std::map<UShortArray, SDRArrayDeque>& popped_resp_map
      = sd.popped_response_map();
    std::map<UShortArray, SDRArrayDeque>::const_iterator r2_cit;
    std::map<UShortArray, SDRArrayDeque>& new_popped_resp_map
      = sdRep->poppedRespData;
    new_popped_resp_map.clear();
    for (r2_cit = popped_resp_map.begin(); r2_cit != popped_resp_map.end();
	 ++r2_cit) {
      const SDRArrayDeque& sdr_2d_array = r2_cit->second;
      num_sdra = sdr_2d_array.size();
      SDRArrayDeque new_popped_sdr_2d(num_sdra);
      for (i=0; i<num_sdra; ++i) {
	const SDRArray& sdra_i = sdr_2d_array[i];
	num_pts = sdra_i.size();
	new_popped_sdr_2d[i].resize(num_pts);
	for (j=0; j<num_pts; ++j)
	  new_popped_sdr_2d[i][j] = sdra_i[j].copy();
      }
      new_popped_resp_map[r2_cit->first] = new_popped_sdr_2d;
    }
  }
  else { // shallow SDR copies based on operator=
    sdRep->respData       = sd.response_data_map();
    sdRep->poppedRespData = sd.popped_response_map();
  }

  sdRep->failedRespData = sd.failed_response_data_map();
  sdRep->activeKey      = sd.active_key();
}


inline SurrogateData SurrogateData::copy(short sdv_mode, short sdr_mode) const
{
  SurrogateData sd;
  sd.copy(*this, sdv_mode, sdr_mode);
  return sd;
}


inline void SurrogateData::
copy_active(const SurrogateData& sd, short sdv_mode, short sdr_mode) const
{
  const UShortArray&           key = sdRep->activeKey;
  SDVArray&          new_sdv_array = sdRep->varsDataIter->second;
  SDVArrayDeque& new_pop_sdv_array = sdRep->poppedVarsData[key];
  if (sdv_mode == DEEP_COPY) {
    const SDVArray&          sdv_array = sd.variables_data();
    const SDVArrayDeque& pop_sdv_array = sd.popped_variables();
    size_t i, j, num_pts = sdv_array.size();
    new_sdv_array.resize(num_pts);
    for (i=0; i<num_pts; ++i)
      new_sdv_array[i] = sdv_array[i].copy();

    size_t num_sdva = pop_sdv_array.size();
    new_pop_sdv_array.resize(num_sdva);
    for (i=0; i<num_sdva; ++i) {
      const SDVArray& sdva_i = pop_sdv_array[i];
      num_pts = sdva_i.size();
      new_pop_sdv_array[i].resize(num_pts);
      for (j=0; j<num_pts; ++j)
	new_pop_sdv_array[i][j] = sdva_i[j].copy();
    }
  }
  else { // shallow SDV copies based on operator=
    new_sdv_array     = sd.variables_data();
    new_pop_sdv_array = sd.popped_variables();
  }

  SDRArray&          new_sdr_array = sdRep->respDataIter->second;
  SDRArrayDeque& new_pop_sdr_array = sdRep->poppedRespData[key];
  if (sdr_mode == DEEP_COPY) {
    const SDRArray&          sdr_array = sd.response_data();
    const SDRArrayDeque& pop_sdr_array = sd.popped_response();
    size_t i, j, num_pts = sdr_array.size();
    new_sdr_array.resize(num_pts);
    for (i=0; i<num_pts; ++i)
      new_sdr_array[i] = sdr_array[i].copy();

    size_t num_sdra = pop_sdr_array.size();
    new_pop_sdr_array.resize(num_sdra);
    for (i=0; i<num_sdra; ++i) {
      const SDRArray& sdra_i = pop_sdr_array[i];
      num_pts = sdra_i.size();
      new_pop_sdr_array[i].resize(num_pts);
      for (j=0; j<num_pts; ++j)
	new_pop_sdr_array[i][j] = sdra_i[j].copy();
    }
  }
  else { // shallow SDR copies based on operator=
    new_sdr_array     = sd.response_data();
    new_pop_sdr_array = sd.popped_response();
  }

  sdRep->failedRespData[key] = sd.failed_response_data();
}


inline void SurrogateData::clear_active()
{
  /*
  // Too aggressive due to DataFitSurrModel::build_approximation() call to
  // approxInterface.clear_current();
  const UShortArray& key = sdRep->activeKey;
  sdRep->varsData.erase(key);
  sdRep->varsDataIter = sdRep->varsData.end();
  sdRep->respData.erase(key);
  sdRep->respDataIter = sdRep->respData.end();
  sdRep->anchorIndex.erase(key);
  sdRep->failedRespData.erase(key);
  */

  // Retain valid {vars,Resp}DataIter when clearing {vars,resp}Data:
  sdRep->varsDataIter->second.clear();
  sdRep->respDataIter->second.clear();
  // anchorIndex and failedRespData can be pruned:
  const UShortArray& key = sdRep->activeKey;
  sdRep->anchorIndex.erase(key);   //sdRep->anchorIndex[key] = _NPOS;
  sdRep->failedRespData.erase(key);//sdRep->failedRespData[key].clear();
}


inline void SurrogateData::clear_inactive()
{
  std::map<UShortArray, SDVArray>::iterator vd_it = sdRep->varsData.begin();
  std::map<UShortArray, SDRArray>::iterator rd_it = sdRep->respData.begin();
  while (vd_it != sdRep->varsData.end())
    if (vd_it == sdRep->varsDataIter) // preserve active
      { ++vd_it; ++rd_it; }
    else {                            //  clear inactive
      const UShortArray& key = vd_it->first;
      // anchorIndex and failedRespData can be pruned:
      sdRep->anchorIndex.erase(key);    // if it exists
      sdRep->failedRespData.erase(key); // if it exists
      // Be more conservative with clearing {vars,resp}Data:
      vd_it->second.clear(); ++vd_it;
      rd_it->second.clear(); ++rd_it;
      // Too aggressive. Note: postfix increments manage iterator invalidations
      //sdRep->varsData.erase(vd_it++); sdRep->respData.erase(rd_it++);
    }
}


inline void SurrogateData::clear_active_popped()
{
  const UShortArray& key = sdRep->activeKey;
  // Ok to prune as will be recreated in pop() if needed:
  sdRep->poppedVarsData.erase(key);//sdRep->poppedVarsData[key].clear();
  sdRep->poppedRespData.erase(key);//sdRep->poppedRespData[key].clear();
  sdRep->popCountStack.erase(key); //sdRep->popCountStack[key].clear();
}


inline void SurrogateData::clear_popped()
{
  sdRep->poppedVarsData.clear();
  sdRep->poppedRespData.clear();
  sdRep->popCountStack.clear();
}


inline void SurrogateData::clear_keys()
{
  sdRep->activeKey.clear();
  sdRep->varsData.clear(); sdRep->varsDataIter = sdRep->varsData.end();
  sdRep->respData.clear(); sdRep->respDataIter = sdRep->respData.end();
  sdRep->anchorIndex.clear();
  sdRep->failedRespData.clear();

  clear_popped();
}


inline SurrogateDataRep* SurrogateData::data_rep() const
{ return sdRep; }


inline bool SurrogateData::is_null() const
{ return (sdRep) ? false : true; }

} // namespace Pecos

#endif
