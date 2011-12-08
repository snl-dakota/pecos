/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef SURROGATE_DATA_HPP
#define SURROGATE_DATA_HPP

#include "pecos_data_types.hpp"
#include "pecos_stat_util.hpp" // for boost::math::isfinite


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
  SurrogateDataVarsRep(const RealVector& c_vars, short mode);
  /// alternate constructor (data sizing only)
  SurrogateDataVarsRep(size_t num_vars);
  /// destructor
  ~SurrogateDataVarsRep();

  //
  //- Heading: Private data members
  //

  RealVector continuousVars; ///< continuous variables
  //IntVector discreteIntVars;
  //RealVector discreteRealVars;
  int referenceCount;        ///< number of handle objects sharing sdvRep
};


inline SurrogateDataVarsRep::
SurrogateDataVarsRep(const RealVector& c_vars, short mode): referenceCount(1)
{
  // Note: provided a way to query DataAccess mode for c_vars, could make
  // greater use of operator= for {DEEP,SHALLOW}_COPY modes
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(c_vars, continuousVars);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    continuousVars = RealVector(Teuchos::View,c_vars.values(),c_vars.length());
  else                           // default: assume existing Copy/View state
    continuousVars = c_vars;
}


inline SurrogateDataVarsRep::
SurrogateDataVarsRep(size_t num_vars): referenceCount(1)
{ continuousVars.sizeUninitialized(num_vars); }


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
  SurrogateDataVars(const RealVector& c_vars, short mode = DEFAULT_COPY);
  /// alternate constructor (data sizing only)
  SurrogateDataVars(size_t num_vars);
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

  /// set i^{th} entry within continuousVars
  void continuous_variable(const Real& c_var, size_t i);
  /// set continuousVars
  void continuous_variables(const RealVector& c_vars,
			    short mode = DEFAULT_COPY);
  /// get continuousVars
  const RealVector& continuous_variables() const;
  /// get view of continuousVars for updating in place
  RealVector continuous_variables_view();

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
SurrogateDataVars(const RealVector& c_vars, short mode):
  sdvRep(new SurrogateDataVarsRep(c_vars, mode))
{ }


inline SurrogateDataVars::SurrogateDataVars(size_t num_vars):
  sdvRep(new SurrogateDataVarsRep(num_vars))
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
  // Decrement old
  if (sdvRep) // Check for NULL
    if ( --sdvRep->referenceCount == 0 ) 
      delete sdvRep;
  // Increment new
  sdvRep = sdv.sdvRep;
  if (sdvRep) // Check for an assignment of NULL
    ++sdvRep->referenceCount;
  return *this;
}


//inline bool SurrogateDataVars::operator==(const SurrogateDataVars& sdv) const
//{ return (sdvRep->continuousVars==sdv.sdvRep->continuousVars) ? true : false;}


inline void SurrogateDataVars::continuous_variable(const Real& c_var, size_t i)
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
  SurrogateDataRespRep(const Real& fn_val, const RealVector& fn_grad,
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
SurrogateDataRespRep(const Real& fn_val, const RealVector& fn_grad,
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
  SurrogateDataResp(const Real& fn_val, const RealVector& fn_grad,
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

  /// set responseFn
  void response_function(const Real& fn);
  /// get responseFn
  const Real& response_function() const;
  /// get "view" of responseFn for updating in place
  Real& response_function_view();

  /// set i^{th} entry within responseGrad
  void response_gradient(const Real& grad_i, size_t i);
  /// set responseGrad
  void response_gradient(const RealVector& grad, short mode = DEFAULT_COPY);
  /// get responseGrad
  const RealVector& response_gradient() const;
  /// get view of responseGrad for updating in place
  RealVector response_gradient_view();

  /// set i-j^{th} entry within responseHess
  void response_hessian(const Real& hess_ij, size_t i, size_t j);
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
SurrogateDataResp(const Real& fn_val, const RealVector& fn_grad,
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
  // Decrement old
  if (sdrRep) // Check for NULL
    if ( --sdrRep->referenceCount == 0 ) 
      delete sdrRep;
  // Increment new
  sdrRep = sdr.sdrRep;
  if (sdrRep) // Check for an assignment of NULL
    ++sdrRep->referenceCount;
  return *this;
}


//inline bool SurrogateDataResp::operator==(const SurrogateDataResp& sdr) const
//{
//  return ( sdrRep->responseFn   == sdr.sdrRep->responseFn   &&
//	     sdrRep->responseGrad == sdr.sdrRep->responseGrad &&
//	     sdrRep->responseHess == sdr.sdrRep->responseHess ) ? true : false;
//}


inline void SurrogateDataResp::response_function(const Real& fn)
{ sdrRep->responseFn = fn; }


inline const Real& SurrogateDataResp::response_function() const
{ return sdrRep->responseFn; }


inline Real& SurrogateDataResp::response_function_view()
{ return sdrRep->responseFn; }


inline void SurrogateDataResp::
response_gradient(const Real& grad_i, size_t i)
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
response_hessian(const Real& hess_ij, size_t i, size_t j)
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
    s << "function value = " << sdrRep->responseFn << '\n';
  if (sdrRep->activeBits & 2)
    { s << "function gradient =\n"; write_data(s, sdrRep->responseGrad); }
  if (sdrRep->activeBits & 4) {
    s << "function Hessian =\n";
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

  /// a special variables sample (often at the center of the
  /// approximation region) for which exact matching is enforced
  /// (e.g., using equality-constrained least squares regression).
  SurrogateDataVars anchorVars;
  /// a set of variables samples used to build the approximation.
  /// These sample points are fit approximately (e.g., using least
  /// squares regression); exact matching is not enforced.
  SDVArray varsData;
  /// sets of variables samples that have been popped off varsData but
  /// which are available for future restoration.  The granularity of
  /// this 2D array corresponds to multiple trial sets, each
  /// contributing a SDVArray.
  SDV2DArray savedVarsTrials;
  /// a set of variables samples that have been stored for future
  /// restoration.  The granularity of this array corresponds to a
  /// wholesale storage of the current varsData state.
  SDVArray storedVarsData;

  /// a special response sample (often at the center of the
  /// approximation region) for which exact matching is enforced
  /// (e.g., using equality-constrained least squares regression).
  SurrogateDataResp anchorResp;
  /// a set of response samples used to build the approximation.  These
  /// sample points are fit approximately (e.g., using least squares
  /// regression); exact matching is not enforced.
  SDRArray respData;
  /// sets of response samples that have been popped off respData but
  /// which are available for future restoration.  The granularity of
  /// this 2D array corresponds to multiple trial sets, each
  /// contributing a SDRArray.
  SDR2DArray savedRespTrials;
  /// a set of response samples that have been stored for future
  /// restoration.  The granularity of this array corresponds to a
  /// wholesale storage of the current respData state.
  SDRArray storedRespData;

  /// failed anchor data bits; defined in sample_checks() and used for
  /// fault tolerance in regression() and expectation()
  short failedAnchorData;
  /// map from failed respData indices to failed data bits; defined
  /// in sample_checks() and used for fault tolerance
  SizetShortMap failedRespData;

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

  /// get anchorVars.continuous_variables()
  const RealVector& anchor_continuous_variables() const;
  /// get anchorResp.response_function()
  const Real& anchor_function() const;
  /// get anchorResp.response_gradient()
  const RealVector& anchor_gradient() const;
  /// get anchorResp.response_hessian()
  const RealSymMatrix& anchor_hessian() const;

  /// get varsData[i].continuous_variables()
  const RealVector& continuous_variables(size_t i) const;
  /// get respData[i].response_function()
  const Real& response_function(size_t i) const;
  /// get respData[i].response_gradient()
  const RealVector& response_gradient(size_t i) const;
  /// get respData[i].response_hessian()
  const RealSymMatrix& response_hessian(size_t i) const;

  /// push sdv onto end of varsData
  void push_back(const SurrogateDataVars& sdv);
  /// push sdr onto end of respData
  void push_back(const SurrogateDataResp& sdr);
  /// push {sdv,sdr} onto ends of {vars,resp}Data
  void push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr);
  /// remove num_pop_pts entries from ends of {vars,resp}Data
  void pop(size_t num_pop_pts, bool save_data = true);
  /// move all entries from {vars,resp}Data to saved{Vars,Resp}Data
  void store();
  /// move all entries from saved{Vars,Resp}Data to {vars,resp}Data
  void restore();
  /// remove num_pop_pts entries from ends of {vars,resp}Data
  size_t restore(size_t index, bool erase_saved = true);

  /// query presence of anchor{Vars,Resp}
  bool anchor() const;
  /// return size of {vars,resp}Data arrays (neglecting anchor point)
  size_t /*data_*/size() const;
  /// return number of 1D arrays within saved{Vars,Resp}Trials 2D arrays
  size_t saved_trials() const;
  /// return number of derivative variables as indicated by size of
  /// gradient/Hessian arrays
  size_t num_derivative_variables() const;

  /// screen data sets for samples with Inf/Nan that should be excluded;
  /// defines failedAnchorData and failedRespData
  void data_checks();
  /// return failedAnchorData
  short failed_anchor_data();
  /// return failedRespData
  const SizetShortMap& failed_response_data();

  /// clear anchor{Vars,Resp}
  void clear_anchor();
  /// clear {vars,resp}Data
  void clear_data();
  /// clear saved{Vars,Resp}Trials
  void clear_saved();
  /// clear stored{Vars,Resp}Data
  void clear_stored();

private:

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
  // Decrement old
  if (sdRep) // Check for NULL
    if ( --sdRep->referenceCount == 0 ) 
      delete sdRep;
  // Increment new
  sdRep = sd.sdRep;
  if (sdRep) // Check for an assignment of NULL
    ++sdRep->referenceCount;
  return *this;
}


inline void SurrogateData::
anchor_point(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr)
{ sdRep->anchorVars = sdv; sdRep->anchorResp = sdr; }


inline void SurrogateData::
data_points(const SDVArray& sdv_array, const SDRArray& sdr_array)
{ sdRep->varsData = sdv_array; sdRep->respData = sdr_array; }


inline void SurrogateData::anchor_variables(const SurrogateDataVars& sdv)
{ sdRep->anchorVars = sdv; }


inline const SurrogateDataVars& SurrogateData::anchor_variables() const
{ return sdRep->anchorVars; }


inline void SurrogateData::anchor_response(const SurrogateDataResp& sdr)
{ sdRep->anchorResp = sdr; }


inline const SurrogateDataResp& SurrogateData::anchor_response() const
{ return sdRep->anchorResp; }


inline void SurrogateData::variables_data(const SDVArray& sdv_array)
{ sdRep->varsData = sdv_array; }


inline const SDVArray& SurrogateData::variables_data() const
{ return sdRep->varsData; }


inline void SurrogateData::response_data(const SDRArray& sdr_array)
{ sdRep->respData = sdr_array; }


inline const SDRArray& SurrogateData::response_data() const
{ return sdRep->respData; }


inline const RealVector& SurrogateData::anchor_continuous_variables() const
{ return sdRep->anchorVars.continuous_variables(); }


inline const RealVector& SurrogateData::continuous_variables(size_t i) const
{ return sdRep->varsData[i].continuous_variables(); }


inline const Real& SurrogateData::anchor_function() const
{ return sdRep->anchorResp.response_function(); }


inline const RealVector& SurrogateData::anchor_gradient() const
{ return sdRep->anchorResp.response_gradient(); }


inline const RealSymMatrix& SurrogateData::anchor_hessian() const
{ return sdRep->anchorResp.response_hessian(); }


inline const Real& SurrogateData::response_function(size_t i) const
{ return sdRep->respData[i].response_function(); }


inline const RealVector& SurrogateData::response_gradient(size_t i) const
{ return sdRep->respData[i].response_gradient(); }


inline const RealSymMatrix& SurrogateData::response_hessian(size_t i) const
{ return sdRep->respData[i].response_hessian(); }


inline void SurrogateData::push_back(const SurrogateDataVars& sdv)
{ sdRep->varsData.push_back(sdv); }


inline void SurrogateData::push_back(const SurrogateDataResp& sdr)
{ sdRep->respData.push_back(sdr); }


inline void SurrogateData::
push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr)
{ sdRep->varsData.push_back(sdv); sdRep->respData.push_back(sdr); }


inline void SurrogateData::pop(size_t num_pop_pts, bool save_data)
{
  if (num_pop_pts) {
    size_t data_size = size();
    if (data_size >= num_pop_pts) {
      if (save_data) {
	// append empty arrays and then update them in place
	SDVArray sdv_array; SDRArray sdr_array;
	sdRep->savedVarsTrials.push_back(sdv_array);
	sdRep->savedRespTrials.push_back(sdr_array);
	SDVArray& last_sdv_array = sdRep->savedVarsTrials.back();
	SDRArray& last_sdr_array = sdRep->savedRespTrials.back();
	/*
	// prevent underflow portability issue w/ compiler coercion of -num_pop
	SDVArray::difference_type reverse_adv_vars = num_pop_pts;
	SDRArray::difference_type reverse_adv_resp = num_pop_pts;
	SDVIter vit = varsData.end(); std::advance(vit, -reverse_adv_vars);
	SDRIter rit = respData.end(); std::advance(rit, -reverse_adv_resp);
	*/
	last_sdv_array.insert(last_sdv_array.begin(), //vit,
			      sdRep->varsData.end() - num_pop_pts,
			      sdRep->varsData.end());
	last_sdr_array.insert(last_sdr_array.begin(), //rit,
			      sdRep->respData.end() - num_pop_pts,
			      sdRep->respData.end());
      }
      size_t new_size = data_size - num_pop_pts;
      sdRep->varsData.resize(new_size); sdRep->respData.resize(new_size);
    }
    else {
      PCerr << "Error: pop count (" << num_pop_pts << ") exceeds data size ("
	    << data_size << ") in SurrogateData::pop(size_t)." << std::endl;
      abort_handler(-1);
    }
  }
}


inline size_t SurrogateData::restore(size_t index, bool erase_saved)
{
  SDV2DArray::iterator vit = sdRep->savedVarsTrials.begin();
  SDR2DArray::iterator rit = sdRep->savedRespTrials.begin();
  std::advance(vit, index); std::advance(rit, index);
  size_t num_pts = std::min(vit->size(), rit->size());
  sdRep->varsData.insert(sdRep->varsData.end(), vit->begin(), vit->end());
  sdRep->respData.insert(sdRep->respData.end(), rit->begin(), rit->end());
  if (erase_saved)
    { sdRep->savedVarsTrials.erase(vit); sdRep->savedRespTrials.erase(rit); }
  return num_pts;
}


inline void SurrogateData::store()
{
  sdRep->storedVarsData = sdRep->varsData; // shallow copies
  sdRep->storedRespData = sdRep->respData; // shallow copies
  clear_data();
}


inline void SurrogateData::restore()
{
  sdRep->varsData.insert(sdRep->varsData.end(), sdRep->storedVarsData.begin(),
			 sdRep->storedVarsData.end());
  sdRep->respData.insert(sdRep->respData.end(), sdRep->storedRespData.begin(),
			 sdRep->storedRespData.end());
  clear_stored();
}


inline bool SurrogateData::anchor() const
{ return (!sdRep->anchorVars.is_null() && !sdRep->anchorResp.is_null()); }


inline size_t SurrogateData::/*data_*/size() const
{ return std::min(sdRep->varsData.size(), sdRep->respData.size()); }


inline size_t SurrogateData::saved_trials() const
{
  return std::min(sdRep->savedVarsTrials.size(), sdRep->savedRespTrials.size());
}


inline size_t SurrogateData::num_derivative_variables() const
{
  return (sdRep->anchorResp.is_null()) ?
    sdRep->respData[0].response_gradient().length() :
    sdRep->anchorResp.response_gradient().length();
}


inline void SurrogateData::data_checks()
{
  // We take a conservative approach of rejecting all data of derivative
  // order greater than or equal to a detected failure:

  size_t i, j, num_resp = sdRep->respData.size(),
    num_deriv_vars = num_derivative_variables();

  sdRep->failedAnchorData = 0;
  if (anchor()) {
    if (!bmth::isfinite(anchor_function()))
      sdRep->failedAnchorData = (num_deriv_vars) ? 3 : 1; // conservative
    else if (num_deriv_vars) {
      const RealVector& grad = anchor_gradient();
      for (j=0; j<num_deriv_vars; ++j)
	if (!bmth::isfinite(grad[j]))
	  { sdRep->failedAnchorData = 2; break;}
    }
  }

  sdRep->failedRespData.clear();
  for (i=0; i<num_resp; ++i) {
    if (!bmth::isfinite(response_function(i)))
      sdRep->failedRespData[i] = (num_deriv_vars) ? 3 : 1; // conservative
    else if (num_deriv_vars) {
      const RealVector& grad_i = response_gradient(i);
      for (j=0; j<num_deriv_vars; ++j)
	if (!bmth::isfinite(grad_i[j]))
	  { sdRep->failedRespData[i] = 2; break;}
    }
  }

#ifdef DEBUG
  if (sdRep->failedAnchorData) {
    PCout << "failedAnchorData = " << sdRep->failedAnchorData << '\n';
  if (!sdRep->failedRespData.empty()) {
    PCout << "failedRespData:\n";
    for (SizetShortMap::iterator it=sdRep->failedRespData.begin();
	 it!=sdRep->failedRespData.end(); ++it)
      PCout << "index: " << std::setw(6) << it->first
	    << " data: " << it->second << '\n';
  }
#endif // DEBUG
}


inline short SurrogateData::failed_anchor_data()
{ return sdRep->failedAnchorData; }


inline const SizetShortMap& SurrogateData::failed_response_data()
{ return sdRep->failedRespData; }


inline void SurrogateData::clear_anchor()
{
  sdRep->anchorVars = SurrogateDataVars();
  sdRep->anchorResp = SurrogateDataResp();
}


inline void SurrogateData::clear_data()
{ sdRep->varsData.clear(); sdRep->respData.clear(); }


inline void SurrogateData::clear_saved()
{ sdRep->savedVarsTrials.clear(); sdRep->savedRespTrials.clear(); }


inline void SurrogateData::clear_stored()
{ sdRep->storedVarsData.clear(); sdRep->storedRespData.clear(); }

} // namespace Pecos

#endif
