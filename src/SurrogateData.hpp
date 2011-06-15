/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef SURROGATE_DATA_HPP
#define SURROGATE_DATA_HPP

#include "pecos_data_types.hpp"


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
  SurrogateDataVarsRep(const RealVector& x, short mode);
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
SurrogateDataVarsRep(const RealVector& x, short mode):
  referenceCount(1)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(x, continuousVars);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    continuousVars = RealVector(Teuchos::View, x.values(), x.length());
  else                           // default: assume existing Copy/View state
    continuousVars = x;
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
  SurrogateDataVars(const RealVector& x, short mode = DEFAULT_COPY);
  /// copy constructor
  SurrogateDataVars(const SurrogateDataVars& sdv);
  /// destructor
  ~SurrogateDataVars();

  /// assignment operator
  SurrogateDataVars& operator=(const SurrogateDataVars& sdv);
  /// equality operator
  bool operator==(const SurrogateDataVars& sdv) const;

  //
  //- Heading: member functions
  //

  const RealVector& continuous_variables() const; ///< return continuousVars

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


inline SurrogateDataVars::SurrogateDataVars(const RealVector& x, short mode):
  sdvRep(new SurrogateDataVarsRep(x, mode))
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


inline bool SurrogateDataVars::operator==(const SurrogateDataVars& sdv) const
{ return (sdvRep->continuousVars == sdv.sdvRep->continuousVars) ? true : false;}


inline const RealVector& SurrogateDataVars::continuous_variables() const
{ return sdvRep->continuousVars; }


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
		       const RealSymMatrix& fn_hess, short mode);
  /// destructor
  ~SurrogateDataRespRep();

  //
  //- Heading: Private data members
  //

  Real          responseFn;   ///< truth response function value
  RealVector    responseGrad; ///< truth response function gradient
  RealSymMatrix responseHess; ///< truth response function Hessian

  /// number of handle objects sharing sdrRep
  int referenceCount;
};


inline SurrogateDataRespRep::
SurrogateDataRespRep(const Real& fn_val, const RealVector& fn_grad,
		     const RealSymMatrix& fn_hess, short mode):
  referenceCount(1)
{
  responseFn = fn_val;              // deep copy for scalars
  if (mode == DEEP_COPY)            // enforce vector/matrix deep copy
    { copy_data(fn_grad, responseGrad); copy_data(fn_hess, responseHess); }
  else if (mode == SHALLOW_COPY) {  // enforce vector/matrix shallow copy
    responseGrad = RealVector(Teuchos::View,fn_grad.values(),fn_grad.length());
    responseHess = RealSymMatrix(Teuchos::View, fn_hess, fn_hess.numRows());
  }
  else                              // default: assume existing Copy/View state
    { responseGrad = fn_grad; responseHess = fn_hess; }
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
		    const RealSymMatrix& fn_hess, short mode = DEFAULT_COPY);
  /// copy constructor
  SurrogateDataResp(const SurrogateDataResp& sdr);
  /// destructor
  ~SurrogateDataResp();

  /// assignment operator
  SurrogateDataResp& operator=(const SurrogateDataResp& sdr);
  /// equality operator
  bool operator==(const SurrogateDataResp& sdr) const;

  //
  //- Heading: member functions
  //

  const Real&          response_function() const; ///< return responseFn
  const RealVector&    response_gradient() const; ///< return responseGrad
  const RealSymMatrix& response_hessian()  const; ///< return responseHess

  /// function to check sdrRep (does this handle contain a body)
  bool is_null() const;

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
		  const RealSymMatrix& fn_hess, short mode):
  sdrRep(new SurrogateDataRespRep(fn_val, fn_grad, fn_hess, mode))
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


inline bool SurrogateDataResp::operator==(const SurrogateDataResp& sdr) const
{
  return ( sdrRep->responseFn   == sdr.sdrRep->responseFn   &&
	   sdrRep->responseGrad == sdr.sdrRep->responseGrad &&
	   sdrRep->responseHess == sdr.sdrRep->responseHess ) ? true : false;
}


inline const Real& SurrogateDataResp::response_function() const
{ return sdrRep->responseFn; }


inline const RealVector& SurrogateDataResp::response_gradient() const
{ return sdrRep->responseGrad; }


inline const RealSymMatrix& SurrogateDataResp::response_hessian() const
{ return sdrRep->responseHess; }


inline bool SurrogateDataResp::is_null() const
{ return (sdrRep) ? false : true; }


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
  /// set of variables samples used to build the approximation.  These
  /// sample points are fit approximately (e.g., using least squares
  /// regression); exact matching is not enforced.
  SDVArray varsData;
  /// set of variables samples that have been popped off varsData but
  /// which are available for future restoration.
  SDV2DArray savedVarsData;

  /// a special response sample (often at the center of the
  /// approximation region) for which exact matching is enforced
  /// (e.g., using equality-constrained least squares regression).
  SurrogateDataResp anchorResp;
  /// set of response samples used to build the approximation.  These
  /// sample points are fit approximately (e.g., using least squares
  /// regression); exact matching is not enforced.
  SDRArray respData;
  /// set of response samples that have been popped off respData but
  /// which are available for future restoration.
  SDR2DArray savedRespData;

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
  /// set varsData
  void variables_data(const SDVArray& sdv_array);
  /// set respData
  void response_data(const SDRArray& sdr_array);

  /// get anchorVars
  const SurrogateDataVars& anchor_variables() const;
  /// get anchorResp
  const SurrogateDataResp& anchor_response() const;
  /// get varsData
  const SDVArray& variables_data() const;
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

  /// push sdv/sdr onto ends of {vars,resp}Data
  void push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr);
  /// remove num_pop_pts entries from ends of {vars,resp}Data
  void pop(size_t num_pop_pts, bool save_data = true);
  /// remove num_pop_pts entries from ends of {vars,resp}Data
  void restore(size_t index, bool erase_saved = true);

  /// query presence of anchor{Vars,Resp}
  bool anchor() const;
  /// return size of {vars,resp}Data arrays (neglecting anchor point)
  size_t /*data_*/size() const;
  /// return size of saved{Vars,Resp}Data arrays
  size_t saved_size() const;
  /// return number of derivative variables as indicated by size of
  /// gradient/Hessian arrays
  size_t num_derivative_variables() const;

  /// clear anchor{Vars,Resp}
  void clear_anchor();
  /// clear {vars,resp}Data
  void clear_data();
  /// clear saved{Vars,Resp}Data
  void clear_saved();

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


inline void SurrogateData::variables_data(const SDVArray& sdv_array)
{ sdRep->varsData = sdv_array; }


inline void SurrogateData::response_data(const SDRArray& sdr_array)
{ sdRep->respData = sdr_array; }


inline const SurrogateDataVars& SurrogateData::anchor_variables() const
{ return sdRep->anchorVars; }


inline const SurrogateDataResp& SurrogateData::anchor_response() const
{ return sdRep->anchorResp; }


inline const SDVArray& SurrogateData::variables_data() const
{ return sdRep->varsData; }


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


inline void SurrogateData::
push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr)
{ sdRep->varsData.push_back(sdv); sdRep->respData.push_back(sdr); }


inline void SurrogateData::pop(size_t num_pop_pts, bool save_data)
{
  if (num_pop_pts) {
    size_t data_size = size();
    if (data_size >= num_pop_pts) {
      if (save_data) {
	SDVArray sdv; sdRep->savedVarsData.push_back(sdv); // append empty array
	SDRArray sdr; sdRep->savedRespData.push_back(sdr); // append empty array
	SDVArray& last_sdv_array = sdRep->savedVarsData.back();//update in place
	SDRArray& last_sdr_array = sdRep->savedRespData.back();//update in place
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


inline void SurrogateData::restore(size_t index, bool erase_saved)
{
  SDV2DArray::iterator vit = sdRep->savedVarsData.begin();
  SDR2DArray::iterator rit = sdRep->savedRespData.begin();
  std::advance(vit, index); std::advance(rit, index);
  sdRep->varsData.insert(sdRep->varsData.end(), vit->begin(), vit->end());
  sdRep->respData.insert(sdRep->respData.end(), rit->begin(), rit->end());
  if (erase_saved)
    { sdRep->savedVarsData.erase(vit); sdRep->savedRespData.erase(rit); }
}


inline bool SurrogateData::anchor() const
{ return (!sdRep->anchorVars.is_null() && !sdRep->anchorResp.is_null()); }


inline size_t SurrogateData::/*data_*/size() const
{ return std::min(sdRep->varsData.size(), sdRep->respData.size()); }


inline size_t SurrogateData::saved_size() const
{ return std::min(sdRep->savedVarsData.size(), sdRep->savedRespData.size()); }


inline size_t SurrogateData::num_derivative_variables() const
{
  return (sdRep->anchorResp.is_null()) ?
    sdRep->respData[0].response_gradient().length() :
    sdRep->anchorResp.response_gradient().length();
}


inline void SurrogateData::clear_anchor()
{
  sdRep->anchorVars = SurrogateDataVars();
  sdRep->anchorResp = SurrogateDataResp();
}


inline void SurrogateData::clear_data()
{ sdRep->varsData.clear(); sdRep->respData.clear(); }


inline void SurrogateData::clear_saved()
{ sdRep->savedVarsData.clear(); sdRep->savedRespData.clear(); }

} // namespace Pecos

#endif
