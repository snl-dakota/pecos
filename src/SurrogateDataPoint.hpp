/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef SURROGATE_DATA_POINT_HPP
#define SURROGATE_DATA_POINT_HPP

#include "pecos_data_types.hpp"


namespace Pecos {


/// The representation of a surrogate data point.  This representation,
/// or body, may be shared by multiple SurrogateDataPoint handle instances.

/** The SurrogateDataPoint/SurrogateDataPointRep pairs utilize a
    handle-body idiom (Coplien, Advanced C++). */

class PECOS_EXPORT SurrogateDataPointRep{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class SurrogateDataPoint;

private:

  //
  //- Heading: Private member functions
  //

  /// constructor
  SurrogateDataPointRep(const RealVector& x, const Real& fn_val,
			const RealVector& fn_grad,
			const RealSymMatrix& fn_hess);
  /// destructor
  ~SurrogateDataPointRep();

  //
  //- Heading: Private data members
  //

  RealVector    continuousVars; ///< continuous variables
  Real          responseFn;     ///< truth response function value
  RealVector    responseGrad;   ///< truth response function gradient
  RealSymMatrix responseHess;   ///< truth response function Hessian

  /// number of handle objects sharing sdpRep
  int referenceCount;
};


inline SurrogateDataPointRep::
SurrogateDataPointRep(const RealVector& x, const Real& fn_val,
		      const RealVector& fn_grad,
		      const RealSymMatrix& fn_hess):
  referenceCount(1)
{
  copy_data(x, continuousVars); // x may be a vector view
  responseFn = fn_val;  responseGrad   = fn_grad; responseHess = fn_hess;
}


inline SurrogateDataPointRep::~SurrogateDataPointRep()
{ }


/// Container class encapsulating basic parameter and response data
/// for defining a "truth" data point.

/** A list of these data points is contained in Approximation instances
    (e.g., Dakota::Approximation::currentPoints) and provides the data
    to build the approximation.  A handle-body idiom is used to avoid
    excessive data copying overhead. */

class PECOS_EXPORT SurrogateDataPoint{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  /// default constructor
  SurrogateDataPoint();
  /// standard constructor
  SurrogateDataPoint(const RealVector& x, const Real& fn_val,
		     const RealVector& fn_grad,
		     const RealSymMatrix& fn_hess);
  /// copy constructor
  SurrogateDataPoint(const SurrogateDataPoint& sdp);
  /// destructor
  ~SurrogateDataPoint();

  /// assignment operator
  SurrogateDataPoint& operator=(const SurrogateDataPoint& sdp);
  /// equality operator
  bool operator==(const SurrogateDataPoint& sdp) const;

  //
  //- Heading: member functions
  //

  const RealVector&    continuous_variables() const; ///< return continuousVars
  const Real&          response_function()    const; ///< return responseFn
  const RealVector&    response_gradient()    const; ///< return responseGrad
  const RealSymMatrix& response_hessian()     const; ///< return responseHess

  /// function to check sdpRep (does this handle contain a body)
  bool is_null() const;

private:

  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  SurrogateDataPointRep* sdpRep;
};


inline SurrogateDataPoint::SurrogateDataPoint(): sdpRep(NULL)
{ }


inline SurrogateDataPoint::
SurrogateDataPoint(const RealVector& x, const Real& fn_val,
		   const RealVector& fn_grad,
		   const RealSymMatrix& fn_hess):
  sdpRep(new SurrogateDataPointRep(x, fn_val, fn_grad, fn_hess))
{ }


inline SurrogateDataPoint::SurrogateDataPoint(const SurrogateDataPoint& sdp)
{
  // Increment new (no old to decrement)
  sdpRep = sdp.sdpRep;
  if (sdpRep) // Check for an assignment of NULL
    sdpRep->referenceCount++;
}


inline SurrogateDataPoint::~SurrogateDataPoint()
{
  if (sdpRep) { // Check for NULL
    --sdpRep->referenceCount; // decrement
    if (sdpRep->referenceCount == 0)
      delete sdpRep;
  }
}


inline SurrogateDataPoint& SurrogateDataPoint::
operator=(const SurrogateDataPoint& sdp)
{
  // Decrement old
  if (sdpRep) // Check for NULL
    if ( --sdpRep->referenceCount == 0 ) 
      delete sdpRep;
  // Increment new
  sdpRep = sdp.sdpRep;
  if (sdpRep) // Check for an assignment of NULL
    sdpRep->referenceCount++;
  return *this;
}


inline bool SurrogateDataPoint::operator==(const SurrogateDataPoint& sdp) const
{
  return ( sdpRep->continuousVars == sdp.sdpRep->continuousVars &&
	   sdpRep->responseFn     == sdp.sdpRep->responseFn     &&
	   sdpRep->responseGrad   == sdp.sdpRep->responseGrad   &&
	   sdpRep->responseHess   == sdp.sdpRep->responseHess ) ? true : false;
}


inline const RealVector& SurrogateDataPoint::continuous_variables() const
{ return sdpRep->continuousVars; }


inline const Real& SurrogateDataPoint::response_function() const
{ return sdpRep->responseFn; }


inline const RealVector& SurrogateDataPoint::response_gradient() const
{ return sdpRep->responseGrad; }


inline const RealSymMatrix& SurrogateDataPoint::response_hessian() const
{ return sdpRep->responseHess; }


inline bool SurrogateDataPoint::is_null() const
{ return (sdpRep) ? false : true; }

} // namespace Pecos

#endif
