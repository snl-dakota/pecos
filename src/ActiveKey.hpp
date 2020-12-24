/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef ACTIVE_KEY_HPP
#define ACTIVE_KEY_HPP

#include "pecos_data_types.hpp"


namespace Pecos {

/// Shared representation for ActiveKey class (body within handle-body idiom).

/** Manages a set of model indices and a set of continuous/discrete
    hyper-parameter (resolution) controls. */

class ActiveKeyRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class ActiveKey;

public:
  /// destructor
  ~ActiveKeyRep();

private:

  //
  //- Heading: Private member functions
  //

  /// default constructor (default empty key)
  ActiveKeyRep();
  /// partial constructor (legacy use case: model indices only)
  ActiveKeyRep(const UShortArray& indices);
  /// full constructor
  ActiveKeyRep(const UShortArray& indices, const RealVector&   c_params,
	       const IntVector& di_params, const SizetVector& ds_params,
	       short mode);

  //
  //- Heading: Private data members
  //

  /// identifies an instance within a multi-dimensional model ensemble
  UShortArray modelIndices;

  // the identified subset of (state) variables that serve as
  // solution control hyper-parameters:
  RealVector  continuousHyperParams;  ///< continuous hyper-parameters
  IntVector   discreteIntHyperParams; ///< discrete int range hyper-parameters
  SizetVector discreteSetHyperParams; ///< discrete set index hyper-parameters
};


inline ActiveKeyRep::~ActiveKeyRep()
{ }


inline ActiveKeyRep::ActiveKeyRep()
{ }


inline ActiveKeyRep::ActiveKeyRep(const UShortArray& indices)
{ modelIndices = indices; }


inline ActiveKeyRep::
ActiveKeyRep(const UShortArray& indices, const RealVector&   c_params,
	     const IntVector& di_params, const SizetVector& ds_params,
	     short mode)
{
  modelIndices = indices;

  // Note: provided a way to query DataAccess mode for c_params, could make
  // greater use of operator= for {DEEP,SHALLOW}_COPY modes
  if (mode == DEEP_COPY) {         // enforce deep vector copy
    if (!c_params.empty())  copy_data( c_params, continuousHyperParams);
    if (!di_params.empty()) copy_data(di_params, discreteIntHyperParams);
    if (!ds_params.empty()) copy_data(ds_params, discreteSetHyperParams);
  }
  else if (mode == SHALLOW_COPY) { // enforce shallow vector copy
    if (!c_params.empty())
      continuousHyperParams
	= RealVector(Teuchos::View,  c_params.values(),  c_params.length());
    if (!di_params.empty())
      discreteIntHyperParams
	= IntVector(Teuchos::View,  di_params.values(), di_params.length());
    if (!ds_params.empty())
      discreteSetHyperParams
	= RealVector(Teuchos::View, ds_params.values(), ds_params.length());
  }
  else {                           // default: assume existing Copy/View state
    if (!c_params.empty())   continuousHyperParams = c_params;
    if (!di_params.empty()) discreteIntHyperParams = di_params;
    if (!ds_params.empty()) discreteSetHyperParams = ds_params;
  }
}


/// Handle class for providing a unique key to a data set comprised of
/// variables and responses.

/** Multiple data sets may co-exist and are keyed based on model
    indices and state hyper-parameters.  A handle-body idiom is used
    to reduce data copying overhead. */

class ActiveKey
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  /// default constructor
  ActiveKey();
  /// partial constructor (legacy use case: model indices only)
  ActiveKey(const UShortArray& indices, short mode = DEFAULT_COPY);
  /// full constructor
  ActiveKey(const UShortArray& indices, const RealVector&   c_params,
	    const IntVector& di_params, const SizetVector& ds_params,
	    short mode = DEFAULT_COPY);

  /// copy constructor
  ActiveKey(const ActiveKey& key);
  /// destructor
  ~ActiveKey();

  /// assignment operator
  ActiveKey& operator=(const ActiveKey& key);
  // equality operator
  //bool operator==(const ActiveKey& key) const;

  //
  //- Heading: member functions
  //

  /// return number of continuous parameters
  size_t chp() const;
  /// return number of discrete integer parameters
  size_t dihp() const;
  /// return number of discrete real parameters
  size_t dshp() const;

  /// return deep copy of ActiveKey instance
  ActiveKey copy() const;

  /// set i^{th} entry within continuousHyperParams
  void continuous_parameter(Real c_param, size_t i);
  /// set continuousHyperParams
  void continuous_parameters(const RealVector& c_params,
			     short mode = DEFAULT_COPY);
  /// get continuousHyperParams
  const RealVector& continuous_parameters() const;
  /// get view of continuousHyperParams for updating in place
  RealVector continuous_parameters_view();

  /// set i^{th} entry within discreteIntHyperParams
  void discrete_int_parameter(int di_param, size_t i);
  /// set discreteIntHyperParams
  void discrete_int_parameters(const IntVector& di_params,
			       short mode = DEFAULT_COPY);
  /// get discreteIntHyperParams
  const IntVector& discrete_int_parameters() const;
  /// get view of discreteIntHyperParams for updating in place
  IntVector discrete_int_parameters_view();

  /// set i^{th} entry within discreteSetHyperParams
  void discrete_set_index(size_t ds_index, size_t i);
  /// set discreteSetHyperParams
  void discrete_set_indices(const SizetVector& ds_indices,
			    short mode = DEFAULT_COPY);
  /// get discreteSetHyperParams
  const SizetVector& discrete_set_indices() const;
  /// get view of discreteSetHyperParams for updating in place
  SizetVector discrete_set_indices_view();

  /// function to check keyRep (does this handle contain a body)
  bool is_null() const;

private:

  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  std::shared_ptr<ActiveKeyRep> keyRep;
};


inline ActiveKey::ActiveKey()
{ } // keyRep is null


// BMA NOTE: The following don't use make_shared<ActiveKeyRep>()
// due to private ctors

inline ActiveKey::ActiveKey(const UShortArray& indices):
  keyRep(new ActiveKeyRep(indices))
{ }


inline ActiveKey::
ActiveKey(const UShortArray& indices, const RealVector&   c_params,
	  const IntVector& di_params, const SizetVector& ds_params, short mode):
  keyRep(new ActiveKeyRep(indices, c_params, di_params, ds_params, mode))
{ }


inline ActiveKey::ActiveKey(const ActiveKey& key):
  keyRep(key.keyRep)  
{ }


inline ActiveKey::~ActiveKey()
{ }


inline ActiveKey& ActiveKey::operator=(const ActiveKey& key)
{
  keyRep = key.keyRep;
  return *this;
}


//inline bool ActiveKey::operator==(const ActiveKey& key) const
//{
//  return (keyRep->modelIndices     == key.keyRep->modelIndices           &&
//    keyRep->continuousHyperParams  == key.keyRep->continuousHyperParams  &&
//    keyRep->discreteIntHyperParams == key.keyRep->discreteIntHyperParams &&
//    keyRep->discreteSetHyperParams == key.keyRep->discreteSetHyperParams);
//}


inline size_t ActiveKey::chp() const
{ return keyRep->continuousHyperParams.length(); }


inline size_t ActiveKey::dihp() const
{ return keyRep->discreteIntHyperParams.length(); }


inline size_t ActiveKey::dshp() const
{ return keyRep->discreteSetHyperParams.length(); }


/// deep copy of ActiveKey instance
inline ActiveKey ActiveKey::copy() const
{
  ActiveKey key(keyRep->modelIndices,           keyRep->continuousHyperParams,
		keyRep->discreteIntHyperParams, keyRep->discreteSetHyperParams,
		DEEP_COPY);
  return key;
}


inline void ActiveKey::continuous_parameter(Real c_param, size_t i)
{ keyRep->continuousHyperParams[i] = c_param; }


inline void ActiveKey::
continuous_parameters(const RealVector& c_params, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(c_params, keyRep->continuousHyperParams);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    keyRep->continuousHyperParams
      = RealVector(Teuchos::View, c_params.values(), c_params.length());
  else                           // default: assume existing Copy/View state
    keyRep->continuousHyperParams = c_params;
}


inline const RealVector& ActiveKey::continuous_parameters() const
{ return keyRep->continuousHyperParams; }


inline RealVector ActiveKey::continuous_parameters_view()
{
  return RealVector(Teuchos::View, keyRep->continuousHyperParams.values(),
		    keyRep->continuousHyperParams.length());
}


inline void ActiveKey::discrete_int_parameter(int di_param, size_t i)
{ keyRep->discreteIntHyperParams[i] = di_param; }


inline void ActiveKey::
discrete_int_parameters(const IntVector& di_params, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(di_params, keyRep->discreteIntHyperParams);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    keyRep->discreteIntHyperParams
      = IntVector(Teuchos::View, di_params.values(), di_params.length());
  else                           // default: assume existing Copy/View state
    keyRep->discreteIntHyperParams = di_params;
}


inline const IntVector& ActiveKey::discrete_int_parameters() const
{ return keyRep->discreteIntHyperParams; }


inline IntVector ActiveKey::discrete_int_parameters_view()
{
  return IntVector(Teuchos::View, keyRep->discreteIntHyperParams.values(),
		   keyRep->discreteIntHyperParams.length());
}


inline void ActiveKey::discrete_set_index(size_t ds_index, size_t i)
{ keyRep->discreteSetHyperParams[i] = ds_index; }


inline void ActiveKey::
discrete_set_indices(const SizetVector& ds_indices, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(ds_indices, keyRep->discreteSetHyperParams);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    keyRep->discreteSetHyperParams
      = SizetVector(Teuchos::View, ds_indices.values(), ds_indices.length());
  else                           // default: assume existing Copy/View state
    keyRep->discreteSetHyperParams = ds_indices;
}


inline const SizetVector& ActiveKey::discrete_set_indices() const
{ return keyRep->discreteSetHyperParams; }


inline SizetVector ActiveKey::discrete_set_indices_view()
{
  return SizetVector(Teuchos::View, keyRep->discreteSetHyperParams.values(),
		    keyRep->discreteSetHyperParams.length());
}


inline bool ActiveKey::is_null() const
{ return (keyRep) ? false : true; }


////////////////////////////////////////////////////////////////////////////////


/// Shared representation for composing a set of active keys plus a
/// group identifier.

/** For example, a model pairing for approximating a discrepancy would
    aggregate a high-fidelity plus a low-fidelity key. */

class AggregateKeyRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class AggregateKey;

public:

  ~AggregateKeyRep(); ///< destructor

private:

  //
  //- Heading: Constructors and destructor
  //

  AggregateKeyRep();  ///< default constructor
  AggregateKeyRep(unsigned short id,
		  const std::vector<ActiveKey>& keys);  ///< constructor

  //
  //- Heading: Member functions
  //

  //
  //- Heading: Private data members
  //

  // Currently there are two types of key aggregations:
  // 1. concatenation: a single UShortArray activeKey like 40302, identifying
  //    data group 4 for a discrepancy comprised of HF model 0 + resolution 3
  //    and LF model 0 + resolution 2
  // 2. SharedApproxData::approxDataKeys is a UShort2Darray that groups related
  //    approx data for Approximation::push/pop/finalize since the approx data
  //    includes raw HF, raw LF, and processed (distinct,recursive discrepancy)
  //    records.  In SharedApproxData::active_model_key(), our 40302 key is
  //    unrolled into keys for 3 datasets: HF = 403, LF = 402, discrep = 40302.
  //
  // Key question is how to redesign this...
  // > ActiveKeyArray --> ActiveKey2DArray ???   Hopefully not.
  // > Note: 3rd key only exists as combination of keys 1,2 and adds no new info
  //   (and it cannot be defined as a single state corresponding to an ActiveKey
  //   instance).  Therefore target the first case above with AggregateKey, and
  //   support the second case through enumeration ops on this aggregation.
  //   >> need a way to mark discrepancy data in the combined SurrogateData
  //      database.  Consider subsetting the DB into map<ActiveKey, ...> raw
  //      data and map<AggregateKey, ...> derived data --> could eliminate some
  //      current dataset filtering but disallows plug-and-play of raw/filtered.
  //      Detail: ActiveKey does not currently contain a group id.
  //      --> could push this down (exclusively) and generalize aggregation
  //          to allow cross group combination
  //      --> or go the other way and use map<AggregateKey, ...> exclusively
  //          where raw data uses an AggregateKey with only the set id and a
  //          single ActiveKey (**** seems preferable to start ***)
  // > Also consider adding an enum data type: RAW, DISTINCT_DISCREPANCY,
  //   SYNTHETIC_DISCREPANCY, ..., to retire some special case logic (USHRT_MAX
  //   for no model, no resolution)

  unsigned short dataSetId;

  std::vector<ActiveKey> activeKeys;

  // Don't need this since the idea is that the ActiveKey instances do not
  // reflect the complete state variables, only the subset identified as part
  // of solution control
  //ActiveKey sharedState;
};


inline AggregateKeyRep::AggregateKeyRep()
{ }


inline AggregateKeyRep::
AggregateKeyRep(unsigned short id, const std::vector<ActiveKey>& keys)
{ dataSetId = id; activeKeys = keys; }


inline AggregateKeyRep::~AggregateKeyRep()
{ }



/// Handle class for managing shared representations for AggregateKeyRep.

/** Provides user APIs for composing a key aggregation. */

class AggregateKey
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  AggregateKey();                            ///< default handle ctor (no body)
  AggregateKey(bool handle);                 ///< minimal handle + body ctor
  AggregateKey(unsigned short id,
	       const std::vector<ActiveKey>& keys); ///< full constructor
  AggregateKey(const AggregateKey& agg_key); ///< copy constructor
  ~AggregateKey();                           ///< destructor

  /// assignment operator
  AggregateKey& operator=(const AggregateKey& agg_key);

  //
  //- Heading: Member functions
  //


private:

  //
  //- Heading: Member functions
  //


  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  std::shared_ptr<AggregateKeyRep> aggKeyRep;
};


inline AggregateKey::AggregateKey()
{ } // aggKeyRep is null


// BMA NOTE: doesn't use make_shared<AggregateKeyRep>() due to private ctors
inline AggregateKey::AggregateKey(bool handle):
  aggKeyRep(new AggregateKeyRep())
{ }


inline AggregateKey::AggregateKey(unsigned short id,
				  const std::vector<ActiveKey>& keys):
  aggKeyRep(new AggregateKeyRep(id, keys))
{ }


inline AggregateKey::AggregateKey(const AggregateKey& agg_key):
  aggKeyRep(agg_key.aggKeyRep)  
{ }


inline AggregateKey::~AggregateKey()
{ }


inline AggregateKey& AggregateKey::operator=(const AggregateKey& agg_key)
{
  aggKeyRep = key.aggKeyRep;
  return *this;
}

} // namespace Pecos

#endif
