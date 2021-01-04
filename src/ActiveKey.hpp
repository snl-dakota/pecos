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

/// Shared representation for ActiveKeyData class (body within
/// handle-body idiom).

/** Manages a set of model indices and a set of continuous/discrete
    hyper-parameter (resolution) controls. */

class ActiveKeyDataRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class ActiveKeyData;

public:

  /// destructor
  ~ActiveKeyDataRep();

private:

  //
  //- Heading: Private member functions
  //

  /// default constructor (default empty key)
  ActiveKeyDataRep();
  /// partial constructor (legacy use case: model indices only)
  ActiveKeyDataRep(const UShortArray& indices);
  /// full constructor
  ActiveKeyDataRep(const UShortArray& indices, const RealVector&   c_params,
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


inline ActiveKeyDataRep::~ActiveKeyDataRep()
{ }


inline ActiveKeyDataRep::ActiveKeyDataRep()
{ }


inline ActiveKeyDataRep::ActiveKeyDataRep(const UShortArray& indices)
{ modelIndices = indices; }


inline ActiveKeyDataRep::
ActiveKeyDataRep(const UShortArray& indices, const RealVector&   c_params,
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

class ActiveKeyData
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  /// default constructor
  ActiveKeyData();
  /// partial constructor (legacy use case: model indices only)
  ActiveKeyData(const UShortArray& indices, short mode = DEFAULT_COPY);
  /// full constructor
  ActiveKeyData(const UShortArray& indices, const RealVector&   c_params,
		const IntVector& di_params, const SizetVector& ds_params,
		short mode = DEFAULT_COPY);

  /// copy constructor
  ActiveKeyData(const ActiveKeyData& key);
  /// destructor
  ~ActiveKeyData();

  /// assignment operator
  ActiveKeyData& operator=(const ActiveKeyData& key);
  // equality operator
  //bool operator==(const ActiveKeyData& key) const;

  //
  //- Heading: member functions
  //

  /// return number of continuous parameters
  size_t chp() const;
  /// return number of discrete integer parameters
  size_t dihp() const;
  /// return number of discrete real parameters
  size_t dshp() const;

  /// return deep copy of ActiveKeyData instance
  ActiveKeyData copy() const;

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

  /// function to check keyDataRep (does this handle contain a body)
  bool is_null() const;

private:

  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  std::shared_ptr<ActiveKeyDataRep> keyDataRep;
};


inline ActiveKeyData::ActiveKeyData()
{ } // keyDataRep is null


// BMA NOTE: The following don't use make_shared<ActiveKeyDataRep>()
// due to private ctors

inline ActiveKeyData::ActiveKeyData(const UShortArray& indices):
  keyDataRep(new ActiveKeyDataRep(indices))
{ }


inline ActiveKeyData::
ActiveKeyData(const UShortArray& indices, const RealVector&   c_params,
	  const IntVector& di_params, const SizetVector& ds_params, short mode):
  keyDataRep(new ActiveKeyDataRep(indices, c_params, di_params, ds_params,mode))
{ }


inline ActiveKeyData::ActiveKeyData(const ActiveKeyData& key):
  keyDataRep(key.keyDataRep)  
{ }


inline ActiveKeyData::~ActiveKeyData()
{ }


inline ActiveKeyData& ActiveKeyData::operator=(const ActiveKeyData& key)
{
  keyDataRep = key.keyDataRep;
  return *this;
}


//inline bool ActiveKeyData::operator==(const ActiveKeyData& key) const
//{
//  return ( keyDataRep->modelIndices    == key.keyDataRep->modelIndices &&
//    keyDataRep->continuousHyperParams  ==
//      key.keyDataRep->continuousHyperParams  &&
//    keyDataRep->discreteIntHyperParams ==
//      key.keyDataRep->discreteIntHyperParams &&
//    keyDataRep->discreteSetHyperParams ==
//      key.keyDataRep->discreteSetHyperParams );
//}


inline size_t ActiveKeyData::chp() const
{ return keyDataRep->continuousHyperParams.length(); }


inline size_t ActiveKeyData::dihp() const
{ return keyDataRep->discreteIntHyperParams.length(); }


inline size_t ActiveKeyData::dshp() const
{ return keyDataRep->discreteSetHyperParams.length(); }


/// deep copy of ActiveKeyData instance
inline ActiveKeyData ActiveKeyData::copy() const
{
  ActiveKeyData data(keyDataRep->modelIndices,
		     keyDataRep->continuousHyperParams,
		     keyDataRep->discreteIntHyperParams,
		     keyDataRep->discreteSetHyperParams, DEEP_COPY);
  return data;
}


inline void ActiveKeyData::continuous_parameter(Real c_param, size_t i)
{ keyDataRep->continuousHyperParams[i] = c_param; }


inline void ActiveKeyData::
continuous_parameters(const RealVector& c_params, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(c_params, keyDataRep->continuousHyperParams);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    keyDataRep->continuousHyperParams
      = RealVector(Teuchos::View, c_params.values(), c_params.length());
  else                           // default: assume existing Copy/View state
    keyDataRep->continuousHyperParams = c_params;
}


inline const RealVector& ActiveKeyData::continuous_parameters() const
{ return keyDataRep->continuousHyperParams; }


inline RealVector ActiveKeyData::continuous_parameters_view()
{
  return RealVector(Teuchos::View, keyDataRep->continuousHyperParams.values(),
		    keyDataRep->continuousHyperParams.length());
}


inline void ActiveKeyData::discrete_int_parameter(int di_param, size_t i)
{ keyDataRep->discreteIntHyperParams[i] = di_param; }


inline void ActiveKeyData::
discrete_int_parameters(const IntVector& di_params, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(di_params, keyDataRep->discreteIntHyperParams);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    keyDataRep->discreteIntHyperParams
      = IntVector(Teuchos::View, di_params.values(), di_params.length());
  else                           // default: assume existing Copy/View state
    keyDataRep->discreteIntHyperParams = di_params;
}


inline const IntVector& ActiveKeyData::discrete_int_parameters() const
{ return keyDataRep->discreteIntHyperParams; }


inline IntVector ActiveKeyData::discrete_int_parameters_view()
{
  return IntVector(Teuchos::View, keyDataRep->discreteIntHyperParams.values(),
		   keyDataRep->discreteIntHyperParams.length());
}


inline void ActiveKeyData::discrete_set_index(size_t ds_index, size_t i)
{ keyDataRep->discreteSetHyperParams[i] = ds_index; }


inline void ActiveKeyData::
discrete_set_indices(const SizetVector& ds_indices, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(ds_indices, keyDataRep->discreteSetHyperParams);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    keyDataRep->discreteSetHyperParams
      = SizetVector(Teuchos::View, ds_indices.values(), ds_indices.length());
  else                           // default: assume existing Copy/View state
    keyDataRep->discreteSetHyperParams = ds_indices;
}


inline const SizetVector& ActiveKeyData::discrete_set_indices() const
{ return keyDataRep->discreteSetHyperParams; }


inline SizetVector ActiveKeyData::discrete_set_indices_view()
{
  return SizetVector(Teuchos::View, keyDataRep->discreteSetHyperParams.values(),
		    keyDataRep->discreteSetHyperParams.length());
}


inline bool ActiveKeyData::is_null() const
{ return (keyDataRep) ? false : true; }


////////////////////////////////////////////////////////////////////////////////


/// Shared representation for composing a set of active key data
/// instances plus a group identifier.

/** For example, a model pairing for approximating a discrepancy would
    aggregate a high-fidelity plus a low-fidelity key. */

class ActiveKeyRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class ActiveKey;

public:

  ~ActiveKeyRep(); ///< destructor

private:

  //
  //- Heading: Constructors and destructor
  //

  ActiveKeyRep(); ///< default constructor
  ActiveKeyRep(unsigned short id,
	       const std::vector<ActiveKeyData>& data); ///< constructor

  //
  //- Heading: Member functions
  //

  //
  //- Heading: Private data members
  //

  // Currently Dakota,Pecos use two types of key aggregations:
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
  // > ActiveKeyDataArray --> ActiveKeyData2DArray ???   Hopefully not.
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

  std::vector<ActiveKeyData> activeKeyDataArray;

  // Don't need this since the idea is that the ActiveKeyData instances do not
  // reflect the complete state variables, only the subset identified as part
  // of solution control
  //ActiveKeyData sharedState;
};


inline ActiveKeyRep::ActiveKeyRep()
{ }


inline ActiveKeyRep::
ActiveKeyRep(unsigned short id, const std::vector<ActiveKeyData>& data)
{ dataSetId = id; activeKeyDataArray = data; }


inline ActiveKeyRep::~ActiveKeyRep()
{ }



/// Handle class for managing shared representations for ActiveKeyRep.

/** Provides user APIs for composing an active key that aggregates a
    group id and one or more key data instances. */

class ActiveKey
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  ActiveKey();                         ///< default handle ctor (no body)
  ActiveKey(bool handle);              ///< minimal handle + body ctor
  ActiveKey(unsigned short id,
	    const std::vector<ActiveKeyData>& data); ///< full constructor
  ActiveKey(const ActiveKey& key); ///< copy constructor
  ~ActiveKey();                        ///< destructor

  /// assignment operator
  ActiveKey& operator=(const ActiveKey& key);

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
  std::shared_ptr<ActiveKeyRep> keyRep;
};


inline ActiveKey::ActiveKey()
{ } // keyRep is null


// BMA NOTE: doesn't use make_shared<ActiveKeyRep>() due to private ctors
inline ActiveKey::ActiveKey(bool handle):
  keyRep(new ActiveKeyRep())
{ }


inline ActiveKey::
ActiveKey(unsigned short id, const std::vector<ActiveKeyData>& data):
  keyRep(new ActiveKeyRep(id, data))
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

} // namespace Pecos

#endif
