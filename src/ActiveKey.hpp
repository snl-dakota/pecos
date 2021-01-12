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
	= SizetVector(Teuchos::View, ds_params.values(), ds_params.length());
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
  ActiveKeyData(const UShortArray& indices);
  /// full constructor
  ActiveKeyData(const UShortArray& indices, const RealVector&   c_params,
		const IntVector& di_params, const SizetVector& ds_params,
		short mode = DEFAULT_COPY);

  /// copy constructor
  ActiveKeyData(const ActiveKeyData& key_data);
  /// destructor
  ~ActiveKeyData();

  /// assignment operator
  ActiveKeyData& operator=(const ActiveKeyData& key_data);
  // equality operator
  bool operator==(const ActiveKeyData& key_data) const;
  // less-than operator
  bool operator<(const ActiveKeyData& key_data) const;

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

  /// set i^{th} entry within modelIndices
  void model_index(unsigned short mi, size_t i);
  /// get i^{th} entry from modelIndices
  unsigned short model_index(size_t i) const;
  /// set modelIndices
  void model_indices(const UShortArray& indices);
  /// get modelIndices
  const UShortArray& model_indices() const;

  /// set i^{th} entry within continuousHyperParams
  void continuous_parameter(Real c_param, size_t i);
  /// get i^{th} entry from continuousHyperParams
  Real continuous_parameter(size_t i) const;
  /// set continuousHyperParams
  void continuous_parameters(const RealVector& c_params,
			     short mode = DEFAULT_COPY);
  /// get continuousHyperParams
  const RealVector& continuous_parameters() const;
  /// get view of continuousHyperParams for updating in place
  RealVector continuous_parameters_view();

  /// set i^{th} entry within discreteIntHyperParams
  void discrete_int_parameter(int di_param, size_t i);
  /// get i^{th} entry from discreteIntHyperParams
  int discrete_int_parameter(size_t i) const;
  /// set discreteIntHyperParams
  void discrete_int_parameters(const IntVector& di_params,
			       short mode = DEFAULT_COPY);
  /// get discreteIntHyperParams
  const IntVector& discrete_int_parameters() const;
  /// get view of discreteIntHyperParams for updating in place
  IntVector discrete_int_parameters_view();

  /// set i^{th} entry within discreteSetHyperParams
  void discrete_set_index(size_t ds_index, size_t i);
  /// get i^{th} entry from discreteSetHyperParams
  size_t discrete_set_index(size_t i) const;
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
	      const IntVector& di_params, const SizetVector& ds_params,
	      short mode):
  keyDataRep(new ActiveKeyDataRep(indices,c_params,di_params,ds_params,mode))
{ }


inline ActiveKeyData::ActiveKeyData(const ActiveKeyData& key_data):
  keyDataRep(key_data.keyDataRep)  
{ }


inline ActiveKeyData::~ActiveKeyData()
{ }


inline ActiveKeyData& ActiveKeyData::operator=(const ActiveKeyData& key_data)
{
  keyDataRep = key_data.keyDataRep;
  return *this;
}


inline bool ActiveKeyData::operator==(const ActiveKeyData& key_data) const
{
  return ( keyDataRep->modelIndices == key_data.keyDataRep->modelIndices &&
    keyDataRep->continuousHyperParams  ==
      key_data.keyDataRep->continuousHyperParams  &&
    keyDataRep->discreteIntHyperParams ==
      key_data.keyDataRep->discreteIntHyperParams &&
    keyDataRep->discreteSetHyperParams ==
      key_data.keyDataRep->discreteSetHyperParams );
}


inline bool ActiveKeyData::operator<(const ActiveKeyData& key_data) const
{
  std::shared_ptr<ActiveKeyDataRep> kdr = key_data.keyDataRep;
  if      (keyDataRep->modelIndices < kdr->modelIndices)
    return true;
  else if (kdr->modelIndices < keyDataRep->modelIndices)
    return false;
  // else equal -> continue to next array

  if      (keyDataRep->continuousHyperParams < kdr->continuousHyperParams)
    return true;
  else if (kdr->continuousHyperParams < keyDataRep->continuousHyperParams)
    return false;
  // else equal -> continue to next array

  if      (keyDataRep->discreteIntHyperParams < kdr->discreteIntHyperParams)
    return true;
  else if (kdr->discreteIntHyperParams < keyDataRep->discreteIntHyperParams)
    return false;
  // else equal -> continue to final array

  if      (keyDataRep->discreteSetHyperParams < kdr->discreteSetHyperParams)
    return true;

  return false;
}


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


inline void ActiveKeyData::model_index(unsigned short mi, size_t i)
{
  UShortArray& indices = keyDataRep->modelIndices;
  size_t len = indices.size();
  if      (i <  len)  indices[i] = mi;
  else if (i == len)  indices.push_back(mi); // allow appends
  else {
    PCerr << "Error: index " << i << " out of bounds in ActiveKeyData::"
	  << "model_indices(unsigned short)" << std::endl;
    abort_handler(-1);
  }
}


inline unsigned short ActiveKeyData::model_index(size_t i) const
{
  const UShortArray& indices = keyDataRep->modelIndices;
  if (i < indices.size()) return indices[i];
  else {
    PCerr << "Error: index " << i << " out of bounds in ActiveKeyData::"
	  << "model_indices()" << std::endl;
    abort_handler(-1);  return 0;
  }
}


inline void ActiveKeyData::model_indices(const UShortArray& indices)
{ keyDataRep->modelIndices = indices; }


inline const UShortArray& ActiveKeyData::model_indices() const
{ return keyDataRep->modelIndices; }


inline void ActiveKeyData::continuous_parameter(Real c_param, size_t i)
{
  RealVector& chp = keyDataRep->continuousHyperParams;
  size_t len = chp.length();
  if (i == len) chp.resize(len+1); // allow appends (resize since no push_back)
  if (i <= len) chp[i] = c_param;
  else {
    PCerr << "Error: index " << i << " out of bounds in ActiveKeyData::"
	  << "continuous_parameter(Real)" << std::endl;
    abort_handler(-1);
  }
}


inline Real ActiveKeyData::continuous_parameter(size_t i) const
{
  const RealVector& chp = keyDataRep->continuousHyperParams;
  size_t len = chp.length();
  if (i < len) return chp[i];
  else {
    PCerr << "Error: index " << i << " out of bounds in ActiveKeyData::"
	  << "continuous_parameter()" << std::endl;
    abort_handler(-1);  return 0.;
  }
}


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
{
  IntVector& dihp = keyDataRep->discreteIntHyperParams;
  size_t len = dihp.length();
  if (i == len) dihp.resize(len+1); // allow appends (resize since no push_back)
  if (i <= len) dihp[i] = di_param;
  else {
    PCerr << "Error: index " << i << " out of bounds in ActiveKeyData::"
	  << "discrete_int_parameter(int)" << std::endl;
    abort_handler(-1);
  }
}


inline int ActiveKeyData::discrete_int_parameter(size_t i) const
{
  const IntVector& dihp = keyDataRep->discreteIntHyperParams;
  if (i < dihp.length()) return dihp[i];
  else {
    PCerr << "Error: index " << i << " out of bounds in ActiveKeyData::"
	  << "discrete_int_parameter()" << std::endl;
    abort_handler(-1);  return 0;
  }
}


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
{
  SizetVector& dshp = keyDataRep->discreteSetHyperParams;
  size_t len = dshp.length();
  if (i == len) dshp.resize(len+1); // allow appends (resize since no push_back)
  if (i <= len) dshp[i] = ds_index;
  else {
    PCerr << "Error: index " << i << " out of bounds in ActiveKeyData::"
	  << "discrete_set_index(size_t)" << std::endl;
    abort_handler(-1);
  }
}


inline size_t ActiveKeyData::discrete_set_index(size_t i) const
{
  const SizetVector& dshp = keyDataRep->discreteSetHyperParams;
  if (i < dshp.length()) return dshp[i];
  else {
    PCerr << "Error: index " << i << " out of bounds in ActiveKeyData::"
	  << "discrete_set_index()" << std::endl;
    abort_handler(-1);  return 0;
  }
}


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
  ActiveKeyRep(unsigned short id); ///< minimal constructor
  ActiveKeyRep(unsigned short id, const std::vector<ActiveKeyData>& data,
	       short mode); ///< constructor
  ActiveKeyRep(unsigned short id, const ActiveKeyData& data,
	       short mode); ///< constructor

  //
  //- Heading: Member functions
  //

  /// assign dataSetId and activeKeyDataArray using shallow/deep/default copy
  void assign(const std::vector<ActiveKeyData>& data_vec, short mode);
  /// assign dataSetId and activeKeyDataArray using shallow/deep/default copy
  void assign(const ActiveKeyData& data, short mode);

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
  // > Note: 3rd key only exists as combination of keys 1,2 and adds no new
  //   info (and it cannot be defined as a single state corresponding to an
  //   ActiveKeyData instance).  Therefore target 1st case above with ActiveKey,
  //   and support 2nd case through enumeration ops on this aggregation.
  //   >> need a way to mark discrepancy data in the combined SurrogateData
  //      database.  Consider subsetting the DB into map<ActiveKeyData, ...> raw
  //      data and map<ActiveKey, ...> derived data --> could eliminate some
  //      current dataset filtering but disallows plug-and-play of raw/filtered.
  //      Detail: ActiveKeyData does not currently contain a group id.
  //      --> could push this down (exclusively) and generalize aggregation
  //          to allow cross group combination
  //      --> or go the other way and use map<ActiveKey, ...> exclusively
  //          where raw data uses an ActiveKey with only the set id and a
  //          single ActiveKeyData (*** seems preferable to start ***)
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


inline ActiveKeyRep::ActiveKeyRep(unsigned short id):
  dataSetId(id)
{ }


inline ActiveKeyRep::
ActiveKeyRep(unsigned short id, const std::vector<ActiveKeyData>& key_data_vec,
	     short mode):
  dataSetId(id)
{ assign(key_data_vec, mode); }


inline ActiveKeyRep::
ActiveKeyRep(unsigned short id, const ActiveKeyData& key_data, short mode):
  dataSetId(id)
{ assign(key_data, mode); }


inline ActiveKeyRep::~ActiveKeyRep()
{ }


inline void ActiveKeyRep::
assign(const std::vector<ActiveKeyData>& key_data_vec, short mode)
{
  if (mode == DEEP_COPY) { // enforce deep copy for each activeKeyData
    size_t i, num_data = key_data_vec.size();
    activeKeyDataArray.resize(num_data);
    for (i=0; i<num_data; ++i)
      activeKeyDataArray[i] = key_data_vec[i].copy();
  }
  else // each activeKeyData shares rep
    activeKeyDataArray = key_data_vec;
}


inline void ActiveKeyRep::assign(const ActiveKeyData& key_data, short mode)
{
  activeKeyDataArray.clear();
  if (mode == DEEP_COPY) activeKeyDataArray.push_back(key_data.copy());
  else                   activeKeyDataArray.push_back(key_data);
}


/// Handle class for managing shared representations for ActiveKeyRep.

/** Provides user APIs for composing an active key that aggregates a
    group id and one or more key data instances. */

class ActiveKey
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  ActiveKey();                     ///< default handle ctor (no body)
  ActiveKey(unsigned short id);    ///< minimal handle + body ctor
  ActiveKey(unsigned short id, const std::vector<ActiveKeyData>& key_data_vec,
	    short mode);           ///< full constructor for multiple data
  ActiveKey(unsigned short id, const ActiveKeyData& key_data,
	    short mode);           ///< full constructor for 1 data
  ActiveKey(const ActiveKey& key); ///< copy constructor
  ~ActiveKey();                    ///< destructor

  /// assignment operator
  ActiveKey& operator=(const ActiveKey& key);
  // equality operator
  bool operator==(const ActiveKey& key) const;
  // inequality operator
  bool operator!=(const ActiveKey& key) const;
  // less-than operator
  bool operator<(const ActiveKey& key) const;

  //
  //- Heading: Member functions
  //

  /// get dataSetId
  unsigned short id() const;
  /// set dataSetId
  void id(unsigned short set_id);

  /// get activeKeyDataArray
  const std::vector<ActiveKeyData>& data() const;
  /// set activeKeyDataArray
  void data(const std::vector<ActiveKeyData>& key_data_vec, short mode);
  /// set activeKeyDataArray
  void data(const ActiveKeyData& key_data, short mode);

  /// assign data to ActiveKey
  void assign(unsigned short id, const std::vector<ActiveKeyData>& key_data_vec,
	      short mode);
  /// assign data to ActiveKey
  void assign(unsigned short id, const ActiveKeyData& key_data, short mode);

  /// assign data to ActiveKey
  void append(const std::vector<ActiveKeyData>& key_data_vec, short mode);
  /// assign data to ActiveKey
  void append(const ActiveKeyData& key_data, short mode);
  /// clear activeKeyDataArray
  void clear_data();

  /// return deep copy of ActiveKey instance
  ActiveKey copy() const;

  /// function to check keyRep (does this handle contain a body)
  bool is_null() const;

  /*
  /// define a model key including data group, model form, and resolution
  /// level indices
  static void form_key(unsigned short group, unsigned short form,
		       unsigned short lev, UShortArray& key);
  /// define an aggregate model key including data group and two sets of
  /// model form and resolution level indices
  static void form_key(unsigned short group, unsigned short form1,
		       unsigned short lev1,  unsigned short form2,
		       unsigned short lev2,  UShortArray& key);
  /// decrement an incoming model key to correspond to the next lower
  /// resolution or fidelity within a model sequence
  static bool decrement_key(UShortArray& key, size_t index);
  */
  /// aggregate two model keys to indicate a data combination
  /// (e.g., a discrepancy)
  static void aggregate_keys(const ActiveKey& key1, const ActiveKey& key2,
			     ActiveKey& aggregate_key);
  /// aggregate first_key and remaining_keys to indicate a data combination
  /// (e.g., a model ensemble)
  static void aggregate_keys(const ActiveKey& key1,
			     const std::vector<ActiveKey>& other_keys,
			     ActiveKey& aggregate_key);
  /*
  /// extract two constituent keys from an aggregated key
  static void extract_keys(const UShortArray& aggregate_key, UShortArray& key1,
			   UShortArray& key2);
  /// extract one or more constituent keys from an aggregated key
  static void extract_keys(const UShortArray& aggregate_key,
			   UShortArray&   first_key,
			   UShort2DArray& remaining_keys);
  /// extract a particular constituent key from an aggregated key
  static void extract_key(const UShortArray& aggregate_key, UShortArray& key,
			  size_t key_index);
  */

  /// test whether key is an aggregated key (for discrepancy or surplus)
  static bool aggregated_key(const ActiveKey& key);
  /*
  /// test whether key is used for synthetic data (e.g., from an interpolant
  /// using a previous level's raw data, preceding a surplus estimation)
  static bool synthetic_key(const UShortArray& key);
  /// test whether key is a raw data key (response data from a single model)
  static bool raw_data_key(const UShortArray& key);
  */

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
inline ActiveKey::ActiveKey(unsigned short id):
  keyRep(new ActiveKeyRep(id))
{ }


inline ActiveKey::
ActiveKey(unsigned short id, const std::vector<ActiveKeyData>& key_data_vec,
	  short mode):
  keyRep(new ActiveKeyRep(id, key_data_vec, mode))
{ }


inline ActiveKey::
ActiveKey(unsigned short id, const ActiveKeyData& key_data, short mode):
  keyRep(new ActiveKeyRep(id, key_data, mode))
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


inline bool ActiveKey::operator==(const ActiveKey& key) const
{
  return ( keyRep->dataSetId          == key.keyRep->dataSetId &&
	   keyRep->activeKeyDataArray == key.keyRep->activeKeyDataArray );
}


inline bool ActiveKey::operator!=(const ActiveKey& key) const
{
  return ( keyRep->dataSetId          != key.keyRep->dataSetId ||
	   keyRep->activeKeyDataArray != key.keyRep->activeKeyDataArray );
}


inline bool ActiveKey::operator<(const ActiveKey& key) const
{
  std::shared_ptr<ActiveKeyRep> kr = key.keyRep;
  if (keyRep->dataSetId < kr->dataSetId)
    return true;
  else if (kr->dataSetId < keyRep->dataSetId)
    return false;
  // else equal -> continue to next array

  if (keyRep->activeKeyDataArray < kr->activeKeyDataArray)
    return true;

  return false;
}


inline unsigned short ActiveKey::id() const
{ return keyRep->dataSetId; }


inline void ActiveKey::id(unsigned short set_id)
{ keyRep->dataSetId = set_id; }


inline const std::vector<ActiveKeyData>& ActiveKey::data() const
{ return keyRep->activeKeyDataArray; }


inline void ActiveKey::
data(const std::vector<ActiveKeyData>& key_data_vec, short mode)
{ keyRep->assign(key_data_vec, mode); }


inline void ActiveKey::data(const ActiveKeyData& key_data, short mode)
{ keyRep->assign(key_data, mode); }


inline void ActiveKey::
assign(unsigned short set_id, const std::vector<ActiveKeyData>& key_data_vec,
       short mode)
{
  if  (keyRep) { id(set_id); data(key_data_vec, mode); }
  else keyRep = std::make_shared<ActiveKeyRep>(set_id, key_data_vec, mode);
}


inline void ActiveKey::
assign(unsigned short set_id, const ActiveKeyData& key_data, short mode)
{
  if  (keyRep) { id(set_id); data(key_data, mode); }
  else keyRep = std::make_shared<ActiveKeyRep>(set_id, key_data, mode);
}


inline void ActiveKey::
append(const std::vector<ActiveKeyData>& key_data_vec, short mode)
{
  std::vector<ActiveKeyData>& act_key_data = keyRep->activeKeyDataArray;
  if (mode == DEEP_COPY) {
    size_t i, len = key_data_vec.size();
    for (i=0; i<len; ++i)
      act_key_data.push_back(key_data_vec[i].copy());
  }
  else
    act_key_data.insert(act_key_data.end(), key_data_vec.begin(),
			key_data_vec.end());
}


inline void ActiveKey::append(const ActiveKeyData& key_data, short mode)
{
  if (mode == DEEP_COPY) keyRep->activeKeyDataArray.push_back(key_data.copy());
  else                   keyRep->activeKeyDataArray.push_back(key_data);
}


inline void ActiveKey::clear_data()
{ keyRep->activeKeyDataArray.clear(); }


/// deep copy of ActiveKey instance
inline ActiveKey ActiveKey::copy() const
{
  ActiveKey key(keyRep->dataSetId, keyRep->activeKeyDataArray, DEEP_COPY);
  return key;
}


inline bool ActiveKey::is_null() const
{ return (keyRep) ? false : true; }


////////////////////////////////////////////////////////////////////////////////


/*
inline void ActiveKey::
form_key(unsigned short group, unsigned short form, unsigned short lev,
	 UShortArray& key)
{ key.resize(3);  key[0] = group;  key[1] = form;  key[2] = lev; }


inline void ActiveKey::
form_key(unsigned short group, unsigned short form1, unsigned short lev1,
	 unsigned short form2, unsigned short lev2,  UShortArray& key)
{
  key.resize(5);
  key[0] = group;                  // data group
  key[1] = form1;  key[2] = lev1;  // HF model form, soln level
  key[3] = form2;  key[4] = lev2;  // LF model form, soln level
}


inline bool ActiveKey::decrement_key(UShortArray& key, size_t index)
{
  // decrement the active index, if present, to create a key within the same
  // group id but with the next lower resolution in the sequence

  //if (key.size() != 3) { // don't allow aggregated keys
  //  PCerr << "Error: wrong size for {group,form,lev} format in Discrepancy"
  //	    << "Calculator::decrement_key()" << std::endl;
  //  abort_handler(-1);    
  //}
  if (index >= key.size())
    return false;

  unsigned short &key_i = key[index];
  if (key_i && key_i != USHRT_MAX)
    { --key_i; return true; }
  else// decrement undefined (e.g., already at coarsest resolution / lowest fid)
    return false;
}


inline bool ActiveKey::decrement_key(UShortArray& key)
{
  // decrement the active index, if present, to create a key within the same
  // group id but with the next lower resolution in the sequence

  if (key.size() != 3) { // don't allow aggregated keys
    PCerr << "Error: wrong size for {group,form,lev} format in Discrepancy"
	  << "Calculator::decrement_key()" << std::endl;
    abort_handler(-1);    
  }

  // Logic is fragile in that it fails if a fixed model index (index that is
  // not part of the sequence) is assigned a value other than 0 or USHRT_MAX
  // > precedence given to lev for this reason, as form is more likely to have
  //   a non-zero/inf fixed value (see NonDExpansion::configure_sequence())
  // > more robust approach would be to pass in a multilev boolean
  unsigned short &form = key[1], &lev = key[2];
  if      (lev  && lev  != USHRT_MAX)
    { --lev;  return true; }
  else if (form && form != USHRT_MAX)
    { --form; return true; }
  //else no op (already at coarsest resolution / lowest fidelity)
  return false;

  // Old logic for {form} | {form,lev} format was simply --key.back();
}


inline bool ActiveKey::synthetic_key(const UShortArray& key)
{ return (key.size() == 3 && key[1] == USHRT_MAX); } // model form is undefined


inline bool ActiveKey::raw_data_key(const UShortArray& key)
{ return (key.size() == 3 && key[1] != USHRT_MAX); } // model form is defined
*/


inline bool ActiveKey::aggregated_key(const ActiveKey& key)
{ return (key.data().size() > 1); }


inline void ActiveKey::
aggregate_keys(const ActiveKey& key1, const ActiveKey& key2,
	       ActiveKey& aggregate_key)
{
  // extract and verify consistency in id number
  unsigned short id = USHRT_MAX;
  bool empty1 = key1.is_null(), empty2 = key2.is_null();
  if (!empty1 && !empty2) {
    id = key1.id();
    if (id != key2.id()) {
      PCerr << "Error: mismatch in group ids in ActiveKey::aggregate_keys()"
	    << std::endl;
      abort_handler(-1);
    }
  }
  else if (!empty1) id = key1.id();
  else if (!empty2) id = key2.id();
  else {
    PCerr << "Error: neither key defined in ActiveKey::aggregate_keys"
	  << "(key1, key2)" << std::endl;
    abort_handler(-1);    
  }

  aggregate_key.id(id);
  aggregate_key.clear_data();
  if (!empty1) aggregate_key.append(key1.data(), DEEP_COPY);
  if (!empty2) aggregate_key.append(key2.data(), DEEP_COPY);
}


inline void ActiveKey::
aggregate_keys(const ActiveKey& key1, const std::vector<ActiveKey>& other_keys,
	       ActiveKey& aggregate_key)
{
  // extract and verify consistency in group number
  unsigned short id = USHRT_MAX;
  bool empty1 = key1.is_null();
  size_t i, num_other_k = other_keys.size();
  if (!empty1) {
    id = key1.id();
    for (i=0; i<num_other_k; ++i)
      if (id != other_keys[i].id()) {
	PCerr << "Error: mismatch in group ids in ActiveKey::aggregate_keys()"
	      << std::endl;
	abort_handler(-1);
      }
  }
  else if (num_other_k) {
    id = other_keys[0].front();
    for (i=1; i<num_other_k; ++i)
      if (id != other_keys[i].id()) {
	PCerr << "Error: mismatch in group ids in ActiveKey::aggregate_keys()"
	      << std::endl;
	abort_handler(-1);
      }
  }
  else {
    PCerr << "Error: neither key set defined in ActiveKey::aggregate_keys"
	  << "(key1, other_keys)" << std::endl;
    abort_handler(-1);
  }

  // form aggregate of group + HF form/lev + LF form/lev
  aggregate_key.id(id);
  aggregate_key.clear_data();
  if (!empty1) aggregate_key.append(key1.data(), DEEP_COPY);
  for (i=0; i<num_other_k; ++i)
    aggregate_key.append(other_keys[i].data(), DEEP_COPY);
}


/*
inline void ActiveKey::
extract_keys(const UShortArray& aggregate_key, UShortArray& key1,
	     UShortArray& key2)
{
  if (aggregate_key.is_null())
    { key1.clear(); key2.clear(); return; }

  // extract one or two aggregated keys
  unsigned short group = aggregate_key.front(),
    len = aggregate_key.size() - 1, num_keys = len / 2;
  switch (num_keys) {
  case 1:
    key1 = aggregate_key;   key2.clear();  break;
  case 2: { // normal case
    UShortArray::const_iterator start1 = aggregate_key.begin() + 1,
      end1 = start1 + 2, end2 = end1 + 2;
    key1.assign(1, group);  key1.insert(key1.end(), start1, end1);
    key2.assign(1, group);  key2.insert(key2.end(),   end1, end2);
    break;
  }
  default:
    PCerr << "Error: bad aggregate key size in ActiveKey::"
	  << "extract_keys()" << std::endl;
    abort_handler(-1);
    break;
  }
}


inline void ActiveKey::
extract_keys(const UShortArray& aggregate_key, UShortArray& key1,
	     UShort2DArray& other_keys)
{
  if (aggregate_key.is_null())
    { key1.clear(); other_keys.clear(); return; }

  // extract one or more aggregated keys
  unsigned short group = aggregate_key.front(),
    len = aggregate_key.size() - 1, num_keys = len / 2;
  switch (num_keys) {
  case 0: case 1:
    key1 = aggregate_key;  other_keys.clear();  break;
  default: {
    UShortArray::const_iterator start = aggregate_key.begin() + 1,
      end = start + 2;
    key1.assign(1, group);  key1.insert(key1.end(), start, end);
    size_t k, num_rem_keys = num_keys - 1;
    other_keys.resize(num_rem_keys);
    for (k=0; k<num_rem_keys; ++k) {
      start += 2;  end +=2;
      UShortArray& key_k = other_keys[k];
      key_k.assign(1, group);  key_k.insert(key_k.end(), start, end);
    }
    break;
  }
  }
}


inline void ActiveKey::
extract_key(const UShortArray& aggregate_key, UShortArray& key,
	    size_t key_index)
{
  if (aggregate_key.is_null())
    { key.clear(); return; }

  // extract one or two aggregated keys
  unsigned short group = aggregate_key.front(),
    len = aggregate_key.size() - 1, num_keys = len / 2;
  if (key_index >= num_keys)
    { key.clear(); return; }

  key.assign(1, group);
  UShortArray::const_iterator start = aggregate_key.begin() + 1 + 2*key_index,
    end = start + 2;
  key.insert(key.end(), start, end);
}
*/

} // namespace Pecos

#endif
