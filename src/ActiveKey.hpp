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

  /// default constructor (default empty key)
  ActiveKeyDataRep();
  /// partial constructor (legacy use case: model indices only)
  ActiveKeyDataRep(const UShortArray& indices);
  /// full constructor
  ActiveKeyDataRep(const UShortArray& indices, const RealVector&   c_params,
		   const IntVector& di_params, const SizetVector& ds_params,
		   short copy_mode);
  /// destructor
  ~ActiveKeyDataRep();

private:

  //
  //- Heading: Private member functions
  //

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


inline ActiveKeyDataRep::ActiveKeyDataRep()
{ }


inline ActiveKeyDataRep::ActiveKeyDataRep(const UShortArray& indices)
{ modelIndices = indices; }


inline ActiveKeyDataRep::
ActiveKeyDataRep(const UShortArray& indices, const RealVector&   c_params,
		 const IntVector& di_params, const SizetVector& ds_params,
		 short copy_mode)
{
  modelIndices = indices;

  // Note: provided a way to query DataAccess mode for c_params, could make
  // greater use of operator= for {DEEP,SHALLOW}_COPY modes
  if (copy_mode == DEEP_COPY) {         // enforce deep vector copy
    if (!c_params.empty())  copy_data( c_params, continuousHyperParams);
    if (!di_params.empty()) copy_data(di_params, discreteIntHyperParams);
    if (!ds_params.empty()) copy_data(ds_params, discreteSetHyperParams);
  }
  else if (copy_mode == SHALLOW_COPY) { // enforce shallow vector copy
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


inline ActiveKeyDataRep::~ActiveKeyDataRep()
{ }


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

  /// default constructor (no keyDataRep instantiation)
  ActiveKeyData();
  /// minimal constructor (keyDataRep instantiated if handle is true)
  ActiveKeyData(bool handle);
  /// partial constructor (legacy use case: model indices only)
  ActiveKeyData(const UShortArray& indices);
  /// full constructor
  ActiveKeyData(const UShortArray& indices, const RealVector&   c_params,
		const IntVector& di_params, const SizetVector& ds_params,
		short copy_mode = DEFAULT_COPY);

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
			     short copy_mode = DEFAULT_COPY);
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
			       short copy_mode = DEFAULT_COPY);
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
			    short copy_mode = DEFAULT_COPY);
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


inline ActiveKeyData::ActiveKeyData(bool handle)
{
  if (handle)
    keyDataRep = std::make_shared<ActiveKeyDataRep>();
}


inline ActiveKeyData::ActiveKeyData(const UShortArray& indices):
  keyDataRep(std::make_shared<ActiveKeyDataRep>(indices))
{ }


inline ActiveKeyData::
ActiveKeyData(const UShortArray& indices, const RealVector&   c_params,
	      const IntVector& di_params, const SizetVector& ds_params,
	      short copy_mode):
  keyDataRep(std::make_shared<ActiveKeyDataRep>(indices, c_params, di_params,
						ds_params, copy_mode))
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
continuous_parameters(const RealVector& c_params, short copy_mode)
{
  if (copy_mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(c_params, keyDataRep->continuousHyperParams);
  else if (copy_mode == SHALLOW_COPY) // enforce shallow vector copy
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
discrete_int_parameters(const IntVector& di_params, short copy_mode)
{
  if (copy_mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(di_params, keyDataRep->discreteIntHyperParams);
  else if (copy_mode == SHALLOW_COPY) // enforce shallow vector copy
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
discrete_set_indices(const SizetVector& ds_indices, short copy_mode)
{
  if (copy_mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(ds_indices, keyDataRep->discreteSetHyperParams);
  else if (copy_mode == SHALLOW_COPY) // enforce shallow vector copy
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

  /// minimal constructor
  ActiveKeyRep(unsigned short set_id, unsigned short r_type);
  /// constructor for aggregated key data
  ActiveKeyRep(unsigned short set_id, unsigned short r_type,
	       const std::vector<ActiveKeyData>& data, short copy_mode);
  /// constructor for a single key data
  ActiveKeyRep(unsigned short set_id, unsigned short r_type,
	       const ActiveKeyData& data, short copy_mode);
  /// destructor
  ~ActiveKeyRep();

private:

  //
  //- Heading: Constructors and destructor
  //

  ActiveKeyRep(); ///< default constructor

  //
  //- Heading: Member functions
  //

  /// assign activeKeyDataArray using shallow/deep/default copy
  void assign(const std::vector<ActiveKeyData>& data_vec, short copy_mode);
  /// assign activeKeyDataArray using shallow/deep/default copy
  void assign(const ActiveKeyData& data, short copy_mode);

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

  unsigned short dataSetId;

  // NO_AGGREGATION = 0 (raw only), AGGREGATION_AUGMENTS (enumerate raw + agg),
  // AGGREGATION_REPLACES (no raw, agg only)
  //unsigned short aggregationMode;
  // RAW, DISTINCT_DISCREPANCY, RECURSIVE_DISCREPANCY, ...;
  // retire some special case logic (USHRT_MAX for no model, no resolution)
  //unsigned short aggregationType;
  // SHOULD BE ABLE TO COMBINE / CONDENSE THESE...

  /// indicates type of data reduction, if present, that augments raw data:
  /// > NO_REDUCTION: new data is NOT generated from an aggregation, which may
  ///   or may not be present --> one or more embedded KeyDatas correspond to
  ///   unique data sets, but ActiveKey's array only serves to group them by id
  /// > DISTINCT_DISCREPANCY: new data is generated from a difference of
  ///   model-based data --> both embedded keys and original aggregate key
  ///   manage unique data sets
  /// > RECURSIVE_DISCREPANCY: new data is generated from a difference of
  ///   model-based and surrogate-based (synthetic) data --> both embedded keys
  ///   and original aggregate key manage unique data sets
  /// > Note: if a reduction replaces the raw data (e.g., Dakota responseMode
  ///   of MODEL_DISCREPANCY instead of AGGREGATED_MODELS), then Pecos will
  ///   view this as raw data with no additional reduction.
  unsigned short reductionType;

  std::vector<ActiveKeyData> activeKeyDataArray;

  // Don't need this since the idea is that the ActiveKeyData instances do not
  // reflect the complete inactive state, only the subset identified as part of
  // solution control
  //ActiveKeyData sharedState;
};


inline ActiveKeyRep::ActiveKeyRep()
{ }


inline ActiveKeyRep::
ActiveKeyRep(unsigned short set_id, unsigned short r_type):
  dataSetId(set_id), reductionType(r_type)
{ }


inline ActiveKeyRep::
ActiveKeyRep(unsigned short set_id, unsigned short r_type,
	     const std::vector<ActiveKeyData>& key_data_vec, short copy_mode):
  dataSetId(set_id), reductionType(r_type)
{ assign(key_data_vec, copy_mode); }


inline ActiveKeyRep::
ActiveKeyRep(unsigned short set_id, unsigned short r_type,
	     const ActiveKeyData& key_data, short copy_mode):
  dataSetId(set_id), reductionType(r_type)
{ assign(key_data, copy_mode); }


inline ActiveKeyRep::~ActiveKeyRep()
{ }


inline void ActiveKeyRep::
assign(const std::vector<ActiveKeyData>& key_data_vec, short copy_mode)
{
  if (copy_mode == DEEP_COPY) { // enforce deep copy for each activeKeyData
    size_t i, num_data = key_data_vec.size();
    activeKeyDataArray.resize(num_data);
    for (i=0; i<num_data; ++i)
      activeKeyDataArray[i] = key_data_vec[i].copy();
  }
  else // each activeKeyData shares rep
    activeKeyDataArray = key_data_vec;
}


inline void ActiveKeyRep::assign(const ActiveKeyData& key_data, short copy_mode)
{
  activeKeyDataArray.clear();
  if (copy_mode == DEEP_COPY) activeKeyDataArray.push_back(key_data.copy());
  else                        activeKeyDataArray.push_back(key_data);
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

  /// default handle ctor (no body)
  ActiveKey();
  /// minimal handle + body ctor
  ActiveKey(unsigned short id, unsigned short type);
  /// constructor for aggregated key data
  ActiveKey(unsigned short id, unsigned short type,
	    const std::vector<ActiveKeyData>& key_data_vec, short copy_mode);
  /// constructor for a single key data
  ActiveKey(unsigned short id, unsigned short type,
	    const ActiveKeyData& key_data, short copy_mode);
  /// copy constructor
  ActiveKey(const ActiveKey& key);
  /// destructor
  ~ActiveKey();

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

  /// get reductionType
  unsigned short type() const;
  /// set reductionType
  void type(unsigned short r_type);

  /// get activeKeyDataArray.size()
  size_t data_size() const;

  /// get activeKeyDataArray
  const std::vector<ActiveKeyData>& data() const;
  /// get activeKeyDataArray[i]
  const ActiveKeyData& data(size_t i) const;
  /// get activeKeyDataArray[i]
  ActiveKeyData& data(size_t i);
  /// set activeKeyDataArray
  void data(const std::vector<ActiveKeyData>& key_data_vec, short copy_mode);
  /// set activeKeyDataArray
  void data(const ActiveKeyData& key_data, short copy_mode);

  /// assign data to ActiveKey
  void assign(unsigned short set_id, unsigned short r_type,
	      const std::vector<ActiveKeyData>& key_data_vec, short copy_mode);
  /// assign data to ActiveKey
  void assign(unsigned short set_id, unsigned short r_type,
	      const ActiveKeyData& key_data, short copy_mode);

  /// assign data to ActiveKey
  void append(const std::vector<ActiveKeyData>& key_data_vec, short copy_mode);
  /// assign data to ActiveKey
  void append(const ActiveKeyData& key_data, short copy_mode);
  /// clear activeKeyDataArray
  void clear_data();
  /// clear activeKeyDataArray and reset dataSetId and reductionType
  void clear();

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
  static void aggregate_keys(const std::vector<ActiveKey>& keys,
			     ActiveKey& aggregate_key);

  /// extract two constituent keys from an aggregated key
  void extract_keys(ActiveKey& key1, ActiveKey& key2) const;
  /// extract one or more constituent keys from an aggregated key
  std::vector<ActiveKey> extract_keys() const;
  /// extract a particular constituent key from an aggregated key
  ActiveKey extract_key(size_t key_index) const;

  /// test whether key is an aggregated key (for discrepancy or surplus)
  bool aggregated() const;
  /// test whether key is a singleton key (response data from a single model)
  bool singleton() const;

  /// test whether this key manages a reduction of data from multiple
  /// embedded keys
  bool reduction() const;
  /// test whether this key manages a reduction of data from multiple embedded
  /// keys, where one or more of these keys manage synthetic data (e.g., from
  /// an interpolant using a previous level's model-based data)
  bool recursive_reduction() const;
  /// test whether this key manages a reduction of data from multiple embedded
  /// keys, where all keys manage model-based data
  bool distinct_reduction() const;

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


inline ActiveKey::ActiveKey(unsigned short set_id, unsigned short r_type):
  keyRep(std::make_shared<ActiveKeyRep>(set_id, r_type))
{ }


inline ActiveKey::
ActiveKey(unsigned short set_id, unsigned short r_type,
	  const std::vector<ActiveKeyData>& key_data_vec, short copy_mode):
  keyRep(std::make_shared<ActiveKeyRep>(set_id, r_type, key_data_vec,copy_mode))
{ }


inline ActiveKey::
ActiveKey(unsigned short set_id, unsigned short r_type,
	  const ActiveKeyData& key_data, short copy_mode):
  keyRep(std::make_shared<ActiveKeyRep>(set_id, r_type, key_data, copy_mode))
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
	   keyRep->reductionType      == key.keyRep->reductionType &&
	   keyRep->activeKeyDataArray == key.keyRep->activeKeyDataArray );
}


inline bool ActiveKey::operator!=(const ActiveKey& key) const
{
  return ( keyRep->dataSetId          != key.keyRep->dataSetId ||
	   keyRep->reductionType      != key.keyRep->reductionType ||
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

  if (keyRep->reductionType < kr->reductionType)
    return true;
  else if (kr->reductionType < keyRep->reductionType)
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


inline unsigned short ActiveKey::type() const
{ return keyRep->reductionType; }


inline void ActiveKey::type(unsigned short r_type)
{ keyRep->reductionType = r_type; }


inline size_t ActiveKey::data_size() const
{ return keyRep->activeKeyDataArray.size(); }


inline const std::vector<ActiveKeyData>& ActiveKey::data() const
{ return keyRep->activeKeyDataArray; }


inline const ActiveKeyData& ActiveKey::data(size_t i) const
{ return keyRep->activeKeyDataArray[i]; }


inline ActiveKeyData& ActiveKey::data(size_t i)
{ return keyRep->activeKeyDataArray[i]; }


inline void ActiveKey::
data(const std::vector<ActiveKeyData>& key_data_vec, short copy_mode)
{ keyRep->assign(key_data_vec, copy_mode); }


inline void ActiveKey::data(const ActiveKeyData& key_data, short copy_mode)
{ keyRep->assign(key_data, copy_mode); }


inline void ActiveKey::
assign(unsigned short set_id, unsigned short r_type,
       const std::vector<ActiveKeyData>& key_data_vec, short copy_mode)
{
  if  (keyRep)
    { id(set_id); type(r_type); data(key_data_vec, copy_mode); }
  else
    keyRep = std::make_shared<ActiveKeyRep>(set_id, r_type,
					    key_data_vec, copy_mode);
}


inline void ActiveKey::
assign(unsigned short set_id, unsigned short r_type,
       const ActiveKeyData& key_data, short copy_mode)
{
  if  (keyRep)
    { id(set_id); type(r_type); data(key_data, copy_mode); }
  else
    keyRep = std::make_shared<ActiveKeyRep>(set_id, r_type,
					    key_data, copy_mode);
}


inline void ActiveKey::
append(const std::vector<ActiveKeyData>& key_data_vec, short copy_mode)
{
  std::vector<ActiveKeyData>& act_key_data = keyRep->activeKeyDataArray;
  if (copy_mode == DEEP_COPY) {
    size_t i, len = key_data_vec.size();
    for (i=0; i<len; ++i)
      act_key_data.push_back(key_data_vec[i].copy());
  }
  else
    act_key_data.insert(act_key_data.end(), key_data_vec.begin(),
			key_data_vec.end());
}


inline void ActiveKey::append(const ActiveKeyData& key_data, short copy_mode)
{
  if (copy_mode == DEEP_COPY)
    keyRep->activeKeyDataArray.push_back(key_data.copy());
  else
    keyRep->activeKeyDataArray.push_back(key_data);
}


inline void ActiveKey::clear_data()
{ keyRep->activeKeyDataArray.clear(); }


inline void ActiveKey::clear()
{
  // only clear keyRep's data:
  //clear_data();
  //keyRep->dataSetId     = USHRT_MAX;
  //keyRep->reductionType = NO_REDUCTION;

  // unbind from this keyRep
  keyRep.reset(); // decrements ref count by 1 and deletes if count becomes 0
}


/// deep copy of ActiveKey instance
inline ActiveKey ActiveKey::copy() const
{
  ActiveKey key(keyRep->dataSetId, keyRep->reductionType,
		keyRep->activeKeyDataArray, DEEP_COPY);
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
*/


/** if a key has one data instance, then it is a singleton key; if it has
    more than one, then it represents an aggregation (e.g., a discrepancy
    between two consecutive singleton keys). */
inline bool ActiveKey::aggregated() const
{ return (keyRep->activeKeyDataArray.size() >  1); } // opposite of singleton


inline bool ActiveKey::singleton() const
{ return (keyRep->activeKeyDataArray.size() <= 1); } // opposite of aggregated


inline bool ActiveKey::reduction() const
{ return (keyRep->reductionType != NO_REDUCTION); } //&& aggregated()


inline bool ActiveKey::recursive_reduction() const
{ return (keyRep->reductionType == RECURSIVE_DISCREPANCY); } //&& aggregated()


inline bool ActiveKey::distinct_reduction() const
{ return (keyRep->reductionType == DISTINCT_DISCREPANCY); } //&& aggregated()


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
  // Note: aggregated() check will be correct but rely on calling context
  //       to assign a reduction type
  //aggregate_key.type(...);
  aggregate_key.clear_data();
  if (!empty1) aggregate_key.append(key1.data(), DEEP_COPY);
  if (!empty2) aggregate_key.append(key2.data(), DEEP_COPY);
}


inline void ActiveKey::
aggregate_keys(const std::vector<ActiveKey>& keys,
	       ActiveKey& aggregate_key)
{
  if (keys.empty() && !aggregate_key.is_null()) // leave rep defined
    { aggregate_key.clear(); return; }

  // extract and verify consistency in group number
  size_t k, num_k = keys.size();
  unsigned short id = keys[0].id();
  for (k=1; k<num_k; ++k)
    if (id != keys[k].id()) {
      PCerr << "Error: mismatch in group ids in ActiveKey::aggregate_keys()."
	    << std::endl;
      abort_handler(-1);
    }

  // form aggregate of group + HF form/lev + LF form/lev
  aggregate_key.id(id);
  // Note: aggregated() check will be correct but rely on calling context
  //       to assign a reduction type
  //aggregate_key.type(...);
  aggregate_key.clear_data();
  for (k=0; k<num_k; ++k)
    aggregate_key.append(keys[k].data(), DEEP_COPY);
}


inline void ActiveKey::
extract_keys(ActiveKey& key1, ActiveKey& key2) const
{
  const std::vector<ActiveKeyData>& key_data = data();
  size_t data_size = key_data.size();
  if (data_size < 1 || data_size > 2) { // allow 1 or 2 key data instances
    PCerr << "Error: wrong key data size in ActiveKey::extract_keys(key1, key2)"
	  << std::endl;
    abort_handler(-1);
  }

  unsigned short key_id = id();
  key1 = (data_size == 1) ? *this :
    ActiveKey(key_id, NO_REDUCTION, key_data[0], SHALLOW_COPY);
  if (data_size > 1)
    key2 = ActiveKey(key_id, NO_REDUCTION, key_data[1], SHALLOW_COPY);
  else
    key2.clear();
}


inline std::vector<ActiveKey> ActiveKey::extract_keys() const
{
  const std::vector<ActiveKeyData>& key_data = data();
  size_t k, data_size = key_data.size();
  std::vector<ActiveKey> embedded_keys(data_size);
  if (data_size > 1) { // not a singleton key
    // create new singleton keys from embedded key data
    unsigned short key_id = id();
    for (k=0; k<data_size; ++k)
      embedded_keys[k]
	= ActiveKey(key_id, NO_REDUCTION, key_data[k], SHALLOW_COPY);
  }
  else if (data_size == 1)
    embedded_keys[0] = *this;
  return embedded_keys;
}


inline ActiveKey ActiveKey::extract_key(size_t index) const
{
  const std::vector<ActiveKeyData>& key_data = data();
  size_t data_size = key_data.size();
  if (index == _NPOS) // special value bypasses indexed extraction
    return *this;
  else if (index < data_size)
    return (data_size == 1) ?
      *this :                            // no extraction, already a single key
      ActiveKey(id(), type(), key_data[index], SHALLOW_COPY); // extract single
  else {
    PCerr << "Error: index out of range in ActiveKey::extract_key(index)."
	  << std::endl;
    abort_handler(-1);
    return ActiveKey();
  }
}

} // namespace Pecos

#endif
