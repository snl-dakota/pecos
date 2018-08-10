/**
 * \file IndexHashSet.hpp
 * \author John D. Jakeman
 * \date 31 December 2011
 * \brief Decleration of the hash storage class. Warning do not modify 
 * data_ using SetIndex at any stage during looping over the elements of
 * idSet_ as this will change idSet_.end() and loop will continue for ever.
 */

#ifndef INDEX_HASH_SET_HPP
#define INDEX_HASH_SET_HPP

#include <ext/hash_set>
#include <ext/hash_map>
#include <queue>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <set>
 
// Sturcture used in containers to determine if two indices are equal.
template<typename T>
struct IndexEqual
{
  /*
   * Determine if index1 == index2
   * \param index1 a multi-dimensional index
   * \param index1 a multi-dimensional index
   * \return true if index1 == index2, false otherwise
   */
  bool operator()(boost::shared_ptr<T> index1, boost::shared_ptr<T> index2) const
  {
    return (*index1)==(*index2);
  }
};

template<typename T>
struct IndexHash
{
  int operator()(boost::shared_ptr<T> index) const
  {
    return index->hash();
  }
};

template<typename T>
struct IndexHashSet
{
  typedef __gnu_cxx::hash_set<boost::shared_ptr<T>, IndexHash<T>, IndexEqual<T> >  hashSet;
  typedef typename __gnu_cxx::hash_set<boost::shared_ptr<T>, IndexHash<T>, IndexEqual<T> >::iterator iterator;
  typedef typename __gnu_cxx::hash_set<boost::shared_ptr<T>, IndexHash<T>, IndexEqual<T> >::const_iterator const_iterator;
};

template<typename T>
void copy_hash_set( typename IndexHashSet<T>::hashSet& source, 
		    typename IndexHashSet<T>::hashSet& destination )
{
  typename IndexHashSet<T>::iterator iter;
  for( iter = source.begin(); iter != source.end(); iter++ )
    {
      boost::shared_ptr<T> index_copy( new T( **iter ) );
      destination.insert( index_copy );
    }
}


template<typename T, typename C>
void copy_multiset( std::multiset<boost::shared_ptr<T>,C>& source, 
		    std::multiset<boost::shared_ptr<T>,C>& destination )
{
  typename std::multiset<boost::shared_ptr<T>,C>::iterator iter;
  for( iter = source.begin(); iter != source.end(); iter++ )
    {
      boost::shared_ptr<T> index_copy( new T( **iter ) );
      destination.insert( index_copy );
    }
}

template<typename T>
void copy_vector( std::vector< boost::shared_ptr<T> >& source, 
		  std::vector< boost::shared_ptr<T> >& destination )
{
  typename std::vector< boost::shared_ptr<T> >::iterator iter;
  for( iter = source.begin(); iter != source.end(); iter++ )
    {
      boost::shared_ptr<T> index_copy( new T( **iter ) );
      destination.push_back( index_copy );
    }
}

template<typename T>
struct IndexHashMap
{
  typedef __gnu_cxx::hash_map<boost::shared_ptr<T>, int, IndexHash<T>, IndexEqual<T> >  hashMap;
  typedef typename __gnu_cxx::hash_map<boost::shared_ptr<T>, int, IndexHash<T>, IndexEqual<T> >::iterator iterator;
  typedef typename __gnu_cxx::hash_map<boost::shared_ptr<T>, int, IndexHash<T>, IndexEqual<T> >::const_iterator const_iterator;
};

template<typename T>
struct CompareSharedPtr
{
  bool operator()( const boost::shared_ptr<T>& s1, const boost::shared_ptr<T>& s2 ) const
  {
    return (*s1) < (*s2);
  };
};


template<typename T>
void copy_priority_queue( std::priority_queue< boost::shared_ptr<T>, std::vector<boost::shared_ptr<T> >,  CompareSharedPtr<T> >& source, 
			  std::priority_queue< boost::shared_ptr<T>, std::vector<boost::shared_ptr<T> >, CompareSharedPtr<T> >& destination )
{
  boost::shared_ptr<T> pointer( new T() );
  // extract all the pointers from the queue. May be better to migrate 
  // to multiset instead of queue because this would not be necessary
  // in the former
  std::vector< boost::shared_ptr<T>  > list( source.size() ) ;
  int num_items = 0;
  while ( !source.empty() )
    {
      pointer = source.top();
      list[num_items] = pointer;
      source.pop();
      num_items++;
    }
  // reconstuct source and destination from list
  for ( int i = 0; i < num_items; i++ )
    {
      source.push( list[i] );
      boost::shared_ptr<T> item_copy ( new T( *list[i] ) );
      destination.push( item_copy );
    }
}

#ifdef ENABLE_SERIALIZATION
namespace boost { 
  namespace serialization {
    template<class Archive, typename T >
    void save( Archive & ar, const __gnu_cxx::hash_set<boost::shared_ptr<T>, IndexHash<T>, IndexEqual<T> >& set, 
	       const unsigned int version )
    {
      int len = set.size();
      ar & len;
      typename __gnu_cxx::hash_set<boost::shared_ptr<T>, IndexHash<T>, IndexEqual<T> >::iterator iter;
      for( iter = set.begin(); iter != set.end(); iter++ )
	ar & (*iter);
    } 

    template<class Archive, typename T >
    void load( Archive & ar, 
	       __gnu_cxx::hash_set<boost::shared_ptr<T>, IndexHash<T>, IndexEqual<T> >& set, 
	       const unsigned int version )
    {
      int len = 0;
      ar & len;
      for ( int m = 0 ; m < len; m++ )
	{
	  boost::shared_ptr<T> tmp( new T() );
	  ar & tmp;
	  set.insert( tmp );
	}
    } 

    template<class Archive, typename T>
    inline void serialize( Archive & ar, 
			   __gnu_cxx::hash_set<boost::shared_ptr<T>, IndexHash<T>, IndexEqual<T> >& set,
			   const unsigned int version)
    {
      split_free( ar, set, version ); 
    }

    template<class Archive, typename T >
    void save( Archive & ar, const std::priority_queue< boost::shared_ptr<T>, std::vector<boost::shared_ptr<T> >, CompareSharedPtr<T> >& queue, 
	       const unsigned int version )
    {
      int len = queue.size();
      ar & len;
      boost::shared_ptr<T> pointer( new T() );
      // extract all the pointers from the queue. May be better to migrate 
      // to multiset instead of queue because this would not be necessary
      // in the former
      std::vector< boost::shared_ptr<T>  > list( len );
      int num_items = 0;
      while ( !queue.empty() )
	{
	  pointer = queue.top();
	  list[num_items] = pointer;
	  ar & pointer; // todo need to make serialization of shared_ptr
	  const_cast< std::priority_queue< boost::shared_ptr<T>, std::vector<boost::shared_ptr<T> >, CompareSharedPtr<T> >& >(queue).pop();
	  num_items++;
	}
      // reconstuct source and destination from list
      for ( int i = 0; i < num_items; i++ )
	{
	  const_cast< std::priority_queue< boost::shared_ptr<T>, std::vector<boost::shared_ptr<T> >, CompareSharedPtr<T> >& >(queue).push( list[i] );
	}
    } 

    template<class Archive, typename T >
    void load( Archive & ar, std::priority_queue< boost::shared_ptr<T>, std::vector<boost::shared_ptr<T> >, CompareSharedPtr<T> >& queue, 
	       const unsigned int version )
    {
      int len = 0;
      ar & len;
      for ( int m = 0 ; m < len; m++ )
	{
	  boost::shared_ptr<T> pointer( new T() );
	  ar & pointer;
	  queue.push( pointer );
	}
    } 

    template<class Archive, typename T>
    inline void serialize( Archive & ar, 
			   std::priority_queue< boost::shared_ptr<T>, std::vector<boost::shared_ptr<T> >, CompareSharedPtr<T> >& queue,
			   const unsigned int version )
    {
      split_free( ar, queue, version ); 
    }

    template<class Archive, typename T >
    void save( Archive & ar, const std::vector< boost::shared_ptr<T> >& list, 
	       const unsigned int version )
    {
      int len = list.size();
      ar & len;
      for( int m = 0 ; m < len; m++ )
	ar & list[m];
    } 

    template<class Archive, typename T >
    void load( Archive & ar, std::vector< boost::shared_ptr<T> >& list, 
	       const unsigned int version )
    {
      int len = 0;
      ar & len;
      list.resize( len );
      for ( int m = 0 ; m < len; m++ )
	{
	  boost::shared_ptr<T> tmp( new T() );
	  ar & tmp;
	  list[m] = tmp;
	}
    } 

    template<class Archive, typename T>
    inline void serialize( Archive & ar, 
			   std::vector< boost::shared_ptr<T> >& list,
			   const unsigned int version)
    {
      split_free( ar, list, version ); 
    }
  }
}
#endif // ENABLE_SERIALIZATION

#endif // INDEX_HASH_SET_HPP
