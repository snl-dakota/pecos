/**
 * @file Printing.hpp
 * @author John D. Jakeman
 * @date 31 October 2011
 * @brief Functions to print various data structures to standard I/O
 */

#ifndef PRINTING_HPP
#define PRINTING_HPP

// print a variable name and its value to std::cout
#define disp( a ) std::cout << #a << ": " << ( a ) << std::endl

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdarg.h>
#include <stdexcept>
#include <set>

#include "linear_algebra.hpp"

template<typename T> 
inline std::string toString( const T& item )
{
  std::string s;
  std::stringstream ss;
  ss.precision( std::numeric_limits<T>::digits10 );
  ss.setf( std::ios::scientific );
  ss << item;
  s=ss.str();
  return s;
};

template <typename T> 
T fromString( const std::string &s )
{
  T t;
  std::istringstream iss( s );
  iss >> t;
  return t;
}

template<typename T> 
void display( const T& item, std::ostream& os = std::cout )
{
  os << item << std::endl;
}

template<typename T> 
void display( const std::vector<T>& v, 
		std::ostream& os = std::cout,
		char delimiter = ',' )
{
  os.setf( std::ios::scientific );
  os.precision( std::numeric_limits<T>::digits10 );
  if ( v.size() == 0 ) display( "The vector is empty" );
  else
    {
      int i;
      for ( i = 0; i < (int)v.size()-1; i++ )
	os << v[i] << delimiter;
      os << v[i] << std::endl;
    }
};

template<typename T> 
void display( const std::vector< std::vector<T> >& v, 
	      std::ostream& os = std::cout,
	      char delimiter = ',')
{
  os.setf( std::ios::scientific );
  os.precision( std::numeric_limits<T>::digits10 );
  int i;
  for ( i = 0; i < v.numRows(); i++ )
    {
      int j;
      for ( j = 0; j < v.numCols()-1; j++ )
	{
	  os << v(i,j) << delimiter;
	}
      os << v(i,j) << std::endl;
    }
  if ( i == 0 ) display( "The 2D vector is empty" );
}

template<typename O, typename T> 
void display( const Teuchos::SerialDenseMatrix<O,T>& v, 
	      std::ostream& os = std::cout,
	      char delimiter = ',')
{
  os.setf( std::ios::scientific );
  os.precision( std::numeric_limits<T>::digits10 );
  int i;
  for ( i = 0; i < v.numRows(); i++ )
    {
      int j;
      for ( j = 0; j < v.numCols()-1; j++ )
	{
	  os << v(i,j) << delimiter;
	}
      os << v(i,j) << std::endl;
    }
  if ( i == 0 ) display( "The matrix is empty" );
}

template<typename O, typename T> 
void display( const Teuchos::SerialDenseVector<O,T>& v, 
	      std::ostream& os = std::cout,
	      char delimiter = ',')
{
  os.setf( std::ios::scientific );
  os.precision( std::numeric_limits<T>::digits10 );
  int i;
  if ( v.length() == 0 ) display( "The vector is empty" );
  else
    {
      for ( i = 0; i < v.length()-1; i++ )
	os << v[i] << delimiter;
      os << v[i] << std::endl;
    }
}

template<typename O,typename T> 
void display( const std::vector< Teuchos::SerialDenseVector<O,T> >& v, 
	      std::ostream& os = std::cout, char delimiter = ',' )
{
  os.setf( std::ios::scientific );
  os.precision( std::numeric_limits<T>::digits10 );
  for  (int i = 0; i < (int)v.size(); i++ )
    {
      display( v[i], os, delimiter );
    }
  if ( v.size() == 0 ) display( "The vector is empty" );
}

template<typename T> 
void display( const std::set<T>& s, std::ostream& os = std::cout )
{
  os.setf( std::ios::scientific );
  os.precision( std::numeric_limits<T>::digits10 );
  typename std::set<T>::iterator it;
  for ( it=s.begin(); it!=s.end(); it++ )
    {
      display( *it, os );
    }  
  if ( s.empty() ) display( "The set is empty" );
};

template<typename T> 
void display( const T* data, int n,
	      std::ostream& os = std::cout, char delimiter = ',' )
{
  os.precision( std::numeric_limits<T>::digits10  );
  os.setf( std::ios::scientific );
  int i;
  for ( i = 0; i < n-1; i++ )
      os << data[i] << delimiter;
  os << data[i] << std::endl; 
};

#endif
