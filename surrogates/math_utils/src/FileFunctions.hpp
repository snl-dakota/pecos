/**
 * @file FileFunctions.hpp
 * @author John D. Jakeman
 * @date 31 October 2011
 * @brief File I/O macros.
 */

#ifndef FILE_FUNCTIONS_HPP
#define FILE_FUNCTIONS_HPP

#include "Printing.hpp"
#include <glob.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <iterator>

int num_rows_in_file( const std::string &filename );

int num_columns_in_file( const std::string &filename, char delimiter = ',' );

/**
 * Determine if a file exists.
 * Note: This will not work if you do not have read permissions
 * to the file, but if you don''t have read, it
 * may as well not exist to begin with.
 */
bool file_exists( const std::string &filename );

void copy_file ( std::string filename_src, std::string filename_dest );

void remove_file( std::string filename );

void write( std::string filename, std::string s );

std::vector< std::string> glob_directory( const std::string &pattern );

std::vector<int> extract_numbers( const std::string &s );

bool mkpath( const std::string &path );

template<typename T> 
void write( const std::vector<T> &data, 
	    const std::string &filename,
	    char delimiter = ',')
{
  if ( data.size() == 0 )
    {
      std::string msg;
      msg = "write( std::vector, filename ) vector is empty. File not written.";
      display( msg );
    }
  else
    {
      std::ofstream file;
      file.open( filename.c_str() );
      if ( file.is_open() )
	{
	  display( data, file, delimiter );
	  file.close();
	}
      else
	{
	  std::string msg;
	  msg = "write( std::vector, filename ) Unable to write the file ";
	  msg += filename;
	  throw( std::runtime_error( msg ) );
	}
    }
};

template<typename O, typename T> 
void write_csv_matrix( const Teuchos::SerialDenseMatrix<O,T> &M, 
		   std::string &filename,
		   char delimiter = ',', bool overwrite = false  )
{
  if ( ( M.numRows() == 0 ) || ( M.numRows() == 0 ) )
    {
      std::string msg;
      msg = "write( matrix, filename ) Matrix is empty. File ";
      msg += filename + " was not written.";
      display( msg );
    }
  else if ( overwrite || !file_exists( filename ) )
    {
      std::ofstream fid;
      fid.open( filename.c_str() );
      if (fid.is_open())
	{
	  display( M, fid, delimiter );
	  fid.close();
	}
      else
	{
	  std::string msg;
	  msg = "write( matrix, filename ) Unable to write the file " + filename;
	  throw( std::runtime_error( msg ) );
	}
    }
  else
    {
      std::string msg;
      msg = "write( matrix, filename ) File " + filename + " already exists.";
      throw( std::runtime_error( msg ) );
    }
}

template<typename O, typename T>
void read_csv_matrix( std::string filename,
	   Teuchos::SerialDenseMatrix<O,T> &M,
	   char delimiter = ',' ){
  int i,j;
  int numRows = num_rows_in_file( filename );
  int numCols = num_columns_in_file( filename, delimiter );
  //Assumes the file is terminated with a "\n" as the last row
  M.reshape( numRows-1, numCols );
  std::ifstream file;
  file.open( filename.c_str()) ;
  if ( file.is_open() ){
    std::string line;
    i = 0;
    bool loop=true;
    while( loop )
      {
	loop = std::getline(file, line);
	j = 0;
	std::stringstream  lineStream( line );
	std::string        cell;
	while( std::getline( lineStream, cell, delimiter ) ){
	  M(i,j) = fromString<T>( cell );
	  j++;
	}
	i++;
      }
    file.close();
  }else {
    std::string msg = "read( filename, SerialDenseMatrix ) Unable to open ";
    msg += "file " + filename;
    throw( std::runtime_error( msg ) );
  }
};

#endif
