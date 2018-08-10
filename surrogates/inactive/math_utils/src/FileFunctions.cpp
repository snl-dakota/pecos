#include "FileFunctions.hpp"

int num_rows_in_file(const std::string& filename)
{
  int numLines = 0;
  std::string line;
  std::ifstream file;
  file.open( filename.c_str() );
  if ( file.is_open() )
    {
      while( !file.eof() )
	{
	  getline( file, line );
	  numLines++;
	}
      file.close();
      return numLines;
    }
  else
    {
      std::string msg = "num_rows_in_file() Unable to open file: ";
      msg += filename;
      throw( std::runtime_error( msg ) );
    }
}

int num_columns_in_file( const std::string& filename, char delimiter )
{
  std::string line;
  std::ifstream file;
  file.open( filename.c_str() );
  if ( file.is_open() )
    {
      getline(file, line);
      std::stringstream  lineStream( line );
      std::string        cell;
      int numCols = 0;
      while(std::getline( lineStream,cell, delimiter ) )
        {
	  numCols++;
        }
      file.close();
      return numCols;
    }
  else
    {
      std::string msg = "num_columns_in_file() Unable to open file: ";
      msg += filename;
      throw(std::runtime_error(msg));
    }
};

/**
 * Determine if a file exists.
 * Note: This will not work if you do not have read permissions
 * to the file, but if you don''t have read, it
 * may as well not exist to begin with.
 */
bool file_exists( const std::string& filename )
{
    FILE* fp = NULL;

    fp = fopen( filename.c_str(), "rb" );
    if( fp != NULL )
    {
        fclose( fp );
        return true;
    }

    return false;
};

void write( std::string filename, std::string s )
{
  std::ofstream file;
  if ( !file_exists( filename ) )
    {
      file.open ( filename.c_str(), std::ios::out ); 
      file << s;
      file.close();
    }
  else
    {
      throw( std::runtime_error(  "The file " + filename + " already exists" ) );
    }
};

void copy_file ( std::string filename_src, std::string filename_dest )
{
  std::ifstream file_src ( filename_src.c_str() );
 
  std::ofstream file_dest ( filename_dest.c_str());

  if ( !file_src.is_open() )
    throw( std::runtime_error( "CopyFile() Source file could not be opened" ) );

  if ( !file_dest.is_open() )
    throw( std::runtime_error( "CopyFile() Destination file could not be opened" ) );

  file_dest << file_src.rdbuf();

  file_src.close();
  file_dest.close();
};

void remove_file( std::string filename )
{
  if( std::remove( filename.c_str() ) != 0 )
    throw( std::runtime_error("RemoveFile() File could not be deleted") );
};

std::vector< std::string> glob_directory( const std::string& pattern )
{
  glob_t glob_result;
  glob( pattern.c_str(), GLOB_TILDE, NULL, &glob_result );
  std::vector<std::string> ret;
  for( unsigned int i = 0; i < glob_result.gl_pathc; ++i )
    {
      ret.push_back( std::string( glob_result.gl_pathv[i] ) );
    }
  globfree( &glob_result );
  return ret;
};

std::vector<int> extract_numbers( const std::string &s )
{
  std::vector<int> numbers;
  const char* digits = {"0123456789"};
  std::string::size_type pos = s.find_first_of( digits );
  std::string workingString = s, num;
  while ( pos != std::string::npos )
    { 
      workingString = workingString.substr( pos );
      pos = workingString.find_first_not_of( digits );
      if ( pos == std::string::npos ) 
	{
	  numbers.push_back( atoi( workingString.c_str() ) );
	  break;
	}
      num = workingString.substr( 0, pos );
      numbers.push_back( atoi( num.c_str() ) );
      workingString = workingString.substr( pos );
      pos = workingString.find_first_of( digits );
    }
  if ( numbers.size() == 0 ) display( "extractNumbers() No Numbers Found" );
  return numbers;
};

bool mkpath( const std::string &path )
{
    bool bSuccess = false;
    int nRC = mkdir( path.c_str(), 0775 );
    if( nRC == -1 )
    {
        switch( errno )
        {
            case ENOENT:
                //parent didn't exist, try to create it
                if( mkpath( path.substr(0, path.find_last_of('/')) ) )
                    //Now, try to create again.
                    bSuccess = 0 == mkdir( path.c_str(), 0775 );
                else
                    bSuccess = false;
                break;
            case EEXIST:
                //Done!
                bSuccess = true;
                break;
            default:
                bSuccess = false;
                break;
        }
    }
    else
        bSuccess = true;
    return bSuccess;
};
