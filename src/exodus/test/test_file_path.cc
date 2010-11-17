// -------------------------------------------------------------
/**
 * @file   test_file_path.cc
 * @author William A. Perkins
 * @date Wed Nov 17 09:06:08 2010
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 17, 2010 by William A. Perkins
// Last Change: Wed Nov 17 09:06:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef EXODUS_FILE_DIR
// #error Path to test ExodusII files is not specified
#define EXODUS_FILE_DIR "test_files"
#endif

#include <boost/filesystem.hpp>
namespace bfs = boost::filesystem;

#include <string>


// -------------------------------------------------------------
// test_file_path
// -------------------------------------------------------------
/// A common way to specify a full path to a test Exodus file
std::string
test_file_path(const std::string& fname)
{
  bfs::path thefile(EXODUS_FILE_DIR);
  thefile /= fname;
  return thefile.string();
}


// -------------------------------------------------------------
// split_file_path
// -------------------------------------------------------------
std::string
split_file_path(const std::string& fname)
{
  bfs::path thefile(EXODUS_FILE_DIR);
  thefile /= "split1";
  thefile /= fname;
  return thefile.string();
}
