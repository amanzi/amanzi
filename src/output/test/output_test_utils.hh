/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! test utilities used in a few tests

#include <string>
#include <fstream>
#include <sstream>

//
// Reads a file and returns a string
//
std::string
read_text_file(const std::string& filename)
{
  std::ifstream t(filename);
  std::stringstream buffer;
  buffer << t.rdbuf();
  return std::move(buffer.str());
}

//
// Compares (text) two files
//
bool
compareTextFiles(const std::string& p1, const std::string& p2)
{
  auto s1 = read_text_file(p1);
  auto s2 = read_text_file(p2);
  return s1 == s2;
}

//
// Compares (bitwise) two arbitrary binary files.
//
bool
compareBinaryFiles(const std::string& p1, const std::string& p2)
{
  std::ifstream f1(p1, std::ifstream::binary | std::ifstream::ate);
  std::ifstream f2(p2, std::ifstream::binary | std::ifstream::ate);

  if (f1.fail() || f2.fail()) {
    return false; // file problem
  }

  if (f1.tellg() != f2.tellg()) {
    return false; // size mismatch
  }

  // seek back to beginning and use std::equal to compare contents
  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
}
