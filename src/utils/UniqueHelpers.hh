/*
  Utils

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  This functionality is available in C++14 and is added here for use in C++11.
  FIXME: Remove me when we update!
*/

#ifndef AMANZI_UNIQUE_HELPERS_HH_
#define AMANZI_UNIQUE_HELPERS_HH_

#include <memory>
namespace std {

template<typename T, typename ...Args>
unique_ptr<T> make_unique(Args&& ...args)
{
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

}  // namespace std

#endif
