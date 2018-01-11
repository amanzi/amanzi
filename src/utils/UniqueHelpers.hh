// Amanzi
//
// Author: Ethan Coon
//
// This functionality is available in C++14 and is added here for use in C++11.
// Remove me when we update!
//

#ifndef AMANZI_UNIQUE_HELPERS_HH_
#define AMANZI_UNIQUE_HELPERS_HH_


#include <memory>
namespace std {

template<typename T, typename ...Args>
unique_ptr<T> make_unique( Args&& ...args )
{
  return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

} // namespace std

#endif
