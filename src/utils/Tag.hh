/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  State

*/

#pragma once

#include <iostream>
#include "errors.hh"

namespace Amanzi {

class Tag {
 public:
  Tag()
    : tag_("")
  {}
  explicit Tag(const std::string& tag)
    : tag_(tag)
  {}
  Tag(const Tag& other) = default;

  void set(const std::string& key) { tag_ = key; }
  std::string get() const { return tag_; }

  // support of hash
  bool operator==(const Tag& other) const { return (tag_ == other.tag_); }
  bool operator!=(const Tag& other) const { return !(*this == other); }

  bool operator<(const Tag& other) const { return (tag_ < other.tag_); }
  bool operator>(const Tag& other) const { return (tag_ > other.tag_); }
  bool operator<=(const Tag& other) const { return !(*this > other); }
  bool operator>=(const Tag& other) const { return !(*this < other); }

  friend std::ostream& operator<<(std::ostream& os, const Tag& t)
  {
    os << t.get();
    return os;
  }

 private:
  std::string tag_;
};


inline Errors::Message&
operator<<(Errors::Message& message, const Tag& tag)
{
  message.add_data(tag.get());
  return message;
}


// non-member function
inline Tag
make_tag(const std::string& key)
{
  Tag tag;
  tag.set(key);
  return tag;
}


} // namespace Amanzi


namespace std {

template<>
struct hash<Amanzi::Tag> {
  std::size_t operator()(const Amanzi::Tag& tag) const
  {
    return std::hash<std::string>()(tag.get());
  }
};

} // namespace std
