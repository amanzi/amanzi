/*
  State 

  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_STATE_TAG_HH_
#define AMANZI_STATE_TAG_HH_

#include "Key.hh"

namespace Amanzi {

class Tag {
 public:
  Tag() : tag_("") {};
  explicit Tag(const std::string& tag) : tag_(tag) {};

  void set(const std::string& tag) { tag_ = tag; }
  std::string get() const { return tag_; }

  // support of hash
  bool operator==(const Tag& other) const { return (tag_ == other.tag_); }

  bool operator<(const Tag& other) const { return (tag_ < other.tag_); }

 private:
  std::string tag_;
};

// non-member function
inline
Tag make_tag(const Key& key) {
  Tag tag;
  tag.set(key);
  return tag;
}  

}  // namespace Amanzi


namespace std {

template <>
struct hash<Amanzi::Tag> {
  std::size_t operator()(const Amanzi::Tag& tag) const {
    return std::hash<std::string>()(tag.get());
  }
};

}  // namespace

# endif
