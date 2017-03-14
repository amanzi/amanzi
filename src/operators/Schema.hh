/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_SCHEMA_HH_
#define AMANZI_SCHEMA_HH_

#include <string>
#include <vector>

#include "Mesh.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace Operators {

struct SchemaItem {
 public:
  SchemaItem() : type(0), num(0) {};
  SchemaItem(AmanziMesh::Entity_kind kind_, int type_, int num_) :
      kind(kind_),
      type(type_),
      num(num_) {};

  void set(AmanziMesh::Entity_kind kind_, int type_, int num_) {
    kind = kind_;
    type = type_;
    num = num_;
  }

 public:
  AmanziMesh::Entity_kind kind;  // geometric location of DOF
  int type;  // scalar, vector component, derivative, etc.
  int num;  // number of times it is repeated.
};


class Schema {
 public:
  // default and code compatibility constructors
  Schema() {};
  Schema(AmanziMesh::Entity_kind kind, int nvec) { Init(kind, nvec); }
  Schema(int schema_old) { Init(schema_old); }  // old schema must go away FIXME
  ~Schema() {};

  // member functions
  void Init(int schema_old);
  void Init(AmanziMesh::Entity_kind kind, int nvec);
  void Init(Teuchos::ParameterList& plist,
            Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  void SetBase(AmanziMesh::Entity_kind base) { base_ = base; }
  void AddItem(AmanziMesh::Entity_kind kind, int type, int num) {
    SchemaItem item(kind, type, num);
    items_.push_back(item);
  }
 
  void Finalize(Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  void ComputeOffset(int c, Teuchos::RCP<const AmanziMesh::Mesh> mesh, std::vector<int>& offset) const;

  // local converters operators/strings/mesh
  int OldSchema() const;

  std::string KindToString(AmanziMesh::Entity_kind kind) const;
  AmanziMesh::Entity_kind StringToKind(std::string& name) const;

  // fancy io
  std::string CreateUniqueName() const;

  // access
  AmanziMesh::Entity_kind base() const { return base_; }
  const std::vector<SchemaItem>& items() const { return items_; } 
  std::vector<SchemaItem>::const_iterator begin() const { return items_.begin(); }
  std::vector<SchemaItem>::const_iterator end() const { return items_.end(); }

  // output 
  friend std::ostream& operator << (std::ostream& os, const Schema& s) {
    os << "base=" << s.KindToString(s.base()) << "\n";
    for (auto it = s.begin(); it != s.end(); ++it) {
      os << " item: kind=" << s.KindToString(it->kind) << ", num=" << it->num << "\n";
    }
    return os;
  }

 private:
  AmanziMesh::Entity_kind base_;
  std::vector<SchemaItem> items_; 
  std::vector<int> offset_;  // starting position of DOF ids

 private:
  explicit Schema(AmanziMesh::Entity_kind kind);
};


// non-member functions
// -- comparison operators
inline bool operator==(const Schema& s1, const Schema& s2) {
  if (s1.base() != s2.base()) return false;
  if (s1.items().size() != s2.items().size()) return false;

  const std::vector<SchemaItem>& it1 = s1.items();
  const std::vector<SchemaItem>& it2 = s2.items();

  for (int i = 0; i < it1.size(); ++i) {
    if (it1[i].kind != it2[i].kind) return false; 
    if (it1[i].type != it2[i].type) return false; 
    if (it1[i].num != it2[i].num) return false; 
  }
  return true;
}

inline bool operator!=(const Schema& s1, const Schema& s2) {
  return !(s1 == s2);
}

}  // namespace Operators
}  // namespace Amanzi

#endif

