/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_SCHEMA_HH_
#define AMANZI_SCHEMA_HH_

#include <string>
#include <vector>

#include "CompositeVectorSpace.hh"
#include "Mesh.hh"
#include "MeshDefs.hh"

#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

struct SchemaItem {
 public:
  SchemaItem() : type(DOF_Type::SCALAR), num(0){};
  SchemaItem(AmanziMesh::Entity_kind kind_, DOF_Type type_, int num_)
    : kind(kind_), type(type_), num(num_){};

  void set(AmanziMesh::Entity_kind kind_, DOF_Type type_, int num_)
  {
    kind = kind_;
    type = type_;
    num = num_;
  }

 public:
  AmanziMesh::Entity_kind kind; // geometric location of DOF
  DOF_Type type; // scalar, vector component, derivative, moment, etc.
  int num;       // number of times it is repeated.
};


class Schema {
 public:
  // default and code compatibility constructors
  Schema(){};
  Schema(AmanziMesh::Entity_kind kind, int nvec) { Init(kind, nvec); }
  Schema(int schema_old) { Init(schema_old); } // old schema must go away FIXME
  ~Schema(){};

  // member functions
  void Init(int schema_old);
  void Init(AmanziMesh::Entity_kind kind, int nvec);
  void Init(Teuchos::ParameterList& plist,
            Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  void SetBase(AmanziMesh::Entity_kind base) { base_ = base; }
  void AddItem(AmanziMesh::Entity_kind kind, DOF_Type type, int num)
  {
    SchemaItem item(kind, type, num);
    items_.push_back(item);
  }

  void Finalize(Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  void ComputeOffset(int c, Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                     std::vector<int>& offset) const;

  // local converters operators/strings/mesh
  int OldSchema() const;

  std::string KindToString(AmanziMesh::Entity_kind kind) const;
  AmanziMesh::Entity_kind StringToKind(std::string& name) const;
  DOF_Type StringToType(std::string& name) const;

  // fancy io
  std::string CreateUniqueName() const;

  // access
  AmanziMesh::Entity_kind base() const { return base_; }
  std::vector<SchemaItem>::const_iterator begin() const
  {
    return items_.begin();
  }
  std::vector<SchemaItem>::const_iterator end() const { return items_.end(); }
  int size() const { return items_.size(); }

  // output
  friend std::ostream& operator<<(std::ostream& os, const Schema& s)
  {
    os << "base=" << s.KindToString(s.base()) << "\n";
    for (auto it = s.begin(); it != s.end(); ++it) {
      os << " item: kind=" << s.KindToString(it->kind) << ", num=" << it->num
         << ", type=" << (int)it->type << "\n";
    }
    return os;
  }

 private:
  AmanziMesh::Entity_kind base_;
  std::vector<SchemaItem> items_;
  std::vector<int> offset_; // starting position of DOF ids

 private:
  explicit Schema(AmanziMesh::Entity_kind kind);
};


// non-member functions
// -- comparison operators
inline bool
operator==(const Schema& s1, const Schema& s2)
{
  if (s1.base() != s2.base()) return false;
  if (s1.size() != s2.size()) return false;

  for (auto it1 = s1.begin(), it2 = s2.begin(); it1 != s1.end(); ++it1, ++it2) {
    if (it1->kind != it2->kind) return false;
    if (it1->type != it2->type) return false;
    if (it1->num != it2->num) return false;
  }
  return true;
}

inline bool
operator!=(const Schema& s1, const Schema& s2)
{
  return !(s1 == s2);
}

inline CompositeVectorSpace
cvsFromSchema(const Schema& schema,
              const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh);
  for (const auto& item : schema) {
    cvs.AddComponent(
      AmanziMesh::entity_kind_string(item.kind), item.kind, item.num);
  }
  return cvs;
}

} // namespace Operators
} // namespace Amanzi

#endif
