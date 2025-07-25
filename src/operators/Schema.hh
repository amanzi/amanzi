/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_SCHEMA_HH_
#define AMANZI_SCHEMA_HH_

#include <string>
#include <tuple>
#include <vector>

#include "BilinearForm.hh"
#include "BilinearFormFactory.hh"
#include "CompositeVectorSpace.hh"
#include "Mesh.hh"
#include "MeshDefs.hh"

#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

class Schema {
 public:
  // default and code compatibility constructors
  Schema() {};
  Schema(AmanziMesh::Entity_kind kind, int nvec) { Init(kind, nvec); }
  Schema(int schema_old) { Init(schema_old); } // old schema must go away FIXME

  ~Schema() {};

  // member functions
  void Init(int schema_old);
  void Init(AmanziMesh::Entity_kind kind, int nvec);
  void Init(Teuchos::RCP<const WhetStone::BilinearForm> form,
            Teuchos::RCP<const AmanziMesh::Mesh> mesh,
            AmanziMesh::Entity_kind base);

  void AddItem(AmanziMesh::Entity_kind kind, WhetStone::DOF_Type type, int num)
  {
    WhetStone::SchemaItem item(kind, type, num);
    items_.push_back(item);
  }

  void Finalize(Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  void ComputeOffset(int c,
                     Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                     std::vector<int>& offset) const;

  // local converters operators/strings/mesh
  int OldSchema() const;

  AmanziMesh::Entity_kind StringToKind(const std::string& name) const;
  WhetStone::DOF_Type StringToType(const std::string& name) const;

  // fancy io
  std::string CreateUniqueName() const;

  // accessers and modifiers
  void set_base(AmanziMesh::Entity_kind base) { base_ = base; }
  AmanziMesh::Entity_kind get_base() const { return base_; }
  const std::vector<WhetStone::SchemaItem>& get_items() { return items_; }

  std::vector<WhetStone::SchemaItem>::const_iterator begin() const { return items_.begin(); }
  std::vector<WhetStone::SchemaItem>::const_iterator end() const { return items_.end(); }
  int size() const { return items_.size(); }

  // only constant access should be used in apps
  const WhetStone::SchemaItem& operator[](int i) const { return items_[i]; }

  // output
  friend std::ostream& operator<<(std::ostream& os, const Schema& s)
  {
    os << "base=" << AmanziMesh::to_string(s.get_base()) << "\n";
    for (auto it = s.begin(); it != s.end(); ++it) {
      os << " item: kind=" << AmanziMesh::to_string(std::get<0>(*it))
         << ", num=" << std::get<2>(*it) << ", type=" << (int)std::get<1>(*it) << "\n";
    }
    return os;
  }

 private:
  AmanziMesh::Entity_kind base_;
  std::vector<WhetStone::SchemaItem> items_;
  std::vector<int> offset_; // starting position of DOF ids

 private:
  explicit Schema(AmanziMesh::Entity_kind kind);
};


// non-member functions
// -- comparison operators
inline bool
operator==(const Schema& s1, const Schema& s2)
{
  if (s1.get_base() != s2.get_base()) return false;
  if (s1.size() != s2.size()) return false;

  for (auto it1 = s1.begin(), it2 = s2.begin(); it1 != s1.end(); ++it1, ++it2) {
    if (std::get<0>(*it1) != std::get<0>(*it2)) return false;
    if (std::get<1>(*it1) != std::get<1>(*it2)) return false;
    if (std::get<2>(*it1) != std::get<2>(*it2)) return false;
  }
  return true;
}

inline bool
operator!=(const Schema& s1, const Schema& s2)
{
  return !(s1 == s2);
}

inline CompositeVectorSpace
cvsFromSchema(const Schema& schema, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, bool ghosted)
{
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh);
  cvs.SetGhosted(ghosted);
  for (const auto& item : schema) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = item;
    cvs.AddComponent(AmanziMesh::to_string(kind), kind, num);
  }
  return cvs;
}


inline Schema
schemaFromPList(const Teuchos::ParameterList& plist,
                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  Schema schema;

  auto& list = plist.sublist("schema");
  auto form = WhetStone::BilinearFormFactory::Create(list, mesh);
  auto base = schema.StringToKind(list.get<std::string>("base"));
  schema.Init(form, mesh, base);

  return schema;
}

} // namespace Operators
} // namespace Amanzi

#endif
