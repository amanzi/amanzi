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
  Schema() : base_(0) {};
  Schema(int schema_old) { Init(schema_old); }
  ~Schema() {};

  // member functions
  void Init(int schema_old);

  void SetBase(int base) { base_ = base; }
  void AddItem(AmanziMesh::Entity_kind kind, int type, int num) {
    SchemaItem item(kind, type, num);
    items_.push_back(item);
  }
 
  void Finalize(Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  void ComputeOffset(int c, Teuchos::RCP<const AmanziMesh::Mesh> mesh, std::vector<int>& offset);

  // local converters operators/strings/mesh
  int OldSchema() const;

  std::string KindToString(AmanziMesh::Entity_kind kind) const;
  AmanziMesh::Entity_kind StringToKind(std::string& name) const;

  // fancy io
  std::string CreateUniqueName() const;

  // access
  std::vector<SchemaItem>& items() { return items_; } 
  std::vector<SchemaItem>::const_iterator begin() const { return items_.begin(); }
  std::vector<SchemaItem>::const_iterator end() const { return items_.end(); }

 private:
  int base_;
  std::vector<SchemaItem> items_; 
  std::vector<int> offset_;  // starting position of DOF ids
};

}  // namespace Operators
}  // namespace Amanzi

#endif

