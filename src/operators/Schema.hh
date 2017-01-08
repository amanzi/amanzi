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
  SchemaItem() : location(0), type(0), num(0) {};
  SchemaItem(int location_, int type_, int num_) :
      location(location_),
      type(type_),
      num(num_) {};

  void set(int location_, int type_, int num_) {
    location = location_;
    type = type_;
    num = num_;
  }

 public:
  int location;  // geometric location of DOF
  int type;  // scalar, vector component, derivative, etc.
  int num;  // how many time it is repeated.
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
  void AddItem(int location, int type, int num) {
    SchemaItem item(location, type, num);
    items_.push_back(item);
  }
 
  void Finalize(Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  // local converters operators/strings/mesh
  int OldSchema() const;

  std::string LocationToString(int loc) const;
  AmanziMesh::Entity_kind LocationToMeshID(int loc) const;
  int StringToLocation(std::string& loc) const;
  AmanziMesh::Entity_kind StringToMeshID(std::string& loc) const {
    int id = StringToLocation(loc);
    return LocationToMeshID(id);
  }

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

