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

#include "MeshDefs.hh"

namespace Amanzi {
namespace Operators {

struct SchemaItem {
 public:
  SchemaItem() : location(0), type(0), num(0) {};

  void set(int location_, int type_, int num_) {
    location = location_;
    type = type_;
    num = num_;
  }

 public:
  int location;
  int type;
  int num; 
};


class Schema {
 public:
  // default and code compatibility constructors
  Schema() : base_(0) {};
  Schema(int schema_old) { Init(schema_old); }
  ~Schema() {};

  // member functions
  void Init(int schema_old);
  int OldSchema() const;

  std::string LocationName(int loc) const;
  AmanziMesh::Entity_kind LocationMeshID(int loc) const;

  std::string CreateUniqueName() const;

  // access
  std::vector<SchemaItem>& items() { return items_; } 

 private:
  int base_;
  std::vector<SchemaItem> items_; 
};

}  // namespace Operators
}  // namespace Amanzi

#endif

