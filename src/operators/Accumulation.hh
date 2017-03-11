/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_ACCUMULATION_HH_
#define AMANZI_OPERATOR_ACCUMULATION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "CompositeVector.hh"

#include "BCs.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Accumulation {
 public:
  Accumulation(AmanziMesh::Entity_kind entity, Teuchos::RCP<Operator> global_op)
    : global_op_(global_op),
      mesh_(Teuchos::null)
  {
    Schema schema(entity);
    InitAccumulation_(schema);
  }

  Accumulation(AmanziMesh::Entity_kind entity, Teuchos::RCP<AmanziMesh::Mesh> mesh)
    : global_op_(Teuchos::null),
      mesh_(mesh)
  {
    Schema schema(entity);
    InitAccumulation_(schema);
  }

  Accumulation(AmanziMesh::Entity_kind entity, Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : global_op_(Teuchos::null),
      mesh_(mesh)
  {
    Schema schema(entity);
    InitAccumulation_(schema);
  }

  Accumulation(Teuchos::ParameterList& plist, Teuchos::RCP<Operator> global_op)
    : global_op_(global_op),
      mesh_(Teuchos::null)
  {
    Schema schema;
    std::string name = plist.get<std::string>("entity kind");

    schema.Init(schema.StringToKind(name));
    InitAccumulation_(schema, plist.get<bool>("surface operator", false));
  }

  Accumulation(Teuchos::ParameterList& plist, Teuchos::RCP<AmanziMesh::Mesh> mesh)
    : global_op_(Teuchos::null),
      mesh_(mesh)
  {
    Schema schema;
    std::string name = plist.get<std::string>("entity kind");

    schema.Init(schema.StringToKind(name));
    InitAccumulation_(schema, plist.get<bool>("surface operator", false));
  }

  Accumulation(Teuchos::ParameterList& plist, Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : global_op_(Teuchos::null),
      mesh_(mesh)
  {
    Schema schema;
    std::string name = plist.get<std::string>("entity kind");
    bool surface = plist.get<bool>("surface operator", false);

    schema.Init(schema.StringToKind(name));
    InitAccumulation_(schema, surface);
  }
  
  Accumulation(const Schema& schema, Teuchos::RCP<Operator> global_op)
    : global_op_(global_op),
      mesh_(Teuchos::null)
  {
    InitAccumulation_(schema, false);
  }

  // update methods
  // -- modifiers for diagonal operators
  void AddAccumulationTerm(const CompositeVector& du,
                           double dT, const std::string& name);

  // -- linearized update methods with storage terms
  void AddAccumulationDelta(const CompositeVector& u0,
                            const CompositeVector& s0, const CompositeVector& ss,
                            double dT, const std::string& name);
  void AddAccumulationDelta(const CompositeVector& u0,
                            double dT, const std::string& name);
  void AddAccumulationDeltaNoVolume(const CompositeVector& u0, const CompositeVector& ss,
                                    const std::string& name);

  // -- operator modification
  void ApplyBCs(const Teuchos::RCP<BCs>& bc);

  // access (for developers only)
  int schema_dofs() { return local_op_schema_.OldSchema(); }
  int schema_prec_dofs() { return global_op_schema_.OldSchema(); }

  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }

 protected:
  void CalculateEntitylVolume_(CompositeVector& entity_volume, const std::string& name);
  void InitAccumulation_(const Schema& schema, bool surface=false);
  Teuchos::RCP<Op> FindOp_(const std::string& name) const;

 protected:
  // operator
  Teuchos::RCP<Operator> global_op_;
  std::vector<Teuchos::RCP<Op> > local_ops_;
  Schema global_op_schema_, local_op_schema_;

  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int ncells_owned;
  int nfaces_owned;
  int nnodes_owned;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


