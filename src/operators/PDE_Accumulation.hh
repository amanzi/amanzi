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
#include "PDE_HelperBCsList.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"

/*!
``PDE_Accumulation`` assembles the discrete form of :math:`\frac{\partial A}{\partial t}`.

This class is usually used as part of a preconditioner, providing the linearization:

.. math::
  \frac{\partial}{\partial A} \left[ \frac{\partial A}{\partial t} \right]_{A_0} i
  = \frac{|\Omega_E|}{\Delta t}

for a grid element :math:`\Omega_E`.


.. _pde-accumulation-spec:
.. admonition:: pde-accumulation-spec

  * `"entity kind`" ``[string]`` **optional** Typically set by the PK
  * `"number of vectors`" ``[int]`` **optional** Typically set by the PK

*/

namespace Amanzi {
namespace Operators {

class PDE_Accumulation : public PDE_HelperBCsList {
 public:
  PDE_Accumulation(AmanziMesh::Entity_kind entity, Teuchos::RCP<Operator> global_op)
    : global_op_(global_op),
      mesh_(Teuchos::null)
  {
    Schema schema(entity, 1);
    InitAccumulation_(schema);
  }

  PDE_Accumulation(AmanziMesh::Entity_kind entity, Teuchos::RCP<AmanziMesh::Mesh> mesh)
    : global_op_(Teuchos::null),
      mesh_(mesh)
  {
    Schema schema(entity, 1);
    InitAccumulation_(schema);
  }

  PDE_Accumulation(AmanziMesh::Entity_kind entity, Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : global_op_(Teuchos::null),
      mesh_(mesh)
  {
    Schema schema(entity, 1);
    InitAccumulation_(schema);
  }

  PDE_Accumulation(Teuchos::ParameterList& plist, Teuchos::RCP<Operator> global_op)
    : global_op_(global_op),
      plist_(plist),
      mesh_(Teuchos::null)
  {
    Schema schema;
    std::string name = plist_.get<std::string>("entity kind");
    int n_vecs = plist_.get<int>("number of vectors", 1);

    schema.Init(schema.StringToKind(name), n_vecs);
    InitAccumulation_(schema, plist_.get<bool>("surface operator", false));
  }

  PDE_Accumulation(Teuchos::ParameterList& plist, Teuchos::RCP<AmanziMesh::Mesh> mesh)
    : global_op_(Teuchos::null),
      plist_(plist),
      mesh_(mesh)
  {
    Schema schema;
    std::string name = plist_.get<std::string>("entity kind");
    int n_vecs = plist_.get<int>("number of vectors", 1);

    schema.Init(schema.StringToKind(name), n_vecs);
    InitAccumulation_(schema, plist_.get<bool>("surface operator", false));
  }

  PDE_Accumulation(Teuchos::ParameterList& plist, Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : global_op_(Teuchos::null),
      plist_(plist),
      mesh_(mesh)
  {
    Schema schema;
    std::string name = plist_.get<std::string>("entity kind");
    int n_vecs = plist_.get<int>("number of vectors", 1);

    schema.Init(schema.StringToKind(name), n_vecs);

    bool surface = plist_.get<bool>("surface operator", false);
    InitAccumulation_(schema, surface);
  }
  
  PDE_Accumulation(const Schema& schema, Teuchos::RCP<Operator> global_op)
    : global_op_(global_op),
      mesh_(Teuchos::null)
  {
    InitAccumulation_(schema, false);
  }

  // update methods
  // -- modifiers for diagonal operators
  void AddAccumulationTerm(const CompositeVector& du, const std::string& name);
  void AddAccumulationTerm(const CompositeVector& du, double dT, const std::string& name, bool volume=true);
  // -- modifiers for diagonal operators and rhs
  void AddAccumulationRhs(const CompositeVector& s1,
                          const CompositeVector& s2,
                          double alpha,
                          const std::string& name,
                          bool volume);
  // -- linearized update methods with storage terms
  void AddAccumulationDelta(const CompositeVector& u0,
                            const CompositeVector& s0, const CompositeVector& ss,
                            double dT, const std::string& name);
  void AddAccumulationDelta(const CompositeVector& u0,
                            double dT, const std::string& name);
  void AddAccumulationDeltaNoVolume(const CompositeVector& u0, const CompositeVector& ss,
                                    const std::string& name);

  // -- operator modification
  void ApplyBCs();

  // access (for developers only)
  int schema_dofs() { return local_op_schema_.OldSchema(); }
  int schema_prec_dofs() { return global_op_schema_.OldSchema(); }

  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }

  Teuchos::RCP<const Op> local_op(int i) const { return local_ops_[i]; }
  Teuchos::RCP<Op> local_op(int i) { return local_ops_[i]; }

 protected:
  void CalculateEntityVolume_(CompositeVector& entity_volume, const std::string& name);
  void InitAccumulation_(const Schema& schema, bool surface=false);
  Teuchos::RCP<Op> FindOp_(const std::string& name) const;

 protected:
  // operator
  Teuchos::RCP<Operator> global_op_;
  std::vector<Teuchos::RCP<Op> > local_ops_;
  Schema global_op_schema_, local_op_schema_;

  Teuchos::ParameterList plist_;
  
  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int ncells_owned;
  int nfaces_owned;
  int nnodes_owned;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


