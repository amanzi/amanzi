/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Tests some basic functiomality of Operators around Ops, Apply, etc.

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "UnitTest++.h"

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "AmanziVector.hh"
#include "CompositeVectorSpace.hh"
#include "MeshFactory.hh"

#include "Operator.hh"
#include "Operator_Cell.hh"
#include "Op.hh"
#include "Op_Cell_Cell.hh"

using namespace Amanzi;
using namespace Amanzi::Operators;

struct test {
  test() 
      : comm(Amanzi::getDefaultComm())
  {
    // read parameter list
    AmanziMesh::MeshFactory meshfactory(comm, Teuchos::null);
    mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
  }

  Comm_ptr_type comm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
};

// helper function
template<class Op_type>
void opPutScalar(Op_type& op, double scalar) {
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> range({0,0}, {op.data.extent(0), op.data.extent(1)});
  Kokkos::parallel_for(
      "operator_basic",
      range,
      KOKKOS_LAMBDA(const int i, const int j) {
        op.data(i,j) = scalar;
      });
}

SUITE(OPERATOR)
{

  TEST_FIXTURE(test, MAT_VEC_CELL_CELL)
  {
    Teuchos::ParameterList plist;

    // create the space
    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh)
        ->SetGhosted()->SetComponents({"cell"}, {AmanziMesh::CELL}, {1});

    // create the operator
    Operator_Cell oper(cvs.CreateSpace(), plist, OPERATOR_SCHEMA_DOFS_CELL);

    // push back an op
    auto op = Teuchos::rcp(new Op_Cell_Cell("diag", mesh));
    oper.OpPushBack(op);

    // make it a diagonal operator
    opPutScalar(*op, 2.0);

    // create a domain and range
    CompositeVector X(oper.getDomainMap());
    CompositeVector Y(oper.getRangeMap());

    // Check spaces
    CHECK(X.HasComponent("cell"));
    CHECK(Y.HasComponent("cell"));
    CHECK(X.Mesh() == mesh);
    CHECK(Y.Mesh() == mesh);
    X.random();

    // check apply does something reasonable
    oper.apply(X,Y);

    Y.update(-2., X, 1.);
    CHECK_CLOSE(0., Y.norm2(), 1.e-10);
    
  }
}
    
        
      
