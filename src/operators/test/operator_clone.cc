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

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "Tensor.hh"

// Operators
#include "Op_Cell_Face.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionMFD.hh"


/* *****************************************************************
* Clone an operator and replace matrices
***************************************************************** */
void Print(const Amanzi::Operators::Operator& global_op)
{
  int nops = global_op.size();
  std::cout << "================\n" << "number of ops: " << nops << std::endl;

  std::set<std::vector<Amanzi::WhetStone::DenseMatrix>*> counts;

  int i(0);
  for (auto it : global_op) {
    std::cout << "Op #" << ++i << std::endl;
    std::cout << " name=" << it->schema_string << std::endl;
    std::cout << " schema: " << it->schema_row();
    std::cout << " matrices=" << &it->matrices << std::endl;
    std::cout << " " << it << std::endl << std::endl;

    counts.insert(&it->matrices);
  }

  CHECK(counts.size() == nops);
}

TEST(OPERATOR_CLONE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: clone operator" << std::endl;

  // create parameter list
  Teuchos::ParameterList op_list;
  op_list.set<std::string>("discretization primary", "mfd: default")
         .set<Teuchos::Array<std::string>>("schema", { "face", "cell" });

  // create a mesh framework
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);

  int ncells = mesh->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int nfaces_wghost = mesh->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);

  // TEST 1
  // create diffusion operator
  auto op = Teuchos::rcp(new PDE_DiffusionMFD(op_list, mesh));
  op->Init(op_list);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto global_op = op->global_operator();
  int nops = global_op->size();
  auto opt = *global_op->begin();

  CHECK_EQUAL(nops, 1);
  Print(*global_op);

  // TEST 2
  // create new op
  std::string name = "new Op";
  auto op2 = Teuchos::rcp(new Op_Cell_Face(name, mesh));
  global_op->OpReplace(op2, 0);
  nops = global_op->size();
  auto ops = *global_op->begin();

  CHECK_EQUAL(nops, 1);
  CHECK(opt.get() != ops.get());
  CHECK(&opt->matrices != &ops->matrices);
  CHECK(opt->diag.get() == nullptr);
  Print(*global_op);

  // TEST 3
  // add new op
  global_op->OpPushBack(opt);
  nops = global_op->size();
  CHECK_EQUAL(nops, 2);
  Print(*global_op);

  auto it0 = global_op->begin();
  auto it1 = it0++;;
  CHECK(&(*it0)->matrices == &opt->matrices);
  CHECK(&(*it1)->matrices == &ops->matrices);

  // TEST 4
  // clone both ops and add
  global_op->OpPushBack(opt->Clone());
  global_op->OpPushBack(ops->Clone());
  Print(*global_op);
}
