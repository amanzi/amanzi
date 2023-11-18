/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Operators

*/

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "MeshFactory.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_Accumulation.hh"
#include "PDE_CouplingFlux.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector.hh"
#include "TreeVector_Utils.hh"
#include "SuperMap.hh"
#include "TreeOperator.hh"
#include "Operator_Schema.hh"
#include "Op_MeshInjection.hh"

SUITE(SURFACE_SUBSURFACE)
{
  TEST(SURFACE_SUBSURFACE_SUPERMAP)
  {
    // test two different meshes
    using namespace Amanzi;
    auto comm = Amanzi::getDefaultComm();

    Teuchos::ParameterList regions("regions");
    regions.sublist("surface")
      .sublist("region: plane")
      .set<Teuchos::Array<double>>("point", std::vector<double>{ 0., 0., 0. });
    regions.sublist("surface")
      .sublist("region: plane")
      .set<Teuchos::Array<double>>("normal", std::vector<double>{ 0., 0., 1. });

    // create meshes
    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, regions, *comm));
    AmanziMesh::MeshFactory meshfactory(comm, gm);
    auto mesh_subsurf = meshfactory.create(-10, -10, -10, 10, 10, 0, 4, 4, 2);
    int ncells_subsurf =
      mesh_subsurf->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    std::vector<std::string> setnames({ "surface" });
    auto mesh_surf =
      meshfactory.create(mesh_subsurf, setnames, AmanziMesh::Entity_kind::FACE, true);
    int ncells_surf =
      mesh_surf->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    // cvs
    auto cvs_surf = Teuchos::rcp(new CompositeVectorSpace());
    cvs_surf->SetMesh(mesh_surf)->SetGhosted(true)->SetComponent(
      "cell", AmanziMesh::Entity_kind::CELL, 1);
    auto tvs_surf =
      Teuchos::rcp(new TreeVectorSpace((Teuchos::RCP<const CompositeVectorSpace>)cvs_surf));

    auto cvs_subsurf = Teuchos::rcp(new CompositeVectorSpace());
    cvs_subsurf->SetMesh(mesh_subsurf)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    auto tvs_subsurf =
      Teuchos::rcp(new TreeVectorSpace((Teuchos::RCP<const CompositeVectorSpace>)cvs_subsurf));

    TreeVectorSpace tvs_global;
    tvs_global.PushBack(tvs_subsurf);
    tvs_global.PushBack(tvs_surf);

    auto supermap = Operators::createSuperMap(tvs_global);

    CHECK_EQUAL(ncells_subsurf, supermap->Indices(0, "cell", 0).size());
    CHECK_EQUAL(0, supermap->Indices(0, "cell", 0).front());
    CHECK_EQUAL(ncells_subsurf - 1, supermap->Indices(0, "cell", 0).back());

    CHECK_EQUAL(ncells_surf, supermap->Indices(1, "cell", 0).size());
    CHECK_EQUAL(ncells_subsurf, supermap->Indices(1, "cell", 0).front());
    CHECK_EQUAL(ncells_surf + ncells_subsurf - 1, supermap->Indices(1, "cell", 0).back());

    Epetra_Vector supervec(*supermap->Map(), 1);
    CHECK_EQUAL(ncells_subsurf + ncells_surf, supervec.MyLength());

    TreeVector tv(tvs_global);
    tv.SubVector(0)->Data()->PutScalar(1.0);
    tv.SubVector(1)->Data()->PutScalar(2.0);
    Operators::copyToSuperVector(*supermap, tv, supervec);
    CHECK_EQUAL(1, supervec[0]);
    CHECK_EQUAL(1, supervec[ncells_subsurf - 1]);
    CHECK_EQUAL(2, supervec[ncells_subsurf]);
    CHECK_EQUAL(2, supervec[ncells_subsurf + ncells_surf - 1]);
  }


  TEST(SURFACE_SUBSURFACE_DIFFUSION)
  {
    // This tests coupling two operators that are not the same map (and therefore
    // have a non-square off-diagonal block).
    using namespace Amanzi;
    auto comm = Amanzi::getDefaultComm();

    Teuchos::ParameterList regions("regions");
    regions.sublist("surface")
      .sublist("region: plane")
      .set<Teuchos::Array<double>>("point", std::vector<double>{ 0., 0., 0. });
    regions.sublist("surface")
      .sublist("region: plane")
      .set<Teuchos::Array<double>>("normal", std::vector<double>{ 0., 0., 1. });

    // create meshes
    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, regions, *comm));
    AmanziMesh::MeshFactory meshfactory(comm, gm);
    auto mesh_subsurf = meshfactory.create(-10, -10, -10, 10, 10, 0, 3, 3, 3);

    std::vector<std::string> setnames({ "surface" });
    auto mesh_surf =
      meshfactory.create(mesh_subsurf, setnames, AmanziMesh::Entity_kind::FACE, true);
    int ncells_surf =
      mesh_surf->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    // create diffusion operators
    Teuchos::ParameterList diff_list_subsurf;
    Operators::PDE_DiffusionFV diff_subsurf(diff_list_subsurf, mesh_subsurf);
    Teuchos::RCP<Operators::BCs> bc_subsurf = Teuchos::rcp(
      new Operators::BCs(mesh_subsurf, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    diff_subsurf.SetBCs(bc_subsurf, bc_subsurf);
    diff_subsurf.Setup(Teuchos::null, Teuchos::null, Teuchos::null);

    Teuchos::ParameterList diff_list_surf;
    Operators::PDE_DiffusionFV diff_surf(diff_list_surf, mesh_surf);
    Teuchos::RCP<Operators::BCs> bc_surf = Teuchos::rcp(
      new Operators::BCs(mesh_surf, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    diff_surf.SetBCs(bc_surf, bc_surf);
    diff_surf.Setup(Teuchos::null, Teuchos::null, Teuchos::null);

    // create the TreeOperator
    // -- note these operators are square
    auto tv_subsurf =
      Teuchos::rcp(new TreeVectorSpace(diff_subsurf.global_operator()->get_row_map()));
    auto tv_surf = Teuchos::rcp(new TreeVectorSpace(diff_surf.global_operator()->get_row_map()));

    // -- row and col maps for the coupled operator
    auto tv_global = Teuchos::rcp(new TreeVectorSpace(comm));
    tv_global->PushBack(tv_subsurf);
    tv_global->PushBack(tv_surf);

    // -- now the operator
    Operators::TreeOperator global_op(tv_global, tv_global);
    global_op.set_operator_block(0, 0, diff_subsurf.global_operator());
    global_op.set_operator_block(1, 1, diff_surf.global_operator());

    // now we need to do something with off-diagonal blocks, which are not square
    Teuchos::ParameterList subsurf_surf_list;
    auto gop_subsurf_surf =
      Teuchos::rcp(new Operators::Operator_Schema(diff_subsurf.global_operator()->get_row_map(),
                                                  diff_surf.global_operator()->get_col_map(),
                                                  subsurf_surf_list,
                                                  diff_subsurf.global_operator()->schema_row(),
                                                  diff_surf.global_operator()->schema_col()));

    // assume this is a coupling block for q_exchange_surf_to_subsurf = alpha * (p_surf - p_subsurf)
    // Create an injection map
    std::vector<int> parents(
      mesh_surf->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
    const Epetra_Map& cell_map = mesh_subsurf->getMap(AmanziMesh::Entity_kind::CELL, false);
    for (int sc = 0; sc != ncells_surf; ++sc) {
      int f = mesh_surf->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
      auto cells = mesh_subsurf->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);
      parents[sc] = cell_map.GID(cells[0]);
    }
    auto injection = Teuchos::rcp(new Epetra_Map(-1, parents.size(), parents.data(), 0, *comm));

    auto lop_subsurf_surf =
      Teuchos::rcp(new Operators::Op_MeshInjection(diff_subsurf.global_operator()->schema_row(),
                                                   mesh_subsurf,
                                                   diff_surf.global_operator()->schema_col(),
                                                   mesh_surf,
                                                   injection,
                                                   true));
    gop_subsurf_surf->OpPushBack(lop_subsurf_surf);
    global_op.set_operator_block(0, 1, gop_subsurf_surf);


    // and the other coupling block
    Teuchos::ParameterList surf_subsurf_list;
    auto gop_surf_subsurf =
      Teuchos::rcp(new Operators::Operator_Schema(diff_surf.global_operator()->get_row_map(),
                                                  diff_subsurf.global_operator()->get_col_map(),
                                                  subsurf_surf_list,
                                                  diff_surf.global_operator()->schema_row(),
                                                  diff_subsurf.global_operator()->schema_col()));

    auto lop_surf_subsurf =
      Teuchos::rcp(new Operators::Op_MeshInjection(diff_surf.global_operator()->schema_row(),
                                                   mesh_surf,
                                                   diff_subsurf.global_operator()->schema_col(),
                                                   mesh_subsurf,
                                                   injection,
                                                   false));
    gop_surf_subsurf->OpPushBack(lop_surf_subsurf);
    global_op.set_operator_block(1, 0, gop_surf_subsurf);

    // now, can we do something with it?
    global_op.SymbolicAssembleMatrix();


    // ok, fill it...
    diff_subsurf.global_operator()->Init();
    diff_subsurf.UpdateMatrices(Teuchos::null, Teuchos::null);

    diff_surf.global_operator()->Init();
    diff_surf.UpdateMatrices(Teuchos::null, Teuchos::null);

    gop_surf_subsurf->Init();
    gop_subsurf_surf->Init();
    lop_surf_subsurf->diag->PutScalar(0.1);
    lop_subsurf_surf->diag->PutScalar(0.2);

    global_op.AssembleMatrix();

    global_op.A()->Print(std::cout);
  }


  TEST(SURFACE_SUBSURFACE_LOTS_OF_DIAGONALS)
  {
    // This tests coupling two operators that are not the same map (and therefore
    // have a non-square off-diagonal block).
    //
    // We simply use a bunch of diagonal entries to make life easier to track and
    // see if we get the cross-block dependencies we expect.
    using namespace Amanzi;
    auto comm = Amanzi::getDefaultComm();

    Teuchos::ParameterList regions("regions");
    regions.sublist("surface")
      .sublist("region: plane")
      .set<Teuchos::Array<double>>("point", std::vector<double>{ 0., 0., 0. });
    regions.sublist("surface")
      .sublist("region: plane")
      .set<Teuchos::Array<double>>("normal", std::vector<double>{ 0., 0., 1. });

    // create meshes
    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, regions, *comm));
    AmanziMesh::MeshFactory meshfactory(comm, gm);
    auto mesh_subsurf = meshfactory.create(-10, -10, -10, 10, 10, 0, 1, 1, 2);

    std::vector<std::string> setnames({ "surface" });
    auto mesh_surf =
      meshfactory.create(mesh_subsurf, setnames, AmanziMesh::Entity_kind::FACE, true);
    int ncells_surf =
      mesh_surf->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    // create primary (diagonal) entries
    Operators::PDE_Accumulation diff1_subsurf(AmanziMesh::Entity_kind::CELL, mesh_subsurf);
    Operators::PDE_Accumulation diff2_subsurf(AmanziMesh::Entity_kind::CELL, mesh_subsurf);
    Operators::PDE_Accumulation diff1_surf(AmanziMesh::Entity_kind::CELL, mesh_surf);
    Operators::PDE_Accumulation diff2_surf(AmanziMesh::Entity_kind::CELL, mesh_surf);

    // create the TreeOperator
    // -- note these operators are square
    auto tv1_subsurf =
      Teuchos::rcp(new TreeVectorSpace(diff1_subsurf.global_operator()->get_row_map()));
    auto tv2_subsurf =
      Teuchos::rcp(new TreeVectorSpace(diff2_subsurf.global_operator()->get_row_map()));
    auto tv1_surf = Teuchos::rcp(new TreeVectorSpace(diff1_surf.global_operator()->get_row_map()));
    auto tv2_surf = Teuchos::rcp(new TreeVectorSpace(diff2_surf.global_operator()->get_row_map()));

    // -- row and col maps for the coupled operator
    auto tv_subsurf = Teuchos::rcp(new TreeVectorSpace(comm));
    tv_subsurf->PushBack(tv1_subsurf);
    tv_subsurf->PushBack(tv2_subsurf);

    auto tv_surf = Teuchos::rcp(new TreeVectorSpace(comm));
    tv_surf->PushBack(tv1_surf);
    tv_surf->PushBack(tv2_surf);

    auto tv_global = Teuchos::rcp(new TreeVectorSpace(comm));
    tv_global->PushBack(tv_subsurf);
    tv_global->PushBack(tv_surf);

    // -- now the operator
    auto op_subsurf = Teuchos::rcp(new Operators::TreeOperator(tv_subsurf, tv_subsurf));
    op_subsurf->set_operator_block(0, 0, diff1_subsurf.global_operator());
    op_subsurf->set_operator_block(1, 1, diff2_subsurf.global_operator());

    auto op_surf = Teuchos::rcp(new Operators::TreeOperator(tv_surf, tv_surf));
    op_surf->set_operator_block(0, 0, diff1_surf.global_operator());
    op_surf->set_operator_block(1, 1, diff2_surf.global_operator());

    // note these are calls to set_block() to set TreeOperators, not Operators
    Operators::TreeOperator op_global(tv_global, tv_global);
    op_global.set_block(0, 0, op_subsurf);
    op_global.set_block(1, 1, op_surf);


    // first pretend we are the flow + energy coupler in the subsurface, and we need to add a diagonal block
    Operators::PDE_Accumulation dWC_dT(AmanziMesh::Entity_kind::CELL, mesh_subsurf);
    op_subsurf->set_operator_block(0, 1, dWC_dT.global_operator());

    // now maybe the surface needs dE/dp
    Operators::PDE_Accumulation dE_dp_surf(AmanziMesh::Entity_kind::CELL, mesh_surf);
    op_surf->set_operator_block(1, 0, dE_dp_surf.global_operator());

    // Now we need to add things to the coupling block, where surface flow is coupled to subsurface flow
    Teuchos::ParameterList subsurf_surf_list;
    auto gop_subsurf_surf_flow =
      Teuchos::rcp(new Operators::Operator_Schema(diff1_subsurf.global_operator()->get_row_map(),
                                                  diff1_surf.global_operator()->get_col_map(),
                                                  subsurf_surf_list,
                                                  diff1_subsurf.global_operator()->schema_row(),
                                                  diff1_surf.global_operator()->schema_col()));

    // assume this is a coupling block for q_exchange_surf_to_subsurf = alpha * (p_surf - p_subsurf)
    // Create an injection map
    std::vector<int> parents(
      mesh_surf->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
    const Epetra_Map& cell_map = mesh_subsurf->getMap(AmanziMesh::Entity_kind::CELL, false);
    for (int sc = 0; sc != ncells_surf; ++sc) {
      int f = mesh_surf->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
      auto cells = mesh_subsurf->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);
      parents[sc] = cell_map.GID(cells[0]);
    }
    auto injection = Teuchos::rcp(new Epetra_Map(-1, parents.size(), parents.data(), 0, *comm));

    auto lop_subsurf_surf_flow =
      Teuchos::rcp(new Operators::Op_MeshInjection(diff1_subsurf.global_operator()->schema_row(),
                                                   mesh_subsurf,
                                                   diff1_surf.global_operator()->schema_col(),
                                                   mesh_surf,
                                                   injection,
                                                   true));
    gop_subsurf_surf_flow->OpPushBack(lop_subsurf_surf_flow);

    // this is actually one block deeper, so we must create the block to sit in the middle level
    auto dsubsurf_dsurf =
      Teuchos::rcp(new Operators::TreeOperator(op_subsurf->get_row_map(), op_surf->get_col_map()));
    dsubsurf_dsurf->set_operator_block(0, 0, gop_subsurf_surf_flow);

    // and then set THAT block into the highest level
    op_global.set_block(0, 1, dsubsurf_dsurf);

    // and the other coupling block (note in reality there would be at least
    // four, a diagonal block for derivative of the coupling term with respect to
    // the primary pressure, and the corresponding block for derivative with
    // respect to the other domain's primary variable.  Here we only do the
    // latter)
    Teuchos::ParameterList surf_subsurf_list;
    auto gop_surf_subsurf_flow =
      Teuchos::rcp(new Operators::Operator_Schema(diff1_surf.global_operator()->get_row_map(),
                                                  diff1_subsurf.global_operator()->get_col_map(),
                                                  subsurf_surf_list,
                                                  diff1_surf.global_operator()->schema_row(),
                                                  diff1_subsurf.global_operator()->schema_col()));

    auto lop_surf_subsurf_flow =
      Teuchos::rcp(new Operators::Op_MeshInjection(diff1_surf.global_operator()->schema_row(),
                                                   mesh_surf,
                                                   diff1_subsurf.global_operator()->schema_col(),
                                                   mesh_subsurf,
                                                   injection,
                                                   false));
    gop_surf_subsurf_flow->OpPushBack(lop_surf_subsurf_flow);
    auto dsurf_dsubsurf =
      Teuchos::rcp(new Operators::TreeOperator(op_surf->get_row_map(), op_subsurf->get_col_map()));
    dsurf_dsubsurf->set_operator_block(0, 0, gop_surf_subsurf_flow);

    op_global.set_block(1, 0, dsurf_dsubsurf);

    // check the block structure
    std::vector<std::vector<Teuchos::RCP<Operators::Operator>>> leaves;
    std::size_t n_row_leaves = getNumTreeVectorLeaves(*op_global.get_row_map());
    std::size_t n_col_leaves = getNumTreeVectorLeaves(*op_global.get_col_map());
    CHECK_EQUAL(4, n_row_leaves);
    CHECK_EQUAL(4, n_col_leaves);
    leaves.resize(n_row_leaves,
                  std::vector<Teuchos::RCP<Operators::Operator>>(n_col_leaves, Teuchos::null));
    auto shape = Operators::Impl::collectTreeOperatorLeaves(op_global, leaves, 0, 0);
    CHECK_EQUAL(4, shape.first);
    CHECK_EQUAL(4, shape.second);
    for (int i = 0; i != 4; ++i) {
      for (int j = 0; j != 4; ++j) {
        if ((i == j) ||           // diagonal
            (i == 0 && j == 1) || // dWC/dT
            (i == 0 && j == 2) || // subsurf_surf
            (i == 2 && j == 0) || // surf_subsurf
            (i == 3 && j == 2)) { // dE/dp surf
          CHECK(leaves[i][j] != Teuchos::null);
        } else {
          CHECK(leaves[i][j] == Teuchos::null);
        }
      }
    }


    // print capabilities
    std::cout << "Tree Operator Structure:\n" << op_global.PrintDiagnostics() << std::endl;

    // now, can we do something with it?
    op_global.SymbolicAssembleMatrix();

    // ok, fill it...
    diff1_subsurf.global_operator()->Init();
    CompositeVector vals1(diff1_subsurf.global_operator()->DomainMap());
    vals1.PutScalar(1.);
    diff1_subsurf.AddAccumulationTerm(vals1, "cell");

    diff2_subsurf.global_operator()->Init();
    CompositeVector vals2(diff2_subsurf.global_operator()->DomainMap());
    vals2.PutScalar(2.);
    diff2_subsurf.AddAccumulationTerm(vals2, "cell");

    diff1_surf.global_operator()->Init();
    CompositeVector vals3(diff1_surf.global_operator()->DomainMap());
    vals3.PutScalar(3.);
    diff1_surf.AddAccumulationTerm(vals3, "cell");

    diff2_surf.global_operator()->Init();
    CompositeVector vals4(diff2_surf.global_operator()->DomainMap());
    vals4.PutScalar(4.);
    diff2_surf.AddAccumulationTerm(vals4, "cell");

    // fill the coupling flux terms
    gop_surf_subsurf_flow->Init();
    gop_subsurf_surf_flow->Init();
    lop_surf_subsurf_flow->diag->PutScalar(0.1);
    lop_subsurf_surf_flow->diag->PutScalar(0.2);

    // fill the tightly coupled flow/energy terms
    dWC_dT.global_operator()->Init();
    CompositeVector dWC_dT_values(dWC_dT.global_operator()->DomainMap());
    dWC_dT_values.PutScalar(1.e-2);
    dWC_dT.AddAccumulationTerm(dWC_dT_values, "cell");

    dE_dp_surf.global_operator()->Init();
    CompositeVector dE_dp_surf_values(dE_dp_surf.global_operator()->DomainMap());
    dE_dp_surf_values.PutScalar(2.e-2);
    dE_dp_surf.AddAccumulationTerm(dE_dp_surf_values, "cell");

    op_global.AssembleMatrix();
    op_global.A()->Print(std::cout);

    // apply it to the ones vector and make sure we get the right answer...

    // apply in assembled form
    Epetra_Vector ones(*op_global.get_supermap()->Map());
    ones.PutScalar(1.0);
    ones[0] = 1.0001;
    Epetra_Vector result(ones);
    result.PutScalar(0.);
    op_global.A()->Apply(ones, result);
    result.Print(std::cout);

    CHECK_EQUAL(6, result.MyLength());      // size
    CHECK_CLOSE(1.0101, result[0], 1.e-10); // subsurface flow + dWC_dT
    CHECK_CLOSE(1.21, result[2], 1.e-10);   // subsurface flow + dwc_dt + coupling to surface

    CHECK_CLOSE(2.0, result[1], 1.e-10); // subsurface energy
    CHECK_CLOSE(2.0, result[3], 1.e-10); // subsurface energy

    CHECK_CLOSE(3.1, result[4], 1.e-10);  // surface flow + coupling to subsurface
    CHECK_CLOSE(4.02, result[5], 1.e-10); // surface energy + dE_dp_surf


    // apply in unassembled form
    TreeVector ones_t(op_global.DomainMap());
    ones_t.PutScalar(1.0);
    // make sure we get the right subsurface gid in our restriction operation
    (*ones_t.SubVector(0)->SubVector(0)->Data()->ViewComponent("cell", false))[0][0] = 1.0001;
    TreeVector result_t(op_global.RangeMap());
    result_t.PutScalar(0.0);
    op_global.Apply(ones_t, result_t);


    CHECK_CLOSE(1.0101,
                (*result_t.SubVector(0)->SubVector(0)->Data()->ViewComponent("cell", false))[0][0],
                1.e-10); // subsurface flow + dWC_dT
    CHECK_CLOSE(1.21,
                (*result_t.SubVector(0)->SubVector(0)->Data()->ViewComponent("cell", false))[0][1],
                1.e-10); // subsurface flow + dWC_dT + coupling to surface

    CHECK_CLOSE(2.0,
                (*result_t.SubVector(0)->SubVector(1)->Data()->ViewComponent("cell", false))[0][0],
                1.e-10); // subsurface energy
    CHECK_CLOSE(2.0,
                (*result_t.SubVector(0)->SubVector(1)->Data()->ViewComponent("cell", false))[0][1],
                1.e-10); // subsurface energy

    CHECK_CLOSE(3.1,
                (*result_t.SubVector(1)->SubVector(0)->Data()->ViewComponent("cell", false))[0][0],
                1.e-10); // surface flow + coupling to subsurface
    CHECK_CLOSE(4.02,
                (*result_t.SubVector(1)->SubVector(1)->Data()->ViewComponent("cell", false))[0][0],
                1.e-10); // surface energy + dE_dp_surf


    Epetra_Vector result2(result);
    Operators::copyToSuperVector(*op_global.get_supermap(), result_t, result2);
    result2.Update(1, result, -1);
    double norm2(0.);
    result2.Norm2(&norm2);
    CHECK_CLOSE(0.0, norm2, 1.e-10);
  }


  TEST(THREE_LEVEL_HIERARCHY)
  {
    // This tests coupling three block with 1x1, 2x2 and 3x3 diagonal operators.
    using namespace Amanzi;
    auto comm = Amanzi::getDefaultComm();

    // create meshes
    AmanziMesh::MeshFactory meshfactory(comm, Teuchos::null);
    auto mesh1 = meshfactory.create(0.0, 0.0, 1.0, 1.0, 2, 2);
    auto mesh2 = meshfactory.create(0.0, 0.0, 1.0, 1.0, 3, 3);
    auto mesh3 = meshfactory.create(0.0, 0.0, 1.0, 1.0, 4, 4);

    // create primary (diagonal) entries
    Operators::PDE_Accumulation diff1(AmanziMesh::Entity_kind::CELL, mesh1);

    Operators::PDE_Accumulation diff2a(AmanziMesh::Entity_kind::CELL, mesh2);
    Operators::PDE_Accumulation diff2b(AmanziMesh::Entity_kind::CELL, mesh2);

    Operators::PDE_Accumulation diff3a(AmanziMesh::Entity_kind::CELL, mesh3);
    Operators::PDE_Accumulation diff3b(AmanziMesh::Entity_kind::CELL, mesh3);
    Operators::PDE_Accumulation diff3c(AmanziMesh::Entity_kind::CELL, mesh3);

    // create the TreeOperator
    // -- note these operators are square
    auto tv1 = Teuchos::rcp(new TreeVectorSpace(diff1.global_operator()->get_row_map()));
    auto tv2 = Teuchos::rcp(new TreeVectorSpace(diff2a.global_operator()->get_row_map()));
    auto tv3 = Teuchos::rcp(new TreeVectorSpace(diff3a.global_operator()->get_row_map()));

    // -- row and col maps for the coupled operator
    auto tv1all = Teuchos::rcp(new TreeVectorSpace(comm));
    tv1all->PushBack(tv1);

    auto tv2all = Teuchos::rcp(new TreeVectorSpace(comm));
    tv2all->PushBack(tv2);
    tv2all->PushBack(tv2);

    auto tv3all = Teuchos::rcp(new TreeVectorSpace(comm));
    tv3all->PushBack(tv3);
    tv3all->PushBack(tv3);
    tv3all->PushBack(tv3);

    auto tv23 = Teuchos::rcp(new TreeVectorSpace(comm));
    tv23->PushBack(tv2all);
    tv23->PushBack(tv3all);

    auto tv123 = Teuchos::rcp(new TreeVectorSpace(comm));
    tv123->PushBack(tv1all);
    tv123->PushBack(tv23);

    // -- now the operator
    auto op1 = Teuchos::rcp(new Operators::TreeOperator(tv1all, tv1all));
    op1->set_operator_block(0, 0, diff1.global_operator());

    auto op2 = Teuchos::rcp(new Operators::TreeOperator(tv2all, tv2all));
    op2->set_operator_block(0, 0, diff2a.global_operator());
    op2->set_operator_block(1, 1, diff2b.global_operator());

    auto op3 = Teuchos::rcp(new Operators::TreeOperator(tv3all, tv3all));
    op3->set_operator_block(0, 0, diff3a.global_operator());
    op3->set_operator_block(1, 1, diff3b.global_operator());
    op3->set_operator_block(2, 2, diff3c.global_operator());

    auto op23 = Teuchos::rcp(new Operators::TreeOperator(tv23, tv23));
    op23->set_block(0, 0, op2);
    op23->set_block(1, 1, op3);

    auto op123 = Teuchos::rcp(new Operators::TreeOperator(tv123, tv123));
    op123->set_block(0, 0, op1);
    op123->set_block(1, 1, op23);

    // -- add coupling of operators across multiple levels (1 -> 2a)
    //    curently we cannot couple between two levels
    /*
  auto coup12 = Teuchos::rcp(new Operators::TreeOperator(tv1all, tv2all));
  op123->set_block(0, 1, coup12);

  auto cvs0 = op1->get_row_map()->SubVector(0)->Data();
  auto cvs1 = op2->get_col_map()->SubVector(0)->Data();

  auto inds_row0 = std::make_shared<std::vector<std::vector<int> > >(4);
  auto inds_col1 = std::make_shared<std::vector<std::vector<int> > >(4);

  for (int k = 0; k < 4; ++k) {
    (*inds_row0)[k].resize(1);
    (*inds_col1)[k].resize(1);

    (*inds_row0)[k][0] = k;  // local row id for op1
    (*inds_col1)[k][0] = k;  // local column id for op23
  }

  auto op_coupling12 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      plist, cvs0, cvs1, inds_row0, inds_col1));
  coup12->set_operator_block(0, 0, op_coupling12->global_operator());
  */

    // -- add coupling of operators on the same level (2a -> 3a)
    //    we add one off-diagonal block to 2x2 tree-operator op23.
    //    currently it is not possible to select components to couple.
    auto coup23 = Teuchos::rcp(new Operators::TreeOperator(tv2all, tv3all));
    op23->set_block(0, 1, coup23);

    auto cvs1 = op23->get_row_map()->SubVector(0)->SubVector(0)->Data();
    auto cvs2 = op23->get_col_map()->SubVector(1)->SubVector(0)->Data();

    auto inds_row1 = std::make_shared<std::vector<std::vector<int>>>(9);
    auto inds_col2 = std::make_shared<std::vector<std::vector<int>>>(9);

    for (int k = 0; k < 9; ++k) {
      (*inds_row1)[k].resize(1);
      (*inds_col2)[k].resize(1);

      (*inds_row1)[k][0] = k; // local row id for op23
      (*inds_col2)[k][0] = k; // local column id for op23
    }

    Teuchos::ParameterList plist;
    auto op_coupling23 =
      Teuchos::rcp(new Operators::PDE_CouplingFlux(plist, cvs1, cvs2, inds_row1, inds_col2));
    coup23->set_operator_block(0, 1, op_coupling23->global_operator());

    // print capabilities
    std::cout << "Tree Operator Structure:\n" << op123->PrintDiagnostics() << std::endl;

    // add data
    diff1.global_operator()->Init();
    CompositeVector val1(diff1.global_operator()->DomainMap());
    val1.PutScalar(1.0);
    diff1.AddAccumulationTerm(val1, "cell");

    diff2a.global_operator()->Init();
    CompositeVector val2a(diff2a.global_operator()->DomainMap());
    val2a.PutScalar(2.1);
    diff2a.AddAccumulationTerm(val2a, "cell");

    diff2b.global_operator()->Init();
    CompositeVector val2b(diff2b.global_operator()->DomainMap());
    val2b.PutScalar(2.2);
    diff2b.AddAccumulationTerm(val2b, "cell");

    diff3a.global_operator()->Init();
    CompositeVector val3a(diff3a.global_operator()->DomainMap());
    val3a.PutScalar(3.1);
    diff3a.AddAccumulationTerm(val3a, "cell");

    diff3b.global_operator()->Init();
    CompositeVector val3b(diff3b.global_operator()->DomainMap());
    val3b.PutScalar(3.2);
    diff3b.AddAccumulationTerm(val3b, "cell");

    diff3c.global_operator()->Init();
    CompositeVector val3c(diff3c.global_operator()->DomainMap());
    val3c.PutScalar(3.3);
    diff3c.AddAccumulationTerm(val3c, "cell");

    // auto values12 = std::make_shared<std::vector<double> >(4);
    // for (int k = 0; k < 4; ++k) (*values12)[k] = 4.1;
    // op_coupling12->Setup(values12, -1.0);
    // op_coupling12->UpdateMatrices(Teuchos::null, Teuchos::null);

    auto values23 = std::make_shared<std::vector<double>>(9);
    for (int k = 0; k < 9; ++k) { (*values23)[k] = 4.2; }
    op_coupling23->Setup(values23, -1.0);
    op_coupling23->UpdateMatrices(Teuchos::null, Teuchos::null);

    // assemble procedure
    op123->SymbolicAssembleMatrix();
    op123->AssembleMatrix();
    // std::cout << *op123->A() << std::endl;

    // verify matrix
    Epetra_Vector ones(*op123->get_supermap()->Map());
    Epetra_Vector result(ones);

    ones.PutScalar(1.0);
    op123->A()->Apply(ones, result);
    // CHECK_CLOSE(result[0], -3.1, 1e-12);
    CHECK_CLOSE(result[4], -2.1, 1e-12);
  }
}
