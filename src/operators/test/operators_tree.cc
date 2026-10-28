/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Operators

  Exercises a 2x2 TreeOperator with cross-mesh Op_MeshInjection /
  Operator_Schema off-diagonal coupling blocks -- the same machinery used by
  MPCCoupledWaterNonInjected
  (src/physics/ats/src/pks/mpc/mpc_coupled_water_noninjected.cc).
*/

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "AmanziTypes.hh"
#include "MeshFactory.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_Accumulation.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector.hh"
#include "DataStructuresHelpers.hh"
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
    CompositeVectorSpace cvs_surf;
    cvs_surf.SetMesh(mesh_surf)->SetGhosted(true)->SetComponent(
      "cell", AmanziMesh::Entity_kind::CELL, 1);
    auto tvs_surf = Teuchos::rcp(new TreeVectorSpace(cvs_surf.CreateSpace()));

    CompositeVectorSpace cvs_subsurf;
    cvs_subsurf.SetMesh(mesh_subsurf)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    auto tvs_subsurf = Teuchos::rcp(new TreeVectorSpace(cvs_subsurf.CreateSpace()));

    auto tvs_global = Teuchos::rcp(new TreeVectorSpace(comm));
    tvs_global->PushBack(tvs_subsurf);
    tvs_global->PushBack(tvs_surf);

    auto supermap = createSuperMap(*tvs_global);

    auto subsurf_inds = supermap->viewIndices<MemSpace_kind::HOST>(0, "cell", 0);
    CHECK_EQUAL(ncells_subsurf, (int)subsurf_inds.extent(0));
    CHECK_EQUAL(0, subsurf_inds(0));
    CHECK_EQUAL(ncells_subsurf - 1, subsurf_inds(ncells_subsurf - 1));

    auto surf_inds = supermap->viewIndices<MemSpace_kind::HOST>(1, "cell", 0);
    CHECK_EQUAL(ncells_surf, (int)surf_inds.extent(0));
    CHECK_EQUAL(ncells_subsurf, surf_inds(0));
    CHECK_EQUAL(ncells_surf + ncells_subsurf - 1, surf_inds(ncells_surf - 1));

    Vector_type supervec(supermap->getMap());
    CHECK_EQUAL(ncells_subsurf + ncells_surf, (int)supervec.getLocalLength());

    TreeVector tv(tvs_global);
    tv.getSubVector(0)->getData()->putScalar(1.0);
    tv.getSubVector(1)->getData()->putScalar(2.0);
    copyToSuperVector(*supermap, tv, supervec);

    auto supervec_v = supervec.getLocalViewHost(Tpetra::Access::ReadOnly);
    CHECK_EQUAL(1, supervec_v(0, 0));
    CHECK_EQUAL(1, supervec_v(ncells_subsurf - 1, 0));
    CHECK_EQUAL(2, supervec_v(ncells_subsurf, 0));
    CHECK_EQUAL(2, supervec_v(ncells_subsurf + ncells_surf - 1, 0));
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
    Operators::PDE_DiffusionFV diff_subsurf(
      diff_list_subsurf, Teuchos::RCP<const AmanziMesh::Mesh>(mesh_subsurf));
    diff_subsurf.Init();
    Teuchos::RCP<Operators::BCs> bc_subsurf = Teuchos::rcp(
      new Operators::BCs(mesh_subsurf, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    diff_subsurf.SetBCs(bc_subsurf, bc_subsurf);
    diff_subsurf.Setup(Teuchos::null, Teuchos::null, Teuchos::null);

    Teuchos::ParameterList diff_list_surf;
    Operators::PDE_DiffusionFV diff_surf(
      diff_list_surf, Teuchos::RCP<const AmanziMesh::Mesh>(mesh_surf));
    diff_surf.Init();
    Teuchos::RCP<Operators::BCs> bc_surf = Teuchos::rcp(
      new Operators::BCs(mesh_surf, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    diff_surf.SetBCs(bc_surf, bc_surf);
    diff_surf.Setup(Teuchos::null, Teuchos::null, Teuchos::null);

    // create the TreeOperator
    // -- note these operators are square
    auto tv_subsurf = Teuchos::rcp(new TreeVectorSpace(diff_subsurf.global_operator()->getRowMap()));
    auto tv_surf = Teuchos::rcp(new TreeVectorSpace(diff_surf.global_operator()->getRowMap()));

    // -- row and col maps for the coupled operator
    auto tv_global = Teuchos::rcp(new TreeVectorSpace(comm));
    tv_global->PushBack(tv_subsurf);
    tv_global->PushBack(tv_surf);

    // -- now the operator
    Operators::TreeOperator global_op(tv_global, tv_global);
    global_op.set_operator_block(0, 0, diff_subsurf.global_operator());
    global_op.set_operator_block(1, 1, diff_surf.global_operator());

    // now we need to do something with off-diagonal blocks, which are not square.
    // These are purely cell-to-cell, regardless of the diagonal blocks' own
    // (FV, so also cell-only here, but mirroring the general pattern used in
    // mpc_coupled_water_noninjected.cc) schema.
    Operators::Schema cell_schema(AmanziMesh::Entity_kind::CELL, 1);

    Teuchos::ParameterList subsurf_surf_list;
    auto gop_subsurf_surf =
      Teuchos::rcp(new Operators::Operator_Schema(diff_subsurf.global_operator()->getRowMap(),
                                                  diff_surf.global_operator()->getDomainMap(),
                                                  subsurf_surf_list,
                                                  cell_schema,
                                                  cell_schema));

    // assume this is a coupling block for q_exchange_surf_to_subsurf = alpha * (p_surf - p_subsurf)
    // Create an injection map: for surface-cell-local-index sc, the GID of its
    // parent subsurface top cell.
    int ncells_surf_owned = ncells_surf;
    Kokkos::View<GO*, Kokkos::HostSpace> parent_gids_host("coupling injection GIDs",
                                                          ncells_surf_owned);
    for (int sc = 0; sc != ncells_surf_owned; ++sc) {
      AmanziMesh::Entity_ID f = mesh_surf->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
      auto cells = mesh_subsurf->getFaceCells(f);
      AMANZI_ASSERT(cells.size() == 1);
      parent_gids_host(sc) = mesh_subsurf->getEntityGID(AmanziMesh::Entity_kind::CELL, cells(0));
    }
    auto injection = Teuchos::rcp(new Map_type(-1, parent_gids_host, 0, comm));

    auto lop_subsurf_surf =
      Teuchos::rcp(new Operators::Op_MeshInjection(cell_schema,
                                                   mesh_subsurf,
                                                   cell_schema,
                                                   mesh_surf,
                                                   injection,
                                                   true));
    gop_subsurf_surf->OpPushBack(lop_subsurf_surf);
    global_op.set_operator_block(0, 1, gop_subsurf_surf);


    // and the other coupling block
    Teuchos::ParameterList surf_subsurf_list;
    auto gop_surf_subsurf =
      Teuchos::rcp(new Operators::Operator_Schema(diff_surf.global_operator()->getRowMap(),
                                                  diff_subsurf.global_operator()->getDomainMap(),
                                                  surf_subsurf_list,
                                                  cell_schema,
                                                  cell_schema));

    auto lop_surf_subsurf =
      Teuchos::rcp(new Operators::Op_MeshInjection(cell_schema,
                                                   mesh_surf,
                                                   cell_schema,
                                                   mesh_subsurf,
                                                   injection,
                                                   false));
    gop_surf_subsurf->OpPushBack(lop_surf_subsurf);
    global_op.set_operator_block(1, 0, gop_surf_subsurf);

    // now, can we do something with it?
    global_op.SymbolicAssembleMatrix();

    // ok, fill it...
    diff_subsurf.global_operator()->Zero();
    diff_subsurf.UpdateMatrices(Teuchos::null, Teuchos::null);

    diff_surf.global_operator()->Zero();
    diff_surf.UpdateMatrices(Teuchos::null, Teuchos::null);

    gop_surf_subsurf->Zero();
    gop_subsurf_surf->Zero();
    lop_surf_subsurf->diag->putScalar(0.1);
    lop_subsurf_surf->diag->putScalar(0.2);

    global_op.AssembleMatrix();

    // Applying the assembled matrix and the matrix-free (unassembled)
    // operator to the same vector must agree exactly.
    TreeVector x(tv_global);
    x.randomize();

    TreeVector y_mf(tv_global);
    y_mf.putScalar(0.0);
    global_op.apply(x, y_mf, 0.0);

    TreeVector y_asm(tv_global);
    y_asm.putScalar(0.0);
    global_op.ApplyAssembled(x, y_asm);

    auto ymf_sub = y_mf.getSubVector(0)->getData()->viewComponent<MemSpace_kind::HOST>("cell", false);
    auto yasm_sub =
      y_asm.getSubVector(0)->getData()->viewComponent<MemSpace_kind::HOST>("cell", false);
    for (int c = 0; c != (int)ymf_sub.extent(0); ++c) {
      CHECK_CLOSE(yasm_sub(c, 0), ymf_sub(c, 0), 1.e-10);
    }

    auto ymf_surf = y_mf.getSubVector(1)->getData()->viewComponent<MemSpace_kind::HOST>("cell", false);
    auto yasm_surf =
      y_asm.getSubVector(1)->getData()->viewComponent<MemSpace_kind::HOST>("cell", false);
    for (int c = 0; c != (int)ymf_surf.extent(0); ++c) {
      CHECK_CLOSE(yasm_surf(c, 0), ymf_surf(c, 0), 1.e-10);
    }
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
      Teuchos::rcp(new TreeVectorSpace(diff1_subsurf.global_operator()->getRowMap()));
    auto tv2_subsurf =
      Teuchos::rcp(new TreeVectorSpace(diff2_subsurf.global_operator()->getRowMap()));
    auto tv1_surf = Teuchos::rcp(new TreeVectorSpace(diff1_surf.global_operator()->getRowMap()));
    auto tv2_surf = Teuchos::rcp(new TreeVectorSpace(diff2_surf.global_operator()->getRowMap()));

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
    Operators::Schema cell_schema(AmanziMesh::Entity_kind::CELL, 1);

    Teuchos::ParameterList subsurf_surf_list;
    auto gop_subsurf_surf_flow =
      Teuchos::rcp(new Operators::Operator_Schema(diff1_subsurf.global_operator()->getRowMap(),
                                                  diff1_surf.global_operator()->getDomainMap(),
                                                  subsurf_surf_list,
                                                  cell_schema,
                                                  cell_schema));

    // assume this is a coupling block for q_exchange_surf_to_subsurf = alpha * (p_surf - p_subsurf)
    // Create an injection map
    Kokkos::View<GO*, Kokkos::HostSpace> parent_gids_host("coupling injection GIDs", ncells_surf);
    for (int sc = 0; sc != ncells_surf; ++sc) {
      AmanziMesh::Entity_ID f = mesh_surf->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
      auto cells = mesh_subsurf->getFaceCells(f);
      AMANZI_ASSERT(cells.size() == 1);
      parent_gids_host(sc) = mesh_subsurf->getEntityGID(AmanziMesh::Entity_kind::CELL, cells(0));
    }
    auto injection = Teuchos::rcp(new Map_type(-1, parent_gids_host, 0, comm));

    auto lop_subsurf_surf_flow =
      Teuchos::rcp(new Operators::Op_MeshInjection(cell_schema,
                                                   mesh_subsurf,
                                                   cell_schema,
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
      Teuchos::rcp(new Operators::Operator_Schema(diff1_surf.global_operator()->getRowMap(),
                                                  diff1_subsurf.global_operator()->getDomainMap(),
                                                  surf_subsurf_list,
                                                  cell_schema,
                                                  cell_schema));

    auto lop_surf_subsurf_flow =
      Teuchos::rcp(new Operators::Op_MeshInjection(cell_schema,
                                                   mesh_surf,
                                                   cell_schema,
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
    diff1_subsurf.global_operator()->Zero();
    CompositeVector vals1(diff1_subsurf.global_operator()->getDomainMap());
    vals1.putScalar(1.);
    diff1_subsurf.AddAccumulationTerm(vals1, "cell");

    diff2_subsurf.global_operator()->Zero();
    CompositeVector vals2(diff2_subsurf.global_operator()->getDomainMap());
    vals2.putScalar(2.);
    diff2_subsurf.AddAccumulationTerm(vals2, "cell");

    diff1_surf.global_operator()->Zero();
    CompositeVector vals3(diff1_surf.global_operator()->getDomainMap());
    vals3.putScalar(3.);
    diff1_surf.AddAccumulationTerm(vals3, "cell");

    diff2_surf.global_operator()->Zero();
    CompositeVector vals4(diff2_surf.global_operator()->getDomainMap());
    vals4.putScalar(4.);
    diff2_surf.AddAccumulationTerm(vals4, "cell");

    // fill the coupling flux terms
    gop_surf_subsurf_flow->Zero();
    gop_subsurf_surf_flow->Zero();
    lop_surf_subsurf_flow->diag->putScalar(0.1);
    lop_subsurf_surf_flow->diag->putScalar(0.2);

    // fill the tightly coupled flow/energy terms
    dWC_dT.global_operator()->Zero();
    CompositeVector dWC_dT_values(dWC_dT.global_operator()->getDomainMap());
    dWC_dT_values.putScalar(1.e-2);
    dWC_dT.AddAccumulationTerm(dWC_dT_values, "cell");

    dE_dp_surf.global_operator()->Zero();
    CompositeVector dE_dp_surf_values(dE_dp_surf.global_operator()->getDomainMap());
    dE_dp_surf_values.putScalar(2.e-2);
    dE_dp_surf.AddAccumulationTerm(dE_dp_surf_values, "cell");

    op_global.AssembleMatrix();

    // apply it to the ones vector and make sure we get the right answer...

    // apply in assembled form
    auto supermap = op_global.getSuperMap();
    Vector_type ones(supermap->getMap());
    ones.putScalar(1.0);
    {
      auto ones_v = ones.getLocalViewHost(Tpetra::Access::ReadWrite);
      ones_v(0, 0) = 1.0001;
    }
    Vector_type result(supermap->getMap());
    result.putScalar(0.);
    op_global.A()->apply(ones, result);

    CHECK_EQUAL(6, (int)result.getLocalLength()); // size
    {
      auto result_v = result.getLocalViewHost(Tpetra::Access::ReadOnly);
      CHECK_CLOSE(1.0101, result_v(0, 0), 1.e-10); // subsurface flow + dWC_dT
      CHECK_CLOSE(1.21, result_v(2, 0), 1.e-10);   // subsurface flow + dwc_dt + coupling to surface

      CHECK_CLOSE(2.0, result_v(1, 0), 1.e-10); // subsurface energy
      CHECK_CLOSE(2.0, result_v(3, 0), 1.e-10); // subsurface energy

      CHECK_CLOSE(3.1, result_v(4, 0), 1.e-10);  // surface flow + coupling to subsurface
      CHECK_CLOSE(4.02, result_v(5, 0), 1.e-10); // surface energy + dE_dp_surf
    }

    // apply in unassembled form
    TreeVector ones_t(op_global.getDomainMap());
    ones_t.putScalar(1.0);
    // make sure we get the right subsurface gid in our restriction operation
    ones_t.getSubVector(0)
      ->getSubVector(0)
      ->getData()
      ->viewComponent<MemSpace_kind::HOST>("cell", false)(0, 0) = 1.0001;
    TreeVector result_t(op_global.getRangeMap());
    result_t.putScalar(0.0);
    op_global.apply(ones_t, result_t);

    CHECK_CLOSE(1.0101,
               result_t.getSubVector(0)->getSubVector(0)->getData()->viewComponent<MemSpace_kind::HOST>(
                 "cell", false)(0, 0),
               1.e-10); // subsurface flow + dWC_dT
    CHECK_CLOSE(1.21,
               result_t.getSubVector(0)->getSubVector(0)->getData()->viewComponent<MemSpace_kind::HOST>(
                 "cell", false)(1, 0),
               1.e-10); // subsurface flow + dWC_dT + coupling to surface

    CHECK_CLOSE(2.0,
               result_t.getSubVector(0)->getSubVector(1)->getData()->viewComponent<MemSpace_kind::HOST>(
                 "cell", false)(0, 0),
               1.e-10); // subsurface energy
    CHECK_CLOSE(2.0,
               result_t.getSubVector(0)->getSubVector(1)->getData()->viewComponent<MemSpace_kind::HOST>(
                 "cell", false)(1, 0),
               1.e-10); // subsurface energy

    CHECK_CLOSE(3.1,
               result_t.getSubVector(1)->getSubVector(0)->getData()->viewComponent<MemSpace_kind::HOST>(
                 "cell", false)(0, 0),
               1.e-10); // surface flow + coupling to subsurface
    CHECK_CLOSE(4.02,
               result_t.getSubVector(1)->getSubVector(1)->getData()->viewComponent<MemSpace_kind::HOST>(
                 "cell", false)(0, 0),
               1.e-10); // surface energy + dE_dp_surf

    // matrix-free (unassembled) apply and assembled apply must agree exactly.
    Vector_type result2(supermap->getMap());
    copyToSuperVector(*supermap, result_t, result2);
    // result2 <- result2 - result  (update's convention is this = beta*this + alpha*A)
    result2.update(-1., result, 1.);
    double norm2 = result2.norm2();
    CHECK_CLOSE(0.0, norm2, 1.e-10);
  }
}
