/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  DG methods for linear advection equations.
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "Explicit_TI_RK.hh"
#include "GMVMesh.hh"
#include "CompositeVector.hh"
#include "DG_Modal.hh"
#include "LinearOperatorGMRES.hh"
#include "MeshFactory.hh"
#include "MeshMapsFactory.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"
#include "VectorPolynomial.hh"

// Operators
#include "AnalyticDG06.hh"
#include "AnalyticDG06b.hh"

#include "OperatorAudit.hh"
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"


namespace Amanzi {

template <class AnalyticDG>
class AdvectionFn : public Explicit_TI::fnBase<CompositeVector> {
 public:
  AdvectionFn(Teuchos::ParameterList& plist,
              int nx, double dt, 
              const std::string& filename,
              const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
              int order)
      : plist_(plist), nx_(nx), dt_(dt), 
        filename_(filename), mesh_(mesh),
        ana_(mesh, order) {

    // create global operator 
    // -- upwind flux term
    Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist("flux operator");
    op_flux = Teuchos::rcp(new Operators::PDE_AdvectionRiemann(op_list, mesh));
    global_op_ = op_flux->global_operator();

    // -- volumetric advection term
    op_list = plist.sublist("PK operator").sublist("advection operator");
    op_adv = Teuchos::rcp(new Operators::PDE_Abstract(op_list, global_op_));

    // -- reaction term
    op_list = plist.sublist("PK operator").sublist("inverse mass operator");
    op_mass = Teuchos::rcp(new Operators::PDE_Abstract(op_list, mesh_));
  }

  void Functional(double t, const CompositeVector& u, CompositeVector& f) override {
    int d = mesh_->space_dimension();
    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

    // update velocity coefficient
    WhetStone::VectorPolynomial v;
    auto velc = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost));
    auto velf = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost));

    for (int c = 0; c < ncells_wghost; ++c) {
      ana_.VelocityTaylor(mesh_->cell_centroid(c), t, v); 
      (*velc)[c] = v;
    }

    for (int f = 0; f < nfaces_wghost; ++f) {
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      ana_.VelocityTaylor(mesh_->face_centroid(f), t, v); 
      (*velf)[f] = v * (-normal);
    }

    // modify analytic Taylor expansions
    // CalculateApproximateVelocity(t, dt_, velc, velf);

    // update problem coefficients
    // -- accumulation
    auto K = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_wghost));
    WhetStone::Polynomial Kc(d, 0);
    Kc(0, 0) = 1.0;
    for (int c = 0; c < ncells_wghost; c++) {
      (*K)[c] = Kc;
    }

    // -- source term
    int order = plist_.sublist("PK operator").sublist("flux operator")
                      .template get<int>("method order");
    int nk = (order + 1) * (order + 2) / 2;

    WhetStone::Polynomial sol, src, pc(2, order);
    WhetStone::NumericalIntegration numi(mesh_);

    CompositeVector& rhs = *global_op_->rhs();
    Epetra_MultiVector& rhs_c = *rhs.ViewComponent("cell");
    rhs_c.PutScalar(0.0);

    for (int c = 0; c < ncells; ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
      double volume = mesh_->cell_volume(c);

      ana_.SourceTaylor(xc, t, src);

      for (auto it = pc.begin(); it.end() <= pc.end(); ++it) {
        int n = it.PolynomialPosition();
        int k = it.MonomialOrder();

        double factor = numi.MonomialNaturalScale(k, volume);
        WhetStone::Polynomial cmono(2, it.multi_index(), factor);
        cmono.set_origin(xc);      

        WhetStone::Polynomial tmp = src * cmono;      

        rhs_c[n][c] = numi.IntegratePolynomialCell(c, tmp);
      }
    }

    // -- boundary data
    auto bc = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, Operators::DOF_Type::VECTOR));
    std::vector<int>& bc_model = bc->bc_model();
    std::vector<std::vector<double> >& bc_value = bc->bc_value_vector(nk);

    WhetStone::Polynomial coefs;
    WhetStone::DenseVector data;

    for (int f = 0; f < nfaces_wghost; f++) {
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
          fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
        AmanziGeometry::Point vp = ana_.VelocityExact(xf, t); 

        if (vp * normal < -1e-12) {
          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;

          ana_.SolutionTaylor(xf, t, coefs);
          coefs.GetPolynomialCoefficients(data);

          for (int i = 0; i < nk; ++i) {
            bc_value[f][i] = data(i);
          }
        }
      }
    }

    // populate the global operator
    op_flux->Setup(velc, velf);
    op_flux->UpdateMatrices(velf.ptr());
    op_flux->ApplyBCs(bc, true);

    op_adv->SetupPolyVector(velc);
    op_adv->UpdateMatrices();

    op_mass->SetupPoly(K);
    op_mass->UpdateMatrices(Teuchos::null);

    // calculate functional
    CompositeVector u1(u);
    global_op_->Apply(u, u1);
    u1.Update(1.0, rhs, 1.0);

    // invert vector
    op_mass->global_operator()->Apply(u1, f);

    // statistics: l2 norm of solution
    u.Dot(u, &l2norm);
  }

  // Modify cell and face velocities using a mesh map
  void CalculateApproximateVelocity(
      double t, double dt,
      const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& velc,
      const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& velf) {
    // create new mesh for time interval (t, t + dt)
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    AmanziMesh::MeshFactory factory(&comm);
    factory.set_partitioner(AmanziMesh::Partitioner_type::ZOLTAN_RCB);
    factory.preference(AmanziMesh::FrameworkPreference({AmanziMesh::MSTK}));
   
    Teuchos::RCP<AmanziMesh::Mesh> mesh_new;
    if (filename_ == "square")
      mesh_new = factory(0.0, 0.0, 1.0, 1.0, nx_, nx_, mesh_->geometric_model());
    else 
      mesh_new = factory(filename_, mesh_->geometric_model(), true, true);

    int dim = mesh_->space_dimension();
    AmanziGeometry::Point xv(dim), yv(dim);
    AmanziMesh::Entity_ID_List faces, nodeids;
    AmanziGeometry::Point_List new_positions, final_positions;

    int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    int nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

    for (int v = 0; v < nnodes_wghost; ++v) {
      mesh_->node_get_coordinates(v, &xv);
      yv = ana_.VelocityExact(xv, t); 
      xv += yv * dt;

      nodeids.push_back(v);
      new_positions.push_back(xv);
    }
    mesh_new->deform(nodeids, new_positions, false, &final_positions);

    // create a mesh map at time t
    Teuchos::ParameterList map_list;
    map_list.set<std::string>("method", "CrouzeixRaviart")
            .set<int>("method order", 2)
            .set<std::string>("projector", "H1 harmonic")
            .set<std::string>("map name", "VEM");
  
    WhetStone::MeshMapsFactory maps_factory;
    auto maps = maps_factory.Create(map_list, mesh_, mesh_new);

    // calculate approximate velocities
    double dtfac = 1.0 / dt;
    for (int c = 0; c < ncells_wghost; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      WhetStone::VectorPolynomial v;
      std::vector<WhetStone::VectorPolynomial> vvf;

      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        maps->VelocityFace(f, v);
        vvf.push_back(v);
        (*velf)[f] = (v * mesh_->face_normal(f)) * (-dtfac);
      }

      maps->VelocityCell(c, vvf, (*velc)[c]);
      (*velc)[c] *= dtfac;
    }
  }

 public:
  double l2norm;

  Teuchos::RCP<Operators::PDE_AdvectionRiemann> op_flux;
  Teuchos::RCP<Operators::PDE_Abstract> op_adv;
  Teuchos::RCP<Operators::PDE_Abstract> op_mass;

 private:
  AnalyticDG ana_;
  Teuchos::RCP<Operators::Operator> global_op_;

  const Teuchos::ParameterList plist_;
  int nx_;
  double dt_;
  std::string filename_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

}  // namespace Amanzi


/* *****************************************************************
* This tests exactness of the transient advection scheme for
* dp/dt + div(v p) = f.
* **************************************************************** */
template <class AnalyticDG>
void AdvectionTransient(std::string filename, int nx, int ny, double dt,
                        const Amanzi::Explicit_TI::method_t& rk_method)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  // read parameter list
  std::string xmlFileName = "test/operator_advection_dg_transient.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  int order = plist.sublist("PK operator")
                   .sublist("flux operator").get<int>("method order");

  if (MyPID == 0) std::cout << "\nTest: 2D dG transient advection problem, " << filename 
                            << ", order=" << order << std::endl;

  // create a mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

  MeshFactory meshfactory(&comm);
  meshfactory.set_partitioner(AmanziMesh::Partitioner_type::ZOLTAN_RCB);
  meshfactory.preference(FrameworkPreference({MSTK,STKMESH}));
  RCP<const Mesh> mesh;
  if (nx == 0 || ny == 0)
    mesh = meshfactory(filename, gm, true, true);
  else
    mesh = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny, gm);

  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // create main advection class
  AdvectionFn<AnalyticDG> fn(plist, nx, dt, filename, mesh, order);

  // create initial guess
  CompositeVector& rhs = *fn.op_flux->global_operator()->rhs();
  CompositeVector sol(rhs), sol_next(rhs);
  sol.PutScalar(0.0);

  AnalyticDG ana(mesh, order);
  WhetStone::DG_Modal dg(order, mesh);
  dg.set_basis(WhetStone::TAYLOR_BASIS_NATURAL);

  Epetra_MultiVector& sol_c = *sol.ViewComponent("cell");
  ana.InitialGuess(dg, sol_c, 0.0);

  int nstep(0);
  double t(0.0), tend(1.0), tprint(0.0);
  Explicit_TI::RK<CompositeVector> rk(fn, rk_method, sol);

  while(t < tend - dt/2) {
    rk.TimeStep(t, dt, sol, sol_next);

    sol = sol_next;

    if (MyPID == 0 && std::fabs(t - tprint) < dt/4) {
      tprint += 0.1; 
      printf("t=%9.6f  |p|=%9.6g\n", t, fn.l2norm);

      // visualization
      const Epetra_MultiVector& p = *sol.ViewComponent("cell");
      GMV::open_data_file(*mesh, (std::string)"operators.gmv");
      GMV::start_data();
      GMV::write_cell_data(p, 0, "solution");
      GMV::write_cell_data(p, 1, "gradx");
      GMV::write_cell_data(p, 2, "grady");
      if (order > 1) {
        GMV::write_cell_data(p, 3, "hesxx");
        GMV::write_cell_data(p, 4, "hesxy");
        GMV::write_cell_data(p, 5, "hesyy");
        GMV::close_data_file();
      }
    }

    t += dt;
    nstep++;
  }

  // compute solution error
  sol.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *sol.ViewComponent("cell", false);

  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, tend, pnorm, pl2_err, pinf_err);

  if (MyPID == 0) {
    printf("nx=%3d  L2(p)=%9.6g  Inf(p)=%9.6g\n", nx, pl2_err, pinf_err);
    CHECK(pl2_err < 0.1 / nx);
  }
}


TEST(OPERATOR_ADVECTION_TRANSIENT_DG) {
  AdvectionTransient<AnalyticDG06>("square",  4,  4, 0.1, Amanzi::Explicit_TI::tvd_3rd_order);

  /*
  AdvectionTransient<AnalyticDG06>("square",  20,  20, 0.005, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("square",  40,  40, 0.005 / 2, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("square",  80,  80, 0.005 / 4, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("square", 160, 160, 0.005 / 8, Amanzi::Explicit_TI::tvd_3rd_order);

  AdvectionTransient<AnalyticDG06>("test/triangular8.exo",    8, 0, 0.01, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("test/triangular16.exo",  16, 0, 0.01 / 2, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("test/triangular32.exo",  32, 0, 0.01 / 4, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("test/triangular64.exo",  64, 0, 0.01 / 8, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("test/triangular128.exo",128, 0, 0.01 / 16,Amanzi::Explicit_TI::tvd_3rd_order);

  AdvectionTransient("test/median15x16.exo", 16, 16, 0.05 / 2, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient("test/median32x33.exo", 32, 0, 0.05 / 4, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient("test/median63x64.exo", 64, 0, 0.05 / 8, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient("test/median127x128.exo", 128, 0, 0.05 / 16, Amanzi::Explicit_TI::tvd_3rd_order);
  */ 
}


