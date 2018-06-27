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
#include "AnalyticDG07.hh"

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
              const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
              Teuchos::RCP<WhetStone::DG_Modal> dg,  
              bool conservative_form, std::string weak_form);

  // functional in dy/dt = F(y)
  void FunctionalTimeDerivative(double t, const CompositeVector& u, CompositeVector& f) override;

  // modify cell and face velocities using a mesh map
  void ApproximateVelocity_Projection(
      double t, double dt,
      const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& velc,
      const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& velf,
      const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& divc);

  void ApproximateVelocity_LevelSet(
      const CompositeVector& u,
      const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& velc,
      const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& velf,
      const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& divc);

 public:
  double l2norm;

  Teuchos::RCP<Operators::PDE_AdvectionRiemann> op_flux;
  Teuchos::RCP<Operators::PDE_Abstract> op_adv;
  Teuchos::RCP<Operators::PDE_Abstract> op_mass;
  Teuchos::RCP<Operators::PDE_Abstract> op_reac;

 private:
  AnalyticDG ana_;
  Teuchos::RCP<Operators::Operator> global_op_;
  const Teuchos::ParameterList plist_;

  int nx_;
  double dt_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_new_;

  int order_;
  Teuchos::RCP<WhetStone::DG_Modal> dg_;  

  double weak_sign_;
  bool conservative_form_;
  bool setup_, high_order_velf_, level_set_velf_;
  std::string pk_name_;
};


/* *****************************************************************
* Constructor
***************************************************************** */
template <class AnalyticDG>
AdvectionFn<AnalyticDG>::AdvectionFn(
    Teuchos::ParameterList& plist, int nx, double dt, 
    const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
    Teuchos::RCP<WhetStone::DG_Modal> dg,  
    bool conservative_form, std::string weak_form)
    : plist_(plist), nx_(nx), dt_(dt), 
      mesh_(mesh),
      dg_(dg), 
      ana_(mesh, dg->order()),
      conservative_form_(conservative_form)
{
  if (weak_form == "dual") {
    weak_sign_ = 1.0;
    pk_name_ = "PK operator";
  } else {
    weak_sign_ = -1.0;
    pk_name_ = "PK operator: primal";
  }

  // create global operator 
  // -- upwind flux term
  Teuchos::ParameterList op_list = plist.sublist(pk_name_).sublist("flux operator");
  op_flux = Teuchos::rcp(new Operators::PDE_AdvectionRiemann(op_list, mesh));
  global_op_ = op_flux->global_operator();

  // -- volumetric advection term
  op_list = plist.sublist(pk_name_).sublist("advection operator");
  op_adv = Teuchos::rcp(new Operators::PDE_Abstract(op_list, global_op_));

  // -- accumulation term
  op_list = plist.sublist(pk_name_).sublist("inverse mass operator");
  op_mass = Teuchos::rcp(new Operators::PDE_Abstract(op_list, mesh_));

  // -- reaction term
  if (!conservative_form_) {
    op_list = plist.sublist(pk_name_).sublist("reaction operator");
    op_reac = Teuchos::rcp(new Operators::PDE_Abstract(op_list, global_op_));
  }

  order_ = dg_->order();

  // create auxiliary mesh
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  AmanziMesh::MeshFactory factory(&comm);
  factory.set_partitioner(AmanziMesh::Partitioner_type::ZOLTAN_RCB);
  factory.preference(AmanziMesh::FrameworkPreference({AmanziMesh::MSTK}));
   
  std::string name = plist.get<std::string>("file name");
  if (name == "square")
    mesh_new_ = factory(0.0, 0.0, 1.0, 1.0, nx_, nx_, mesh_->geometric_model());
  else 
    mesh_new_ = factory(name, mesh_->geometric_model(), true, true);

  // cotrol variables
  name = plist.get<std::string>("face velocity method");
  high_order_velf_ = (name == "high order");
  level_set_velf_ = (name == "level set");

  setup_ = true;
}


/* *****************************************************************
* Functional F in dy/dt = F(y)
***************************************************************** */
template <class AnalyticDG>
void AdvectionFn<AnalyticDG>::FunctionalTimeDerivative(
    double t, const CompositeVector& u, CompositeVector& f)
{
  int d = mesh_->space_dimension();
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // update velocity coefficient
  WhetStone::VectorPolynomial v;
  auto velc = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost));
  auto velf = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost));
  auto divc = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_wghost));

  for (int c = 0; c < ncells_wghost; ++c) {
    ana_.VelocityTaylor(mesh_->cell_centroid(c), t, v); 
    (*velc)[c] = v;
    (*velc)[c] *= -weak_sign_;
  }

  for (int f = 0; f < nfaces_wghost; ++f) {
    ana_.VelocityTaylor(mesh_->face_centroid(f), t, v); 
    (*velf)[f] = v * (mesh_->face_normal(f) * weak_sign_);
  }

  if (! conservative_form_) {
    for (int c = 0; c < ncells_wghost; ++c) {
      (*divc)[c] = Divergence((*velc)[c]);
    }
  }

  // modify analytic Taylor expansions
  if (level_set_velf_) {
    ApproximateVelocity_LevelSet(u, velc, velf, divc);
  } else {
    ApproximateVelocity_Projection(t, dt_, velc, velf, divc);
  }

  // update problem coefficients
  // -- accumulation
  auto K = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_wghost));
  WhetStone::Polynomial Kc(d, 0);
  Kc(0, 0) = 1.0;
  for (int c = 0; c < ncells_wghost; c++) {
    (*K)[c] = Kc;
  }

  // -- source term
  int nk = (order_ + 1) * (order_ + 2) / 2;

  WhetStone::Polynomial sol, src, pc(2, order_);
  WhetStone::NumericalIntegration numi(mesh_);

  CompositeVector& rhs = *global_op_->rhs();
  Epetra_MultiVector& rhs_c = *rhs.ViewComponent("cell");
  rhs_c.PutScalar(0.0);

  /*
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
  */

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
        if (weak_sign_ < 0.0) bc_model[f] = Operators::OPERATOR_BC_DIRICHLET_TYPE2;

        ana_.SolutionTaylor(xf, t, coefs);
        coefs.GetPolynomialCoefficients(data);

        for (int i = 0; i < nk; ++i) {
          bc_value[f][i] = data(i);
        }
      } else if (weak_sign_ < 0.0) {
        bc_model[f] = Operators::OPERATOR_BC_REMOVE;
      }
    }
  }

  // populate the global operator
  op_flux->SetBCs(bc, bc);
  op_flux->Setup(velc, velf);
  op_flux->UpdateMatrices(velf.ptr());
  op_flux->ApplyBCs(true, true, true);

  op_adv->SetupPolyVector(velc);
  op_adv->UpdateMatrices();

  if (!conservative_form_) {
    op_reac->SetupPoly(divc);
    op_reac->UpdateMatrices();
  }

  if (setup_) {
    op_mass->SetupPoly(K);
    op_mass->UpdateMatrices(Teuchos::null);
    setup_ = false;
  }

  // calculate functional
  CompositeVector u1(u);
  global_op_->ComputeResidual(u, u1);

  // invert vector
  op_mass->global_operator()->Apply(u1, f);

  // statistics: l2 norm of solution
  u.Dot(u, &l2norm);
}


/* *****************************************************************
* Change original definitions of velocities: projection algorithm
***************************************************************** */
template <class AnalyticDG>
void AdvectionFn<AnalyticDG>::ApproximateVelocity_Projection(
    double t, double dt,
    const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& velc,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& velf,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& divc)
{
  // deform auxiliary mesh for time interval (t, t + dt)
  int dim = mesh_->space_dimension();
  AmanziGeometry::Point xv(dim), yv(dim), sv(dim);
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
  mesh_new_->deform(nodeids, new_positions, false, &final_positions);

  // create a mesh map at time t
  Teuchos::ParameterList map_list;
  map_list.set<std::string>("method", "Lagrange serendipity")
          .set<int>("method order", 2)
          .set<std::string>("projector", "L2")
          .set<std::string>("map name", "VEM");
  
  WhetStone::MeshMapsFactory maps_factory;
  auto maps = maps_factory.Create(map_list, mesh_, mesh_new_);

  // calculate approximate velocities
  double dtfac = 1.0 / dt;
  for (int c = 0; c < ncells_wghost; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::VectorPolynomial v;
    std::vector<WhetStone::VectorPolynomial> vvf;

    if (high_order_velf_) { // use high-order face velocity 
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        ana_.VelocityTaylor(mesh_->face_centroid(f), t, v);
        vvf.push_back(v * dt);
        (*velf)[f] = v * (mesh_->face_normal(f) * weak_sign_);
      }
    } else {
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        maps->VelocityFace(f, v);
        vvf.push_back(v);
        (*velf)[f] = (v * mesh_->face_normal(f)) * dtfac * weak_sign_;
      }
    }

    maps->VelocityCell(c, vvf, (*velc)[c]);
    (*velc)[c] *= -dtfac * weak_sign_;

    if (! conservative_form_) (*divc)[c] = Divergence((*velc)[c]);
  }
}

/* *****************************************************************
* Change original definitions of velocities: level set algorithm
***************************************************************** */
template <class AnalyticDG>
void AdvectionFn<AnalyticDG>::ApproximateVelocity_LevelSet(
    const CompositeVector& u,
    const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& velc,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& velf,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& divc)
{
  int dim = mesh_->space_dimension();
  const Epetra_MultiVector& u_c = *u.ViewComponent("cell", true);
  int nk = u_c.NumVectors();

  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // cell-based velocity is constant for dGP1
  // we approximate it with a linear function for dGP2
  double norm, norm2, norm3;
  AmanziMesh::Entity_ID_List cells;
  AmanziGeometry::Point zero(dim);

  for (int c = 0; c < ncells_wghost; ++c) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    // polynomial representation of solution
    /*
    Amanzi::WhetStone::DenseVector data(nk);
    Amanzi::WhetStone::Polynomial poly(dim, order_); 

    for (int i = 0; i < nk; ++i) {
      data(i) = u_c[i][c];
    }
    poly.SetPolynomialCoefficients(data);
    poly.set_origin(xc);
    */

    norm2 = 0.0;
    for (int i = 0; i < dim; ++i) {
      norm2 += u_c[i + 1][c] * u_c[i + 1][c];
    }
    norm = std::pow(norm2, 0.5);

    (*velc)[c].resize(dim);
    for (int i = 0; i < dim; ++i) {
      (*velc)[c][i].Reshape(dim, order_ - 1);
      (*velc)[c][i](0, 0) = u_c[i + 1][c] * weak_sign_ / norm;
    }
    (*velc)[c].set_origin(xc);

    if (! conservative_form_) (*divc)[c] = Divergence((*velc)[c]);
  }

  for (int f = 0; f < nfaces_wghost; ++f) {
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    (*velf)[f].Reshape(dim, 0, true);

    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    WhetStone::VectorPolynomial tmp(dim, dim);
    tmp.set_origin(zero);

    for (int n = 0; n < ncells; ++n) {
      WhetStone::VectorPolynomial velc_tmp = (*velc)[cells[n]];
      velc_tmp.ChangeOrigin(zero);
      tmp -= velc_tmp;
    }
    
    tmp.Value(xf).Norm2(&norm);
    (*velf)[f] = tmp * (mesh_->face_normal(f) / norm);
  }
}

}  // namespace Amanzi


/* *****************************************************************
* This tests exactness of the transient advection scheme for
* dp/dt + div(v p) = f.
***************************************************************** */
bool inside(const Amanzi::AmanziGeometry::Point& p) {
  Amanzi::AmanziGeometry::Point c(0.5, 0.5);
  return (norm(p - c) < 0.06); 
}

template <class AnalyticDG>
void AdvectionTransient(std::string filename, int nx, int ny, double dt,
                        const Amanzi::Explicit_TI::method_t& rk_method,
                        bool conservative_form = true, 
                        std::string weak_form = "dual",
                        std::string face_velocity_method = "high order")                      
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

  std::string pk_name = "PK operator";
  if (weak_form == "primal") pk_name = "PK operator: primal";

  int order = plist.sublist(pk_name)
                   .sublist("flux operator").get<int>("method order");

  std::string problem = (conservative_form) ? ", conservative formulation" : "";
  if (MyPID == 0) std::cout << "\nTest: 2D dG transient advection problem, " << filename 
                            << ", order=" << order << problem 
                            << ", weak formulation=" << weak_form << std::endl;

  // create a mesh framework
  ParameterList region_list = plist.sublist("regions");
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
  plist.set<std::string>("file name", filename);
  plist.set<std::string>("face velocity method", face_velocity_method);

  auto dg = Teuchos::rcp(new WhetStone::DG_Modal(order, mesh, "regularized"));
  AdvectionFn<AnalyticDG> fn(plist, nx, dt, mesh, dg, conservative_form, weak_form);

  // create initial guess
  AnalyticDG ana(mesh, order);

  CompositeVector& rhs = *fn.op_flux->global_operator()->rhs();
  CompositeVector sol(rhs), sol_next(rhs);
  sol.PutScalar(0.0);

  Epetra_MultiVector& sol_c = *sol.ViewComponent("cell");
  ana.InitialGuess(*dg, sol_c, 0.0);

  int nstep(0);
  double t(0.0), tend(1.0), tprint(0.0);
  Explicit_TI::RK<CompositeVector> rk(fn, rk_method, sol);

  while(t < tend - dt/2) {
    rk.TimeStep(t, dt, sol, sol_next);

    sol = sol_next;

    // visualization
    if (MyPID == 0 && std::fabs(t - tprint) < dt/4) {
      tprint += 0.1; 
      printf("t=%9.6f  |p|=%10.8g\n", t, fn.l2norm);

      const Epetra_MultiVector& p = *sol.ViewComponent("cell");
      GMV::open_data_file(*mesh, (std::string)"operators.gmv");
      GMV::start_data();
      GMV::write_cell_data(p, 0, "solution");
      if (order > 0) {
        GMV::write_cell_data(p, 1, "gradx");
        GMV::write_cell_data(p, 2, "grady");
      }
      if (order > 1) {
        GMV::write_cell_data(p, 3, "hesxx");
        GMV::write_cell_data(p, 4, "hesxy");
        GMV::write_cell_data(p, 5, "hesyy");
      }
      GMV::close_data_file();
    }

    t += dt;
    nstep++;

    // modify solution at the origin
    if (face_velocity_method == "level set") {
      ana.InitialGuess(*dg, sol_c, t, inside);
    }
  }

  // compute solution error
  sol.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *sol.ViewComponent("cell", false);

  double pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean;
  ana.ComputeCellError(*dg, p, tend, pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean);

  if (MyPID == 0) {
    printf("nx=%3d (mean) L2(p)=%9.6g  Inf(p)=%9.6g\n", nx, pl2_mean, pinf_mean);
    printf("      (total) L2(p)=%9.6g  Inf(p)=%9.6g\n", pl2_err, pinf_err);
    CHECK(pl2_mean < 0.1 / nx);
  }
}


TEST(OPERATOR_ADVECTION_TRANSIENT_DG) {
  // AdvectionTransient<AnalyticDG07>("square", 50, 50, 0.001, Amanzi::Explicit_TI::tvd_3rd_order, false, "primal", "level set");

  AdvectionTransient<AnalyticDG06b>("square",  4,  4, 0.1, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("square",  4,  4, 0.1, Amanzi::Explicit_TI::tvd_3rd_order, false);
  AdvectionTransient<AnalyticDG06>("square",  4,  4, 0.1, Amanzi::Explicit_TI::tvd_3rd_order, false, "primal");

  /*
  AdvectionTransient<AnalyticDG06>("square",  20,  20, 0.01, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("square",  40,  40, 0.01 / 2, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("square",  80,  80, 0.01 / 4, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("square", 160, 160, 0.01 / 8, Amanzi::Explicit_TI::tvd_3rd_order);

  AdvectionTransient<AnalyticDG06>("test/triangular8.exo",    8, 0, 0.01, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("test/triangular16.exo",  16, 0, 0.01 / 2, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("test/triangular32.exo",  32, 0, 0.01 / 4, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("test/triangular64.exo",  64, 0, 0.01 / 8, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("test/triangular128.exo",128, 0, 0.01 / 16,Amanzi::Explicit_TI::tvd_3rd_order);

  double dT0 = 0.01;
  AdvectionTransient<AnalyticDG06>("test/median15x16.exo",   16, 0, dT0, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("test/median32x33.exo",   32, 0, dT0 / 2, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("test/median63x64.exo",   64, 0, dT0 / 4, Amanzi::Explicit_TI::tvd_3rd_order);
  AdvectionTransient<AnalyticDG06>("test/median127x128.exo",128, 0, dT0 / 8, Amanzi::Explicit_TI::tvd_3rd_order);
  */ 
}


