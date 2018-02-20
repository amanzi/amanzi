/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
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
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "LinearOperatorGMRES.hh"
#include "Tensor.hh"

// Operators
#include "Analytic00.hh"
//#include "Analytic02.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_DiffusionMFD.hh"
#include "PDE_DiffusionMFD_Tracer.hh"


using namespace Teuchos;
using namespace Amanzi;

void DefineSurfacePlane(int surf_id,
                        double *surf_parm,
                        Teuchos::Ptr<const AmanziMesh::Mesh> mesh,
                        const Teuchos::Ptr<CompositeVector>& surf_presence,
                        const Teuchos::Ptr<CompositeVector>& surface_param);

/* *****************************************************************
* This test diffusion solver with full tensor and source term.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_TRACER) {
  using namespace Teuchos;
  using namespace Amanzi;  
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: 3D elliptic solver, nodal tracer discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_3d.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, &comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  //RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 0.0,  1.0, 1.0, 1.0, 4, 4, 4, gm, true, true);
  RCP<const Mesh> mesh = meshfactory("test/GMV_F/dbls_10.exo", gm, true, true);

  // modify diffusion coefficient
  // -- since rho=mu=1.0, we do not need to scale the diffusion tensor.
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  int num_surfaces = 1;
  CompositeVectorSpace cvs1, cvs4;
  cvs1.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, num_surfaces);
  cvs4.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 4*num_surfaces);

  Teuchos::RCP<CompositeVector> surf_presence = Teuchos::rcp(new CompositeVector(cvs1));
  Teuchos::RCP<CompositeVector> surface_param = Teuchos::rcp(new CompositeVector(cvs4));

  double s_param[4] = {-1,  -1 , 0, 1};
  // double s_param[4] = {1.3,  1.2 , -5, 1.3};
  // double s_param[4] = {1.,  0.5 , 0.1, -1.3};
  //double s_param[4] = {0.8,  1.2 , -5, 0.3};
  //double s_param[4] = {0.1,  0 , 1, -0.24}; 
  DefineSurfacePlane(0, s_param, mesh.ptr(), surf_presence.ptr(), surface_param.ptr());

  Analytic003D ana(mesh, 0., 3., 1.);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K->push_back(Kc);
  }
  double rho(1.0), mu(1.0);
  AmanziGeometry::Point g(0.0, -1.0);

  // create boundary data
  Point xv(2);
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::NODE, SCHEMA_DOFS_SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int v = 0; v < nnodes_wghost; v++) {
    mesh->node_get_coordinates(v, &xv);
    if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 ||
        fabs(xv[1]) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6 ||
        fabs(xv[2]) < 1e-6 || fabs(xv[2] - 1.0) < 1e-6) {
      bc_model[v] = OPERATOR_BC_DIRICHLET;
      bc_value[v] = ana.pressure_exact(xv, 0.0);
    }
  }

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator nodal");
  Teuchos::RCP<PDE_DiffusionMFD_Tracer> op = Teuchos::rcp(new PDE_DiffusionMFD_Tracer(op_list, mesh));

  op->UpdateMatrices(Teuchos::null, Teuchos::null, surf_presence.ptr(), surface_param.ptr());
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // create source and add it to the operator
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("node", true);
  src.PutScalar(0.0);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    AmanziMesh::Entity_ID_List nodes;
    mesh->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int k = 0; k < nnodes; k++) {
      int v = nodes[k];
      src[0][v] += ana.source_exact(xc, 0.0) * volume / nnodes;
    }
  }
  source.GatherGhostedToMaster();

  
  // populate the diffusion operator
  op->Setup(K, Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null, surf_presence.ptr(), surface_param.ptr());

  // update the source term
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->UpdateRHS(source, true);

  // apply BCs (primary=true, eliminate=true) and assemble
  op->ApplyBCs(*surf_presence->ViewComponent("cell", true), true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // create preconditoner using the base operator class
  ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
  global_op->InitPreconditioner("ILU", slist);

  // Test SPD properties of the preconditioner.
  CompositeVector a(cvs), ha(cvs), b(cvs), hb(cvs);
  a.Random();
  b.Random();
  global_op->ApplyInverse(a, ha);
  global_op->ApplyInverse(b, hb);

  double ahb, bha, aha, bhb;
  a.Dot(hb, &ahb);
  b.Dot(ha, &bha);
  a.Dot(ha, &aha);
  b.Dot(hb, &bhb);

  if (MyPID == 0) {
    std::cout << "Preconditioner:\n"
              << "  Symmetry test: " << ahb << " = " << bha << std::endl;
    std::cout << "  Positivity test: " << aha << " " << bhb << std::endl;
  }
  CHECK_CLOSE(ahb, bha, 1e-12 * fabs(ahb));

  
  
  CHECK(aha > 0.0);
  CHECK(bhb > 0.0);

  // solve the problem
  ParameterList lop_list = plist.sublist("solvers")
                                .sublist("GMRES").sublist("gmres parameters");
  AmanziSolvers::LinearOperatorGMRES<Operator, CompositeVector, CompositeVectorSpace>
      solver(global_op, global_op);
  
  solver.Init(lop_list);

  CompositeVector rhs = *global_op->rhs();
  CompositeVector solution(rhs), solution_ex(rhs);
  Epetra_MultiVector& sol_ex = *solution_ex.ViewComponent("node", true);
  solution.PutScalar(0.0);
  for (int v=0; v<nnodes_wghost;v++){
    Point xv;
    mesh->node_get_coordinates(v, &xv);
    sol_ex[0][v] = ana.pressure_exact(xv,0.);
  }

  double r_norm, b_norm;
  //std::cout<<"rhs\n"<<*rhs.ViewComponent("node", true)<<"\n";

  b.PutScalar(0.);
  global_op->Apply(solution_ex, b);

  //std::cout<<"b\n"<<*b.ViewComponent("node", true)<<"\n";
  rhs.ViewComponent("node", true)->Norm2(&r_norm);
  b.ViewComponent("node", true)->Norm2(&b_norm);

  std::cout<<std::setprecision(10)<<"r "<<r_norm<<" b "<<b_norm<<"\n";

  //  exit(0);
  
  int ierr = solver.ApplyInverse(rhs, solution);

  double err = 0;
  int num_active = 0;
  Epetra_MultiVector& surf_v = *surf_presence->ViewComponent("cell", true);
  for (int c=0; c<ncells; c++){
    if (surf_v[0][c]>0){
      num_active++;
      AmanziMesh::Entity_ID_List nodes;
      mesh->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();
      for (int k = 0; k < nnodes; k++) {
        int v = nodes[k];
        // std::cout<<"cell "<<c<<" node "<<v<<" : solution exact "<<sol_ex[0][v]<<" "<<(*solution.ViewComponent("node", true))[0][v]<<"    "<<
        //   sol_ex[0][v]- (*solution.ViewComponent("node", true))[0][v] <<"\n";
        err += std::abs(sol_ex[0][v]- (*solution.ViewComponent("node", true))[0][v] );
      }
    }
  }

  err *= 1./num_active;

  op->OutputGMV_Surface(0, *solution.ViewComponent("node", true), "surface.gmv");   

  if (MyPID == 0) {
    std::cout << "Number of active cells "<<num_active<<"\n";
    std::cout << "pressure solver (pcg): ||r||=" << solver.residual() 
              << " itr=" << solver.num_itrs()
              << " code=" << solver.returned_code() 
              << " error= "<<err<< std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("node");
    GMV::open_data_file(*mesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_node_data(p, 0, "solution");
    GMV::close_data_file();
  }

  // CHECK(solver.num_itrs() < 10);

  // compute pressure error
  solution.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *solution.ViewComponent("node", false);

  double pnorm, pl2_err, pinf_err, hnorm, ph1_err;
  ana.ComputeNodeError(p, 0.0, pnorm, pl2_err, pinf_err, hnorm, ph1_err);

  if (MyPID == 0) {
    pl2_err /= pnorm;
    ph1_err /= hnorm;
    printf("L2(p)=%9.6f  H1(p)=%9.6f  itr=%3d\n", pl2_err, ph1_err, solver.num_itrs());

    // CHECK(pl2_err < 3e-3 && ph1_err < 2e-2);
    // CHECK(solver.num_itrs() < 10);    
  }
}


void DefineSurfacePlane(int surf_id,
                        double *surf_parm,
                        Teuchos::Ptr<const AmanziMesh::Mesh> mesh,
                        const Teuchos::Ptr<CompositeVector>& surf_presence,
                        const Teuchos::Ptr<CompositeVector>& surface_param){

  Epetra_MultiVector& surf_presence_vec = *surf_presence->ViewComponent("cell", true);
  Epetra_MultiVector& surface_param_vec = *surface_param->ViewComponent("cell", true);

  int ncells = surf_presence_vec.MyLength();
  AmanziMesh::Entity_ID_List cells, faces, edges, vertices;
  AmanziGeometry::Point p0(3), l0, l1, n(surf_parm[0], surf_parm[1], surf_parm[2]);


  for (int i=0;i<3;i++){
    if (std::abs(surf_parm[i]) > 1e-10){
      p0[i]= -surf_parm[3]/surf_parm[i];
      p0[(i+1)%3]=0.;
      p0[(i+2)%3]=0.;
      break;
    }
  }

  //p0[0]=1.3;  p0[1]=0.;   p0[2]=0.;
  // std::cout<<"p0 "<<p0<<"\n";
  // std::cout<<"n "<<n<<"\n";

  bool intersection;
  for (int c=0;c<ncells;c++){
    intersection = false;
    mesh->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    for (int i=0;i<nfaces;i++){
      int f = faces[i];
      mesh->face_get_nodes(f, &vertices);
      int nv = vertices.size();
      for (int j=0;j<nv;j++){
        mesh->node_get_coordinates(vertices[j], &l0);
        mesh->node_get_coordinates(vertices[(j+1)%nv], &l1);
        l1 -= l0;
        double val = l1*n;
        if (std::abs(val) > 1e-14){
          double d = ((p0-l0)*n)/val;
          // std::cout<<"d="<<d<<"\n";
          if (d>0 && d<1){
            intersection=true;
            break;
          }
        }
      }
      if (intersection) break;
    }
    if (intersection){
      // std::cout<<"Cell cnt: "<<mesh->cell_centroid(c)<<"\n";
      surf_presence_vec[surf_id][c] = surf_id+1;
      for (int i=0;i<4;i++) surface_param_vec[surf_id+i][c] = surf_parm[i];
    }
  }

  std::cout<<surf_presence_vec<<"\n";
  //std::cout<<surface_param_vec<<"\n";

  //exit(0);
}

// /* *****************************************************************
// * Exactness test for mixed diffusion solver.
// * NOTE. Mixed boundary condition requires to use mass matrix. We
// *       lump it which leads to a small error.
// * **************************************************************** */
// TEST(OPERATOR_DIFFUSION_NODAL_EXACTNESS) {
//   using namespace Teuchos;
//   using namespace Amanzi;
//   using namespace Amanzi::AmanziMesh;
//   using namespace Amanzi::AmanziGeometry;
//   using namespace Amanzi::Operators;

//   Epetra_MpiComm comm(MPI_COMM_WORLD);
//   int MyPID = comm.MyPID();
//   if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, exactness" 
//                             << " test for nodal discretization" << std::endl;

//   // read parameter list
//   std::string xmlFileName = "test/operator_diffusion.xml";
//   ParameterXMLFileReader xmlreader(xmlFileName);
//   ParameterList plist = xmlreader.getParameters();

//   // create an SIMPLE mesh framework
//   ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
//   Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

//   FrameworkPreference pref;
//   pref.clear();
//   pref.push_back(MSTK);
//   pref.push_back(STKMESH);

//   MeshFactory meshfactory(&comm);
//   meshfactory.preference(pref);
//   // RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 4, 4, gm);
//   RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);

//   // modify diffusion coefficient
//   // -- since rho=mu=1.0, we do not need to scale the diffusion coefficient.
//   Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
//   int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
//   int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
//   int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

//   Analytic02 ana(mesh);

//   for (int c = 0; c < ncells_wghost; c++) {
//     const Point& xc = mesh->cell_centroid(c);
//     const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
//     K->push_back(Kc);
//   }
//   double rho(1.0), mu(1.0);
//   AmanziGeometry::Point g(0.0, -1.0);

//   // create boundary data (no mixed bc)
//   Point xv(2);
//   Teuchos::RCP<BCs> bc_v = Teuchos::rcp(new BCs(mesh, AmanziMesh::NODE, SCHEMA_DOFS_SCALAR));
//   std::vector<int>& bc_model_v = bc_v->bc_model();
//   std::vector<double>& bc_value_v = bc_v->bc_value();

//   for (int v = 0; v < nnodes_wghost; v++) {
//     mesh->node_get_coordinates(v, &xv);
//     if(fabs(xv[0] - 1.0) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
//       bc_model_v[v] = OPERATOR_BC_DIRICHLET;
//       bc_value_v[v] = ana.pressure_exact(xv, 0.0);
//     }
//   }

//   Teuchos::RCP<BCs> bc_f = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, SCHEMA_DOFS_SCALAR));
//   std::vector<int>& bc_model_f = bc_f->bc_model();
//   std::vector<double>& bc_value_f = bc_f->bc_value();
//   std::vector<double>& bc_mixed_f = bc_f->bc_mixed();

//   int nn=0; int nm=0;
//   for (int f = 0; f < nfaces_wghost; f++) {
//     const Point& xf = mesh->face_centroid(f);
//     if (fabs(xf[0]) < 1e-6) {
//       nn++;
//       bc_model_f[f] = OPERATOR_BC_NEUMANN;
//       bc_value_f[f] = -(ana.velocity_exact(xf, 0.0))[0];  // We assume exterior normal.
//     } else if (fabs(xf[1]) < 1e-6) {
//       nm++;
//       bc_model_f[f] = OPERATOR_BC_MIXED;
//       bc_value_f[f] = -(ana.velocity_exact(xf, 0.0))[1];  // We assume exterior normal.

//       double tmp = ana.pressure_exact(xf, 0.0);
//       bc_mixed_f[f] = 1.0;
//       bc_value_f[f] -= bc_mixed_f[f] * tmp;
//     }
//   }

//   // create diffusion operator 
//   ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator nodal");
//   PDE_DiffusionFactory opfactory;
//   Teuchos::RCP<PDE_Diffusion> op = opfactory.Create(op_list, mesh, bc_v, rho, g);
//   op->AddBCs(bc_f, bc_f);
  
//   // populate the diffusion operator
//   Teuchos::RCP<Operator> global_op = op->global_operator();
//   global_op->Init();
//   op->Setup(K, Teuchos::null, Teuchos::null);
//   op->UpdateMatrices(Teuchos::null, Teuchos::null);
//   op->ApplyBCs(true, true);

//   global_op->SymbolicAssembleMatrix();
//   global_op->AssembleMatrix();

//   // create preconditoner using the base operator class
//   ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
//   global_op->InitPreconditioner("Hypre AMG", slist);

//   // solve the problem
//   ParameterList lop_list = plist.sublist("solvers")
//                                 .sublist("AztecOO CG").sublist("pcg parameters");
//   AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
//       solver(global_op, global_op);
//   solver.Init(lop_list);

//   CompositeVector rhs = *global_op->rhs();
//   CompositeVector solution(rhs);
//   solution.PutScalar(0.0);

//   int ierr = solver.ApplyInverse(rhs, solution);

//   if (MyPID == 0) {
//     std::cout << "pressure solver (pcg): ||r||=" << solver.residual() 
//               << " itr=" << solver.num_itrs()
//               << " code=" << solver.returned_code() << std::endl;
//   }

//   // compute pressure error
//   solution.ScatterMasterToGhosted();
//   Epetra_MultiVector& p = *solution.ViewComponent("node", false);

//   double pnorm, pl2_err, pinf_err, hnorm, ph1_err;
//   ana.ComputeNodeError(p, 0.0, pnorm, pl2_err, pinf_err, hnorm, ph1_err);

//   if (MyPID == 0) {
//     pl2_err /= pnorm;
//     ph1_err /= hnorm;
//     printf("L2(p)=%9.6f  H1(p)=%9.6f  itr=%3d\n", pl2_err, ph1_err, solver.num_itrs());

//     CHECK(pl2_err < 1e-5 && ph1_err < 2e-5);
//     CHECK(solver.num_itrs() < 10);
//   }
// }

