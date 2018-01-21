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
#include "EpetraExt_RowMatrixOut.h"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "Tensor.hh"

// Operators
#include "AnalyticMultiMat.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionMFD.hh"
#include "PDE_DiffusionMFD_XMOF.hh"
// #include "UpwindSecondOrder.hh"


int BoundaryFaceGetCell(const Amanzi::AmanziMesh::Mesh& mesh, int f);
// {
//   Amanzi::AmanziMesh::Entity_ID_List cells;
//   mesh.face_get_cells(f, Amanzi::AmanziMesh::USED, &cells);
//   return cells[0];
// }

XMOF2D::CellsMatData read_mat_data(const std::string& matdata_fname);


/* *****************************************************************
* Exactness test for mixed diffusion solver.
***************************************************************** */
void RunTestDiffusionMixedXMOF(double gravity) {
  using namespace Teuchos; 
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, exactness test" 
                            << " for mixed XMOF discretization" << std::endl;

  // read parameter list
  // -- it specifies details of the mesh, diffusion operator, and solver
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create the MSTK mesh framework 
  // -- geometric model is defined in the region sublist of XML list
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

  // -- provide at lest one framework to the mesh factory. The first available
  // -- framework will be used
  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({MSTK, STKMESH}));
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 8, 8, gm);
  // RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  AnalyticMultiMat ana(mesh, gravity);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, SCHEMA_DOFS_SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();
  std::vector<double>& bc_mixed = bc->bc_mixed();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    double area = mesh->face_area(f);
    int dir, c = BoundaryFaceGetCell(*mesh, f);
    const Point& normal = mesh->face_normal(f, false, c, &dir);

    if (fabs(xf[0]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;
    } else if(fabs(xf[1]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_MIXED;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;

      double tmp = ana.pressure_exact(xf, 0.0);
      bc_mixed[f] = 1.0;
      bc_value[f] -= bc_mixed[f] * tmp;
    } else if(fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }

  std::string matdata_fname =  "test/line_8x8_matdata.dat";
  XMOF2D::CellsMatData mat_data = read_mat_data(matdata_fname);
   
  int num_mat = 2;
  CompositeVectorSpace cvs1, cvs2;
  cvs1.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, num_mat);
  cvs2.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 2*num_mat);

  Teuchos::RCP<CompositeVector> vol_frac = Teuchos::rcp(new CompositeVector(cvs1)); 
  Epetra_MultiVector& vol_frac_vec = *vol_frac->ViewComponent("cell", true);
  Teuchos::RCP<CompositeVector> centroids = Teuchos::rcp(new CompositeVector(cvs2));
  Epetra_MultiVector& centroids_vec = *centroids->ViewComponent("cell", true);

  
  for (int i=0; i!= mat_data.cells_materials.size(); ++i){
    // std::cout<<i<<": ";
    // for (int j=0; j!= mat_data.cells_materials[i].size(); ++j) 
    //   std::cout<<mat_data.cells_materials[i][j]<<" ";
    // std::cout<<"\n";
    // std::cout<<i<<": ";
    // for (int j=0; j!= mat_data.cells_materials[i].size(); ++j) 
    //   std::cout<<mat_data.cells_vfracs[i][j]<<" ";
    // std::cout<<"\n";
    for (int j=0; j!= mat_data.cells_materials[i].size(); ++j){
      int mat_id = int(mat_data.cells_materials[i][j]);
      vol_frac_vec[mat_id][i] = mat_data.cells_vfracs[i][j];
      centroids_vec[mat_id*2][i] = mat_data.cells_centroids[i][j].x;
      centroids_vec[mat_id*2+1][i] = mat_data.cells_centroids[i][j].y;
    }      

  }
  std::cout<<"\n\n";
  // create diffusion operator 
  double rho(1.0);
  AmanziGeometry::Point g(0.0, -gravity);
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator mixed");
  Teuchos::RCP<PDE_DiffusionMFD_XMOF> op = Teuchos::rcp(new PDE_DiffusionMFD_XMOF(op_list, mesh));
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // set up the diffusion operator
  op->Setup(K, Teuchos::null, Teuchos::null);
  op->ConstructMiniMesh(vol_frac.ptr(), centroids.ptr());
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // // get and assmeble the global operator
  // Teuchos::RCP<Operator> global_op = op->global_operator();
  // op->ApplyBCs(true, true);
  // global_op->SymbolicAssembleMatrix();
  // global_op->AssembleMatrix();

  // // create preconditoner using the base operator class
  // ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
  // global_op->InitPreconditioner("Hypre AMG", slist);

  // // Test SPD properties of the preconditioner.
  // CompositeVector a(cvs), ha(cvs), b(cvs), hb(cvs);
  // a.Random();
  // b.Random();
  // global_op->ApplyInverse(a, ha);
  // global_op->ApplyInverse(b, hb);

  // double ahb, bha, aha, bhb;
  // a.Dot(hb, &ahb);
  // b.Dot(ha, &bha);
  // a.Dot(ha, &aha);
  // b.Dot(hb, &bhb);

  // if (MyPID == 0) {
  //   std::cout << "Preconditioner:\n"
  //             << "  Symmetry test: " << ahb << " = " << bha << std::endl;
  //   std::cout << "  Positivity test: " << aha << " " << bhb << std::endl;
  // }
  // CHECK_CLOSE(ahb, bha, 1e-12 * fabs(ahb));
  // CHECK(aha > 0.0);
  // CHECK(bhb > 0.0);

  // // solve the problem
  // ParameterList lop_list = plist.sublist("solvers")
  //                               .sublist("AztecOO CG").sublist("pcg parameters");
  // AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
  //     solver(global_op, global_op);
  // solver.Init(lop_list);

  // CompositeVector rhs = *global_op->rhs();
  // Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(rhs));
  // Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(rhs));
  // solution->PutScalar(0.0);

  // int ierr = solver.ApplyInverse(rhs, *solution);

  // if (MyPID == 0) {
  //   std::cout << "pressure solver (pcg): ||r||=" << solver.residual() 
  //             << " itr=" << solver.num_itrs()
  //             << " code=" << solver.returned_code() << std::endl;
  // }

  // // compute pressure error
  // Epetra_MultiVector& p = *solution->ViewComponent("cell", false);
  // double pnorm, pl2_err, pinf_err;
  // ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  // // calculate flux error
  // Epetra_MultiVector& flx = *flux->ViewComponent("face", true);
  // double unorm, ul2_err, uinf_err;

  // op->UpdateFlux(solution.ptr(), flux.ptr());
  // ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  // if (MyPID == 0) {
  //   pl2_err /= pnorm; 
  //   ul2_err /= unorm;
  //   printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f  itr=%3d\n",
  //       pl2_err, pinf_err, ul2_err, uinf_err, solver.num_itrs());

  //   CHECK(pl2_err < 1e-12 && ul2_err < 1e-12);
  //   CHECK(solver.num_itrs() < 10);
  // }
}


TEST(OPERATOR_DIFFUSION_MIXED_XMOF) {
  RunTestDiffusionMixedXMOF(0.0);
}

// TEST(OPERATOR_DIFFUSION_MIXED_wGRAVITY) {
//   RunTestDiffusionMixed(0.1);
// }


// /* *****************************************************************
// * Exactness test for cell-based diffusion solver.
// ***************************************************************** */
// TEST(OPERATOR_DIFFUSION_CELL_EXACTNESS) {
//   using namespace Teuchos;
//   using namespace Amanzi;
//   using namespace Amanzi::AmanziMesh;
//   using namespace Amanzi::AmanziGeometry;
//   using namespace Amanzi::Operators;

//   Epetra_MpiComm comm(MPI_COMM_WORLD);
//   int MyPID = comm.MyPID();

//   if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, exactness" 
//                             << " test for cell-based discretization" << std::endl;

//   // read parameter list
//   std::string xmlFileName = "test/operator_diffusion.xml";
//   ParameterXMLFileReader xmlreader(xmlFileName);
//   ParameterList plist = xmlreader.getParameters();

//   // create an SIMPLE mesh framework
//   ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
//   Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

//   MeshFactory meshfactory(&comm);
//   meshfactory.preference(FrameworkPreference({MSTK, STKMESH}));
//   RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 15, 8, gm);

//   // modify diffusion coefficient
//   Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
//   int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
//   int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

//   Analytic02 ana(mesh);

//   for (int c = 0; c < ncells; c++) {
//     const Point& xc = mesh->cell_centroid(c);
//     WhetStone::Tensor Kc(2, 2);

//     Kc(0, 0) = 3.0;
//     Kc(1, 1) = 1.0;
//     Kc(0, 1) = 0.0;
//     Kc(1, 0) = 0.0;

//     K->push_back(Kc);
//   }
//   AmanziGeometry::Point g(0.0, -1.0);

//   // create boundary data.
//   Teuchos::RCP<BCs> bc_f = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, SCHEMA_DOFS_SCALAR));
//   std::vector<int>& bc_model = bc_f->bc_model();
//   std::vector<double>& bc_value = bc_f->bc_value();

//   for (int f = 0; f < nfaces_wghost; f++) {
//     const Point& xf = mesh->face_centroid(f);
//     double area = mesh->face_area(f);

//     if (fabs(xf[0]) < 1e-6) {
//       bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
//       bc_value[f] = 3.0;
//     /*
//     } else if(fabs(xf[1]) < 1e-6) {
//       bc_model[f] = Operators::OPERATOR_BC_MIXED;
//       bc_value[f] = 2.0;

//       double tmp = ana.pressure_exact(xf, 0.0);
//       bc_mixed[f] = 1.0;
//       bc_value[f] -= bc_mixed[f] * tmp;
//     */
//     } else if(fabs(xf[0] - 1.0) < 1e-6 || 
//               fabs(xf[1] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6) {
//       bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
//       bc_value[f] = ana.pressure_exact(xf, 0.0);
//     }
//   }

//   // create diffusion operator 
//   ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator cell");
//   Teuchos::RCP<PDE_Diffusion> op = Teuchos::rcp(new PDE_DiffusionFV(op_list, mesh));
//   op->SetBCs(bc_f, bc_f);
//   const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

//   // set up the diffusion operator
//   op->Setup(K, Teuchos::null, Teuchos::null);
//   op->UpdateMatrices(Teuchos::null, Teuchos::null);

//   // get and assmeble the global operator
//   Teuchos::RCP<Operator> global_op = op->global_operator();
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
//   Epetra_MultiVector& p = *solution.ViewComponent("cell", false);
//   double pnorm, pl2_err, pinf_err;
//   ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

//   if (MyPID == 0) {
//     pl2_err /= pnorm; 
//     printf("L2(p)=%9.6f  Inf(p)=%9.6f  itr=%3d\n", pl2_err, pinf_err, solver.num_itrs());

//     CHECK(pl2_err < 1e-5);
//     CHECK(solver.num_itrs() < 10);
//   }
// }

XMOF2D::CellsMatData read_mat_data(const std::string& matdata_fname) {
  std::ifstream os(matdata_fname.c_str(), std::ifstream::binary);
  if (!os.good()) {
    std::ostringstream os;
    os << std::endl <<
      "Cannot open " << matdata_fname << " for binary input" << std::endl;
    throw XMOF2D::Exception(os.str());
  }
  
  XMOF2D::CellsMatData mat_data;
  
  int data_dim;
  os.read(reinterpret_cast<char *>(&data_dim), sizeof(int));
  if(data_dim != 2) {
    std::ostringstream os;
    os << std::endl << "Material data should be for a 2D mesh!" << std::endl;
    throw XMOF2D::Exception(os.str());
  }

  int ncells;
  os.read(reinterpret_cast<char *>(&ncells), sizeof(int));

  mat_data.cells_materials.resize(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    int nmats;
    os.read(reinterpret_cast<char *>(&nmats), sizeof(int));
    mat_data.cells_materials[icell].resize(nmats);
    for (int im = 0; im < nmats; im++)
      os.read(reinterpret_cast<char *>(&mat_data.cells_materials[icell][im]), sizeof(int));
  }
  mat_data.cells_vfracs.resize(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    int nmats = mat_data.cells_materials[icell].size();
    mat_data.cells_vfracs[icell].resize(nmats);
    if (nmats == 1) {
      mat_data.cells_vfracs[icell][0] = 1.0;
      continue;
    }
    for (int im = 0; im < nmats; im++)
      os.read(reinterpret_cast<char *>(&mat_data.cells_vfracs[icell][im]), sizeof(double));
  }
  mat_data.cells_centroids.resize(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    int nmats = mat_data.cells_materials[icell].size();
    mat_data.cells_centroids[icell].resize(nmats);
    if (nmats == 1) {
      mat_data.cells_centroids[icell][0] = XMOF2D::BAD_POINT;
      continue;
    }
    for (int im = 0; im < nmats; im++) {
      double cen_x, cen_y;
      os.read(reinterpret_cast<char *>(&cen_x), sizeof(double));
      os.read(reinterpret_cast<char *>(&cen_y), sizeof(double));
      mat_data.cells_centroids[icell][im] = XMOF2D::Point2D(cen_x, cen_y);
    }
  }

  os.close();
  
  return mat_data;
}
