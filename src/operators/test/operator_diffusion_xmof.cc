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
#include "RegionBoxVolumeFractions.hh"
#include "LinearOperatorPCG.hh"
#include "Tensor.hh"

// Operators

//#include "Analytic00.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionMFD.hh"
#include "PDE_DiffusionMFD_XMOF.hh"
// #include "UpwindSecondOrder.hh"

//Analytic
#include "AnalyticMultiMat00.hh"
#include "AnalyticMultiMat01.hh"

using namespace Teuchos; 
using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Operators;

void write_data_example_simple(const Amanzi::AmanziMesh::Mesh& mesh, const std::string& matdata_fname);
void write_data_example_poly(const Amanzi::AmanziMesh::Mesh& mesh, const std::string& matdata_fname);
void SplitPolyCell(const std::vector<Amanzi::AmanziGeometry::Point>& xy1, 
                   const Amanzi::AmanziGeometry::Point& pnt, 
                   const Amanzi::AmanziGeometry::Point& vec, 
                   std::vector<Amanzi::AmanziGeometry::Point>& xy3);
void Intersection_Centroid_Area(const std::vector<Amanzi::AmanziGeometry::Point>& xy1,
                                Amanzi::AmanziGeometry::Point& cntr, double *inter_area);
void TestSPD(Teuchos::RCP<Operator> global_op);

int BoundaryFaceGetCell(const Amanzi::AmanziMesh::Mesh& mesh, int f)
{
  Amanzi::AmanziMesh::Entity_ID_List cells;
  mesh.face_get_cells(f, Amanzi::AmanziMesh::USED, &cells);
  return cells[0];
}

XMOF2D::CellsMatData read_mat_data(const std::string& matdata_fname);


/* *****************************************************************
* Exactness test for mixed diffusion solver.
***************************************************************** */
void RunTestDiffusionMixedXMOF(double gravity) {

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, transient sandwich" 
                            << " for mixed XMOF discretization" << std::endl;

  // read parameter list
  // -- it specifies details of the mesh, diffusion operator, and solver
  std::string xmlFileName = "test/operator_diffusion_xmof.xml";
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
  // generate orthogonal mesh
  //RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 8, 8, gm);
  // read mesh from file
  RCP<const Mesh> mesh = meshfactory("test/median15x16.exo", gm);

  // modify diffusion coefficient for multimaterial diffusion
  Teuchos::RCP<std::vector<std::vector<WhetStone::Tensor> > >KMulti = Teuchos::rcp(new std::vector<std::vector<WhetStone::Tensor> >());

  //number of cells in the mesh
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  // number of faces in the mesh
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  // number of faces with ghost faces 
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  
  //Analytical model
  AnalyticMultiMat01 ana(mesh, 0., -1. ,1);

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, SCHEMA_DOFS_SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();
  std::vector<double>& bc_mixed = bc->bc_mixed();

  //fill boundary data
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    double area = mesh->face_area(f);
    int dir, c = BoundaryFaceGetCell(*mesh, f);
    const Point& normal = mesh->face_normal(f, false, c, &dir);

    if ((fabs(xf[1]) < 1e-6)){
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = 1.;
    }else if (fabs(xf[1]) > 1 - 1e-6){
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = 0;
    }
  }

  //Create multimaterial data
  std::string matdata_fname2 =  "test/sandwich_15x16_matdata.dat";
  write_data_example_poly(*mesh, matdata_fname2);
  //Read multimaterial data  
  XMOF2D::CellsMatData mat_data = read_mat_data(matdata_fname2);
   
  int num_mat = 2;
  CompositeVectorSpace cvs1, cvs2, cvs_pl;
  // composite vector space for volume fractions
  cvs1.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, num_mat);
  // composite vector space for centroids  
  cvs2.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 2*num_mat);
  // composite vector space for solution  
  cvs_pl.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1 + num_mat);
  cvs_pl.AddComponent("face", AmanziMesh::FACE, 1);
  

  Teuchos::RCP<CompositeVector> vol_frac = Teuchos::rcp(new CompositeVector(cvs1)); 
  Epetra_MultiVector& vol_frac_vec = *vol_frac->ViewComponent("cell", true);
  Teuchos::RCP<CompositeVector> centroids = Teuchos::rcp(new CompositeVector(cvs2));
  Epetra_MultiVector& centroids_vec = *centroids->ViewComponent("cell", true);

  // Fill data for volume fraction, centroids, 
  for (int i=0; i!= mat_data.cells_materials.size(); ++i){   
    const Point& xc = mesh->cell_centroid(i);
    std::vector<WhetStone::Tensor> Ks_cell;

    if (mat_data.cells_materials[i].size() > 1) {  //mixed cell
      Ks_cell.resize(num_mat);
      
      for (int j=0; j!= mat_data.cells_materials[i].size(); ++j){
        int mat_id = int(mat_data.cells_materials[i][j]);
        const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0, mat_id);
        Ks_cell[mat_id] = Kc;
        vol_frac_vec[mat_id][i] = mat_data.cells_vfracs[i][j];// + 0.2*sin(1.*j);
        centroids_vec[mat_id*2][i] = mat_data.cells_centroids[i][j].x;
        centroids_vec[mat_id*2+1][i] = mat_data.cells_centroids[i][j].y;
      }
    }else{
      int mat_id = int(mat_data.cells_materials[i][0]);
      const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0, mat_id);
      Ks_cell.push_back(Kc);
      vol_frac_vec[mat_id][i] = 1.0;
    }      
    KMulti->push_back(Ks_cell);    
  }

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator mixed xmof");
  Teuchos::RCP<PDE_DiffusionMFD_XMOF> op = Teuchos::rcp(new PDE_DiffusionMFD_XMOF(op_list, mesh));
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // set up the diffusion tensor for diffusion operator
  op->SetupMultiMatK(KMulti);
  // Create minimesh for each cell
  op->ConstructMiniMesh(vol_frac.ptr(), centroids.ptr());
  Teuchos::RCP<Operator> global_op = op->global_operator();

  //Create local matrices
  op->UpdateMatrices(Teuchos::null, Teuchos::null, 0.);

  // get and assmeble the global operator 
  op->ApplyBCs(true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();
 

  // create preconditoner using the base operator class
  ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
  global_op->InitPreconditioner("Hypre AMG", slist);

  //Test SPD properties of the preconditioner.
  //TestSPD(global_op);

  // solve the problem
  ParameterList lop_list = plist.sublist("solvers")
                                .sublist("AztecOO CG").sublist("pcg parameters");
  AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
      solver(global_op, global_op);
  solver.Init(lop_list);

  //RHS vector
  CompositeVector rhs = *global_op->rhs();
  // Solution vector lambdas
  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(rhs));
  // Flux solution 
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(rhs));
  solution->PutScalar(0.0);

  //
  Teuchos::RCP<CompositeVector> sol_pl = Teuchos::rcp(new CompositeVector(cvs_pl));
  *sol_pl->ViewComponent("face", false) = *solution->ViewComponent("face", false);
  op->UpdateFlux(sol_pl.ptr(), flux.ptr(), 0);  

  
  
  double dt = 1e-2;
  double time_total = 0.;
  int time_step = 0;

  Epetra_MultiVector& lmd = *sol_pl->ViewComponent("face", false);
  Epetra_MultiVector& p = *sol_pl->ViewComponent("cell", false);
  //const Epetra_MultiVector& vol_frac_vec = *vol_frac->ViewComponent("cell", true);
  
  while (time_total < 1e+1){
    op->UpdateMatrices(Teuchos::null, sol_pl.ptr(), dt);
    // get and assmeble the global operator 
    op->ApplyBCs(true, true);
    global_op->SymbolicAssembleMatrix();
    global_op->AssembleMatrix();
  
    int ierr = solver.ApplyInverse(rhs, *solution);
       
    *sol_pl->ViewComponent("face", false) = *solution->ViewComponent("face", false);
    op->UpdateFlux(sol_pl.ptr(), flux.ptr(), dt);

    // // compute pressure error
    // double pmin(0.), pmax(0.);
    // p.MinValue(&pmin);
    // p.MaxValue(&pmax);
    // std::cout<<" Min= "<<pmin<<" Max= "<<pmax<<"\n";

    time_total += dt;

    if (time_step == 50) dt *=10.;
    
    time_step++;

    if ((MyPID == 0)&&( (time_step-1)%10==0)) {
      
      std::string output_file, fileNum;	
      std::ostringstream ss;
      ss<<time_step;
      if (time_step<10) fileNum = "000" + ss.str();
      else if (time_step<100) fileNum = "00" + ss.str();
      else if (time_step<1000) fileNum = "0" + ss.str();
      else fileNum = ss.str();
      output_file ="gmv/operators_xmof.gmv" + fileNum ;
      
      // GMV::open_data_file(*mesh, output_file );
      // GMV::start_data();
      // GMV::write_cell_data(p, 0, "solution");
      // GMV::close_data_file();
      op->WriteSpecialGMV(output_file, vol_frac_vec, *sol_pl->ViewComponent("cell", false));
    }
      
    //op->WriteSpecialGMV((std::string)"special.gmv", *sol_pl->ViewComponent("cell", false));
    //if (time_step == 100)  
    //break;
  }

 
  // double pnorm(0.), pl2_err(0.), pinf_err(0.);
  // double lnorm(0.), ll2_err(0.), linf_err(0.);
  // ana.ComputeLambdaError(lmd, 0.0, lnorm, ll2_err, linf_err);
  // ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);


}


void RunTestDiffusionMixedXMOF_Linear() {
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
  std::string xmlFileName = "test/operator_diffusion_xmof.xml";
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

  // modify diffusion coefficient for multimaterial diffusion
  Teuchos::RCP<std::vector<std::vector<WhetStone::Tensor> > >KMulti = Teuchos::rcp(new std::vector<std::vector<WhetStone::Tensor> >());
  //number of cells in the mesh
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  // number of faces in the mesh
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  // number of faces with ghost faces
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  AnalyticMultiMat00 ana(mesh, -1., 1.);
  //AnalyticMultiMat01 ana(mesh, -1., 1.);

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

    if ((fabs(xf[0]) < 1e-6)||
        (fabs(xf[1]) < 1e-6)||
        (fabs(xf[1]) > 1 - 1e-6)||
        (fabs(xf[0]) > 1 - 1e-6)){
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  
  }

  std::string matdata_fname = "test/simple_matdata.dat";
  write_data_example_simple(*mesh, matdata_fname);
  // std::string matdata_fname1 =  "test/line_4x4_matdata.dat";
  // std::string matdata_fname2 =  "test/sandwich_15x16_matdata.dat";
  // write_data_example_poly(*mesh, matdata_fname2); 
  XMOF2D::CellsMatData mat_data = read_mat_data(matdata_fname);
   
  int num_mat = 2;
  CompositeVectorSpace cvs1, cvs2, cvs_pl;
  cvs1.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, num_mat);
  cvs2.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 2*num_mat);
  cvs_pl.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1 + num_mat);
  cvs_pl.AddComponent("face", AmanziMesh::FACE, 1);
  

  Teuchos::RCP<CompositeVector> vol_frac = Teuchos::rcp(new CompositeVector(cvs1)); 
  Epetra_MultiVector& vol_frac_vec = *vol_frac->ViewComponent("cell", true);
  Teuchos::RCP<CompositeVector> centroids = Teuchos::rcp(new CompositeVector(cvs2));
  Epetra_MultiVector& centroids_vec = *centroids->ViewComponent("cell", true);

  
  for (int i=0; i!= mat_data.cells_materials.size(); ++i){
   
    const Point& xc = mesh->cell_centroid(i);
    std::vector<WhetStone::Tensor> Ks_cell;

    if (mat_data.cells_materials[i].size() > 1) {  //mixed cell
      Ks_cell.resize(num_mat);
      
      for (int j=0; j!= mat_data.cells_materials[i].size(); ++j){
        int mat_id = int(mat_data.cells_materials[i][j]);
        const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0, mat_id);
        Ks_cell[mat_id] = Kc;
        vol_frac_vec[mat_id][i] = mat_data.cells_vfracs[i][j];
        centroids_vec[mat_id*2][i] = mat_data.cells_centroids[i][j].x;
        centroids_vec[mat_id*2+1][i] = mat_data.cells_centroids[i][j].y;
      }
    }else{
      int mat_id = int(mat_data.cells_materials[i][0]);
      const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0, mat_id);
      Ks_cell.push_back(Kc);
      vol_frac_vec[mat_id][i] = 1.0;
    }      
    KMulti->push_back(Ks_cell);       
  }

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator mixed xmof");
  Teuchos::RCP<PDE_DiffusionMFD_XMOF> op = Teuchos::rcp(new PDE_DiffusionMFD_XMOF(op_list, mesh));
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // set up the diffusion operator
  op->SetupMultiMatK(KMulti);
  op->ConstructMiniMesh(vol_frac.ptr(), centroids.ptr());
  Teuchos::RCP<Operator> global_op = op->global_operator();

  op->UpdateMatrices(Teuchos::null, Teuchos::null, 0.);
  // get and assmeble the global operator 
  op->ApplyBCs(true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // create preconditoner using the base operator class
  ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
  global_op->InitPreconditioner("Hypre AMG", slist);

  TestSPD(global_op);

  // solve the problem
  ParameterList lop_list = plist.sublist("solvers")
                                .sublist("AztecOO CG").sublist("pcg parameters");
  AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
      solver(global_op, global_op);
  solver.Init(lop_list);

  CompositeVector rhs = *global_op->rhs();
  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(rhs));
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(rhs));
  solution->PutScalar(0.0);
  Teuchos::RCP<CompositeVector> sol_pl = Teuchos::rcp(new CompositeVector(cvs_pl));
  *sol_pl->ViewComponent("face", false) = *solution->ViewComponent("face", false);
  op->UpdateFlux(sol_pl.ptr(), flux.ptr(), 0.);  

  int ierr = solver.ApplyInverse(rhs, *solution);
  *sol_pl->ViewComponent("face", false) = *solution->ViewComponent("face", false);
  op->UpdateFlux(sol_pl.ptr(), flux.ptr(), 0.);


  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << solver.residual() 
              << " itr=" << solver.num_itrs()
              << " code=" << solver.returned_code() << std::endl;
  }

  

  // compute pressure error
  Epetra_MultiVector& lmd = *sol_pl->ViewComponent("face", false);
  Epetra_MultiVector& p = *sol_pl->ViewComponent("cell", false);
  //const Epetra_MultiVector& vol_frac_vec = *vol_frac->ViewComponent("cell", true);
  double pnorm, pl2_err, pinf_err;
  double lnorm, ll2_err, linf_err;
  ana.ComputeLambdaError(lmd, 0.0, lnorm, ll2_err, linf_err);
  
  // // calculate flux error
  // Epetra_MultiVector& flx = *flux->ViewComponent("face", true);
  // double unorm, ul2_err, uinf_err;

  // op->UpdateFlux(solution.ptr(), flux.ptr());
  // ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm; 
    // ul2_err /= unorm;
    // printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f  itr=%3d\n",
    //     pl2_err, pinf_err, ul2_err, uinf_err, solver.num_itrs());
    // printf("Pressure L2(p)=%9.6f  Inf(p)=%9.6f  itr=%3d\n",
    //     pl2_err, pinf_err, solver.num_itrs());    

    ll2_err /= lnorm;
    printf("Lambda L2(p)=%15.9f  Inf(p)=%15.9f  itr=%3d\n",
        ll2_err, linf_err, solver.num_itrs());    
    CHECK(ll2_err < 1e-6 && linf_err < 1e-6);
    CHECK(solver.num_itrs() < 10);
    // GMV::open_data_file(*mesh, (std::string)"operators_xmof.gmv");
    // GMV::start_data();
    // GMV::write_cell_data(p, 0, "solution");
    // GMV::close_data_file();
    op->WriteSpecialGMV( (std::string)"linear_xmof.gmv", vol_frac_vec, *sol_pl->ViewComponent("cell", false));
    
  }
}


TEST(OPERATOR_DIFFUSION_MIXED_XMOF_LINEAR) {
  RunTestDiffusionMixedXMOF_Linear();
}

TEST(OPERATOR_DIFFUSION_MIXED_XMOF) {
  RunTestDiffusionMixedXMOF(0.0);
}




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
      //std::cout<<cen_x<<" "<< cen_y<<"\n";
      mat_data.cells_centroids[icell][im] = XMOF2D::Point2D(cen_x, cen_y);
    }
  }

  os.close();
  
  return mat_data;
}


void write_data_example_simple(const Amanzi::AmanziMesh::Mesh& mesh, const std::string& matdata_fname){

  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  

  std::ofstream os(matdata_fname.c_str(), std::ofstream::binary);

    // if (!os.good())
    //   THROW_EXCEPTION("Cannot open " << cfg.bin_mesh_data_fname << " for binary output");

    int out_int = 2;  //Dimension: 2 for 2D, 3 for 3D

    int ncells = mesh.num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

    os.write((char*) &out_int, sizeof(int));
    os.write((char*) &ncells, sizeof(int));
    AmanziGeometry::Point cell_cntr, a;

    for(int icell = 0; icell < ncells; icell++){
      cell_cntr = mesh.cell_centroid(icell);
      int nmats = 1, mat_id[2];
      if (fabs(cell_cntr[0] - cell_cntr[1]) < 1e-3){
        nmats = 2.;
        mat_id[0] = 0; mat_id[1] = 1;
      }else if (cell_cntr[0] - cell_cntr[1] > 1e-3){
        mat_id[0] = 0; nmats = 1;
      }else if (cell_cntr[1] - cell_cntr[0] > 1e-3){
        mat_id[0] = 1; nmats = 1;
      }
      os.write((char*) &nmats, sizeof(int));
      for (int imat = 0; imat < nmats; imat++){
        os.write((char*) &mat_id[imat], sizeof(int));
      }
    }

    for(int icell = 0; icell < ncells; icell++){
      cell_cntr = mesh.cell_centroid(icell);
      int nmats = 1;
      double vol_frac = 1.;
      if (fabs(cell_cntr[0] - cell_cntr[1]) < 1e-3){
        nmats = 2; vol_frac = 0.5;
      }
      if (nmats > 1) {
        for (int imat = 0; imat < nmats; imat++){
          vol_frac = 0.5;// + 1e-7*(1 - 2*imat);
          os.write((char*) &vol_frac, sizeof(double));
        }
      }
    }

    for(int icell = 0; icell < ncells; icell++){
      cell_cntr = mesh.cell_centroid(icell);
      double h = sqrt(mesh.cell_volume(icell));
      int nmats = 1;
      if (fabs(cell_cntr[0] - cell_cntr[1]) < 1e-3){
        nmats = 2; 
      }
      if (nmats > 1) {
        a[0] = cell_cntr[0] + h/6.;
        a[1] = cell_cntr[0] - h/6.;
        os.write((char*) &a[0], sizeof(double));
        os.write((char*) &a[1], sizeof(double));
        a[0] = cell_cntr[0] - h/6.;
        a[1] = cell_cntr[0] + h/6.;
        os.write((char*) &a[0], sizeof(double));
        os.write((char*) &a[1], sizeof(double));
      }
    }

    // for(int icell = 0; icell < ncells; icell++) {

    //   int nmats = (int) cfg.cells_mat[icell].size();
    //   os.write((char*) &nmats, sizeof(int));
    //   for (int imat = 0; imat < nmats; imat++)
    //     os.write((char*) &cfg.cells_mat[icell][imat], sizeof(int));
    // }

    // for(int icell = 0; icell < ncells; icell++) {
    //   int nmats = (int) cfg.cells_mat[icell].size();
    //   if (nmats > 1) {
    //     for (int imat = 0; imat < nmats; imat++)
    //       os.write((char*) &cfg.cells_vol_frac[icell][imat], sizeof(double));
    //   }
    // }

    // for(int icell = 0; icell < ncells; icell++) {
    //   int nmats = (int) cfg.cells_mat[icell].size();
    //   if (nmats > 1) {
    //     for (int imat = 0; imat < nmats; imat++) {
    //       os.write((char*) &cfg.cells_centroid[icell][imat].x, sizeof(double));
    //       os.write((char*) &cfg.cells_centroid[icell][imat].y, sizeof(double));
    //     }
    //   }
    // }
    os.close();

}


void write_data_example_poly(const Amanzi::AmanziMesh::Mesh& mesh, const std::string& matdata_fname){

  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;  
  

  std::ofstream os(matdata_fname.c_str(), std::ofstream::binary);

  // if (!os.good())
  //   THROW_EXCEPTION("Cannot open " << cfg.bin_mesh_data_fname << " for binary output");

  int out_int = 2;  //Dimension: 2 for 2D, 3 for 3D

  int ncells = mesh.num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  os.write((char*) &out_int, sizeof(int));
  os.write((char*) &ncells, sizeof(int));
  Point cell_l(2), cell_r(2), p1, p2, pl1(0.35, -1), pl2(0.65,-1), vec(0, 1);
  double area_l, area_r;

  AmanziMesh::Entity_ID_List nodes;

  std::vector<std::vector<int> > mat_ids(ncells);
  std::vector<std::vector<double> > vol_frac(ncells);
  std::vector<std::vector<Point> > part_cntr(ncells);
  

  for (int c=0; c<ncells; c++){
    mesh.cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();
    std::vector<Point> nodes_pnt, intersect_l, intersect_r;

    for (int k = 0; k < nnodes; k++){
      mesh.node_get_coordinates(nodes[k], &p1);
      nodes_pnt.push_back(p1);
    }
    
    SplitPolyCell(nodes_pnt, pl1,  vec, intersect_r);       
    SplitPolyCell(nodes_pnt, pl1, -1*vec, intersect_l);

    if ((intersect_l.size() > 0)&&(intersect_r.size() > 0)){
      mat_ids[c].resize(2);
      mat_ids[c][0] = 0;
      mat_ids[c][1] = 1;
      Intersection_Centroid_Area(intersect_l, cell_l, &area_l);
      Intersection_Centroid_Area(intersect_r, cell_r, &area_r);
      double cell_size = mesh.cell_volume(c);
      vol_frac[c].resize(2);
      vol_frac[c][0] =  area_l/cell_size;
      vol_frac[c][1] =  area_r/cell_size;
      ///Noise
      vol_frac[c][0] *= (1 - 0.2*sin(1.*c));  vol_frac[c][1] = 1 - vol_frac[c][0];
      //
      part_cntr[c].resize(2);
      part_cntr[c][0] = cell_l;
      part_cntr[c][1] = cell_r;
      // std::cout<<"left "<<part_cntr[c][0]<<" right "<< part_cntr[c][1]<<"\n";      
      if (std::abs(area_l + area_r - cell_size) > 1e-10){
        std::cout<< "Wrong area calculation "<<area_l + area_r<<" "<<cell_size<<"\n";
        exit(-1);
      }      
      
    }else if ((intersect_l.size() > 0)&&(intersect_r.size()==0)){
      int id = 0;
      mat_ids[c].push_back(id);
    }           
  }

  for (int c=0; c<ncells; c++){
    mesh.cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();
    std::vector<Point> nodes_pnt, intersect_l, intersect_r;

    for (int k = 0; k < nnodes; k++){
      mesh.node_get_coordinates(nodes[k], &p1);
      nodes_pnt.push_back(p1);
    }
    
    SplitPolyCell(nodes_pnt, pl2,  vec, intersect_r);       
    SplitPolyCell(nodes_pnt, pl2, -1*vec, intersect_l);

    if (mat_ids[c].size()==0){
      if ((intersect_l.size() > 0)&&(intersect_r.size() > 0)){
        mat_ids[c].resize(2);
        mat_ids[c][0] = 1;
        mat_ids[c][1] = 0;
        Intersection_Centroid_Area(intersect_l, cell_l, &area_l);
        Intersection_Centroid_Area(intersect_r, cell_r, &area_r);
        double cell_size = mesh.cell_volume(c);
        vol_frac[c].resize(2);
        vol_frac[c][0] =  area_l/cell_size;
        vol_frac[c][1] =  area_r/cell_size;
        ///Noise
        vol_frac[c][0] *= (1 - 0.2*sin(1.*c));  vol_frac[c][1] = 1 - vol_frac[c][0];
        //
        part_cntr[c].resize(2);
        part_cntr[c][0] = cell_l;
        part_cntr[c][1] = cell_r;
        // std::cout<<"left "<<part_cntr[c][0]<<" right "<< part_cntr[c][1]<<"\n";
      }else if ((intersect_l.size() > 0)&&(intersect_r.size()==0)){
        int id = 1;
        mat_ids[c].push_back(id);
      }else if ((intersect_r.size() > 0)&&(intersect_l.size()==0)){
        int id = 0;
        mat_ids[c].push_back(id);
      }
    }    
  }

  for(int c = 0; c < ncells; c++) {   
    int nmats = mat_ids[c].size();
    os.write((char*) &nmats, sizeof(int));
    for (int imat = 0; imat < nmats; imat++)
        os.write((char*) &mat_ids[c][imat], sizeof(int));
  }

  for(int c = 0; c < ncells; c++) {   
    int nmats = mat_ids[c].size();
    if (nmats > 1) {
      for (int imat = 0; imat < nmats; imat++)
        os.write((char*) &vol_frac[c][imat], sizeof(double));
    }
  }
  for(int c = 0; c < ncells; c++) {   
    int nmats = mat_ids[c].size();
    if (nmats > 1) {
      double val;
      for (int imat = 0; imat < nmats; imat++){
        os.write((char*) &part_cntr[c][imat][0], sizeof(double));
        os.write((char*) &part_cntr[c][imat][1], sizeof(double));
      }
    }
  }
   
  os.close();

  // const Epetra_BlockMap& cmap_owned = mesh.cell_map(false);
  // Epetra_Vector nmat(cmap_owned);

  // for (int c=0; c<ncells; c++){
  //   std::cout<<"Cell cnt: "<<mesh.cell_centroid(c)<<":";
  //   for (int k = 0; k < mat_ids[c].size(); k++) std::cout<<" "<<mat_ids[c][k];
  //   // std::cout<<"\n";
  //   // nmat[c] = mat_ids[c][0] + 2*mat_ids[c].size();
  // }

  // GMV::open_data_file(mesh, (std::string)"nmat_xmof.gmv");
  // GMV::start_data();
  // GMV::write_cell_data(nmat, 0, "nmat");
  // GMV::close_data_file();
  

}

void Intersection_Centroid_Area(const std::vector<Amanzi::AmanziGeometry::Point>& xy1,
                                Amanzi::AmanziGeometry::Point& cntr, double *inter_area)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;

  Point ref_pnt(2);
  int num_pnt = xy1.size();

  ref_pnt.set(0.0);
  if (num_pnt < 1) return;
  
  for (int i=0;i<num_pnt;i++){
    ref_pnt += xy1[i];
  }
  ref_pnt *= 1./num_pnt;

  Point v1(2), v2(2);
  double tri_area;

  cntr.set(0.); 
  (*inter_area) = 0.;
  
  for (int i=0;i<num_pnt;i++){
    v1 = xy1[i] - ref_pnt;
    v2 = xy1[(i+1)%num_pnt] - ref_pnt;
    tri_area = 0.5*norm(v1^v2);
    (*inter_area) += tri_area;
    cntr += (xy1[i] + ref_pnt +  xy1[(i+1)%num_pnt])*(tri_area/3.);    
  }
  cntr *= 1./(*inter_area);
}

                              
                                
void SplitPolyCell(const std::vector<Amanzi::AmanziGeometry::Point>& xy1, 
                   const Amanzi::AmanziGeometry::Point& pnt, 
                   const Amanzi::AmanziGeometry::Point& vec, 
                   std::vector<Amanzi::AmanziGeometry::Point>& xy3){

   using namespace Amanzi;
   using namespace Amanzi::AmanziGeometry;


   std::vector<Point> box;
   Point l_vec = vec;
   Point new_pnt = pnt;

   double a, tmp;

   //std::cout<<"new_pnt "<<new_pnt<<"\n";
   box.push_back(new_pnt);
   for (int i=0; i<4; i++){
     if (i==2) a=2e+6;
     else a=1e+6;
            
     new_pnt = new_pnt + a*l_vec;
     //std::cout<<"new_pnt "<<new_pnt<<"\n";
     box.push_back(new_pnt);
     tmp = l_vec[0];
     l_vec[0] = l_vec[1];
     l_vec[1] = -tmp;
   }

   IntersectConvexPolygons(xy1, box, xy3);

}

void TestSPD(Teuchos::RCP<Operator> global_op){
 
  const CompositeVectorSpace& cvs = global_op->DomainMap();
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  
  //Test SPD properties of the preconditioner.
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


}
