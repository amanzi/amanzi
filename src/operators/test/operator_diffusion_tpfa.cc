/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
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

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "LinearOperatorFactory.hh"

#include "tensor.hh"
#include "mfd3d_diffusion.hh"

#include "OperatorDefs.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorSource.hh"
#include "UpwindSecondOrder.hh"
#include "UpwindStandard.hh"

#include "Analytic04.hh"

namespace Amanzi{

// This class wraps scalar diffusion coefficient.
class HeatConduction {
 public:
  HeatConduction(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh), ana_(mesh) { 
    int dim = mesh_->space_dimension();
    cvs_.SetMesh(mesh_);
    cvs_.SetGhosted(true);
    cvs_.SetComponent("cell", AmanziMesh::CELL, 1);
    cvs_.SetOwned(false);

    values_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));
    derivatives_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));
  }
  ~HeatConduction() {};

  // main members
  void UpdateValues(const CompositeVector& u) { 
    Epetra_MultiVector& vcell = *values_->ViewComponent("cell", true); 
    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
      const WhetStone::Tensor& Kc = ana_.Tensor(xc, 0.0);
      vcell[0][c] = Kc(0, 0);
    }
  }

  void UpdateValuesPostUpwind() { 
    if (!values_->HasComponent("twin")) {
      cvs_.AddComponent("twin", AmanziMesh::FACE, 1);
      Teuchos::RCP<CompositeVector> tmp = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));

      *tmp->ViewComponent("cell") = *values_->ViewComponent("cell"); 
      *tmp->ViewComponent("face") = *values_->ViewComponent("face"); 
      *tmp->ViewComponent("grad") = *values_->ViewComponent("grad"); 
      values_ = tmp;
    }

    AmanziMesh::Entity_ID_List cells;
    Epetra_MultiVector& vcell = *values_->ViewComponent("cell", true); 
    Epetra_MultiVector& vface = *values_->ViewComponent("face", true); 
    Epetra_MultiVector& vtwin = *values_->ViewComponent("twin", true); 

    vtwin = vface;
    int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

    for (int f = 0; f < nfaces; f++) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int ncells = cells.size();
      
      if (ncells == 2) {
        double v1 = vcell[0][cells[0]];
        double v2 = vcell[0][cells[1]];
        if (fabs(v1 - v2) > 2 * std::min(fabs(v1), fabs(v2))) {
          vface[0][f] = v1;
          vtwin[0][f] = v2;
        }  
      } 
    }
  }

  double Conduction(int c, double T) const {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana_.Tensor(xc, 0.0);
    return Kc(0, 0);
  }

  Teuchos::RCP<CompositeVector> values() { return values_; }
  Teuchos::RCP<CompositeVector> derivatives() { return derivatives_; }
   
 private:
  CompositeVectorSpace cvs_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<CompositeVector> values_, derivatives_;
  mutable Analytic04 ana_;
};

typedef double(HeatConduction::*ModelUpwindFn)(int c, double T) const; 
}  // namespace Amanzi


int BoundaryFaceGetCell(const Amanzi::AmanziMesh::Mesh& mesh, int f)
{
  Amanzi::AmanziMesh::Entity_ID_List cells;
  mesh.face_get_cells(f, Amanzi::AmanziMesh::USED, &cells);
  return cells[0];
}


/* *****************************************************************
* This tests diffusion solvers with zero coefficients
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_TPFA_ZEROCOEF) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, nodal discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_strip.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a mesh
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(-4.0, 0.0, 4.0, 1.0, 30, 1, gm);

  // model
  Analytic04 ana(mesh);

  // modify diffusion coefficient
  // -- tensor part
  std::vector<WhetStone::Tensor> K;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  for (int c = 0; c != ncells; ++c) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K.push_back(Kc);
  }

  // -- scalar part
  Teuchos::RCP<CompositeVectorSpace> cell_space = Teuchos::rcp(new CompositeVectorSpace());
  cell_space->SetMesh(mesh)->SetComponent("cell",CELL,1)->SetGhosted();

  Teuchos::RCP<CompositeVectorSpace> face_space = Teuchos::rcp(new CompositeVectorSpace());
  face_space->SetMesh(mesh)->SetComponent("face",FACE,1)->SetGhosted();
  Teuchos::RCP<CompositeVector> coef = Teuchos::rcp(new CompositeVector(*face_space));

  {
    Epetra_MultiVector& coef_faces = *coef->ViewComponent("face",false);

    for (int f = 0; f != nfaces; ++f) {
      const Point& xf = mesh->face_centroid(f);
      coef_faces[0][f] = ana.ScalarCoefficient(xf, 0.0);
    }
  }
  coef->ScatterMasterToGhosted("face");
  double rho(1.0), mu(1.0);
  
  // create boundary data
  Point xv(2);
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);
  std::vector<double> bc_mixed;

  AmanziMesh::Entity_ID_List left;
  mesh->get_set_entities("Left side", AmanziMesh::FACE, AmanziMesh::USED, &left);
  for (int f=0; f!=left.size(); ++f) {
    bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
    mesh->face_centroid(f, &xv);
    bc_value[f] = ana.pressure_exact(xv, 0.);
  }

  AmanziMesh::Entity_ID_List right;
  mesh->get_set_entities("Right side", AmanziMesh::FACE, AmanziMesh::USED, &right);
  for (int f=0; f!=right.size(); ++f) {
    bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
    mesh->face_centroid(f, &xv);
    bc_value[f] = ana.pressure_exact(xv, 0.);
  }
  
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator");
  OperatorDiffusionFactory opfactory;
  Point g(2); g[0] = 0.; g[1] = 0.;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc, op_list, g, 0);
  const CompositeVectorSpace& cvs = op->DomainMap();
  
  // populate the diffusion operator
  op->Setup(K, coef, Teuchos::null, rho, mu);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);
  op->ApplyBCs();
  op->SymbolicAssembleMatrix(op->schema_prec_dofs());
  op->AssembleMatrix(op->schema_prec_dofs());

  // create source and add it to the operator
  CompositeVector source(*cell_space);
  Epetra_MultiVector& src = *source.ViewComponent("cell", false);
  src.PutScalar(0.0);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);
    src[0][c] = ana.source_exact(xc, 0.0) * volume;
  }

  Teuchos::RCP<OperatorSource> op1 = Teuchos::rcp(new OperatorSource(*op));
  op1->UpdateMatrices(source);

  CompositeVector solution(*cell_space);
  {
    Epetra_MultiVector& p_cell = *solution.ViewComponent("cell");
    for (int c = 0; c < ncells; c++) {
      const Point& xc = mesh->cell_centroid(c);
      p_cell[0][c] = ana.pressure_exact(xc, 0.0);
    }
  }

  CompositeVector residual(*cell_space);
  //  int ierr = op1->ComputeNegativeResidual(solution, residual);
  int ierr = op1->Apply(solution, residual);
  CHECK(!ierr);

  solution.Print(std::cout);
  residual.Print(std::cout);
  op1->rhs()->Print(std::cout);
  
  double res_norm(0.);
  ierr |= residual.Norm2(&res_norm);
  CHECK(!ierr);
  CHECK_CLOSE(0., res_norm, 1.e-8);
}


