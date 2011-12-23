/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "Epetra_Vector.h"
#include "Epetra_Import.h"

#include "errors.hh"
#include "exceptions.hh"

#include "mfd3d.hpp"
#include "tensor.hpp"

#include "Flow_State.hpp"
#include "Flow_BC_Factory.hpp"
#include "Darcy_PK.hpp"
#include "Matrix_MFD.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* We set up only default values and call Init() routine to complete
* each variable initialization
****************************************************************** */
Darcy_PK::Darcy_PK(Teuchos::ParameterList& dp_list_, Teuchos::RCP<Flow_State> FS_MPC)
{
  FS = FS_MPC;
  dp_list = dp_list_;

  mesh_ = FS->get_mesh();
  dim = mesh_->space_dimension();

  // Create the combined cell/face DoF map.
  super_map_ = create_super_map();
 
  // Other fundamental physical quantaties
  double rho = *(FS->get_fluid_density());
  double mu = *(FS->get_fluid_viscosity()); 
  gravity.init(dim);
  for (int k=0; k<dim; k++) gravity[k] = (*(FS->get_gravity()))[k];

#ifdef HAVE_MPI
  const  Epetra_Comm & comm = mesh_->cell_map(false).Comm(); 
  MyPID = comm.MyPID();

  const Epetra_Map& source_cmap = mesh_->cell_map(false);
  const Epetra_Map& target_cmap = mesh_->cell_map(true);

  cell_importer_ = Teuchos::rcp(new Epetra_Import(target_cmap, source_cmap));

  const Epetra_Map& source_fmap = mesh_->face_map(false);
  const Epetra_Map& target_fmap = mesh_->face_map(true);

  face_importer_ = Teuchos::rcp(new Epetra_Import(target_fmap, source_fmap));
#endif

  Flow_PK::Init(FS_MPC);
}



/* ******************************************************************
* Extract information from Diffusion Problem parameter list.
****************************************************************** */
void Darcy_PK::Init(Matrix_MFD* matrix_, Matrix_MFD* preconditioner_)
{
  if (matrix_ == NULL) matrix = new Matrix_MFD(FS, *super_map_);
  else matrix = matrix_;

  if (preconditioner_ == NULL) preconditioner = matrix;
  else preconditioner = preconditioner_;

  // Create the solution vectors.
  solution = Teuchos::rcp(new Epetra_Vector(*super_map_));
  solution_cells = Teuchos::rcp(FS->create_cell_view(*solution));
  solution_faces = Teuchos::rcp(FS->create_face_view(*solution));
  rhs = Teuchos::rcp(new Epetra_Vector(*super_map_));

  solver = new AztecOO;
  solver->SetUserOperator(matrix);
  solver->SetPrecOperator(preconditioner);
  solver->SetAztecOption(AZ_solver, AZ_cg);

  // Get some solver parameters from the flow parameter list.
  process_parameter_list();

  // Process boundary data
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  std::vector<int> bc_markers(FLOW_BC_FACE_NULL, ncells);
  std::vector<double> bc_values(0.0, ncells);

  T_physical = FS->get_time();
  double time = (standalone_mode) ? T_internal : T_physical;

  bc_pressure->Compute(time);
  update_data_boundary_faces(bc_pressure, bc_head, bc_flux, bc_markers, bc_values);

  // Process other fundamental structures
  K.resize(number_owned_cells);
  matrix->symbolicAssembleGlobalMatrices(*super_map_);

  // Preconditioner
  Teuchos::ParameterList ML_list = dp_list.sublist("ML Parameters");
  preconditioner->init_ML_preconditioner(ML_list); 
};


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Darcy_PK::process_parameter_list()
{
  Teuchos::ParameterList preconditioner_list;
  preconditioner_list = dp_list.get<Teuchos::ParameterList>("Diffusion Preconditioner");

  max_itrs = dp_list.get<int>("Max Iterations");
  err_tol = dp_list.get<double>("Error Tolerance");

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList> bc_list = Teuchos::rcpFromRef(dp_list.sublist("boundary conditions", true));
  FlowBCFactory bc_factory(mesh_, bc_list);

  bc_pressure = bc_factory.CreatePressure();
  bc_head = bc_factory.CreateStaticHead(0.0, rho, gravity[dim - 1]);
  bc_flux = bc_factory.CreateMassFlux();

  validate_boundary_conditions(bc_pressure, bc_head, bc_flux);  

  double time = (standalone_mode) ? T_internal : T_physical;
  bc_pressure->Compute(time);
  bc_head->Compute(time);
  bc_flux->Compute(time);
}


/* ******************************************************************
*  Calculates steady-state solution assuming that abosolute permeability 
* does not depend on time.                                                    
****************************************************************** */
int Darcy_PK::advance_to_steady_state()
{
  int n = number_owned_cells + number_owned_faces;

  // work-around limited support for tensors
  std::vector<WhetStone::Tensor> K(number_owned_cells);
  populate_absolute_permeability_tensor(K);

  for (int c=0; c<K.size(); c++) K[c] *= rho / mu;

  // calculate and assemble elemental stifness matrices
  matrix->createMFDstiffnessMatrices(K);
  matrix->assembleGlobalMatrices(*rhs);
  addGravityFluxes_MFD(matrix);
  matrix->applyBoundaryConditions(bc_markers, bc_values);
  matrix->computeSchurComplement();

  Epetra_Vector b(*rhs);
  solver->SetRHS(&b);
  solver->SetLHS(&*solution);  // initial solution guess 

  solver->Iterate(max_itrs, err_tol);
  num_itrs = solver->NumIters();
  residual = solver->TrueResidual();

  std::cout << "Darcy solver performed " << solver->NumIters() << " iterations." << std::endl
            << "Norm of true residual = " << solver->TrueResidual() << std::endl;

  Epetra_Vector& darcy_flux = FS->ref_darcy_flux();
  matrix->deriveDarcyFlux(*solution, *rhs, *face_importer_, darcy_flux);
  return 0;
}


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void Darcy_PK::populate_absolute_permeability_tensor(std::vector<WhetStone::Tensor>& K)
{
  const Epetra_Vector& permeability = FS->ref_absolute_permeability();

  for (int c=cmin; c<=cmax; c++) {
    K[c].init(dim, 1);
    K[c](0, 0) = permeability[c];
  }
}


/* ******************************************************************
* .                                               
****************************************************************** */
void Darcy_PK::addGravityFluxes_MFD(Matrix_MFD* matrix)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<double> gravity_flux;
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c=0; c<ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    calculateGravityFluxes(c, K[c], gravity_flux);
    
    Teuchos::SerialDenseMatrix<int, double>& Bff = matrix->get_Aff_cells()[c];
    Epetra_SerialDenseVector& Ff = matrix->get_Ff_cells()[c];

    for (int n=0; n<nfaces; n++) {
      double colsum = 0.0;
      for (int m=0; m<nfaces; m++) colsum += Bff(n, m) * gravity_flux[m];
      Ff[n] = colsum;
    }
  }
}


/* ******************************************************************
*  Printing information about Flow status.                                                     
****************************************************************** */
void Darcy_PK::print_statistics() const
{
  if (!MyPID && verbosity_level > 0) {
    cout << "Flow PK:" << endl;
    cout << "    Execution mode = " << (standalone_mode ? "standalone" : "MPC") << endl;
    cout << "    Verbosity level = " << verbosity_level << endl;
    cout << "    Enable internal tests = " << (internal_tests ? "yes" : "no")  << endl;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

