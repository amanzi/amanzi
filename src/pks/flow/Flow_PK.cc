/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1) 
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "GMVMesh.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PK_Utils.hh"
#include "State.hh"
#include "WhetStoneDefs.hh"

#include "Flow_PK.hh"
#include "FracturePermModelPartition.hh"
#include "FracturePermModelEvaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* default constructor that initializes all pointers to NULL
****************************************************************** */
Flow_PK::Flow_PK(Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& soln) :
  PK_PhysicalBDF(pk_tree, glist, S, soln),
  passwd_("flow"),
  peaceman_model_(false)
{
  vo_ = Teuchos::null;
  Teuchos::RCP<Teuchos::ParameterList> units_list = Teuchos::sublist(glist, "units");
  units_.Init(*units_list);
};


Flow_PK::Flow_PK() : passwd_("flow") { vo_ = Teuchos::null; }


/* ******************************************************************
* Setup of static fields common for Darcy and Richards.
****************************************************************** */
void Flow_PK::Setup(const Teuchos::Ptr<State>& S)
{
  // Work flow can be affected by the list of models
  auto physical_models = Teuchos::sublist(fp_list_, "physical models and assumptions");

  // -- type of the flow (in matrix or on manifold)
  flow_on_manifold_ = physical_models->get<bool>("flow in fractures", false);
  flow_on_manifold_ &= (mesh_->manifold_dimension() != mesh_->space_dimension());

  // -- coupling with other PKs
  coupled_to_matrix_ = physical_models->get<std::string>("coupled matrix fracture flow", "") == "fracture";
  coupled_to_fracture_ = physical_models->get<std::string>("coupled matrix fracture flow", "") == "matrix";

  // register fields
  // -- keys
  darcy_flux_key_ = Keys::getKey(domain_, "darcy_flux"); 
  permeability_key_ = Keys::getKey(domain_, "permeability"); 

  // -- constant fields
  if (!S->HasField("const_fluid_density")) {
    S->RequireScalar("const_fluid_density", passwd_);
  }

  if (!S->HasField("atmospheric_pressure")) {
    S->RequireScalar("atmospheric_pressure", passwd_);
  }

  if (!S->HasField("gravity")) {
    S->RequireConstantVector("gravity", passwd_, dim);  // state resets ownership.
  } 

  // -- effective fracture permeability
  if (flow_on_manifold_) {
    if (!S->HasField(permeability_key_)) {
      S->RequireField(permeability_key_, permeability_key_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);

      auto fpm_list = Teuchos::sublist(fp_list_, "fracture permeability models", true);
      Teuchos::RCP<FracturePermModelPartition> fpm = CreateFracturePermModelPartition(mesh_, fpm_list);

      Teuchos::ParameterList elist;
      elist.set<std::string>("permeability key", permeability_key_)
           .set<std::string>("aperture key", Keys::getKey(domain_, "aperture"));
      Teuchos::RCP<FracturePermModelEvaluator> eval = Teuchos::rcp(new FracturePermModelEvaluator(elist, fpm));
      S->SetFieldEvaluator(permeability_key_, eval);
    }
  // -- matrix absolute permeability
  } else {
    if (!S->HasField(permeability_key_)) {
      S->RequireField(permeability_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, dim);
    }
  }

  // -- darcy flux
  if (!S->HasField(darcy_flux_key_)) {
    if (flow_on_manifold_) {
      auto cvs = Operators::CreateNonManifoldCVS(mesh_);
      *S->RequireField(darcy_flux_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true) = *cvs;
    } else {
      S->RequireField(darcy_flux_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("face", AmanziMesh::FACE, 1);
    }
  }

  if (!S->HasFieldEvaluator(darcy_flux_key_)) {
    AddDefaultPrimaryEvaluator(darcy_flux_key_);
  }

  // Wells
  if (!S->HasField("well_index")) {
    if (fp_list_->isSublist("source terms")) {
      Teuchos::ParameterList& src_list = fp_list_->sublist("source terms");
      for (auto it = src_list.begin(); it != src_list.end(); ++it) {
        std::string name = it->first;
        if (src_list.isSublist(name)) {
          Teuchos::ParameterList& spec = src_list.sublist(name);
          if (IsWellIndexRequire(spec)) {
            S->RequireField("well_index", passwd_)->SetMesh(mesh_)->SetGhosted(true)
                ->SetComponent("cell", AmanziMesh::CELL, 1);
            peaceman_model_ = true;
            break;
          }
        }
      }
    }
  }    
}


/* ******************************************************************
* Initiazition of fundamental flow sturctures.
****************************************************************** */
void Flow_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  dt_ = 0.0;

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  nseepage_prev = 0;
  ti_phase_counter = 0;

  // Fundamental physical quantities
  // -- temporarily these quantities are constant
  double* gravity_data;
  S->GetConstantVectorData("gravity")->ExtractView(&gravity_data);
  gravity_.set(dim, &(gravity_data[0]));  // do it in complicated way because we
                                          // are not sure if gravity_data is an
                                          // array or vector
  g_ = fabs(gravity_[dim - 1]);
  rho_ = *S->GetScalarData("const_fluid_density");

  // -- molar rescaling of some quantatities.
  molar_rho_ = rho_ / CommonDefs::MOLAR_MASS_H2O;
  flux_units_ = 0.0;  // scaling from kg to moles

  // parallel execution data
  MyPID = 0;
#ifdef HAVE_MPI
  MyPID = mesh_->cell_map(false).Comm().MyPID();
#endif

  InitializeFields_();
}


/* ****************************************************************
* This completes initialization of common fields that were not 
* initialized by the state.
**************************************************************** */
void Flow_PK::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // set popular default values for missed fields.
  if (S_->GetField("const_fluid_density")->owner() == passwd_) {
    if (!S_->GetField("const_fluid_density", passwd_)->initialized()) {
      *(S_->GetScalarData("const_fluid_density", passwd_)) = 1000.0;
      S_->GetField("const_fluid_density", passwd_)->set_initialized();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized const_fluid_density to default value 1000.0" << std::endl;  
    }
  }

  if (S_->HasField("const_fluid_viscosity")) {
    if (!S_->GetField("const_fluid_viscosity", passwd_)->initialized()) {
      *(S_->GetScalarData("const_fluid_viscosity", passwd_)) = CommonDefs::ISOTHERMAL_VISCOSITY;
      S_->GetField("const_fluid_viscosity", passwd_)->set_initialized();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized const_fluid_viscosity to default value 1.002e-3" << std::endl;  
    }
  }

  if (S_->HasField("atmospheric_pressure")) {
    if (!S_->GetField("atmospheric_pressure", passwd_)->initialized()) {
      *(S_->GetScalarData("atmospheric_pressure", passwd_)) = FLOW_PRESSURE_ATMOSPHERIC;
      S_->GetField("atmospheric_pressure", passwd_)->set_initialized();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized atmospheric_pressure to default value " << FLOW_PRESSURE_ATMOSPHERIC << std::endl;  
    }
  }

  if (!S_->GetField("gravity", "state")->initialized()) {
    Epetra_Vector& gvec = *S_->GetConstantVectorData("gravity", "state");
    gvec.PutScalar(0.0);
    gvec[dim - 1] = -9.80;
    S_->GetField("gravity", "state")->set_initialized();

    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
        *vo_->os() << "initialized gravity to default value -9.8" << std::endl;  
  }

  InitializeField(S_.ptr(), porosity_key_, porosity_key_, 0.2);

  InitializeField(S_.ptr(), passwd_, specific_storage_key_, 0.0);
  InitializeField(S_.ptr(), passwd_, specific_yield_key_, 0.0);

  InitializeField(S_.ptr(), passwd_, pressure_key_, 0.0);
  InitializeField(S_.ptr(), passwd_, hydraulic_head_key_, 0.0);
  InitializeField(S_.ptr(), passwd_, pressure_head_key_, 0.0);

  InitializeField(S_.ptr(), passwd_, darcy_flux_key_, 0.0);
}


/* ****************************************************************
* Hydraulic head support for Flow PKs.
**************************************************************** */
void Flow_PK::UpdateLocalFields_(const Teuchos::Ptr<State>& S) 
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
    *vo_->os() << "Secondary fields: hydraulic head, darcy_velocity, etc." << std::endl;  
  }  

  Epetra_MultiVector& hydraulic_head = *(S->GetFieldData(hydraulic_head_key_, passwd_)->ViewComponent("cell"));
  const Epetra_MultiVector& pressure = *(S->GetFieldData(pressure_key_)->ViewComponent("cell"));
  double rho = *(S->GetScalarData("const_fluid_density"));

  // calculate hydraulic head
  double g = fabs(gravity_[dim - 1]);
 
  for (int c = 0; c != ncells_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c); 
    double z = xc[dim - 1]; 
    hydraulic_head[0][c] = z + (pressure[0][c] - atm_pressure_) / (g * rho);
  }

  // calculate optional fields
  Key optional_key = Keys::getKey(domain_, "pressure_head"); 
  if (S->HasField(optional_key)) {
    auto& field_c = *S->GetFieldData(optional_key, passwd_)->ViewComponent("cell");
    for (int c = 0; c != ncells_owned; ++c) {
      field_c[0][c] = pressure[0][c] / (g * rho);
    }
  }

  // calculate full velocity vector
  darcy_flux_eval_->SetFieldAsChanged(S);
  S->GetFieldEvaluator(darcy_velocity_key_)->HasFieldChanged(S, darcy_velocity_key_);
}


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Flow_PK::InitializeBCsSources_(Teuchos::ParameterList& plist)
{
  // Process main one-line options (not sublists)
  atm_pressure_ = *S_->GetScalarData("atmospheric_pressure");

  // Create BC objects
  // -- memory
  op_bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  Teuchos::RCP<FlowBoundaryFunction> bc;
  auto& bc_list = plist.sublist("boundary conditions");

  bcs_.clear();

  // -- pressure 
  if (bc_list.isSublist("pressure")) {
    PK_DomainFunctionFactory<FlowBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("pressure");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc = bc_factory.Create(spec, "boundary pressure", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("pressure");
        bcs_.push_back(bc);
      }
    }
  }

  // -- hydraulic head
  if (bc_list.isSublist("static head")) {
    PK_DomainFunctionFactory<FlowBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("static head");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc = bc_factory.Create(spec, "static head", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("head");
        bcs_.push_back(bc);
      }
    }
  }

  // -- Darcy velocity
  if (bc_list.isSublist("mass flux")) {
    PK_DomainFunctionFactory<FlowBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("mass flux");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      if (it->second.isList()) {
        Teuchos::ParameterList spec = Teuchos::getValue<Teuchos::ParameterList>(it->second);
        bc = bc_factory.Create(spec, "outward mass flux", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("flux");
        bcs_.push_back(bc);
      }
    }
  }

  // -- seepage face
  if (bc_list.isSublist("seepage face")) {
    PK_DomainFunctionFactory<FlowBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("seepage face");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc = bc_factory.Create(spec, "outward mass flux", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("seepage");
        bcs_.push_back(bc);
      }
    }
  }

  VV_ValidateBCs();

  // Create source objects
  // -- evaluate the well index
  if (S_->HasField("well_index")) {
    if (!S_->GetField("well_index", passwd_)->initialized()) {
      S_->GetFieldData("well_index", passwd_)->PutScalar(0.0);

      Teuchos::ParameterList& src_list = plist.sublist("source terms");
      for (auto it = src_list.begin(); it != src_list.end(); ++it) {
        std::string name = it->first;
        if (src_list.isSublist(name)) {
          Teuchos::ParameterList& spec = src_list.sublist(name);
          if (IsWellIndexRequire(spec)) {
            ComputeWellIndex(spec);
          }
        }
      }
    }
    S_->GetField("well_index", passwd_)->set_initialized();
  }

  // -- wells
  srcs.clear();
  if (plist.isSublist("source terms")) {
    PK_DomainFunctionFactory<PK_DomainFunction> factory(mesh_, S_);
    PKUtils_CalculatePermeabilityFactorInWell(S_.ptr(), Kxy);

    Teuchos::ParameterList& src_list = plist.sublist("source terms");
    for (auto it = src_list.begin(); it != src_list.end(); ++it) {
      std::string name = it->first;
      if (src_list.isSublist(name)) {
        Teuchos::ParameterList& spec = src_list.sublist(name);
        srcs.push_back(factory.Create(spec, "well", AmanziMesh::CELL, Kxy));
      }
    }
  }
}


/* ******************************************************************
* Compute well index
****************************************************************** */
void Flow_PK::ComputeWellIndex(Teuchos::ParameterList& spec)
{
  AmanziMesh::Entity_ID_List cells, faces;
  Epetra_MultiVector& wi = *S_->GetFieldData("well_index", passwd_)->ViewComponent("cell");
  const Epetra_MultiVector& perm = *S_->GetFieldData(permeability_key_)->ViewComponent("cell");

  double kx, ky, dx, dy, h, r0, rw;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  int d = mesh_->space_dimension();

  std::vector<std::string> regions = spec.get<Teuchos::Array<std::string> >("regions").toVector();
  Teuchos::ParameterList well_list = spec.sublist("well");
  rw = well_list.get<double>("well radius");
  
  for (auto it = regions.begin(); it != regions.end(); ++it) {
    mesh_->get_set_entities(*it, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &cells);

    for (int k = 0; k < cells.size(); k++) {
      int c = cells[k];
      mesh_->cell_get_faces(c, &faces);
      xmin = ymin = zmin = 1e+23;
      xmax = ymax = zmax = 1e-23;
      for (int j = 0; j < faces.size(); j++) {
        int f = faces[j];
        const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
        xmax = std::max(xmax, xf[0]);
        xmin = std::min(xmin, xf[0]);
        ymax = std::max(ymax, xf[1]);
        ymin = std::min(ymin, xf[1]);
        if (d > 2) {
          zmax = std::max(zmax, xf[2]);
          zmin = std::min(zmin, xf[2]);
        }        
      }
      dx = xmax - xmin;
      dy = ymax - ymin;
      if (d > 2) h = zmax - zmin;
      else h = 1.0;

      kx = perm[0][c];
      ky = perm[1][c];

      r0 = 0.28 * sqrt(dy*dy * sqrt(kx/ky) + dx*dx * sqrt(ky/kx));
      r0 = r0 / (std::pow(kx/ky, 0.25) + std::pow(ky/kx, 0.25));
                      
      // r0 = 0.28 * sqrt(sqrt(kx/ky)*dx + sqrt(ky/kx)*dy);
      // r0 = r0 / (sqrt(sqrt(kx/ky)) + sqrt(sqrt(ky/kx)));
      
      wi[0][c] = 2 * M_PI * sqrt(kx/ky) * h / log(r0/rw);
    }
  }
}


/* ******************************************************************
* Analyze spec for well signature
****************************************************************** */
bool Flow_PK::IsWellIndexRequire(Teuchos::ParameterList& spec)
{
  if (spec.isParameter("spatial distribution method")) {
    std::string model = spec.get<std::string>("spatial distribution method");

    if (model == "simple well") {
      Teuchos::ParameterList well_list = spec.sublist("well");
      if (well_list.isParameter("submodel")) {
        if (well_list.get<std::string>("submodel") == "bhp") {
          return true;
        }
      }
    }
  }
  return false;
}


/* ******************************************************************
* A wrapper for updating boundary conditions.
****************************************************************** */
void Flow_PK::UpdateSourceBoundaryData(double t_old, double t_new, const CompositeVector& u)
{
  for (int i = 0; i < srcs.size(); ++i) {
    srcs[i]->Compute(t_old, t_new); 
  }

  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(t_old, t_new);
    bcs_[i]->ComputeSubmodel(mesh_);
  }

  ComputeOperatorBCs(u);
}


/* ******************************************************************
* Add a boundary marker to used faces.
* WARNING: we can skip update of ghost boundary faces, b/c they 
* should be always owned. 
****************************************************************** */
void Flow_PK::ComputeOperatorBCs(const CompositeVector& u)
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();
  std::vector<double>& bc_mixed = op_bc_->bc_mixed();

  for (int n = 0; n < bc_model.size(); n++) {
    bc_model[n] = Operators::OPERATOR_BC_NONE;
    bc_value[n] = 0.0;
    bc_mixed[n] = 0.0;
  }

  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "pressure") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = it->second[0];
      }
    }

    if (bcs_[i]->get_bc_name() == "head") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        if (bcs_[i]->no_flow_above_water_table()) {
          if (it->second[0] < atm_pressure_) {
            bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
            bc_value[f] = 0.0;
            continue;
          }
        }
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = it->second[0];
      }
    }

    if (bcs_[i]->get_bc_name() == "flux") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value[f] = it->second[0] * flux_units_;
      }
    }
  }

  // Seepage face BC is implemented for p-lambda discretization only.
  int nseepage_add, nseepage = 0;
  double area_add, area_seepage = 0.0;

  SeepageFacePFloTran(u, &nseepage_add, &area_add);
  nseepage += nseepage_add;
  area_seepage += area_add;

  SeepageFaceFACT(u, &nseepage_add, &area_add);
  nseepage += nseepage_add;
  area_seepage += area_add;

  // mark missing boundary conditions as zero flux conditions
  AmanziMesh::Entity_ID_List cells;
  missed_bc_faces_ = 0;
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] == Operators::OPERATOR_BC_NONE) {
      cells.clear();
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value[f] = 0.0;
        missed_bc_faces_++;
      }
    }
  }

  dirichlet_bc_faces_ = 0;
  for (int f = 0; f < nfaces_owned; ++f) {
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) dirichlet_bc_faces_++;
  }
  int flag_essential_bc = (dirichlet_bc_faces_ > 0) ? 1 : 0;

  // verify that the algebraic problem is consistent
#ifdef HAVE_MPI
  int flag = flag_essential_bc;
  mesh_->get_comm()->MaxAll(&flag, &flag_essential_bc, 1);  // find the global maximum
#endif
  if (! flag_essential_bc &&
      domain_ == "domain" && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "WARNING: no essential boundary conditions, solver may fail" << std::endl;
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
#ifdef HAVE_MPI
    int nseepage_tmp = nseepage;
    double area_tmp = area_seepage;
    mesh_->get_comm()->SumAll(&area_tmp, &area_seepage, 1);
    mesh_->get_comm()->SumAll(&nseepage_tmp, &nseepage, 1);
#endif
    if (MyPID == 0 && nseepage > 0 && nseepage != nseepage_prev) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "seepage face: " << area_seepage << " [m^2], from "
                 << nseepage_prev << " to " << nseepage << " faces" << std::endl;
    }
  }
  nseepage_prev = nseepage;
}


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void Flow_PK::SetAbsolutePermeabilityTensor()
{
  const CompositeVector& cv = *S_->GetFieldData(permeability_key_);
  cv.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& perm = *cv.ViewComponent("cell", true);

  // For permeabilities given in local (layer-based) coordinates
  AmanziGeometry::Point n1(dim), n2(dim), normal(dim), tau(dim);
  WhetStone::Tensor N(dim, 2), Ninv(dim, 2), D(dim, 2);

  K.resize(ncells_owned);
  bool cartesian = (coordinate_system_ == "cartesian");
  bool off_diag = cv.HasComponent("offd");

  // most common cases of diagonal permeability
  if (cartesian && dim == 2) {
    for (int c = 0; c < ncells_owned; c++) {
      if (!off_diag && perm[0][c] == perm[1][c]) {
        K[c].Init(dim, 1);
        K[c](0, 0) = perm[0][c];
      } else {
        K[c].Init(dim, 2);
        K[c](0, 0) = perm[0][c];
        K[c](1, 1) = perm[1][c];
      }
    }    
  } 
  else if (cartesian && dim == 3) {
    for (int c = 0; c < K.size(); c++) {
      if (!off_diag && perm[0][c] == perm[1][c] && perm[0][c] == perm[2][c]) {
        K[c].Init(dim, 1);
        K[c](0, 0) = perm[0][c];
      } else {
        K[c].Init(dim, 2);
        K[c](0, 0) = perm[0][c];
        K[c](1, 1) = perm[1][c];
        K[c](2, 2) = perm[2][c];
      }
    }        
  }

  // special case of layer-oriented permeability
  if (!cartesian && dim == 2) {
    for (int c = 0; c < ncells_owned; c++) {
      VerticalNormals(c, n1, n2);
      normal = (n1 - n2) / 2;
      normal /= norm(normal);

      tau[0] = normal[1];
      tau[1] = -normal[0];
        
      N.SetColumn(0, tau); 
      N.SetColumn(1, normal); 

      Ninv = N;
      Ninv.Inverse();

      D(0, 0) = perm[0][c];
      D(1, 1) = perm[1][c];
      K[c] = N * D * Ninv;
    }    
  } 

  // special case of permeability with off-diagonal components 
  if (cartesian && off_diag) {
    const Epetra_MultiVector& offd = *cv.ViewComponent("offd");

    for (int c = 0; c < ncells_owned; c++) {
      if (dim == 2) {
        K[c](0, 1) = K[c](1, 0) = offd[0][c];
      } 
      else if (dim == 3) {
        K[c](0, 1) = K[c](1, 0) = offd[0][c];
        K[c](0, 2) = K[c](2, 0) = offd[1][c];
        K[c](1, 2) = K[c](2, 1) = offd[2][c];
      }  
    }
  }
}


/* ******************************************************************
* Add source and sink terms.                                   
****************************************************************** */
void Flow_PK::AddSourceTerms(CompositeVector& rhs)
{
  Epetra_MultiVector& rhs_cell = *rhs.ViewComponent("cell");

  for (int i = 0; i < srcs.size(); ++i) {
    for (auto it = srcs[i]->begin(); it != srcs[i]->end(); ++it) {
      int c = it->first;
      rhs_cell[0][c] += mesh_->cell_volume(c) * it->second[0];
    }
  }
}


/* ******************************************************************
* BDF methods need a good initial guess.
* WARNING: Each owned face must have at least one owned cell. 
* Probability that this assumption is violated is close to zero. 
* Even when it happens, the code will not crash.
****************************************************************** */
void Flow_PK::DeriveFaceValuesFromCellValues(
    const Epetra_MultiVector& ucells, Epetra_MultiVector& ufaces)
{
  AmanziMesh::Entity_ID_List cells;
  auto& fmap = ufaces.Map(); 

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  for (int f = 0; f < nfaces; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n = 0; n < ncells; n++) face_value += ucells[0][cells[n]];
    double pmean = face_value / ncells;

    int first = fmap.FirstPointInElement(f);
    int ndofs = fmap.ElementSize(f);
    for (int k = 0; k < ndofs; ++k) ufaces[0][first + k] = pmean;
  }
}


/* ******************************************************************
* Calculate change of water volume per second due to boundary flux.                                          
****************************************************************** */
double Flow_PK::WaterVolumeChangePerSecond(const std::vector<int>& bc_model,
                                           const Epetra_MultiVector& darcy_flux) const
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  double volume = 0.0;
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      if (bc_model[f] != Operators::OPERATOR_BC_NONE && f < nfaces_owned) {
        if (fdirs[i] >= 0) {
          volume -= darcy_flux[0][f];
        } else {
          volume += darcy_flux[0][f];
        }
      }
    }
  }
  return volume;
}


/* ******************************************************************
* Returns the first cell attached to a boundary face.   
****************************************************************** */
int Flow_PK::BoundaryFaceGetCell(int f) const
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  return cells[0];
}


/* ******************************************************************
* Returns approximation of a solution on a boundary face   
****************************************************************** */
double Flow_PK::BoundaryFaceValue(int f, const CompositeVector& u)
{
  double face_value;
  if (u.HasComponent("face")) {
    const Epetra_MultiVector& u_face = *u.ViewComponent("face");
    face_value = u_face[0][f];
  } else {
    const Epetra_MultiVector& u_cell = *u.ViewComponent("cell");
    int c = BoundaryFaceGetCell(f);
    face_value = u_cell[0][c];
  }
  return face_value;
}


/* ******************************************************************
* Find cell normals that have direction close to gravity (n1) and
* anti-gravity (n2).
****************************************************************** */
void Flow_PK::VerticalNormals(int c, AmanziGeometry::Point& n1, AmanziGeometry::Point& n2)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  int i1, i2;
  double amax(-1e+50), amin(1e+50), a;
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    double area = mesh_->face_area(f);
    const AmanziGeometry::Point normal = mesh_->face_normal(f);

    a = normal[dim - 1] * dirs[i] / area;
    if (a > amax) { 
      i1 = i;
      amax = a;
    } 
    if (a < amin) { 
      i2 = i;
      amin = a;
    } 
  }

  n1 = mesh_->face_normal(faces[i1]) * dirs[i1];
  n2 = mesh_->face_normal(faces[i2]) * dirs[i2];
}

}  // namespace Flow
}  // namespace Amanzi

