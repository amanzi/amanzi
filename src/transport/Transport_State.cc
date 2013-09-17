
#include "Transport_State.hh"

namespace Amanzi {
namespace AmanziTransport {

Transport_State::Transport_State(Teuchos::RCP<AmanziMesh::Mesh> mesh, const int ncomp) :
    PK_State(std::string("state"), mesh)
{
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State());
  S_time_ = S;

  // create dummy component names
  for (int i=0; i<ncomp; ++i) {
    std::stringstream ss;
    ss << "Component " << i; 
    comp_names_.push_back(ss.str());
    comp_numbers_[ss.str()] = i;
  }

  Construct_(ncomp);
}


Transport_State::Transport_State(Teuchos::RCP<State> S, const int ncomp) :
  PK_State(std::string("state"), S),
  S_time_(S)
{
  Construct_(ncomp);
}


Transport_State::Transport_State(Teuchos::RCP<State> S, std::vector<std::string> comp_names) :
    PK_State(std::string("state"), S),
    S_time_(S),
    comp_names_(comp_names)
{
  for (int i=0; i<comp_names_.size(); ++i) {
    comp_numbers_[comp_names_[i]] = i;
  }  
 
  Construct_(comp_names_.size());
}


Transport_State::Transport_State(Transport_State& other,
        PKStateConstructMode mode) :
    PK_State(other, STATE_CONSTRUCT_MODE_COPY_POINTERS),
    comp_numbers_(other.comp_numbers_),
    comp_names_(other.comp_names_),
    S_time_(other.S_time_)
{
  if (mode == CONSTRUCT_MODE_VIEW_DATA) {
    ghosted_ = false;  // non-ghosted views
  } else if (mode == CONSTRUCT_MODE_VIEW_DATA_GHOSTED) {
    ghosted_ = true; // no guarantees -- if other is not ghosted, this is not
                     // ghosted either!
  } else if (mode == CONSTRUCT_MODE_COPY_DATA) {
    // Not currently needed by Transport?
    ASSERT(0);
  } else if (mode == CONSTRUCT_MODE_COPY_DATA_GHOSTED) {
    // Pointers for all but a copy on ghosted TCC and flux
    ghosted_ = true;

    CompositeVectorFactory fac_tcc;
    fac_tcc.SetMesh(mesh_);
    fac_tcc.SetComponent("cell", AmanziMesh::CELL, other.total_component_concentration()->NumVectors()   );

    Teuchos::RCP<CompositeVector> tcc = fac_tcc.CreateVector(true);
    tcc->CreateData();
    *tcc->ViewComponent("cell",false) = *other.total_component_concentration();
    tcc->ScatterMasterToGhosted();

    // setting the data also requires making a new field record
    Teuchos::RCP<Field> other_tcc_field =
        S_->GetField("total_component_concentration", name_);
    Teuchos::RCP<Field> new_tcc_field = other_tcc_field->Clone();
    S_->SetField("total_component_concentration", name_, new_tcc_field);
    S_->SetData("total_component_concentration", name_, tcc);

    Teuchos::RCP<Epetra_MultiVector> tcc_now = total_component_concentration();
    Teuchos::RCP<Epetra_MultiVector> tcc_other_now = other.total_component_concentration();

    CompositeVectorFactory fac;
    fac.SetMesh(mesh_);
    fac.SetComponent("face", AmanziMesh::FACE, 1);

    Teuchos::RCP<CompositeVector> flux = fac.CreateVector(true);
    flux->CreateData();
    *flux->ViewComponent("face",false) = *other.darcy_flux();
    flux->ScatterMasterToGhosted();

    // clone the field, assign the vector
    Teuchos::RCP<Field> other_flux_field =
        S_->GetField("darcy_flux", name_);
    Teuchos::RCP<Field> new_flux_field = other_flux_field->Clone();
    S_->SetField("darcy_flux", name_, new_flux_field);
    S_->SetData("darcy_flux", name_, flux);
  }
}


void Transport_State::Construct_(const int ncomp) 
{
  // Require data, all owned by "state" to cheat the system.
  S_->RequireScalar("fluid_density", name_);
  if (!S_->HasField("porosity")) {
    S_->RequireField("porosity", name_)->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("water_saturation")) {
    S_->RequireField("water_saturation", name_)->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("prev_water_saturation")) {
    S_->RequireField("prev_water_saturation", name_)->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("darcy_flux")) {
    S_->RequireField("darcy_flux", name_)->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("face", AmanziMesh::FACE, 1);
  }
  
  if (!S_->HasField("total_component_concentration")) {
    std::vector<std::vector<std::string> > subfield_names(1);
    if (comp_names_.size() == ncomp) {
      for (std::vector<std::string>::const_iterator compname = comp_names_.begin();
	   compname != comp_names_.end(); ++compname) {    
	subfield_names[0].push_back(*compname + std::string(" conc"));
      }
    } else {
      for (int icn = 0; icn != ncomp; ++icn) {
	std::stringstream ss;
	ss << "Component " << icn; 
	subfield_names[0].push_back(ss.str());
      }
    }
    S_->RequireField("total_component_concentration", name_, subfield_names)->SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, ncomp);
  }
}


void Transport_State::Initialize() {
  if (standalone_mode_) {
    S_->Setup();
    S_->GetField("total_component_concentration",name_)->set_initialized();
    S_->GetField("fluid_density",name_)->set_initialized();
    S_->GetField("porosity",name_)->set_initialized();
    S_->GetField("water_saturation",name_)->set_initialized();
    S_->GetField("prev_water_saturation",name_)->set_initialized();
    S_->GetField("darcy_flux",name_)->set_initialized();
    S_->Initialize();
  } else {
    // fields that might be initialized through the input
    // file need to be tested and initialized 'by hand' here

    if (!S_->GetField("darcy_flux",name_)->initialized()) {
      darcy_flux()->PutScalar(0.0);
      S_->GetField("darcy_flux",name_)->set_initialized();
    }
    if (!S_->GetField("water_saturation",name_)->initialized()) {
      water_saturation()->PutScalar(1.0);
      S_->GetField("water_saturation",name_)->set_initialized();
    }
    if (!S_->GetField("prev_water_saturation",name_)->initialized()) {
      prev_water_saturation()->PutScalar(1.0);
      S_->GetField("prev_water_saturation",name_)->set_initialized();
    }

    // // BEGIN REMOVE ME once flow tests pass --etc
    // S_->GetFieldData("porosity", name_)->PutScalar(0.2);
    // S_->GetFieldData("total_component_concentration", name_)->PutScalar(0.0);
    // S_->GetFieldData("water_saturation", name_)->PutScalar(1.0);
    // S_->GetFieldData("prev_water_saturation", name_)->PutScalar(1.0);
    // S_->GetFieldData("darcy_flux", name_)->PutScalar(0.0);
    // // END REMOVE ME
  }
}


int Transport_State::get_component_number(const std::string component_name)
{
  std::map<std::string,int>::const_iterator lb =
    comp_numbers_.lower_bound(component_name);
  if (lb != comp_numbers_.end() && !(comp_numbers_.key_comp()(component_name, lb->first))) {
    return lb->second;
  } else {
    return -1;
  }
}


std::string Transport_State::get_component_name(const int component_number) {
  return comp_names_[component_number];
}


/* *******************************************************************
* Interpolate linearly in time between two values v0 and v1. The time 
* is measuared relative to value v0; so that v1 is at time dT. The
* interpolated data are at time dT_int.            
******************************************************************* */
void Transport_State::InterpolateCellVector(const Epetra_Vector& v0,
        const Epetra_Vector& v1, double dT_int, double dT,
        Epetra_Vector& v_int) {
  double a = dT_int / dT;
  double b = 1.0 - a;
  v_int.Update(b, v0, a, v1, 0.);
}


/* *******************************************************************
* Debug methods
******************************************************************* */
void Transport_State::set_porosity(const double phi) {
  porosity()->PutScalar(phi);
}

void Transport_State::set_water_saturation(const double ws) {
  water_saturation()->PutScalar(ws);
}

void Transport_State::set_prev_water_saturation(const double ws) {
  prev_water_saturation()->PutScalar(ws);
}

void Transport_State::set_water_density(const double rho) {
  ref_water_density() = rho;
}

void Transport_State::set_darcy_flux(f_flux_t func, const double t) {
  const Epetra_BlockMap& fmap = darcy_flux()->Map();

  for (int f = fmap.MinLID(); f <= fmap.MaxLID(); f++) {
    const AmanziGeometry::Point& normal = mesh()->face_normal(f);
    const AmanziGeometry::Point& fc = mesh()->face_centroid(f);
    ref_darcy_flux()[f] = func(fc, t) * normal;
  }
}

void Transport_State::set_darcy_flux(const AmanziGeometry::Point& u) {
  const Epetra_BlockMap& fmap = darcy_flux()->Map();

  for (int f = fmap.MinLID(); f <= fmap.MaxLID(); f++) {
    const AmanziGeometry::Point& normal = mesh()->face_normal(f);
    ref_darcy_flux()[f] = u * normal;
  }
}

void Transport_State::set_total_component_concentration(f_conc_t func, const double t, const int ind) {
  const Epetra_BlockMap& cmap = total_component_concentration()->Map();

  for (int c = cmap.MinLID(); c <= cmap.MaxLID(); c++) {
    const AmanziGeometry::Point& xc = mesh()->cell_centroid(c);
    ref_total_component_concentration()[ind][c] = func(xc, t);
  }
}

void Transport_State::set_total_component_concentration(const double val, const int i) 
{
  const Epetra_BlockMap& cmap = total_component_concentration()->Map();

  for (int c = cmap.MinLID(); c <= cmap.MaxLID(); c++) {
    const AmanziGeometry::Point& xc = mesh()->cell_centroid(c);
    ref_total_component_concentration()[i][c] = val;
  }
}


void Transport_State::error_total_component_concentration(f_conc_t f, double t, double* L1, double* L2)
{
  int i, j, c;
  double d;
  const Epetra_BlockMap& cmap = total_component_concentration()->Map();

  *L1 = *L2 = 0.0;
  for (c = cmap.MinLID(); c <= cmap.MaxLID(); c++) {
    const AmanziGeometry::Point& xc = mesh()->cell_centroid(c);
    d = (*total_component_concentration())[0][c] - f(xc, t);
    
    double volume = mesh()->cell_volume(c);
    *L1 += fabs(d) * volume;
    *L2 += d * d * volume;
  }

  *L2 = sqrt(*L2);
}

}  // namespace AmanziTransport
}  // namespace Amanzi

