#include "Epetra_IntVector.h"

#include <algorithm>

#include "RichardsProblem.hpp"
#include "boundary-function.hh"
#include "vanGenuchtenModel.hpp"
#include "Mesh.hh"
#include "flow-bc-factory.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

namespace Amanzi
{

RichardsProblem::RichardsProblem(const Teuchos::RCP<AmanziMesh::Mesh> &mesh,
			   Teuchos::ParameterList &list) : mesh_(mesh)
{
  // Create the combined cell/face DoF map.
  dof_map_ = create_dof_map_(CellMap(), FaceMap());

  face_importer_ = new Epetra_Import(FaceMap(true),FaceMap(false));
  cell_importer_ = new Epetra_Import(CellMap(true),CellMap(false));

  // Create the MimeticHexLocal objects.
  init_mimetic_disc_(*mesh, MD);
  md_ = new MimeticHex(mesh); // evolving replacement for mimetic_hex
  
  upwind_k_rel_ = list.get<bool>("Upwind relative permeability", true);
  
  p_atm_   = list.get<double>("Atmospheric pressure");
  rho_ = list.get<double>("fluid density");
  mu_ = list.get<double>("fluid viscosity");
  gravity_ = list.get<double>("gravity");
  gvec_[0] = 0.0; gvec_[1] = 0.0; gvec_[2] = -gravity_;
  
  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList> bc_list = Teuchos::rcpFromRef(list.sublist("boundary conditions",true));
  FlowBCFactory bc_factory(mesh, bc_list);
  bc_press_ = bc_factory.CreatePressure();
  bc_head_  = bc_factory.CreateStaticHead(p_atm_, rho_, gravity_);
  bc_flux_  = bc_factory.CreateMassFlux();
  validate_boundary_conditions();

  // Create the diffusion matrix (structure only, no values)
  D_ = Teuchos::rcp<DiffusionMatrix>(create_diff_matrix_(mesh));

  // Create the preconditioner (structure only, no values)
  Teuchos::ParameterList diffprecon_plist = list.sublist("Diffusion Preconditioner");
  precon_ = new DiffusionPrecon(D_, diffprecon_plist, Map());
  
  // resize the vectors for permeability
  k_.resize(CellMap(true).NumMyElements()); 
  k_rl_.resize(CellMap(true).NumMyElements());
  
  if ( ! list.isSublist("Water retention models") )
    {
      Errors::Message m("There is no Water retention models list");
      Exceptions::amanzi_throw(m);
    }
  Teuchos::ParameterList &vGsl = list.sublist("Water retention models");

  // read the water retention model sublist and create the WRM array

  // first figure out how many entries there are
  int nblocks = 0;
  for (Teuchos::ParameterList::ConstIterator i = vGsl.begin(); i != vGsl.end(); i++)
    {
      // only count sublists
      if (vGsl.isSublist(vGsl.name(i))) 
	{
	  nblocks++;
	}
      else
	{
	  // currently we only support van Genuchten, if a user
	  // specifies something else, throw a meaningful error...
	  Errors::Message m("RichardsProblem: the Water retention models sublist contains an entry that is not a sublist!");
	  Exceptions::amanzi_throw(m);
	}
    }

  WRM.resize(nblocks);
  
  int iblock = 0;
  
  for (Teuchos::ParameterList::ConstIterator i = vGsl.begin(); i != vGsl.end(); i++) {
    if (vGsl.isSublist(vGsl.name(i))) {
   	  Teuchos::ParameterList &wrmlist = vGsl.sublist( vGsl.name(i) );
	  
	    // which water retention model are we using? currently we only have van Genuchten
	    if ( wrmlist.get<string>("Water retention model") == "van Genuchten") {
	      // read the mesh block number that this model applies to
	      //int meshblock = wrmlist.get<int>("Region ID");
	      std::string region = wrmlist.get<std::string>("Region");

	      // read values for the van Genuchten model
	      double vG_m_      = wrmlist.get<double>("van Genuchten m");
	      double vG_alpha_  = wrmlist.get<double>("van Genuchten alpha");
	      double vG_sr_     = wrmlist.get<double>("van Genuchten residual saturation");
	      
	      WRM[iblock] = Teuchos::rcp(new vanGenuchtenModel(region,vG_m_,vG_alpha_,
							       vG_sr_,p_atm_));
	    }
      iblock++;
	  }
  }

  // store the cell volumes in a convenient way 
  cell_volumes = new Epetra_Vector(mesh->cell_map(false)); 
  for (int n=0; n<(md_->CellMap()).NumMyElements(); n++) {
    (*cell_volumes)[n] = mesh_->cell_volume(n);
  }  
}


RichardsProblem::~RichardsProblem()
{
  delete dof_map_;
  delete face_importer_;
  delete cell_importer_;
  delete md_;
  delete precon_;

  delete cell_volumes;
  delete bc_press_;
  delete bc_head_;
  delete bc_flux_;
}

  
void RichardsProblem::validate_boundary_conditions() const
{
  // Create sets of the face indices belonging to each BC type.
  std::set<int> press_faces, head_faces, flux_faces;
  BoundaryFunction::Iterator bc;
  for (bc = bc_press_->begin(); bc != bc_press_->end(); ++bc) press_faces.insert(bc->first);
  for (bc = bc_head_->begin();  bc != bc_head_->end();  ++bc) head_faces.insert(bc->first);
  for (bc = bc_flux_->begin();  bc != bc_flux_->end();  ++bc) flux_faces.insert(bc->first);
  
  std::set<int> overlap;
  std::set<int>::iterator overlap_end;
  int local_overlap, global_overlap;
  
  // Check for overlap between pressure and static head BC.
  std::set_intersection(press_faces.begin(), press_faces.end(),
                         head_faces.begin(),  head_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  Comm().SumAll(&local_overlap, &global_overlap, 1); //TODO: this will over count ghost faces
  if (global_overlap != 0) {
    Errors::Message m;
    std::stringstream s;
    s << global_overlap;
    m << "RichardsProblem: \"static head\" BC overlap \"dirichlet\" BC on "
      << s.str().c_str() << " faces";
    Exceptions::amanzi_throw(m);
  }
  
  // Check for overlap between pressure and flux BC.
  overlap.clear();
  std::set_intersection(press_faces.begin(), press_faces.end(),
                         flux_faces.begin(),  flux_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  Comm().SumAll(&local_overlap, &global_overlap, 1); //TODO: this will over count ghost faces
  if (global_overlap != 0) {
    Errors::Message m;
    std::stringstream s;
    s << global_overlap;
    m << "RichardsProblem: \"flux\" BC overlap \"dirichlet\" BC on "
      << s.str().c_str() << " faces";
    Exceptions::amanzi_throw(m);
  }
  
  // Check for overlap between static head and flux BC.
  overlap.clear();
  std::set_intersection(head_faces.begin(), head_faces.end(),
                        flux_faces.begin(), flux_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  Comm().SumAll(&local_overlap, &global_overlap, 1); //TODO: this will over count ghost faces
  if (global_overlap != 0) {
    Errors::Message m;
    std::stringstream s;
    s << global_overlap;
    m << "RichardsProblem: \"flux\" BC overlap \"static head\" BC on "
      << s.str().c_str() << " faces";
    Exceptions::amanzi_throw(m);
  }
  
  //TODO: Verify that a BC has been applied to every boundary face.
  //      Right now faces without BC are considered no-mass-flux.
}


DiffusionMatrix* RichardsProblem::create_diff_matrix_(const Teuchos::RCP<AmanziMesh::Mesh> &mesh) const
{
  // Generate the list of all Dirichlet-type faces.
  // The provided list should include USED BC faces.
  std::vector<int> dir_faces;
  for (BoundaryFunction::Iterator i = bc_press_->begin(); i != bc_press_->end(); ++i)
    dir_faces.push_back(i->first);
  for (BoundaryFunction::Iterator i = bc_head_->begin();  i != bc_head_->end();  ++i)
    dir_faces.push_back(i->first);
  
  return new DiffusionMatrix(mesh, dir_faces);
}


void RichardsProblem::init_mimetic_disc_(AmanziMesh::Mesh &mesh, std::vector<MimeticHexLocal> &MD) const
{
  // Local storage for the 8 vertex coordinates of a hexahedral cell.
  double x[8][3];
  double *xBegin = &x[0][0];  // begin iterator
  double *xEnd = xBegin+24;   // end iterator

  MD.resize(mesh.cell_map(true).NumMyElements());
  for (int j = 0; j < MD.size(); ++j) {
    mesh.cell_to_coordinates((unsigned int) j, xBegin, xEnd);
    MD[j].update(x);
  }
}


void RichardsProblem::SetPermeability(double k)
{
  /// should verify k > 0
  for (int i = 0; i < k_.size(); ++i) k_[i] = k;
}


void RichardsProblem::SetPermeability(const Epetra_Vector &k)
{
  /// should verify k.Map() is CellMap()
  /// should verify k values are all > 0
  Epetra_Vector k_ovl(mesh_->cell_map(true));
  Epetra_Import importer(mesh_->cell_map(true), mesh_->cell_map(false));

  k_ovl.Import(k, importer, Insert);

  for (int i = 0; i < k_.size(); ++i) k_[i] = k_ovl[i];
}

void RichardsProblem::SetFlowState( Teuchos::RCP<const Flow_State> FS_ )
{
  FS = FS_;
}


void RichardsProblem::ComputeRelPerm(const Epetra_Vector &P, Epetra_Vector &k_rel) const
{
  ASSERT(P.Map().SameAs(CellMap(true)));
  ASSERT(k_rel.Map().SameAs(CellMap(true)));
  for (int mb=0; mb<WRM.size(); mb++) {
    // get mesh block cells
    std::string region = WRM[mb]->region();
    unsigned int ncells = mesh_->get_set_size(region,AmanziMesh::CELL,AmanziMesh::USED);
    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region,AmanziMesh::CELL,AmanziMesh::USED,&block);
    std::vector<unsigned int>::iterator j;
    for (j = block.begin(); j!=block.end(); j++) {
      k_rel[*j] = WRM[mb]->k_relative(P[*j]);
    }
  }
}


void RichardsProblem::UpdateVanGenuchtenRelativePermeability(const Epetra_Vector &P)
{
  for (int mb=0; mb<WRM.size(); mb++) 
    {
      // get mesh block cells
      std::string region = WRM[mb]->region();
      unsigned int ncells = mesh_->get_set_size(region,AmanziMesh::CELL,AmanziMesh::USED);
      std::vector<unsigned int> block(ncells);

      mesh_->get_set_entities(region,AmanziMesh::CELL,AmanziMesh::USED,&block);
      
      std::vector<unsigned int>::iterator j;
      for (j = block.begin(); j!=block.end(); j++)
	{
	  k_rl_[*j] = WRM[mb]->k_relative(P[*j]);
	}
    }
}

void RichardsProblem::dSofP(const Epetra_Vector &P, Epetra_Vector &dS)
{
  for (int mb=0; mb<WRM.size(); mb++) 
    {
      // get mesh block cells
      std::string region = WRM[mb]->region();
      unsigned int ncells = mesh_->get_set_size(region,AmanziMesh::CELL,AmanziMesh::OWNED);
      std::vector<unsigned int> block(ncells);

      mesh_->get_set_entities(region,AmanziMesh::CELL,AmanziMesh::OWNED,&block);
      
      std::vector<unsigned int>::iterator j;
      for (j = block.begin(); j!=block.end(); j++)
	{
	  dS[*j] = WRM[mb]->d_saturation(P[*j]);
	}

    }
}


void RichardsProblem::DeriveVanGenuchtenSaturation(const Epetra_Vector &P, Epetra_Vector &S)
{
  for (int mb=0; mb<WRM.size(); mb++) 
    {
      // get mesh block cells
      std::string region = WRM[mb]->region();
      unsigned int ncells = mesh_->get_set_size(region,AmanziMesh::CELL,AmanziMesh::OWNED);
      std::vector<unsigned int> block(ncells);

      mesh_->get_set_entities(region,AmanziMesh::CELL,AmanziMesh::OWNED,&block);
      
      std::vector<unsigned int>::iterator j;
      for (j = block.begin(); j!=block.end(); j++)
	{
	  S[*j] = WRM[mb]->saturation(P[*j]);
	}

    }
}


void RichardsProblem::ComputePrecon(const Epetra_Vector &P)
{
  std::vector<double> K(k_);

  Epetra_Vector &Pcell_own = *CreateCellView(P);
  Epetra_Vector Pcell(CellMap(true));
  Pcell.Import(Pcell_own, *cell_importer_, Insert);

  Epetra_Vector &Pface_own = *CreateFaceView(P);
  Epetra_Vector Pface(FaceMap(true));
  Pface.Import(Pface_own, *face_importer_, Insert);

  // Fill the diffusion matrix with values.
  if (upwind_k_rel_) {
    Epetra_Vector K_upwind(Pface.Map());
    ComputeUpwindRelPerm(Pcell, Pface, K_upwind);
    for (int j = 0; j < K.size(); ++j) K[j] = rho_ * k_[j] / mu_;
    D_->Compute(K, K_upwind);
  } else {
    Epetra_Vector k_rel(Pcell.Map());
    ComputeRelPerm(Pcell, k_rel);
    for (int j = 0; j < K.size(); ++j) K[j] = rho_ * k_[j] * k_rel[j] / mu_;
    D_->Compute(K);
  }

  // Compute the face Schur complement of the diffusion matrix.
  D_->ComputeFaceSchur();

  // Compute the preconditioner from the newly computed diffusion matrix and Schur complement.
  precon_->Compute();

  delete &Pcell_own, &Pface_own;
}



void RichardsProblem::ComputePrecon(const Epetra_Vector &P, const double h)
{
  std::vector<double> K(k_);

  Epetra_Vector &Pcell_own = *CreateCellView(P);
  Epetra_Vector Pcell(CellMap(true));
  Pcell.Import(Pcell_own, *cell_importer_, Insert);

  Epetra_Vector &Pface_own = *CreateFaceView(P);
  Epetra_Vector Pface(FaceMap(true));
  Pface.Import(Pface_own, *face_importer_, Insert);
  
  if (upwind_k_rel_) {
    Epetra_Vector K_upwind(Pface.Map());
    ComputeUpwindRelPerm(Pcell, Pface, K_upwind);
    for (int j = 0; j < K.size(); ++j) K[j] = rho_ * k_[j] / mu_;
    D_->Compute(K, K_upwind);
  } else {
    Epetra_Vector k_rel(Pcell.Map());
    ComputeRelPerm(Pcell, k_rel);
    for (int j = 0; j < K.size(); ++j) K[j] = rho_ * k_[j] * k_rel[j] / mu_;
    D_->Compute(K);
  }

  // add the time derivative to the diagonal
  Epetra_Vector celldiag(CellMap(false));
  dSofP(Pcell_own, celldiag);
 
  // get the porosity
  const Epetra_Vector& phi = FS->porosity();

  celldiag.Multiply(rho_,celldiag,phi,0.0);
 
  celldiag.Multiply(1.0/h,celldiag,*cell_volumes,0.0);
  
  D_->add_to_celldiag(celldiag);

  // Compute the face Schur complement of the diffusion matrix.
  D_->ComputeFaceSchur();

  // Compute the preconditioner from the newly computed diffusion matrix and Schur complement.
  precon_->Compute();

  delete &Pcell_own, &Pface_own;
}

void RichardsProblem::ComputeF(const Epetra_Vector &X, Epetra_Vector &F)
{
  ComputeF(X,F,0.0);
}

void RichardsProblem::ComputeF(const Epetra_Vector &X, Epetra_Vector &F, double time)
{
  // The cell and face-based DoF are packed together into the X and F Epetra
  // vectors: cell-based DoF in the first part, followed by the face-based DoF.
  // In addition, only the owned DoF belong to the vectors.


  // Create views into the cell and face segments of X and F
  Epetra_Vector &Pcell_own = *CreateCellView(X);
  Epetra_Vector &Pface_own = *CreateFaceView(X);

  Epetra_Vector &Fcell_own = *CreateCellView(F);
  Epetra_Vector &Fface_own = *CreateFaceView(F);

  // Create input cell and face pressure vectors that include ghosts.
  Epetra_Vector Pcell(CellMap(true));
  Pcell.Import(Pcell_own, *cell_importer_, Insert);
  Epetra_Vector Pface(FaceMap(true));
  Pface.Import(Pface_own, *face_importer_, Insert);

  // Precompute the Dirichlet-type BC residuals for later use.
  // Impose the Dirichlet boundary data on the face pressure vector.
  bc_press_->Compute(time);
  std::map<int,double> Fpress;
  for (BoundaryFunction::Iterator bc = bc_press_->begin(); bc != bc_press_->end(); ++bc) {
    Fpress[bc->first] = Pface[bc->first] - bc->second;
    Pface[bc->first] = bc->second;
  }
  bc_head_->Compute(time);
  std::map<int,double> Fhead;
  for (BoundaryFunction::Iterator bc = bc_head_->begin(); bc != bc_head_->end(); ++bc) {
    Fhead[bc->first] = Pface[bc->first] - bc->second;
    Pface[bc->first] = bc->second;
  }
  
  // GENERIC COMPUTATION OF THE FUNCTIONAL /////////////////////////////////////

  Epetra_Vector K(CellMap(true)), K_upwind(FaceMap(true));
  if (upwind_k_rel_) {
    ComputeUpwindRelPerm(Pcell, Pface, K_upwind);
    for (int j = 0; j < K.MyLength(); ++j) K[j] = (rho_ * k_[j] / mu_);
  } else {
    ComputeRelPerm(Pcell, K);
    for (int j = 0; j < K.MyLength(); ++j) K[j] = (rho_ * k_[j] * K[j] / mu_);
  }

  // Create cell and face result vectors that include ghosts.
  Epetra_Vector Fcell(CellMap(true));
  Epetra_Vector Fface(FaceMap(true));

  int cface[6];
  double aux1[6], aux2[6], aux3[6], gflux[6];

  Fface.PutScalar(0.0);
  for (int j = 0; j < Pcell.MyLength(); ++j) {
    // Get the list of process-local face indices for this cell.
    mesh_->cell_to_faces((unsigned int) j, (unsigned int*) cface, (unsigned int*) cface+6);
    // Gather the local face pressures int AUX1.
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    // Compute the local value of the diffusion operator.
    if (upwind_k_rel_) {
      for (int k = 0; k < 6; ++k) aux3[k] = K_upwind[cface[k]];
      MD[j].diff_op(K[j], aux3, Pcell[j], aux1, Fcell[j], aux2);
      // Gravity contribution
      MD[j].GravityFlux(gvec_, gflux);
      for (int k = 0; k < 6; ++k) gflux[k] *= rho_ * K[j] * aux3[k];
      for (int k = 0; k < 6; ++k) {
        aux2[k] -= gflux[k];
        Fcell[j] += gflux[k];
      }
    } else {
      MD[j].diff_op(K[j], Pcell[j], aux1, Fcell[j], aux2);
      // Gravity contribution
      MD[j].GravityFlux(gvec_, gflux);
      for (int k = 0; k < 6; ++k) aux2[k] -= rho_ * K[j] * gflux[k];
    }
    // Scatter the local face result into FFACE.
    for (int k = 0; k < 6; ++k) Fface[cface[k]] += aux2[k];
  }

  // Dirichlet-type condition residuals; overwrite with the pre-computed values.
  std::map<int,double>::const_iterator i;
  for (i = Fpress.begin(); i != Fpress.end(); ++i) Fface[i->first] = i->second;
  for (i = Fhead.begin();  i != Fhead.end();  ++i) Fface[i->first] = i->second;
  
  // Mass flux BC contribution.
  bc_flux_->Compute(time);
  for (BoundaryFunction::Iterator bc = bc_flux_->begin(); bc != bc_flux_->end(); ++bc)
    Fface[bc->first] += bc->second * md_->face_area_[bc->first];
  
  // Copy owned part of result into the output vectors.
  for (int j = 0; j < Fcell_own.MyLength(); ++j) Fcell_own[j] = Fcell[j];
  for (int j = 0; j < Fface_own.MyLength(); ++j) Fface_own[j] = Fface[j];
    
  delete &Pcell_own, &Pface_own, &Fcell_own, &Fface_own;
}


void RichardsProblem::ComputeUpwindRelPerm(const Epetra_Vector &Pcell,
    const Epetra_Vector &Pface, Epetra_Vector &k_rel) const
{
  ASSERT(Pcell.Map().SameAs(CellMap(true)));
  ASSERT(Pface.Map().SameAs(FaceMap(true)));
  ASSERT(k_rel.Map().SameAs(FaceMap(true)));

  int fdirs[6];
  unsigned int cface[6];
  double aux1[6], aux2[6], gflux[6], dummy;

  // Calculate the mimetic face 'fluxes' that omit the relative permeability.
  // When P is a converged solution, the following result on interior faces
  // will be twice the true value (double counting).  For a non-converged
  // solution (typical case) we view it as an approximation (twice the
  // average of the adjacent cell fluxes).  Here we are only interested
  // in the sign of the flux, and then only on interior faces, so we omit
  // bothering with the scaling.

  // Looping over all cells gives desired result on *owned* faces.
  Epetra_Vector Fface(FaceMap(true)); // fills with 0, includes ghosts
  for (unsigned int j = 0; j < Pcell.MyLength(); ++j) {
    // Get the list of process-local face indices for this cell.
    mesh_->cell_to_faces(j, cface, cface+6);
    // Gather the local face pressures int AUX1.
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    // Compute the local value of the diffusion operator.
    double K = (rho_ * k_[j] / mu_);
    MD[j].diff_op(K, Pcell[j], aux1, dummy, aux2); // aux2 is inward flux
    // Gravity contribution; aux2 becomes outward flux
    MD[j].GravityFlux(gvec_, gflux);
    for (int k = 0; k < 6; ++k) aux2[k] = rho_ * K * gflux[k] - aux2[k];
    // Scatter the local face result into FFACE.
    mesh_->cell_to_face_dirs(j, fdirs, fdirs+6);
    for (int k = 0; k < 6; ++k) Fface[cface[k]] += fdirs[k] * aux2[k];
  }

  // Generate the face-to-cell data structure.  For each face, the pair of
  // adjacent cell indices is stored, accessed via the members first and
  // second.  The face is oriented outward with respect to the first cell
  // of the pair and inward with respect to the second.  Boundary faces are
  // identified by a -1 value for the second member.

  typedef std::pair<int,int> CellPair;
  int nface = FaceMap(true).NumMyElements();
  int nface_own = FaceMap(false).NumMyElements();
  CellPair *fcell = new CellPair[nface];
  for (int j = 0; j < nface; ++j) {
    fcell[j].first  = -1;
    fcell[j].second = -1;
  }
  // Looping over all cells gives desired result on *owned* faces.
  for (unsigned int j = 0; j < CellMap(true).NumMyElements(); ++j) {
    // Get the list of process-local face indices for this cell.
    mesh_->cell_to_faces(j, cface, cface+6);
    mesh_->cell_to_face_dirs(j, fdirs, fdirs+6);
    for (int k = 0; k < 6; ++k) {
      if (fdirs[k] > 0) {
        ASSERT(fcell[cface[k]].first == -1);
        fcell[cface[k]].first = j;
      } else {
        ASSERT(fcell[cface[k]].second == -1);
        fcell[cface[k]].second = j;
      }
    }
  }
  // Fix-up at boundary faces: move the cell index to the first position.
  for (int j = 0; j < nface_own; ++j) {
    if (fcell[j].first == -1) {
      ASSERT(fcell[j].second != -1);
      fcell[j].first = fcell[j].second;
      fcell[j].second = -1; // marks a boundary face
    }
  }
  
  // Compute the relative permeability on cells and then upwind on owned faces.
  Epetra_Vector k_rel_cell(Pcell.Map());
  ComputeRelPerm(Pcell, k_rel_cell);
  for (int j = 0; j < nface_own; ++j) {
    if (fcell[j].second == -1) // boundary face
      k_rel[j] = k_rel_cell[fcell[j].first];
    else
      k_rel[j] = k_rel_cell[((Fface[j] >= 0.0) ? fcell[j].first : fcell[j].second)];
  }

  // Communicate the values computed above to the ghosts.
  double *k_rel_data;
  k_rel.ExtractView(&k_rel_data);
  Epetra_Vector k_rel_own(View, FaceMap(false), k_rel_data);
  k_rel.Import(k_rel_own, *face_importer_, Insert);
  
  delete [] fcell;
}


Epetra_Map* RichardsProblem::create_dof_map_(const Epetra_Map &cell_map, const Epetra_Map &face_map) const
{
  // Create the combined cell/face DoF map
  int ncell_tot = cell_map.NumGlobalElements();
  int ndof_tot = ncell_tot + face_map.NumGlobalElements();
  int ncell = cell_map.NumMyElements();
  int ndof = ncell + face_map.NumMyElements();
  int *gids = new int[ndof];
  cell_map.MyGlobalElements(&(gids[0]));
  face_map.MyGlobalElements(&(gids[ncell]));
  for (int i = ncell; i < ndof; ++i) gids[i] += ncell_tot;
  Epetra_Map *map = new Epetra_Map(ndof_tot, ndof, gids, 0, cell_map.Comm());
  delete [] gids;
  return map;
}


Epetra_Vector* RichardsProblem::CreateCellView(const Epetra_Vector &X) const
{
  // should verify that X.Map() is the same as Map()
  double *data;
  X.ExtractView(&data);
  return new Epetra_Vector(View, CellMap(), data);
}


Epetra_Vector* RichardsProblem::CreateFaceView(const Epetra_Vector &X) const
{
  // should verify that X.Map() is the same as Map()
  double *data;
  X.ExtractView(&data);
  int ncell = CellMap().NumMyElements();
  return new Epetra_Vector(View, FaceMap(), data+ncell);
}


void RichardsProblem::DeriveDarcyVelocity(const Epetra_Vector &X, Epetra_MultiVector &Q) const
{
  ASSERT(X.Map().SameAs(Map()));
  ASSERT(Q.Map().SameAs(CellMap(false)));
  ASSERT(Q.NumVectors() == 3);

  // Cell pressure vectors without and with ghosts.
  Epetra_Vector &Pcell_own = *CreateCellView(X);
  Epetra_Vector Pcell(CellMap(true));
  Pcell.Import(Pcell_own, *cell_importer_, Insert);

  // Face pressure vectors without and with ghosts.
  Epetra_Vector &Pface_own = *CreateFaceView(X);
  Epetra_Vector Pface(FaceMap(true));
  Pface.Import(Pface_own, *face_importer_, Insert);

  // The pressure gradient coefficient split as K * K_upwind: K gets
  // incorporated into the mimetic discretization and K_upwind is applied
  // afterwards to the computed mimetic flux.  Note that density (rho) is
  // not included here because we want the Darcy velocity and not mass flux.
  Epetra_Vector K(CellMap(true)), K_upwind(FaceMap(true));
  if (upwind_k_rel_) {
    ComputeUpwindRelPerm(Pcell, Pface, K_upwind);
    for (int j = 0; j < K.MyLength(); ++j) K[j] = (k_[j] / mu_);
  } else {
    ComputeRelPerm(Pcell, K);
    for (int j = 0; j < K.MyLength(); ++j) K[j] = (k_[j] / mu_) * K[j];
  }

  int cface[6];
  double aux1[6], aux2[6], aux3[6], aux4[3], gflux[6], dummy;

  for (int j = 0; j < Pcell_own.MyLength(); ++j) {
    mesh_->cell_to_faces((unsigned int) j, (unsigned int*) cface, (unsigned int*) cface+6);
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    if (upwind_k_rel_) {
      for (int k = 0; k < 6; ++k) aux3[k] = K_upwind[cface[k]];
      MD[j].diff_op(K[j], aux3, Pcell_own[j], aux1, dummy, aux2);
      MD[j].GravityFlux(gvec_, gflux);
      for (int k = 0; k < 6; ++k) aux2[k] = aux3[k] * K[j] * rho_ * gflux[k] - aux2[k];
    } else {
      MD[j].diff_op(K[j], Pcell_own[j], aux1, dummy, aux2);
      MD[j].GravityFlux(gvec_, gflux);
      for (int k = 0; k < 6; ++k) aux2[k] = K[j] * rho_ * gflux[k] - aux2[k];
    }
    MD[j].CellFluxVector(aux2, aux4);
    Q[0][j] = aux4[0];
    Q[1][j] = aux4[1];
    Q[2][j] = aux4[2];
  }

  delete &Pcell_own, &Pface_own;
}


void RichardsProblem::DeriveDarcyFlux(const Epetra_Vector &P, Epetra_Vector &F, double &l2_error) const
{
  /// should verify P.Map() is Map()
  /// should verify F.Map() is FaceMap()

  int fdirs[6];
  unsigned int cface[6];
  double aux1[6], aux2[6], aux3[6], gflux[6], dummy;

  // Create a view into the cell pressure segment of P.
  Epetra_Vector &Pcell_own = *CreateCellView(P);
  Epetra_Vector Pcell(CellMap(true));
  Pcell.Import(Pcell_own, *cell_importer_, Insert);

  // Create a copy of the face pressure segment of P that includes ghosts.
  Epetra_Vector &Pface_own = *CreateFaceView(P);
  Epetra_Vector Pface(FaceMap(true));
  Pface.Import(Pface_own, *face_importer_, Insert);

  // Create face flux and face count vectors that include ghosts.
  Epetra_Vector Fface(FaceMap(true)); // fills with 0
  Epetra_Vector error(FaceMap(true)); // fills with 0
  Epetra_IntVector count(FaceMap(true)); // fills with 0

  // The pressure gradient coefficient split as K * K_upwind: K gets
  // incorporated into the mimetic discretization and K_upwind is applied
  // afterwards to the computed mimetic flux.
  Epetra_Vector K(CellMap(true)), K_upwind(FaceMap(true));
  if (upwind_k_rel_) {
    ComputeUpwindRelPerm(Pcell, Pface, K_upwind);
    for (int j = 0; j < K.MyLength(); ++j) K[j] = (rho_ * k_[j] / mu_);
  } else {
    ComputeRelPerm(Pcell, K);
    for (int j = 0; j < K.MyLength(); ++j) K[j] = (rho_ * k_[j] * K[j] / mu_);
    K_upwind.PutScalar(1.0);  // whole coefficient used in mimetic disc
  }

  // Process-local assembly of the mimetic face fluxes.
  for (unsigned int j = 0; j < Pcell_own.MyLength(); ++j) {
    // Get the list of process-local face indices for this cell.
    mesh_->cell_to_faces(j, cface, cface+6);
    // Gather the local face pressures int AUX1.
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    // Compute the local value of the diffusion operator.
    MD[j].diff_op(K[j], Pcell_own[j], aux1, dummy, aux2);
    // Gravity contribution
    MD[j].GravityFlux(gvec_, gflux);
    for (int k = 0; k < 6; ++k) aux2[k] = K[j] * rho_ * gflux[k] - aux2[k];
    mesh_->cell_to_face_dirs(j, fdirs, fdirs+6);
    // Scatter the local face result into FFACE.
    for (int k = 0; k < 6; ++k) {
      Fface[cface[k]] += fdirs[k] * aux2[k];
      error[cface[k]] += aux2[k]; // sums to the flux discrepancy
      count[cface[k]]++;
    }
  }

  // Global assembly of face mass fluxes into the return vector.
  F.Export(Fface, *face_importer_, Add);

  // Create an owned face count vector that overlays the count vector with ghosts.
  int *count_data;
  count.ExtractView(&count_data);
  Epetra_IntVector count_own(View, FaceMap(false), count_data);

  // Final global assembly of the face count vector.
  count_own.Export(count, *face_importer_, Add);

  // Correct the double counting of fluxes on interior faces and convert
  // the mimetic flux to Darcy flux by dividing by the constant fluid density
  // and multiplying by the K_upwind on faces.
  for (int j = 0; j < F.MyLength(); ++j)
    F[j] = K_upwind[j] * F[j] / (rho_ * count[j]);

  // Create an owned face error vector that overlays the error vector with ghosts.
  double *error_data;
  error.ExtractView(&error_data);
  Epetra_Vector error_own(View, FaceMap(false), error_data);

  // Final global assembly of the flux discrepancy vector.
  error_own.Export(error, *face_importer_, Add);

  // Set the flux discrepancy error to 0 on boundary faces where there was only one flux computed.
  for (int j = 0; j < F.MyLength(); ++j) {
    error_own[j] = K_upwind[j] * error_own[j];
    if (count[j] == 1) error_own[j] = 0.0;
  }

  // Compute the norm of the flux discrepancy.
  error_own.Norm2(&l2_error);

  delete &Pcell_own, &Pface_own;
}


void RichardsProblem::Compute_udot(const double t, const Epetra_Vector& u, Epetra_Vector &udot)
{
  ComputeF(u,udot);

  // zero out the face part
  Epetra_Vector *udot_face = CreateFaceView(udot);
  udot_face->PutScalar(0.0);

}



void RichardsProblem::SetInitialPressureProfileCells(double height, Epetra_Vector *pressure)
{
  for (int j=0; j<pressure->MyLength(); j++)
    {
      std::vector<double> coords;
      coords.resize(24);
      FS->mesh()->cell_to_coordinates(j, coords.begin(), coords.end());
      
      // average the x coordinates
      double zavg = 0.0;
      for (int k=2; k<24; k+=3)
	zavg += coords[k];
      zavg /= 8.0;

      (*pressure)[j] = p_atm_ + rho_ * gvec_[2] * ( zavg - height );
    }
}


void RichardsProblem::SetInitialPressureProfileFaces(double height, Epetra_Vector *pressure)
{
  for (int j=0; j<pressure->MyLength(); j++)
    {
      std::vector<double> coords;
      coords.resize(12);
      FS->mesh()->face_to_coordinates(j, coords.begin(), coords.end());

      // average the x coordinates
      double zavg = 0.0;
      for (int k=2; k<12; k+=3)
	zavg += coords[k];
      zavg /= 4.0;

      (*pressure)[j] = p_atm_ + rho_*gvec_[2] * ( zavg - height );
    }
}

void RichardsProblem::SetInitialPressureProfileFromSaturationCells(double saturation, Epetra_Vector *pressure)
{
  for (int mb=0; mb<WRM.size(); mb++) 
    {
      // get mesh block cells
      std::string region = WRM[mb]->region();
      unsigned int ncells = mesh_->get_set_size(region,AmanziMesh::CELL,AmanziMesh::OWNED);
      std::vector<unsigned int> block(ncells);

      mesh_->get_set_entities(region,AmanziMesh::CELL,AmanziMesh::OWNED,&block);
      
      std::vector<unsigned int>::iterator j;
      for (j = block.begin(); j!=block.end(); j++)
	{
	  (*pressure)[*j] = WRM[mb]->pressure(saturation);
	}
    }  
}
void RichardsProblem::SetInitialPressureProfileFromSaturationFaces(double saturation, Epetra_Vector *pressure)
{
  for (int mb=0; mb<WRM.size(); mb++) 
    {
      // get mesh block cells
      std::string region = WRM[mb]->region();

      unsigned int ncells = mesh_->get_set_size(region,AmanziMesh::CELL,AmanziMesh::OWNED);
      std::vector<unsigned int> block(ncells);

      mesh_->get_set_entities(region,AmanziMesh::CELL,AmanziMesh::OWNED,&block);
      
      std::vector<unsigned int>::iterator j;
      for (j = block.begin(); j!=block.end(); j++)
	{
	  (*pressure)[*j] = WRM[mb]->pressure(saturation);
	}
    }
}

} // close namespace Amanzi
