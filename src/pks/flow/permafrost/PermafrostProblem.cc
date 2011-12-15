/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "Epetra_IntVector.h"

#include <algorithm>
#include <iostream>

#include "PermafrostProblem.hh"
#include "boundary-function.hh"
#include "vanGenuchtenModel.hh"
#include "Mesh.hh"
#include "flow-bc-factory.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

namespace Amanzi {

PermafrostProblem::PermafrostProblem(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
        Teuchos::ParameterList& permafrost_plist) : mesh_(mesh) {
  // Create the combined cell/face DoF map.
  dof_map_ = create_dof_map_(CellMap(), FaceMap());

  face_importer_ = Teuchos::rcp(new Epetra_Import(FaceMap(true),FaceMap(false)));
  cell_importer_ = Teuchos::rcp(new Epetra_Import(CellMap(true),CellMap(false)));

  // Create the MimeticHexLocal objects.
  init_mimetic_disc_(mesh_, MD_);
  md_ = new MimeticHex(mesh); // evolving replacement for mimetic_hex

  // work vectors for within the problem
  pressure_gas_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(true)));
  sat_gas_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(true)));
  sat_liquid_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(true)));
  sat_ice_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(true)));
  rho_gas_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(true)));
  rho_liquid_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(true)));
  rho_ice_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(true)));
  mu_liquid_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(true)));
  k_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(true)));
  k_rel_liquid_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(true)));

  // store the cell volumes and mean height in a convenient way
  cell_volumes_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(false)));
  DeriveCellVolumes(cell_volumes_);
  cell_heights_ = Teuchos::rcp(new Epetra_Vector(mesh->cell_map(false)));
  DeriveCellHeights(cell_heights_);
};

void PermafrostProblem::InitializeProblem(Teuchos::ParameterList& permafrost_plist) {
  // get problem-specific parameters
  p_atm_ = permafrost_plist.get<double>("Atmospheric pressure");
  upwind_k_rel_ = permafrost_plist.get<bool>("Upwind relative permeability", true);

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList> bc_list = Teuchos::rcpFromRef(permafrost_plist.sublist("boundary conditions",true));
  FlowBCFactory bc_factory(mesh_, bc_list);
  bc_press_ = Teuchos::rcp(bc_factory.CreatePressure());
  bc_head_  = Teuchos::rcp(bc_factory.CreateStaticHead(p_atm_, rho_, gravity_));
  bc_flux_  = Teuchos::rcp(bc_factory.CreateMassFlux());
  validate_boundary_conditions_();

  // Create the diffusion matrix (structure only, no values)
  D_ = Teuchos::rcp<DiffusionMatrix>(create_diff_matrix_(mesh_));

  // Create the preconditioner (structure only, no values)
  Teuchos::ParameterList diffprecon_plist = permafrost_plist.sublist("Diffusion Preconditioner");
  precon_ = new DiffusionPrecon(D_, diffprecon_plist, Map());

  // parse sublists to create saturations curves
  SaturationCurveFactory sat_curve_factory;
  Teuchos::ParameterList &satcurve_plist = permafrost_plist.sublist("Saturation Curves");
  sat_curve_factory.parse_saturation_curves(satcurve_plist, sat_curves_);

  // create EOS curves
  EOSFactory eos_factory;
  Teuchos::ParameterList &gas_eos_plist = permafrost_plist.sublist("Gas EOS");
  gas_eos_ = Teuchos::rcp(eos_factory.create_eos(gas_eos_plist));
  Teuchos::ParameterList &liquid_eos_plist = permafrost_plist.sublist("Liquid EOS");
  liquid_eos_ = Teuchos::rcp(eos_factory.create_eos(liquid_eos_plist));
  Teuchos::ParameterList &ice_eos_plist = permafrost_plist.sublist("Ice EOS");
  ice_eos_ = Teuchos::rcp(eos_factory.create_eos(ice_eos_plist));
};

PermafrostProblem::~PermafrostProblem() {
  delete md_;
  delete precon_;
};

// private methods
Teuchos::RCP<Epetra_Map> PermafrostProblem::create_dof_map_(const Epetra_Map &cell_map,
        const Epetra_Map &face_map) const {
  // Create the combined cell/face DoF map
  int ncell_tot = cell_map.NumGlobalElements();
  int ndof_tot = ncell_tot + face_map.NumGlobalElements();
  int ncell = cell_map.NumMyElements();
  int ndof = ncell + face_map.NumMyElements();
  int *gids = new int[ndof];
  cell_map.MyGlobalElements(&(gids[0]));
  face_map.MyGlobalElements(&(gids[ncell]));
  for (int i = ncell; i < ndof; ++i) gids[i] += ncell_tot;
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(ndof_tot, ndof, gids, 0, cell_map.Comm()));
  delete [] gids;
  return map;
}

void PermafrostProblem::validate_boundary_conditions_() const {
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
    m << "PermafrostProblem: \"static head\" BC overlap \"dirichlet\" BC on "
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
    m << "PermafrostProblem: \"flux\" BC overlap \"dirichlet\" BC on "
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
    m << "PermafrostProblem: \"flux\" BC overlap \"static head\" BC on "
      << s.str().c_str() << " faces";
    Exceptions::amanzi_throw(m);
  }

  //TODO: Verify that a BC has been applied to every boundary face.
  //      Right now faces without BC are considered no-mass-flux.
};

DiffusionMatrix* PermafrostProblem::create_diff_matrix_(
        Teuchos::RCP<AmanziMesh::Mesh> &mesh) const {
  // Generate the list of all Dirichlet-type faces.
  // The provided list should include USED BC faces.
  std::vector<int> dir_faces;
  for (BoundaryFunction::Iterator i = bc_press_->begin(); i != bc_press_->end(); ++i)
    dir_faces.push_back(i->first);
  for (BoundaryFunction::Iterator i = bc_head_->begin();  i != bc_head_->end();  ++i)
    dir_faces.push_back(i->first);

  return new DiffusionMatrix(mesh, dir_faces);
}

void PermafrostProblem::init_mimetic_disc_(Teuchos::RCP<AmanziMesh::Mesh> &mesh,
        std::vector<MimeticHexLocal> &MD) const {
  // Local storage for the 8 vertex coordinates of a hexahedral cell.
  double x[8][3];
  double *xBegin = &x[0][0];  // begin iterator
  double *xEnd = xBegin+24;   // end iterator

  MD.resize(mesh->cell_map(true).NumMyElements());
  for (int j = 0; j < MD.size(); ++j) {
    mesh->cell_to_coordinates((unsigned int) j, xBegin, xEnd);
    MD[j].update(x);
  }
};

// Set Methods
void PermafrostProblem::SetGravity(const double g[3]) {
  if (g[0] != 0. || g[1] != 0.) {
    Errors::Message message("PermafrostProblem currently requires gravity g = g_z");
    Exceptions::amanzi_throw(message);
  } else {
    (*gvec_)[0] = 0.;
    (*gvec_)[1] = 0.;
    (*gvec_)[2] = g[2];
    gravity_ = -g[2];
  }
};

void PermafrostProblem::SetGravity(double g) {
  (*gvec_)[0] = 0.;
  (*gvec_)[1] = 0.;
  (*gvec_)[2] = -g;
  gravity_ = g;
};

// Derive methods
void PermafrostProblem::DeriveCellVolumes(Teuchos::RCP<Epetra_Vector> &cell_volumes) const {
  ASSERT(cell_volumes->Map().SameAs(md_->CellMap()));
  // cell volumes must be non-ghosted
  for (int n=0; n<(md_->CellMap()).NumMyElements(); n++) {
    (*cell_volumes)[n] = md_->Volume(n);
  }
};

void PermafrostProblem::DeriveCellHeights(Teuchos::RCP<Epetra_Vector> &cell_heights) const {
  ASSERT(cell_heights->Map().SameAs(md_->CellMap()));
  // cell heights must be non-ghosted
  for (int n=0; n<(md_->CellMap()).NumMyElements(); n++) {
    std::vector<double> coords;
    coords.resize(24);
    mesh_->cell_to_coordinates(n, coords.begin(), coords.end());

    double zavg = 0.0;
    for (int i=2; i<24; i+=3) {
      zavg += coords[i];
    }
    (*cell_heights)[n] = zavg / 8.0;
  }
};

void PermafrostProblem::DerivePermeability(const Epetra_Vector &phi,
                                           Teuchos::RCP<Epetra_Vector> &perm) const {
  ASSERT(phi.Map().SameAs(perm->Map()));
  // currently do nothing... this might get set eventually, but for now
  // assuming fixed permeability
};

void PermafrostProblem::DeriveLiquidRelPerm(const Epetra_Vector &saturation,
        Teuchos::RCP<Epetra_Vector> &rel_perm) const {
  ASSERT(saturation.Map().SameAs(rel_perm->Map()));

  for (int mb=0; mb<sat_curves_.size(); ++mb) {
    // get mesh block cells
    unsigned int mb_id = sat_curves_[mb]->mesh_block();
    unsigned int ncells = mesh_->get_set_size(mb_id,AmanziMesh::CELL,AmanziMesh::USED);
    std::vector<unsigned int> block(ncells);

    mesh_->get_set(mb_id,AmanziMesh::CELL,AmanziMesh::USED,block.begin(),block.end());

    std::vector<unsigned int>::iterator j;
    for (j = block.begin(); j!=block.end(); ++j) {
      (*rel_perm)[*j] = sat_curves_[mb]->k_relative(saturation[*j])
	}
  }
};

void PermafrostProblem::DeriveGasDensity(const Epetra_Vector &pressure,
        const Epetra_Vector &temp, Teuchos::RCP<Epetra_Vector> &rho) const {
  ASSERT(pressure.Map().SameAs(temp.Map()));
  ASSERT(pressure.Map().SameAs(rho->Map()));
  for (int i=0; i < (md_->CellMap()).NumMyElements(); ++i) {
    (*rho)[i] = gas_eos_->density(pressure[i], temp[i]);
  }
};

void PermafrostProblem::DeriveLiquidDensity(const Epetra_Vector &pressure,
        const Epetra_Vector &temp, Teuchos::RCP<Epetra_Vector> &rho) const {
  ASSERT(pressure.Map().SameAs(temp.Map()));
  ASSERT(pressure.Map().SameAs(rho->Map()));
  for (int i=0; i < (md_->CellMap()).NumMyElements(); ++i) {
    (*rho)[i] = liquid_eos_->density(pressure[i], temp[i]);
  }
};

void PermafrostProblem::DeriveIceDensity(const Epetra_Vector &pressure,
                                         const Epetra_Vector &temp,
                                         Teuchos::RCP<Epetra_Vector> &rho) const {
  ASSERT(pressure.Map().SameAs(temp.Map()));
  ASSERT(pressure.Map().SameAs(rho->Map()));
  for (int i=0; i < (md_->CellMap()).NumMyElements(); ++i) {
    (*rho)[i] = ice_eos_->density(pressure[i], temp[i]);
  }
};

void PermafrostProblem::DeriveLiquidViscosity(const Epetra_Vector &pressure,
        const Epetra_Vector &temp, Teuchos::RCP<Epetra_Vector> &viscosity)const {
  ASSERT(pressure.Map().SameAs(temp.Map()));
  ASSERT(pressure.Map().SameAs(viscosity->Map()));
  for (int i=0; i < (md_->CellMap()).NumMyElements(); ++i) {
    (*rho)[i] = liquid_eos_->viscosity(pressure[i], temp[i]);
  }
};

void PermafrostProblem::DeriveSaturation(const Epetra_Vector &pressure_gas,
        const Epetra_Vector &pressure_liquid, const Epetra_Vector &temp,
        const Epetra_Vector &density_ice, Teuchos::RCP<Epetra_Vector> &sat_gas,
        Teuchos::RCP<Epetra_Vector> &sat_liquid,
        Teuchos::RCP<Epetra_Vector> &sat_ice) const {
  ASSERT(pressure_gas.Map().SameAs(pressure_liquid.Map()));
  ASSERT(pressure_gas.Map().SameAs(temp.Map()));
  ASSERT(pressure_gas.Map().SameAs(density_ice.Map()));
  ASSERT(pressure_gas.Map().SameAs(sat_gas->Map()));
  ASSERT(pressure_gas.Map().SameAs(sat_liquid->Map()));
  ASSERT(pressure_gas.Map().SameAs(sat_ice->Map()));

  for (int mb=0; mb<sat_curves_.size(); ++mb) {
    // get mesh block cells
    unsigned int mb_id = sat_curves_[mb]->mesh_block();
    unsigned int ncells = mesh_->get_set_size(mb_id,AmanziMesh::CELL,AmanziMesh::USED);
    std::vector<unsigned int> block(ncells);

    mesh_->get_set(mb_id,AmanziMesh::CELL,AmanziMesh::USED,block.begin(),block.end());

    std::vector<unsigned int>::iterator j;
    for (j = block.begin(); j!=block.end(); ++j) {
      // workA and workB correspond to A and B notation in permafrost notes, 1.2.2
      double workA = 1.0/sat_curves_[mb]->s_star(sigma_gl_/sigma_il_*h_iw0 *
                                                density_ice[j]*(T0_-temp[j])/T0_);
      double workB = 1.0/sat_curves_[mb]->s_star(pressure_gas[j] - pressure_liquid[j]);

      double sat_l = 1.0/(workA + workB - 1.0);
      (*sat_liquid)[j] = sat_l;
      (*sat_gas)[j] = sat_l*(workB - 1.0);
      (*sat_ice)[j] = sat_l*(workA - 1.0);
	}
  }
};

void PermafrostProblem::DeriveLiquidFlux(const Epetra_Vector &pressure,
        const Epetra_Vector &rho, const Epetra_Vector &k_rel, const Epetra_Vector &k,
        const Epetra_Vector &mu, Teuchos::RCP<Epetra_Vector> &flux,double &l2_error) const {
  ASSERT(pressure.Map().SameAs(Map())); // cells and faces
  ASSERT(rho.Map().SameAs(CellMap(true))); // ghosted cells
  ASSERT(rho.Map().SameAs(k_rel.Map()));
  ASSERT(rho.Map().SameAs(k.Map()));
  ASSERT(rho.Map().SameAs(mu.Map()));
  ASSERT(flux->Map().SameAs(FaceMap(true)));

  int fdirs[6];
  unsigned int cface[6];
  double aux1[6], aux2[6], aux3[6], gflux[6], dummy;

  // Create a view into the cell pressure segment of P.
  Epetra_Vector &Pcell_own = *CreateCellView(pressure);
  Epetra_Vector Pcell(CellMap(true));
  Pcell.Import(Pcell_own, *cell_importer_, Insert);

  // Create a copy of the face pressure segment of P that includes ghosts.
  Epetra_Vector &Pface_own = *CreateFaceView(pressure);
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
    for (int j = 0; j < K.MyLength(); ++j) {
      K[j] = (rho[j] * k[j] / mu[j]);
    }
    ComputeUpwindRelPerm(Pcell, Pface, k_rel, K, rho, K_upwind);
  } else {
    for (int j = 0; j < K.MyLength(); ++j) {
      K[j] = (rho[j] * k[j] * k_rel[j]) / mu[j];
    }
    K_upwind.PutScalar(1.0);  // whole coefficient used in mimetic disc
  }

  // Process-local assembly of the mimetic face fluxes.
  for (unsigned int j = 0; j < Pcell_own.MyLength(); ++j) {
    // Get the list of process-local face indices for this cell.
    mesh_->cell_to_faces(j, cface, cface+6);
    // Gather the local face pressures int AUX1.
    for (int i = 0; i < 6; ++i) aux1[i] = Pface[cface[i]];
    // Compute the local value of the diffusion operator.
    MD_[j].diff_op(K[j], Pcell_own[j], aux1, dummy, aux2);
    // Gravity contribution
    MD_[j].GravityFlux(*gvec_, gflux);

    for (int i = 0; i < 6; ++i) aux2[i] = K[j] * rho[j] * gflux[i] - aux2[i];
    mesh_->cell_to_face_dirs(j, fdirs, fdirs+6);
    // Scatter the local face result into FFACE.
    for (int i = 0; i < 6; ++i) {
      Fface[cface[i]] += fdirs[i] * aux2[i];
      error[cface[i]] += aux2[i]; // sums to the flux discrepancy
      count[cface[i]]++;
    }
    if (j == 0) std::cout << "Fface = " << Fface[cface[0]] << std::endl;
    if (j == 0) std::cout << "error = " << error[cface[0]] << std::endl;
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
    F[j] = K_upwind[j] * F[j] / (rho[j] * count[j]);

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
};

void PermafrostProblem::DeriveLiquidVelocity(const Epetra_Vector &pressure,
        const Epetra_Vector &rho, const Epetra_Vector &k_rel, const Epetra_Vector &k,
        const Epetra_Vector &mu, Teuchos::RCP<Epetra_MultiVector> &velocity) {
  ASSERT(pressure.Map().SameAs(Map())); // cells and faces
  ASSERT(rho.Map().SameAs(CellMap(true))); // ghosted cells
  ASSERT(rho.Map().SameAs(k_rel.Map()));
  ASSERT(rho.Map().SameAs(k.Map()));
  ASSERT(rho.Map().SameAs(mu.Map()));
  ASSERT(rho.Map().SameAs(velocity->Map()));
  ASSERT(velocity.NumVectors() == 3);

  // Cell pressure vectors without and with ghosts.
  Epetra_Vector &Pcell_own = *CreateCellView(pressure);
  Epetra_Vector Pcell(CellMap(true));
  Pcell.Import(Pcell_own, *cell_importer_, Insert);

  // Face pressure vectors without and with ghosts.
  Epetra_Vector &Pface_own = *CreateFaceView(pressure);
  Epetra_Vector Pface(FaceMap(true));
  Pface.Import(Pface_own, *face_importer_, Insert);

  // The pressure gradient coefficient split as K * K_upwind: K gets
  // incorporated into the mimetic discretization and K_upwind is applied
  // afterwards to the computed mimetic flux.  Note that density (rho) is
  // not included here because we want the Darcy velocity and not mass flux.
  // The pressure gradient coefficient split as K * K_upwind: K gets
  // incorporated into the mimetic discretization and K_upwind is applied
  // afterwards to the computed mimetic flux.
  Epetra_Vector K(CellMap(true)), K_upwind(FaceMap(true));
  if (upwind_k_rel_) {
    for (int j = 0; j < K.MyLength(); ++j) {
      K[j] = (k[j] / mu[j]);
    }
    ComputeUpwindRelPerm(Pcell, Pface, k_rel, K, rho, K_upwind);
  } else {
    for (int j = 0; j < K.MyLength(); ++j) {
      K[j] = (k[j] * k_rel[j]) / mu[j];
    }
    K_upwind.PutScalar(1.0);  // whole coefficient used in mimetic disc
  }

  unsigned int cface[6];
  double aux1[6], aux2[6], aux3[6], aux4[3], gflux[6], dummy;

  for (unsigned int j = 0; j < Pcell_own.MyLength(); ++j) {
    mesh_->cell_to_faces(j, cface, cface+6);
    for (int i = 0; i < 6; ++i) aux1[i] = Pface[cface[i]];

    for (int i = 0; i < 6; ++i) aux3[i] = K_upwind[cface[i]];
    MD_[j].diff_op(K[j], aux3, Pcell_own[j], aux1, dummy, aux2);
    MD_[j].GravityFlux(*gvec_, gflux);
    for (int i = 0; i < 6; ++i) aux2[i] = aux3[i] * K[j] * rho[j] * gflux[i] - aux2[i];

    MD_[j].CellFluxVector(aux2, aux4);
    velocity[0][j] = aux4[0];
    velocity[1][j] = aux4[1];
    velocity[2][j] = aux4[2];
  }

  delete &Pcell_own, &Pface_own;
};

void PermafrostProblem::ComputeUpwindRelPerm(const Epetra_Vector& Pcell,
        const Epetra_Vector& Pface, const Epetra_Vector& k_rel_cell, const Epetra_Vector& k,
        const Epetra_Vector& rho, Epetra_Vector& k_rel_face) const {
  ASSERT(Pcell.Map().SameAs(CellMap(true)));
  ASSERT(k_rel_cell.Map().SameAs(CellMap(true)));
  ASSERT(k.Map().SameAs(CellMap(true)));
  ASSERT(rho.Map().SameAs(CellMap(true)));
  ASSERT(Pface.Map().SameAs(FaceMap(true)));
  ASSERT(k_rel_face.Map().SameAs(FaceMap(true)));

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
    for (int i = 0; i < 6; ++i) aux1[i] = Pface[cface[i]];
    // Compute the local value of the diffusion operator.
    MD_[j].diff_op(k[j], Pcell[j], aux1, dummy, aux2); // aux2 is inward flux
    // Gravity contribution; aux2 becomes outward flux
    MD_[j].GravityFlux(*gvec_, gflux);
    for (int i = 0; i < 6; ++i) aux2[i] = rho[j] * K * gflux[i] - aux2[i];
    // Scatter the local face result into FFACE.
    mesh_->cell_to_face_dirs(j, fdirs, fdirs+6);
    for (int i = 0; i < 6; ++i) Fface[cface[i]] += fdirs[i] * aux2[i];
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
    for (int i = 0; i < 6; ++i) {
      if (fdirs[i] > 0) {
        ASSERT(fcell[cface[i]].first == -1);
        fcell[cface[i]].first = j;
      } else {
        ASSERT(fcell[cface[i]].second == -1);
        fcell[cface[i]].second = j;
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
  for (int j = 0; j < nface_own; ++j) {
    if (fcell[j].second == -1) // boundary face
      k_rel_face[j] = k_rel_cell[fcell[j].first];
    else
      k_rel_face[j] = k_rel_cell[((Fface[j] >= 0.0) ? fcell[j].first : fcell[j].second)];
  }

  // Communicate the values computed above to the ghosts.
  double *k_rel_data;
  k_rel_face.ExtractView(&k_rel_data);
  Epetra_Vector k_rel_own(View, FaceMap(false), k_rel_data);
  k_rel_face.Import(k_rel_own, *face_importer_, Insert);

  delete [] fcell;
};

Epetra_Vector* PermafrostProblem::CreateCellView(const Epetra_Vector &X) const {
  // should verify that X.Map() is the same as Map()
  double *data;
  X.ExtractView(&data);
  return new Epetra_Vector(View, CellMap(ghosted), data);
}


Epetra_Vector* PermafrostProblem::CreateFaceView(const Epetra_Vector &X) const {
  // should verify that X.Map() is the same as Map()
  double *data;
  X.ExtractView(&data);
  int ncell = CellMap().NumMyElements();
  return new Epetra_Vector(View, FaceMap(ghosted), data+ncell);
}

void PermafrostProblem::ComputeUDot(const double t, const Epetra_Vector& u,
        Epetra_Vector &udot) {
  ComputeF(u,udot);

  // zero out the face part
  Epetra_Vector *udot_face = CreateFaceView(udot);
  udot_face->PutScalar(0.0);
};

void PermafrostProblem::SetInitialPressureProfileCells(double ref_height, double rho,
        Epetra_Vector &cell_heights, Teuchos::RCP<Epetra_Vector> &pressure) {
  ASSERT(pressure->Map().SameAs(cell_heights.Map()));
  ASSERT(rho.Map().SameAs(cell_heights.Map()));
  for (int j=0; j<pressure->MyLength(); ++j) {
    (*pressure)[j] = p_atm_ + rho*(*gvec_)[2]*(cell_heights[j]-ref_height);
  }
};

void PermafrostProblem::SetInitialPressureProfileFaces(double ref_height,
        Teuchos::RCP<Epetra_Vector> &pressure) {
  ASSERT(pressure->Map().SameAs(FaceMap()));
  for (int j=0; j<pressure->MyLength(); ++j) {
    std::vector<double> coords;
    coords.resize(12);
    mesh_->face_to_coordinates(j, coords.begin(), coords.end());

    // average the x coordinates
    double zavg = 0.0;
    for (int k=2; k<12; k+=3) {
      zavg += coords[k];
    }
    zavg /= 4.0;

    (*pressure)[j] = p_atm_ + rho_*(*gvec_)[2]*(zavg-ref_height);
  }
};
} // close namespace Amanzi
