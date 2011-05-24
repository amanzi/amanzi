#include <vector>

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"

#include "Mesh.hh"
#include "dbc.hh"

#include "Transport_PK.hpp"
#include "Reconstruction.hpp"

using namespace Teuchos;

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* We set up only default values and call Init() routine to complete
* each variable initialization
****************************************************************** */
Transport_PK::Transport_PK(ParameterList &parameter_list_MPC,
			   RCP<Transport_State> TS_MPC)
{ 
  parameter_list = parameter_list_MPC;
  number_components = TS_MPC->get_total_component_concentration()->NumVectors();

  TS = rcp(new Transport_State(*TS_MPC) );

  dT = dT_debug = 0.0;
  status = TRANSPORT_NULL;
  verbosity_level = 0;
  internal_tests = 0;
  tests_tolerance = TRANSPORT_CONCENTRATION_OVERSHOOT;

  MyPID = 0;

  Init();
}


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
int Transport_PK::Init()
{
  TS_nextBIG = rcp(new Transport_State(*TS, CopyMemory) );  
  TS_nextMPC = rcp(new Transport_State(*TS_nextBIG, ViewMemory) );

  RCP<AmanziMesh::Mesh>  mesh = TS->get_mesh_maps();

  const Epetra_Map& cmap = mesh->cell_map(true);
  const Epetra_Map& fmap = mesh->face_map(true);

  cmin = cmap.MinLID();
  cmax = cmap.MaxLID();

  number_owned_cells = mesh->count_entities( AmanziMesh::CELL, AmanziMesh::OWNED );
  cmax_owned = cmin + number_owned_cells - 1;

  fmin = fmap.MinLID();
  fmax = fmap.MaxLID(); 

  number_owned_faces = mesh->count_entities( AmanziMesh::FACE, AmanziMesh::OWNED );
  fmax_owned = fmin + number_owned_faces - 1;

  number_wghost_cells = cmax + 1;  // assume that enumartion starts with 0 
  number_wghost_faces = fmax + 1;

#ifdef HAVE_MPI
  const  Epetra_Comm & comm = cmap.Comm(); 
  MyPID = comm.MyPID();

  const Epetra_Map&  source_cmap = mesh->cell_map(false);
  const Epetra_Map&  target_cmap = mesh->cell_map(true);

  cell_importer = rcp(new Epetra_Import(target_cmap, source_cmap));

  const Epetra_Map&  source_fmap = mesh->face_map(false);
  const Epetra_Map&  target_fmap = mesh->face_map(true);

  face_importer = rcp(new Epetra_Import(target_fmap, source_fmap));
#endif
 
  process_parameter_list();

  upwind_cell_ = rcp(new Epetra_IntVector(fmap));  // The maps include both owned and ghosts
  downwind_cell_ = rcp(new Epetra_IntVector(fmap));
  component_ = rcp(new Epetra_Vector(cmap));

  return 0;
}



/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Transport_PK::process_parameter_list()
{
  RCP<AmanziMesh::Mesh> mesh = TS->get_mesh_maps();

  // global transport parameters
  cfl = parameter_list.get<double>("CFL", 1.0);

  verbosity_level = parameter_list.get<int>("verbosity level", 0);
  internal_tests = parameter_list.get<string>("enable internal tests", "no") == "yes";
  tests_tolerance = parameter_list.get<double>("internal tests tolerance", TRANSPORT_CONCENTRATION_OVERSHOOT);
  dT_debug = parameter_list.get<double>("maximal time step", TRANSPORT_LARGE_TIME_STEP);
 
  // read number of boundary consitions
  ParameterList  BC_list;
 
  BC_list = parameter_list.get<ParameterList>("Transport BCs");
  int nBCs = BC_list.get<int>("number of BCs");

  // create list of boundary data
  bcs.resize(nBCs);

  int i, k;
  for (i=0; i<nBCs; i++) bcs[i] = Transport_BCs(0, number_components);
  
  for (i=0; i<nBCs; i++) {
    char bc_char_name[10];
    
    sprintf(bc_char_name, "BC %d", i);
    string bc_name(bc_char_name);

    if (!BC_list.isSublist(bc_name)) {
      cout << "MyPID = " << MyPID << endl;
      cout << "Boundary condition with name " << bc_char_name << " does not exist" << endl;
      ASSERT(0);
    }

    ParameterList bc_ss = BC_list.sublist(bc_name);
 
    double value;
    int ssid = bc_ss.get<int>("Side set ID");
    string type = bc_ss.get<string>("Type");

    // check all existing components: right now we check by id 
    // but it is possible to check by name in the future 
    for (k=0; k<number_components; k++) {
      char tcc_char_name[20];

      sprintf( tcc_char_name, "Component %d", k );
      string tcc_name( tcc_char_name );

      if (bc_ss.isParameter(tcc_name)) { 
        value = bc_ss.get<double>(tcc_name); 
      }
      else { 
        value = 0.0; 
      }

      bcs[i].values[k] = value;
    }

    if (type == "Constant") bcs[i].type = TRANSPORT_BC_CONSTANT_INFLUX;

    bcs[i].side_set_id = ssid;
     if ( !mesh->valid_set_id( ssid, AmanziMesh::FACE ) ) {
      cout << "MyPID = " << MyPID << endl;
      cout << "Invalid set of mesh faces with ID " << ssid << endl;
      ASSERT(0);
    }

    // populate list of n boundary faces: it could be empty
    int n = mesh->get_set_size(ssid, AmanziMesh::FACE, AmanziMesh::OWNED);
     n = mesh->get_set_size(ssid, AmanziMesh::FACE, AmanziMesh::OWNED);
    bcs[i].faces.resize(n);

     mesh->get_set(ssid, AmanziMesh::FACE, AmanziMesh::OWNED, bcs[i].faces.begin(), bcs[i].faces.end());

    // allocate memory for influx and outflux vectors
    bcs[i].influx.resize(number_components);
    bcs[i].outflux.resize(number_components);
  }
}




/* ************************************************************* */
/* Printing information about Transport status                   */
/* ************************************************************* */
void Transport_PK::print_statistics() const
{
  if ( !MyPID && verbosity_level > 0 ) {
     cout << "Transport PK: CFL = " << cfl << endl;
     cout << "              Total number of components = " << number_components << endl;
     cout << "              Verbosity level = " << verbosity_level << endl;
     cout << "              Enable internal tests = " << (internal_tests ? "yes" : "no")  << endl;
  }
}




/* *******************************************************************
 * Estimation of the time step based on T.Barth (Lecture Notes   
 * presented at VKI Lecture Series 1994-05, Theorem 4.2.2.       
 * Routine must be called every time we update a flow field      
 ****************************************************************** */
double Transport_PK::calculate_transport_dT()
{
  // flow could not be available at initialization, copy it again 
  if (status == TRANSPORT_NULL) {
    TS->copymemory_multivector( TS->ref_total_component_concentration(), TS_nextBIG->ref_total_component_concentration() );
    TS->copymemory_vector( TS->ref_darcy_flux(), TS_nextBIG->ref_darcy_flux() );

    check_divergence_property();
    identify_upwind_cells();

    status = TRANSPORT_FLOW_AVAILABLE;
  }

  // loop over faces and accumulate upwinding fluxes
  int  i, f, c, c1;

  RCP<AmanziMesh::Mesh> mesh = TS->get_mesh_maps();
  const Epetra_Map& fmap = mesh->face_map( true );
  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();

  std::vector<double> total_influx(number_wghost_cells, 0.0);

  for (f=fmin; f<=fmax; f++) {
    c = (*downwind_cell_)[f];
    if (c >= 0) total_influx[c] += fabs(darcy_flux[f]); 
  }

  // loop over cells and calculate minimal dT
  double influx, dT_cell; 

  const Epetra_Vector& ws  = TS->ref_water_saturation();
  const Epetra_Vector& phi = TS->ref_porosity();

  dT = dT_cell = TRANSPORT_LARGE_TIME_STEP;
  for (c=cmin; c<=cmax_owned; c++) {
    influx = total_influx[c];
    if( influx ) dT_cell = mesh->cell_volume(c) * phi[c] * ws[c] / influx;

    dT = std::min(dT, dT_cell);
  }

#ifdef HAVE_MPI
  double dT_global;
  const  Epetra_Comm & comm = ws.Comm(); 
 
  comm.MinAll(&dT, &dT_global, 1);
  dT = dT_global;
#endif

  // incorporate developers and CFL constraints
  dT = std::min(dT, dT_debug);

  dT *= cfl;
  return dT;
}


/* ******************************************************************* 
 * MPC will call this function to advance the transport state    
 ****************************************************************** */
void Transport_PK::advance(double dT_MPC)
{
  if (TRANSPORT_AMANZI_VERSION == 1) {  // temporary solution (lipnikov@lanl.gov)
    advance_donor_upwind(dT_MPC);
  }
  else if (TRANSPORT_AMANZI_VERSION == 2) {
    advance_second_order_upwind(dT_MPC);
  }
}


/* ******************************************************************* 
 * We have to advance each component independently due to different
 * reconstructions.
 ****************************************************************** */
void Transport_PK::advance_second_order_upwind(double dT_MPC)
{
  status = TRANSPORT_STATE_BEGIN;
  dT = dT_MPC;  // overwrite the transport step

  int i, f, c, c1, c2;
  RCP<AmanziMesh::Mesh> mesh = TS->get_mesh_maps();

  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();
  const Epetra_Vector& ws  = TS_nextBIG->ref_water_saturation();
  const Epetra_Vector& phi = TS_nextBIG->ref_porosity();

  // populating next state of concentrations
  RCP<Epetra_MultiVector> tcc = TS->get_total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_nextBIG->get_total_component_concentration();
  TS_nextBIG->copymemory_multivector(*tcc, *tcc_next);

  // prepare conservative state in master and slave cells 
  double u, vol_phi_ws, tcc_flux;
  int num_components = tcc->NumVectors();
 
  Reconstruction lifting(mesh, component_);

  for (i=0; i<num_components; i++) {
    // calculate conservative quantatity
    for (c=cmin; c<=cmax_owned; c++) {
      vol_phi_ws = mesh->cell_volume(c) * phi[c] * ws[c]; 
      (*tcc_next)[i][c] = (*component_)[c] = (*tcc)[i][c] * vol_phi_ws;
    }

    lifting.reset_field(mesh, component_);
    lifting.calculateCellGradient();

    for (f=fmin; f<=fmax; f++) {  // loop over master and slave faces
      c1 = (*upwind_cell_)[f]; 
      c2 = (*downwind_cell_)[f]; 

      u = fabs(darcy_flux[f]);

      if (c1 >=0 && c1 <= cmax_owned && c2 >= 0 && c2 <= cmax_owned) {
        tcc_flux = dT * u * (*tcc)[i][c1];
        (*component_)[c1] -= tcc_flux;
        (*component_)[c2] += tcc_flux;
      } 
      else if (c1 >=0 && c1 <= cmax_owned && (c2 > cmax_owned || c2 < 0)) {
        tcc_flux = dT * u * (*tcc)[i][c1];
        (*component_)[c1] -= tcc_flux;
      } 
      else if (c1 > cmax_owned && c2 >= 0 && c2 <= cmax_owned) {
        tcc_flux = dT * u * (*component_)[c1];
        (*component_)[c2] += tcc_flux;
      }
    } 

    // loop over exterior boundary sets
    int k, n;
    for (n=0; n<bcs.size(); n++) {
      for (k=0; k<bcs[n].faces.size(); k++) {
        f = bcs[n].faces[k];
        c2 = (*downwind_cell_)[f]; 

        if (c2 >= 0) {
          u = fabs(darcy_flux[f]);

          if (bcs[n].type == TRANSPORT_BC_CONSTANT_INFLUX) {
            tcc_flux = dT * u * bcs[n].values[i];
            (*component_)[c2] += tcc_flux;
            bcs[n].influx[i] += tcc_flux;
          }
        } 
      }
    }

    // recover concentration from new conservative state
    for (c=cmin; c<=cmax_owned; c++) {
      vol_phi_ws = mesh->cell_volume(c) * phi[c] * ws[c]; 
      (*tcc_next)[i][c] = (*component_)[c] / vol_phi_ws;
    }
  }

  if (internal_tests) {
    RCP<Epetra_MultiVector> tcc_nextMPC = TS_nextMPC->get_total_component_concentration();
    check_GEDproperty(*tcc_nextMPC);
  }

  status = TRANSPORT_STATE_COMPLETE;
}



/* ******************************************************************* 
 * A simple first-order transport method 
 ****************************************************************** */
void Transport_PK::advance_donor_upwind(double dT_MPC)
{
  status = TRANSPORT_STATE_BEGIN;
  dT = dT_MPC;  // overwrite the transport step

  int i, f, c, c1, c2;
  RCP<AmanziMesh::Mesh> mesh = TS->get_mesh_maps();

  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();
  const Epetra_Vector& ws  = TS_nextBIG->ref_water_saturation();
  const Epetra_Vector& phi = TS_nextBIG->ref_porosity();

  // populating next state of concentrations
  RCP<Epetra_MultiVector> tcc = TS->get_total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_nextBIG->get_total_component_concentration();
  TS_nextBIG->copymemory_multivector(*tcc, *tcc_next);

  // prepare conservative state in master and slave cells 
  double vol_phi_ws, tcc_flux;
  int num_components = tcc->NumVectors();

  for (c=cmin; c<=cmax_owned; c++) {
    vol_phi_ws = mesh->cell_volume(c) * phi[c] * ws[c]; 

    for (i=0; i<num_components; i++) 
      (*tcc_next)[i][c] = (*tcc)[i][c] * vol_phi_ws;
  }

  // advance all components at once
  double u;
  for (f=fmin; f<=fmax; f++) {  // loop over master and slave faces
    c1 = (*upwind_cell_)[f]; 
    c2 = (*downwind_cell_)[f]; 

    u = fabs(darcy_flux[f]);

    if (c1 >=0 && c1 <= cmax_owned && c2 >= 0 && c2 <= cmax_owned) {
      for( i=0; i<num_components; i++ ) {
        tcc_flux = dT * u * (*tcc)[i][c1];
        (*tcc_next)[i][c1] -= tcc_flux;
        (*tcc_next)[i][c2] += tcc_flux;
      }
    } 
    else if (c1 >=0 && c1 <= cmax_owned && (c2 > cmax_owned || c2 < 0)) {
      for (i=0; i<num_components; i++) {
        tcc_flux = dT * u * (*tcc)[i][c1];
        (*tcc_next)[i][c1] -= tcc_flux;
      }
    } 
    else if (c1 > cmax_owned && c2 >= 0 && c2 <= cmax_owned) {
      for (i=0; i<num_components; i++) {
        tcc_flux = dT * u * (*tcc_next)[i][c1];
        (*tcc_next)[i][c2] += tcc_flux;
      }
    } 
  }

  // loop over exterior boundary sets
  int k, n;
  for (n=0; n<bcs.size(); n++) {
    for (k=0; k<bcs[n].faces.size(); k++) {
      f = bcs[n].faces[k];
      c2 = (*downwind_cell_)[f]; 

      if (c2 >= 0) {
        u = fabs(darcy_flux[f]);

        if (bcs[n].type == TRANSPORT_BC_CONSTANT_INFLUX) {
          for (i=0; i<num_components; i++) {
            tcc_flux = dT * u * bcs[n].values[i];
            (*tcc_next)[i][c2] += tcc_flux;
            bcs[n].influx[i] += tcc_flux;
          }
        } 
      }
    }
  }

  // recover concentration from new conservative state
  for (c=cmin; c<=cmax_owned; c++) {
    vol_phi_ws = mesh->cell_volume(c) * phi[c] * ws[c]; 
    for (i=0; i<num_components; i++) (*tcc_next)[i][c] /= vol_phi_ws;
  }

  if (internal_tests) {
    RCP<Epetra_MultiVector> tcc_nextMPC = TS_nextMPC->get_total_component_concentration();
    check_GEDproperty(*tcc_nextMPC);
  }

  status = TRANSPORT_STATE_COMPLETE;
}



/* ************************************************************* */
/*  MPC will call this function to indicate to the transport PK  */
/*  that it can commit the advanced state it has created.        */
/*  This  call indicates that the MPC has accepted the new state */
/* ************************************************************* */
void Transport_PK::commit_state( RCP<Transport_State> TS )
{
  /* nothing is done her since a pointer to the state is kept */ 
};




/* ************************************************************* */
/* Identify flux direction based on direction of the face normal */
/* andsign of the  Darcy velocity.                               */
/* ************************************************************* */
void Transport_PK::identify_upwind_cells()
{
  RCP<AmanziMesh::Mesh> mesh = TS->get_mesh_maps();
 
  for (int f=fmin; f<=fmax; f++) {
    (*upwind_cell_)[f] = -1;  // negative value is indicator of a boundary
    (*downwind_cell_)[f] = -1;
  }

  AmanziMesh::Entity_ID_List faces; 
  std::vector<int> fdirs;
  Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();

  for (int c=cmin; c<=cmax; c++) {
    mesh->cell_get_faces(c, &faces);
    mesh->cell_get_face_dirs(c, &fdirs);

    for (int i=0; i<6; i++) {
      int f = faces[i];
      if (darcy_flux[f] * fdirs[i] >= 0) { 
        (*upwind_cell_)[f] = c; 
      } else { 
        (*downwind_cell_)[f] = c; 
      }
    }
  }
}



/* ************************************************************* */
/* verifies that the velocity field is divergence free           */
/* Demo I:  divergence must be zero or almost zero               */
/* Demo II: divergence must be non-negative if dS/dP >= 0        */
/* ************************************************************* */
void Transport_PK::check_divergence_property()
{
  int i, c, f;
  double div, u, umax, L8_error;

  RCP<AmanziMesh::Mesh> mesh = TS->get_mesh_maps();
  Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();

  L8_error = 0;

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  for (int c=cmin; c<=cmax_owned; c++) {
    mesh->cell_get_faces(c, &faces);
    mesh->cell_get_face_dirs(c, &fdirs);

    div = umax = 0;
    for (i=0; i<6; i++) {
      f = faces[i];
      u = darcy_flux[f];
      div += u * fdirs[i];
      umax = std::max(umax, fabs(u) / pow(mesh->face_area(f), 0.5));
    }
    div /= mesh->cell_volume(c);

    if (umax) L8_error = std::max(L8_error, fabs(div) / umax);

    /* verify that divergence complies with the flow model  */
    int flag = 0;
    if (TRANSPORT_AMANZI_VERSION == 1 && fabs(div) > tests_tolerance * umax) {
      cout << "TRANSPORT: The flow violates conservation property." << endl;
      cout << "    Modify either flow convergence criteria or transport tolerance." << endl;
      flag = 1;
    }
    if (TRANSPORT_AMANZI_VERSION == 2 && div < -tests_tolerance * umax) {
      cout << "TRANSPORT: The flow has large artificial sinks."<< endl;
      flag = 1;
    }

    if (flag && verbosity_level > 1) {
      cout << "    MyPID = " << MyPID << endl;
      cout << "    cell  = " << c << endl;
      cout << "    divergence = " << div << endl;
      cout << "    maximal velocity = " << umax << endl;
      ASSERT(0);
    }
  }

  if ( verbosity_level > 3 ) {
#ifdef HAVE_MPI
     double L8_global;
     const  Epetra_Comm & comm = darcy_flux.Comm(); 
 
     comm.MinAll( &L8_error, &L8_global, 1 );
     L8_error = L8_global;
#endif
     if ( !MyPID ) cout << "Transport_PK: maximal (divergence / flux) = " << L8_error << endl;
  }
}


/* *******************************************************************
 * Check that global extrema diminished                          
 ****************************************************************** */
void Transport_PK::check_GEDproperty(Epetra_MultiVector& tracer) const
{ 
  int i, num_components = tracer.NumVectors();
  double tol; 

  double tr_min[num_components];
  double tr_max[num_components];

  tracer.MinValue(tr_min);
  tracer.MaxValue(tr_max);

  if (TRANSPORT_AMANZI_VERSION == 1) {
     for (i=0; i<num_components; i++) {
        if (tr_min[i] < 0) {
           cout << "Transport_PK: concentration violated GED property" << endl; 
           cout << "    Make an Amanzi ticket or turn off internal transport tests" << endl;
           cout << "    MyPID = " << MyPID << endl;
           cout << "    component = " << i << endl;
           cout << "    min/max values = " << tr_min[i] << " " << tr_max[i] << endl;

           ASSERT(0); 
        }
     }
  }
}

}  // namespace AmanziTransport
}  // namespace Amanzi
