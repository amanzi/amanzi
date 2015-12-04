/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt  
*/

#include <map>

// Amanzi
#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "PlaneRegion.hh"
#include "PolygonRegion.hh"
#include "ReconstructionCell.hh"
#include "Units.hh"

// MPC
#include "Unstructured_observations.hh"

namespace Amanzi {

/* ******************************************************************
* Constructor.
****************************************************************** */
Unstructured_observations::Unstructured_observations(
    Teuchos::RCP<Teuchos::ParameterList> obs_list,
    Teuchos::RCP<Teuchos::ParameterList> units_list,
    Amanzi::ObservationData& observation_data,
    Epetra_MpiComm* comm)
    : observation_data_(observation_data),
      obs_list_(obs_list)
{
  rank_ = comm->MyPID();

  // initialize units
  units_.Init(*units_list);

  Teuchos::ParameterList tmp_list;
  tmp_list.set<std::string>("Verbosity Level", "high");
  vo_ = new VerboseObject("Observations", tmp_list);

  // loop over the sublists and create an observation for each
  for (Teuchos::ParameterList::ConstIterator i = obs_list_->begin(); i != obs_list_->end(); i++) {

    if (obs_list_->isSublist(obs_list_->name(i))) {
      Teuchos::ParameterList observable_plist = obs_list_->sublist(obs_list_->name(i));

      std::vector<double> times;
      std::vector<std::vector<double> > time_sps;

      // get the observation times
      if (observable_plist.isSublist("time start period stop")) {
        Teuchos::ParameterList& tsps_list = observable_plist.sublist("time start period stop");
        for (Teuchos::ParameterList::ConstIterator it = tsps_list.begin(); it != tsps_list.end(); ++it) {
          std::string name = it->first;
          if (tsps_list.isSublist(name)) {
            Teuchos::ParameterList& itlist = tsps_list.sublist(name);
            if (itlist.isParameter("start period stop")) {
              Teuchos::Array<double> sps = itlist.get<Teuchos::Array<double> >("start period stop");
              time_sps.push_back(sps.toVector());
            }
          }
        }
      }
      if (observable_plist.isParameter("times")) {
        Teuchos::Array<double> vtimes = observable_plist.get<Teuchos::Array<double> >("times");
        times = vtimes.toVector();
      }

      std::vector<int> cycles;
      std::vector<std::vector<int> > cycle_sps;

      // get the observation cycles
      if (observable_plist.isSublist("cycle start period stop")) {
        Teuchos::ParameterList& csps_list = observable_plist.sublist("cycle start period stop");
        for (Teuchos::ParameterList::ConstIterator it = csps_list.begin(); it != csps_list.end(); ++it) {
          std::string name = it->first;
          if (csps_list.isSublist(name)) {
            Teuchos::ParameterList& itlist = csps_list.sublist(name);
            if (itlist.isParameter("start period stop")) {
              Teuchos::Array<int> csps = itlist.get<Teuchos::Array<int> >("start period stop");
              cycle_sps.push_back(csps.toVector());
            }
          }
        }
      }

      if (observable_plist.isParameter("cycles")) {
        Teuchos::Array<int> vcycles = observable_plist.get<Teuchos::Array<int> >("cycles");
        cycles = vcycles.toVector();
      }

      // loop over all variables listed and create an observable for each
      std::string var = observable_plist.get<std::string>("variable");
      observations.insert(std::pair<std::string, Observable>(
          obs_list_->name(i), 
          Observable(var, observable_plist.get<std::string>("region"),
                     observable_plist.get<std::string>("functional"),
                     observable_plist, comm)));
    }
  }
}


/* ******************************************************************
* Process data to extract observations.
****************************************************************** */
int Unstructured_observations::MakeObservations(State& S)
{
  Errors::Message msg;
  int num_obs(0);

  // loop over all observables
  for (std::map<std::string, Observable>::iterator i = observations.begin(); i != observations.end(); i++) {
    if ((i->second).DumpRequested(S.time()) || (i->second).DumpRequested(S.cycle())) {
      num_obs++; 
      
      // for now we can only observe Integrals and Values
      if ((i->second).functional != "Observation Data: Integral"  &&
          (i->second).functional != "Observation Data: Point" )  {
        msg << "Unstructured_observations: can only handle Functional == Observation Data:"
            << " Integral, or Functional == Observation Data: Point";
        Exceptions::amanzi_throw(msg);
      }
      
      std::string label = i->first;
      
      // we need to make an observation for each variable in the observable
      std::string var = (i->second).variable;
      
      // data structure to store the observation
      Amanzi::ObservationData::DataTriple data_triplet;
      
      // build the name of the observation
      std::stringstream ss;
      ss << label << ", " << var;
      
      std::vector<Amanzi::ObservationData::DataTriple>& od = observation_data_[label]; 

      // check if observation is planar
      bool obs_planar(false);

      AmanziGeometry::GeometricModelPtr gm_ptr = S.GetMesh()->geometric_model();
      AmanziGeometry::RegionPtr reg_ptr = gm_ptr->FindRegion((i->second).region);
      AmanziGeometry::Point reg_normal;
      if (reg_ptr->type() == AmanziGeometry::POLYGON) {
        AmanziGeometry::PolygonRegion *poly_reg = dynamic_cast<AmanziGeometry::PolygonRegion*>(reg_ptr);
        reg_normal = poly_reg->normal();
        obs_planar = true;
      } else if (reg_ptr->type() == AmanziGeometry::PLANE) {
        AmanziGeometry::PlaneRegion *plane_reg = dynamic_cast<AmanziGeometry::PlaneRegion*>(reg_ptr);
        reg_normal = plane_reg->normal();
        obs_planar = true;
      }

      // check if observation of solute was requested
      bool obs_solute_liquid(false), obs_solute_gas(false), obs_aqueous(true);     
      int tcc_index(-1);
      for (tcc_index = 0; tcc_index != comp_names_.size(); ++tcc_index) {
        int pos = var.find(comp_names_[tcc_index]);
        if (pos == 0) { 
          (tcc_index < num_liquid_) ? obs_solute_liquid = true : obs_solute_gas = true;
          obs_aqueous = false;
          break;
        }
      }
      bool obs_solute = obs_solute_liquid || obs_solute_gas;

      // check if observation is on faces of cells. 
      bool obs_boundary(false);
      unsigned int mesh_block_size(0);
      AmanziMesh::Entity_ID_List entity_ids;
      std::string solute_var;
      if (obs_solute) solute_var = comp_names_[tcc_index] + " volumetric flow rate";
      if (var == "aqueous mass flow rate" || 
          var == "aqueous volumetric flow rate" ||
          var == solute_var) {  // flux needs faces
        mesh_block_size = S.GetMesh()->get_set_size((i->second).region,
                                                    Amanzi::AmanziMesh::FACE,
                                                    Amanzi::AmanziMesh::OWNED);
        entity_ids.resize(mesh_block_size);
        S.GetMesh()->get_set_entities((i->second).region, 
                                      Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED,
                                      &entity_ids);
        obs_boundary = true;
        for (int i = 0; i != mesh_block_size; ++i) {
          int f = entity_ids[i];
          Amanzi::AmanziMesh::Entity_ID_List cells;
          S.GetMesh()->face_get_cells(f, Amanzi::AmanziMesh::USED, &cells);
          if (cells.size() == 2) {
            obs_boundary = false;
            break;
          }
        }
      } else { // all others need cells
        mesh_block_size = S.GetMesh()->get_set_size((i->second).region,
                                                    Amanzi::AmanziMesh::CELL,
                                                    Amanzi::AmanziMesh::OWNED);    
        entity_ids.resize(mesh_block_size);
        S.GetMesh()->get_set_entities((i->second).region,
                                      Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                                      &entity_ids);
      }
      
      // find global meshblocksize
      int dummy = mesh_block_size; 
      int global_mesh_block_size(0);
      S.GetMesh()->get_comm()->SumAll(&dummy, &global_mesh_block_size, 1);
      
      if (global_mesh_block_size == 0) {  // bail if this region is empty
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "Cannot make observation on region \"" << (i->second).region 
                   << "\", that is empty or has incorrect face/cell type, skipping it." << std::endl;
        continue;
      }

      double value(0.0), volume(0.0);

      // the user is asking for an observation on tcc
      if (obs_solute) { 
        if (!S.HasField("total_component_concentration")) {  // bail out if this field is not yet created
          Teuchos::OSTab tab = vo_->getOSTab();
          *vo_->os() << "Field \"total_component_concentration\" does not exist, skipping it." << std::endl;
          continue;
        }

        const Epetra_MultiVector& ws = *S.GetFieldData("saturation_liquid")->ViewComponent("cell");
        const Epetra_MultiVector& tcc = *S.GetFieldData("total_component_concentration")->ViewComponent("cell");
        const Epetra_MultiVector& porosity = *S.GetFieldData("porosity")->ViewComponent("cell");    

        if (var == comp_names_[tcc_index] + " aqueous concentration") { 
          for (int i = 0; i < mesh_block_size; i++) {
            int c = entity_ids[i];
            double factor = porosity[0][c] * ws[0][c] * S.GetMesh()->cell_volume(c);
            factor *= units_.concentration_factor();

            value += tcc[tcc_index][c] * factor;
            volume += factor;
          }

        } else if (var == comp_names_[tcc_index] + " gaseous concentration") { 
          for (int i = 0; i < mesh_block_size; i++) {
            int c = entity_ids[i];
            double factor = porosity[0][c] * (1.0 - ws[0][c]) * S.GetMesh()->cell_volume(c);
            factor *= units_.concentration_factor();

            value += tcc[tcc_index][c] * factor;
            volume += factor;
          }

        } else if (var == comp_names_[tcc_index] + " volumetric flow rate") {
          const Epetra_MultiVector& darcy_flux = *S.GetFieldData("darcy_flux")->ViewComponent("face");
          Amanzi::AmanziMesh::Entity_ID_List cells;

          if (obs_boundary) { // observation is on a boundary set
            for (int i = 0; i != mesh_block_size; ++i) {
              int f = entity_ids[i];
              S.GetMesh()->face_get_cells(f, Amanzi::AmanziMesh::USED, &cells);

              int sign, c = cells[0];
              const AmanziGeometry::Point& face_normal = S.GetMesh()->face_normal(f, false, c, &sign);
              double area = S.GetMesh()->face_area(f);
              double factor = units_.concentration_factor();

              value += std::max(0.0, sign * darcy_flux[0][f]) * tcc[tcc_index][c] * factor;
              volume += area * factor;
            }

          } else if (obs_planar) {  // observation is on an interior planar set
            for (int i = 0; i != mesh_block_size; ++i) {
              int f = entity_ids[i];
              S.GetMesh()->face_get_cells(f, Amanzi::AmanziMesh::USED, &cells);

              int csign, c = cells[0];
              const AmanziGeometry::Point& face_normal = S.GetMesh()->face_normal(f, false, c, &csign);
              if (darcy_flux[0][f] * csign < 0) c = cells[1];

              double area = S.GetMesh()->face_area(f);
              double sign = (reg_normal * face_normal) * csign / area;
              double factor = units_.concentration_factor();
    
              value += sign * darcy_flux[0][f] * tcc[tcc_index][c] * factor;
              volume += area * factor;
            }

          } else {
            msg << "Observations of \"SOLUTE volumetric flow rate\""
                << " is only possible for Polygon, Plane and Boundary side sets";
            Exceptions::amanzi_throw(msg);
          }
        } else {
          msg << "Cannot make an observation for solute variable \"" << var << "\"";
          Exceptions::amanzi_throw(msg);
        }
      }

      // aqueous observations
      if (obs_aqueous) {
        double rho = *S.GetScalarData("fluid_density");
        const Epetra_MultiVector& porosity = *S.GetFieldData("porosity")->ViewComponent("cell");    
        const Epetra_MultiVector& ws = *S.GetFieldData("saturation_liquid")->ViewComponent("cell");
        const Epetra_MultiVector& pressure = *S.GetFieldData("pressure")->ViewComponent("cell");
  
        if (var == "volumetric water content") {
          for (int i = 0; i < mesh_block_size; i++) {
            int c = entity_ids[i];
            double vol = S.GetMesh()->cell_volume(c);
            volume += vol;
            value += porosity[0][c] * ws[0][c] * vol;
          }
        } else if (var == "gravimetric water content") {
          if (!S.HasField("particle_density")) {
            msg << "Observation \""  << var << "\" requires field \"particle_density\".\n";
            Exceptions::amanzi_throw(msg);
          }
          const Epetra_MultiVector& pd = *S.GetFieldData("particle_density")->ViewComponent("cell");    
  
          for (int i = 0; i < mesh_block_size; i++) {
            int c = entity_ids[i];
            double vol = S.GetMesh()->cell_volume(c);
            volume += vol;
            value += porosity[0][c] * ws[0][c] * rho / (pd[0][c] * (1.0 - porosity[0][c])) * vol;
          }    
        } else if (var == "aqueous pressure") {
          for (int i = 0; i < mesh_block_size; i++) {
            int c = entity_ids[i];
            double vol = S.GetMesh()->cell_volume(c);
            volume += vol;
            value += pressure[0][c] * vol;
          }
        } else if (var == "water table") {
          value = CalculateWaterTable_(S, entity_ids);
          volume = 1.0;
        } else if (var == "aqueous saturation") {
          for (int i = 0; i < mesh_block_size; i++) {
            int c = entity_ids[i];
            double vol = S.GetMesh()->cell_volume(c);
            volume += vol;
            value += ws[0][c] * vol;
          }    
        } else if (var == "hydraulic head") {
          const Epetra_MultiVector& hydraulic_head = *S.GetFieldData("hydraulic_head")->ViewComponent("cell");
  
          for (int i = 0; i < mesh_block_size; ++i) {
            int c = entity_ids[i];
            double vol = S.GetMesh()->cell_volume(c);
            volume += vol;
            value += hydraulic_head[0][c] * vol;
          }
        } else if (var == "drawdown") {
          const Epetra_MultiVector& hydraulic_head = *S.GetFieldData("hydraulic_head")->ViewComponent("cell");
  
          for (int i = 0; i < mesh_block_size; ++i) {
            int c = entity_ids[i];
            double vol = S.GetMesh()->cell_volume(c);
            volume += vol;
            value += hydraulic_head[0][c] * vol;
          }

          // zero drawdown at time = t0 wil be written directly to the file.
          if (od.size() > 0) { 
            value = od.begin()->value - value;
          }
        } else if (var == "aqueous mass flow rate" || 
                   var == "aqueous volumetric flow rate") {
          double density(1.0);
          if (var == "aqueous mass flow rate") density = rho;
          const Epetra_MultiVector& darcy_flux = *S.GetFieldData("darcy_flux")->ViewComponent("face");
  
          if (obs_boundary) { // observation is on a boundary set
            Amanzi::AmanziMesh::Entity_ID_List cells;

            for (int i = 0; i != mesh_block_size; ++i) {
              int f = entity_ids[i];
              S.GetMesh()->face_get_cells(f, Amanzi::AmanziMesh::USED, &cells);

              int sign, c = cells[0];
              const AmanziGeometry::Point& face_normal = S.GetMesh()->face_normal(f, false, c, &sign);
              double area = S.GetMesh()->face_area(f);

              value += sign * darcy_flux[0][f] * density;
              volume += area;
            }
          } else if (obs_planar) {  // observation is on an interior planar set
            for (int i = 0; i != mesh_block_size; ++i) {
              int f = entity_ids[i];
              const AmanziGeometry::Point& face_normal = S.GetMesh()->face_normal(f);
              double area = S.GetMesh()->face_area(f);
              double sign = reg_normal * face_normal / area;
    
              value += sign * darcy_flux[0][f] * density;
              volume += area;
            }
          } else {
            msg << "Observations of \"aqueous mass flow rate\" and \"aqueous volumetric flow rate\""
                << " are only possible for Polygon, Plane and Boundary side sets";
            Exceptions::amanzi_throw(msg);
          }

        } else if (var == "pH") {
          const Epetra_MultiVector& pH = *S.GetFieldData("pH")->ViewComponent("cell");
  
          for (int i = 0; i < mesh_block_size; ++i) {
            int c = entity_ids[i];
            double vol = S.GetMesh()->cell_volume(c);
            volume += vol;
            value += pH[0][c] * vol;
          }
        } else {
          msg << "Cannot make an observation for aqueous variable \"" << var << "\"";
          Exceptions::amanzi_throw(msg);
        }
      }
      
      // syncronize the result across processors
      double result;
      S.GetMesh()->get_comm()->SumAll(&value, &result, 1);
      
      double vresult;
      S.GetMesh()->get_comm()->SumAll(&volume, &vresult, 1);
 
      if ((i->second).functional == "Observation Data: Integral") {  
        data_triplet.value = result;  
      } else if ((i->second).functional == "Observation Data: Point") {
        data_triplet.value = result / vresult;
      }
      
      data_triplet.is_valid = true;
      data_triplet.time = S.time();

      bool time_exist = false;
      for (std::vector<Amanzi::ObservationData::DataTriple>::iterator it = od.begin(); it != od.end(); ++it) {
        if (it->time == data_triplet.time) {
          time_exist = true;
          break;
        }
      }
            
      if (!time_exist) od.push_back(data_triplet);
    }
  }

  FlushObservations();
  return num_obs;
}


/* ******************************************************************
* Auxiliary routine: calculate maximum water table in a region.
****************************************************************** */
double Unstructured_observations::CalculateWaterTable_(
    State& S, AmanziMesh::Entity_ID_List& ids)
{
  Teuchos::RCP<const Epetra_MultiVector> pressure = S.GetFieldData("pressure")->ViewComponent("cell", true);
  double patm = *S.GetScalarData("atmospheric_pressure");

  // initilize and apply the reconstruction operator
  Teuchos::ParameterList plist;
  Operators::ReconstructionCell lifting(S.GetMesh());
  std::vector<AmanziGeometry::Point> gradient; 

  lifting.Init(pressure, plist);
  lifting.ComputeGradient(ids, gradient);

  // set up extreme values for water table
  int dim = S.GetMesh()->space_dimension();
  double zmin(1e+99), zmax(-1e+99), pref(-1e+99), value(-1e+99);

  // estimate water table
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int found(0);
  for (int i = 0; i < ids.size(); i++) {
    int c = ids[i];
    const AmanziGeometry::Point& xc = S.GetMesh()->cell_centroid(c);
    double pf, pc = (*pressure)[0][c];
    pref = pc;

    S.GetMesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
    for (int n = 0; n < faces.size(); ++n) {
      const AmanziGeometry::Point& xf = S.GetMesh()->face_centroid(faces[n]);
      zmin = std::min(zmin, xf[dim - 1]);
      zmax = std::max(zmax, xf[dim - 1]);

      double dp = gradient[i] * (xf - xc);
      pf = pc + dp;

      if ((pf - patm) * (pc - patm) <= 0.0) {
        if (fabs(dp) > 1e-8) {
          double a = (patm - pc) / dp;
          value = xc[dim - 1] + (xf[dim - 1] - xc[dim - 1]) * a;
          found = 1;
          break;
        }
      }
    }
  }

  // parallel update
  double tmp_loc[3] = {value, pref, zmax};
  double tmp_glb[3];
  S.GetMesh()->get_comm()->MaxAll(tmp_loc, tmp_glb, 3);
  value = tmp_glb[0];
  pref = tmp_glb[1];
  zmax = tmp_glb[2];

  double zmin_tmp(zmin);
  S.GetMesh()->get_comm()->MinAll(&zmin_tmp, &zmin, 1);

  int found_tmp = found;
  S.GetMesh()->get_comm()->MaxAll(&found_tmp, &found, 1);

  // process fully saturated and dry cases
  if (found == 0) {
    if (pref < patm) value = zmin;
    if (pref > patm) value = zmax;
  }

  return value;
}


/* ******************************************************************
* Save observation based on time or cycle.
****************************************************************** */
bool Unstructured_observations::DumpRequested(const double time) {
  bool result = false;
  for (std::map<std::string, Observable>::iterator i = observations.begin(); i != observations.end(); i++) {
    result = result || (i->second).DumpRequested(time);
  }  
  return result;
}


bool Unstructured_observations::DumpRequested(const int cycle) {
  bool result = false;
  for (std::map<std::string, Observable>::iterator i = observations.begin(); i != observations.end(); i++) {
    result = result || (i->second).DumpRequested(cycle);
  }  
  return result;  
}


bool Unstructured_observations::DumpRequested(const int cycle, const double time) {
  return DumpRequested(time) || DumpRequested(cycle);
}


/* ******************************************************************
* Loop over all observations and register each of them with the time 
* step manager.
****************************************************************** */
void Unstructured_observations::RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm) {
  for (std::map<std::string, Observable>::iterator i = observations.begin(); i != observations.end(); i++) {
    (i->second).RegisterWithTimeStepManager(tsm);
  }  
}


/* ******************************************************************
* Write observatins to a file. Clsoe the file to flush data.
****************************************************************** */
void Unstructured_observations::FlushObservations()
{
  if (obs_list_->isParameter("observation output filename")) {
    std::string obs_file = obs_list_->get<std::string>("observation output filename");
    int precision = obs_list_->get<int>("precision", 16);

    if (rank_ == 0) {
      std::ofstream out;
      out.open(obs_file.c_str(), std::ios::out);
      
      out.precision(precision);
      out.setf(std::ios::scientific);

      out << "Observation Name, Region, Functional, Variable, Time, Value\n";
      out << "===========================================================\n";

      for (Teuchos::ParameterList::ConstIterator i = obs_list_->begin(); i != obs_list_->end(); ++i) {
        std::string label = obs_list_->name(i);
        const Teuchos::ParameterEntry& entry = obs_list_->getEntry(label);
        if (entry.isList()) {
          const Teuchos::ParameterList& ind_obs_list = obs_list_->sublist(label);
          std::vector<Amanzi::ObservationData::DataTriple>& od = observation_data_[label]; 

          for (int j = 0; j < od.size(); j++) {
            if (od[j].is_valid) {
              if (!out.good()) {
                std::cout << "PROBLEM BEFORE" << std::endl;
              }
              std::string var = ind_obs_list.get<std::string>("variable");
              out << label << ", "
                  << ind_obs_list.get<std::string>("region") << ", "
                  << ind_obs_list.get<std::string>("functional") << ", "
                  << var << ", "
                  << od[j].time << ", "
                  << ((var == "drawdown" && !j) ? 0.0 : od[j].value) << '\n';
              if (!out.good()) {
                std::cout << "PROBLEM AFTER" << std::endl;
              }
            }
          }
        }
      }
      out.close();
    }
  }
}

}  // namespace Amanzi
