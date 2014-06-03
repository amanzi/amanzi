#include "Unstructured_observations.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "PolygonRegion.hh"
#include "PlaneRegion.hh"

#include <map>

namespace Amanzi {

Unstructured_observations::Unstructured_observations(Teuchos::ParameterList obs_list,
                                                    Amanzi::ObservationData& observation_data,
                                                    Epetra_MpiComm* comm)
    : observation_data_(observation_data), obs_list_(obs_list)
{
  rank_ = comm->MyPID();

  Teuchos::ParameterList tmp_list;
  tmp_list.set<std::string>("Verbosity Level", "high");
  vo_ = new VerboseObject("Observations", tmp_list);

  // interpret paramerter list
  // loop over the sublists and create an observation for each
  for (Teuchos::ParameterList::ConstIterator i = obs_list_.begin(); i != obs_list_.end(); i++) {

    if (obs_list_.isSublist(obs_list_.name(i))) {
      Teuchos::ParameterList observable_plist = obs_list_.sublist(obs_list_.name(i));

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
      observations.insert(std::pair
                          <std::string, Observable>(obs_list_.name(i),
                                                    Observable(var,
                                                               observable_plist.get<std::string>("region"),
                                                               observable_plist.get<std::string>("functional"),
                                                               observable_plist,
                                                               comm)));
    }
  }
}


void Unstructured_observations::MakeObservations(State& state)
{
  // loop over all observables
  for (std::map<std::string, Observable>::iterator i = observations.begin(); i != observations.end(); i++) {
    
    if ((i->second).DumpRequested(state.time()) || (i->second).DumpRequested(state.cycle())) {
      
      // for now we can only observe Integrals and Values
      if ( (i->second).functional != "Observation Data: Integral"  &&
	   (i->second).functional != "Observation Data: Point" )  {
	Errors::Message m("Unstructured_observations: can only handle Functional == Observation Data: Integral, or Functional == Observation Data: Point");
	Exceptions::amanzi_throw(m);
      }
      
      std::string label = i->first;
      
      // we need to make an observation for each variable in the observable
      std::string var = (i->second).variable;
      
      // data structure to store the observation
      Amanzi::ObservationData::DataTriple data_triplet;
      
      // build the name of the observation
      std::stringstream ss;
      ss << label << ", " << var;
      
      std::vector<Amanzi::ObservationData::DataTriple> &od = observation_data_[label]; 
      
      double value(0.0);
      double volume(0.0);
      
      // check if Aqueous concentration was requested
      std::string name = var;
      
      int pos = name.find("Aqueous concentration");
      if (pos != std::string::npos) {
	var = name.substr(0, pos-1);
      } else {
	var = name;
      }
      
      unsigned int mesh_block_size(0);
      Amanzi::AmanziMesh::Entity_ID_List entity_ids;
      if ((var == "Aqueous mass flow rate") || (var == "Aqueous volumetric flow rate")) { // for flux we need faces
	mesh_block_size = state.GetMesh()->get_set_size((i->second).region,
							Amanzi::AmanziMesh::FACE,
							Amanzi::AmanziMesh::OWNED);
	entity_ids.resize(mesh_block_size);
	state.GetMesh()->get_set_entities((i->second).region, Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED,
					  &entity_ids);
      } else { // all others need cells
	mesh_block_size = state.GetMesh()->get_set_size((i->second).region,
							Amanzi::AmanziMesh::CELL,
							Amanzi::AmanziMesh::OWNED);    
	entity_ids.resize(mesh_block_size);
	state.GetMesh()->get_set_entities((i->second).region, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
					  &entity_ids);
      }
      
      // find global meshblocksize
      int dummy = mesh_block_size; 
      int global_mesh_block_size(0);
      state.GetMesh()->get_comm()->SumAll(&dummy, &global_mesh_block_size, 1);
      
      if (global_mesh_block_size == 0) {
	// warn that this region is empty and bail
	Teuchos::OSTab tab = vo_->getOSTab();
	*vo_->os() << "Cannot make an observation on an empty region: " 
                   << (i->second).region << ", skipping" << std::endl;
	continue;
      }

      // is the user asking for a component concentration?
      int comp_index(0);
      if (comp_names_.size() > 0) {
	for (comp_index = 0; comp_index != comp_names_.size(); ++comp_index) {
	  if (comp_names_[comp_index] == var) break;
	}
      }
      
      // the user is asking to for an observation on tcc
      if (comp_names_.size() > 0 && comp_index != comp_names_.size() ) { 
	value = 0.0;
	volume = 0.0;
	
	Teuchos::RCP<const Epetra_MultiVector> total_component_concentration = 
	  state.GetFieldData("total_component_concentration")->ViewComponent("cell", false);
	
	for (int i=0; i<mesh_block_size; i++) {
	  int ic = entity_ids[i];
	  value += (*(*total_component_concentration)(comp_index))[ic] * state.GetMesh()->cell_volume(ic);
	  
	  volume += state.GetMesh()->cell_volume(ic);
	}
      } else if (var == "Volumetric water content") {
	value = 0.0;
	volume = 0.0;
	
	Teuchos::RCP<const Epetra_Vector> porosity = 
	  Teuchos::rcpFromRef(*(*state.GetFieldData("porosity")->ViewComponent("cell", false))(0));	  
	Teuchos::RCP<const Epetra_Vector> water_saturation = 
	  Teuchos::rcpFromRef(*(*state.GetFieldData("water_saturation")->ViewComponent("cell", false))(0));
	
	for (int i=0; i<mesh_block_size; i++) {
	  int ic = entity_ids[i];
	  value += (*porosity)[ic] * (*water_saturation)[ic] * state.GetMesh()->cell_volume(ic);
	  volume += state.GetMesh()->cell_volume(ic);
	}
      } else if (var == "Gravimetric water content") {
	value = 0.0;
	volume = 0.0;
	
	Teuchos::RCP<const Epetra_Vector> water_saturation = 
	  Teuchos::rcpFromRef(*(*state.GetFieldData("water_saturation")->ViewComponent("cell", false))(0));
	double water_density =  *state.GetScalarData("fluid_density");
	double particle_density(1.0); // does not exist in new state, yet... TODO
	Teuchos::RCP<const Epetra_Vector> porosity = 
	  Teuchos::rcpFromRef(*(*state.GetFieldData("porosity")->ViewComponent("cell", false))(0));
	
	for (int i=0; i<mesh_block_size; i++) {
	  int ic = entity_ids[i];
	  value += (*porosity)[ic] * (*water_saturation)[ic] * water_density 
	    / ( particle_density * (1.0 - (*porosity)[ic] ) )  * state.GetMesh()->cell_volume(ic);
	  volume += state.GetMesh()->cell_volume(ic);
	}    
      } else if (var == "Aqueous pressure") {
	value = 0.0;
	volume = 0.0;
	
	Teuchos::RCP<const Epetra_Vector> pressure = 
	  Teuchos::rcpFromRef(*(*state.GetFieldData("pressure")->ViewComponent("cell", false))(0));	  
	
	for (int i=0; i<mesh_block_size; i++) {
	  int ic = entity_ids[i];
	  value += (*pressure)[ic] * state.GetMesh()->cell_volume(ic);
	  volume += state.GetMesh()->cell_volume(ic);
	}
      } else if (var == "Aqueous saturation") {
	value = 0.0;
	volume = 0.0;
	
	Teuchos::RCP<const Epetra_Vector> water_saturation = 
	  Teuchos::rcpFromRef(*(*state.GetFieldData("water_saturation")->ViewComponent("cell", false))(0));	  
	
	for (int i=0; i<mesh_block_size; i++) {
	  int ic = entity_ids[i];
	  value += (*water_saturation)[ic] * state.GetMesh()->cell_volume(ic);
	  volume += state.GetMesh()->cell_volume(ic);
	}    
      } else if (var == "Hydraulic Head") {
	value = 0.0;
	volume = 0.0;
	
	Teuchos::RCP<const Epetra_Vector> hydraulic_head = 
	  Teuchos::rcpFromRef(*(*state.GetFieldData("hydraulic_head")->ViewComponent("cell", false))(0));
	
	for (int i=0; i<mesh_block_size; ++i) {
	  int ic = entity_ids[i];
	  Amanzi::AmanziGeometry::Point p = state.GetMesh()->cell_centroid(ic);
	  value += (*hydraulic_head)[ic] * state.GetMesh()->cell_volume(ic);
	  volume += state.GetMesh()->cell_volume(ic);
	}
      } else if (var == "Drawdown") {
	value = 0.0;
	volume = 0.0;
	
	Teuchos::RCP<const Epetra_Vector> hydraulic_head = 
	  Teuchos::rcpFromRef(*(*state.GetFieldData("hydraulic_head")->ViewComponent("cell", false))(0));
	
	for (int i=0; i<mesh_block_size; ++i) {
	  int ic = entity_ids[i];
	  Amanzi::AmanziGeometry::Point p = state.GetMesh()->cell_centroid(ic);
	  value += (*hydraulic_head)[ic] * state.GetMesh()->cell_volume(ic);
	  volume += state.GetMesh()->cell_volume(ic);
	}

        std::map<std::string, double>::iterator it = drawdown_.find(label);
        if (it == drawdown_.end()) { 
          drawdown_[label] = value;
          value = 0.0;
        } else {
          value -= it->second;
        }
      } else if ( (var == "Aqueous mass flow rate") || (var == "Aqueous volumetric flow rate")) {
	value = 0.0;
	volume = 0.0;
	
	// get the region object
	AmanziGeometry::GeometricModelPtr gm_ptr = state.GetMesh()->geometric_model();
	AmanziGeometry::RegionPtr reg_ptr = gm_ptr->FindRegion((i->second).region);
	AmanziGeometry::Point reg_normal;
	if (reg_ptr->type() == AmanziGeometry::POLYGON) {
	  AmanziGeometry::PolygonRegion *poly_reg = dynamic_cast<AmanziGeometry::PolygonRegion*>(reg_ptr);
	  reg_normal = poly_reg->normal();
	} else if (reg_ptr->type() == AmanziGeometry::PLANE) {
	  AmanziGeometry::PlaneRegion *plane_reg = dynamic_cast<AmanziGeometry::PlaneRegion*>(reg_ptr);
	  reg_normal = plane_reg->normal();
	} else {
	  // error
	  Exceptions::amanzi_throw(Errors::Message("Observations of Aqueous mass flow rate and Aqueous volumetric flow rate are only possible for Polygon and Plane regions"));
	}

	Teuchos::RCP<const Epetra_Vector> darcy_flux = 
	  Teuchos::rcpFromRef(*(*state.GetFieldData("darcy_flux")->ViewComponent("face", false))(0));      
	
	double density(1.0);
	if (var == "Aqueous mass flow rate") {
	  density = *state.GetScalarData("fluid_density");
	}
	
	for (int i = 0; i != mesh_block_size; ++i) {
	  int iface = entity_ids[i];
	  Amanzi::AmanziGeometry::Point face_normal = state.GetMesh()->face_normal(iface);
	  double area =  state.GetMesh()->face_area(iface);
	  double sign = reg_normal * face_normal / area;
	  
	  value += sign * (*darcy_flux)[iface] * density;
	  volume += area;
	}
	
      } else {
	std::stringstream ss;
	ss << "Unstructured_observations::make_observations: cannot make an observation for variable " << name;
	Errors::Message m(ss.str().c_str());
	Exceptions::amanzi_throw(m);
      }
      
      // syncronize the result across processors
      double result;
      state.GetMesh()->get_comm()->SumAll(&value,&result,1);
      
      double vresult;
      state.GetMesh()->get_comm()->SumAll(&volume,&vresult,1);
 
      if ((i->second).functional == "Observation Data: Integral") {  
	data_triplet.value = result;	
      } else if ((i->second).functional == "Observation Data: Point") {
	data_triplet.value = result/vresult;
      }
      
      data_triplet.is_valid = true;
      data_triplet.time = state.time();
      
      od.push_back(data_triplet);
    }
  }

  FlushObservations();
}


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


void Unstructured_observations::RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm) {
  // loop over all observations and register each of them with the time step manager...
  for (std::map<std::string, Observable>::iterator i = observations.begin(); i != observations.end(); i++) {
    (i->second).RegisterWithTimeStepManager(tsm);
  }  
}


void Unstructured_observations::FlushObservations()
{
  // print out observation file in ASCII format
  if (obs_list_.isParameter("Observation Output Filename")) {
    std::string obs_file = obs_list_.get<std::string>("Observation Output Filename");

    if (rank_ == 0) {
      std::ofstream out;
      out.open(obs_file.c_str(),std::ios::out);

      out.precision(16);
      out.setf(std::ios::scientific);

      out << "Observation Name, Region, Functional, Variable, Time, Value\n";
      out << "===========================================================\n";

      for (Teuchos::ParameterList::ConstIterator i=obs_list_.begin(); i!=obs_list_.end(); ++i) {
        std::string label = obs_list_.name(i);
        const Teuchos::ParameterEntry& entry = obs_list_.getEntry(label);
        if (entry.isList()) {
          const Teuchos::ParameterList& ind_obs_list = obs_list_.sublist(label);

          for (int j = 0; j < observation_data_[label].size(); j++) {
            if (observation_data_[label][j].is_valid) {
              if (!out.good()) {
                std::cout << "PROBLEM BEFORE" << std::endl;
              }
              out << label << ", "
                  << ind_obs_list.get<std::string>("region") << ", "
                  << ind_obs_list.get<std::string>("functional") << ", "
                  << ind_obs_list.get<std::string>("variable") << ", "
                  << observation_data_[label][j].time << ", "
                  << observation_data_[label][j].value << '\n';
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
