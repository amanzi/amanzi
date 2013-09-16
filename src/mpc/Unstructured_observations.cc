#include "Unstructured_observations.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "PolygonRegion.hh"

#include <map>

namespace Amanzi {

Unstructured_observations::Unstructured_observations(Teuchos::ParameterList observations_plist_,
                                                     Amanzi::ObservationData& observation_data_,
						     Epetra_MpiComm* comm):
    observation_data(observation_data_)
{
  // interpret paramerter list
  // loop over the sublists and create an observation for each
  for (Teuchos::ParameterList::ConstIterator i = observations_plist_.begin(); i != observations_plist_.end(); i++) {

    if (observations_plist_.isSublist(observations_plist_.name(i))) {
      Teuchos::ParameterList observable_plist = observations_plist_.sublist(observations_plist_.name(i));

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
      std::string var = observable_plist.get<string>("Variable");
      observations.insert(std::pair
                          <std::string, Observable>(observations_plist_.name(i),
                                                    Observable(var,
                                                               observable_plist.get<string>("Region"),
                                                               observable_plist.get<string>("Functional"),
							       observable_plist,
							       comm)));
    }
  }
}


void Unstructured_observations::make_observations(State& state)
{
  // loop over all observables
  for (std::map<std::string, Observable>::iterator i = observations.begin(); i != observations.end(); i++) {
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
    
    std::vector<Amanzi::ObservationData::DataTriple> &od = observation_data[label]; 
    
    double value(0.0);
    double volume(0.0);

    // check if Aqueous concentration was requested
    std::string name = var;
    
    int pos = name.find("Aqueous concentration");
    if (pos != string::npos) {
      var = name.substr(0, pos-1);
    } else {
      var = name;
    }
    
    unsigned int mesh_block_size(0);
    Amanzi::AmanziMesh::Entity_ID_List entity_ids;
    if (var == "Aqueous mass flux") { // for flux we need faces
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

    // is the user asking for a component concentration?
    int comp_index(0);
    if (comp_names_.size() > 0) {
      for (comp_index = 0; comp_index != comp_names_.size(); ++comp_index) {
	if (comp_names_[comp_index] == var) break;
      }
    }
    
    if (comp_names_.size() > 0 && comp_index != comp_names_.size() ) { // the user is asking to for an observation on tcc
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
    } else if (var == "Aqueous mass flux") {
      value = 0.0;
      volume = 0.0;
      
      // get the region object
      AmanziGeometry::GeometricModelPtr gm_ptr = state.GetMesh()->geometric_model();
      AmanziGeometry::RegionPtr reg_ptr = gm_ptr->FindRegion((i->second).region);
      if (reg_ptr->type() != AmanziGeometry::POLYGON) {

      }
      AmanziGeometry::PolygonRegion *poly_reg = dynamic_cast<AmanziGeometry::PolygonRegion*>(reg_ptr);
      
      AmanziGeometry::Point reg_normal = poly_reg->normal();

      Teuchos::RCP<const Epetra_Vector> darcy_flux = 
	Teuchos::rcpFromRef(*(*state.GetFieldData("darcy_flux")->ViewComponent("face", false))(0));      
      
      for (int i = 0; i != mesh_block_size; ++i) {
	int iface = entity_ids[i];
	Amanzi::AmanziGeometry::Point face_normal = state.GetMesh()->face_normal(iface);
	double sign = reg_normal * face_normal;
	double area =  state.GetMesh()->face_area(iface);
	
	value += sign * (*darcy_flux)[iface] * area;
	volume += area;
      }
     
    } else {
      std::stringstream ss;
      ss << "State::point_value: cannot make an observation for variable " << name;
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



}  // namespace Amanzi
