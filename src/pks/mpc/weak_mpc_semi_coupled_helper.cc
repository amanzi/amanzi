#include "weak_mpc_semi_coupled_helper.hh"

namespace Amanzi{

void
UpdateIntermediateStateParameters(Teuchos::RCP<Amanzi::State>& S_next_, Teuchos::RCP<Amanzi::State>& S_inter_, int id){

  std::stringstream name, name_ss;
     
  name << "column_" << id <<"_surface";
  name_ss << "column_" << id;
      
	  *S_inter_->GetFieldData(Keys::getKey(name_ss.str(),"pressure"), S_inter_->GetField(Keys::getKey(name_ss.str(),"pressure"))->owner()) = 
	    *S_next_->GetFieldData(Keys::getKey(name_ss.str(),"pressure"));
					     
	  
	  *S_inter_->GetFieldData(Keys::getKey(name_ss.str(),"temperature"), S_inter_->GetField(Keys::getKey(name_ss.str(),"temperature"))->owner()) = 
	    *S_next_->GetFieldData(Keys::getKey(name_ss.str(),"temperature")); 
	  
	  *S_inter_->GetFieldData(Keys::getKey(name.str(),"pressure"), S_inter_->GetField(Keys::getKey(name.str(),"pressure"))->owner()) = 
	    *S_next_->GetFieldData(Keys::getKey(name.str(),"pressure"));
					     
	  
	  *S_inter_->GetFieldData(Keys::getKey(name.str(),"temperature"), S_inter_->GetField(Keys::getKey(name.str(),"temperature"))->owner()) = 
	    *S_next_->GetFieldData(Keys::getKey(name.str(),"temperature"));
	  
	  *S_inter_->GetFieldData(Keys::getKey(name.str(),"snow_depth"), S_inter_->GetField(Keys::getKey(name.str(),"snow_depth"))->owner()) = 
	    *S_next_->GetFieldData(Keys::getKey(name.str(),"snow_depth"));

	   *S_inter_->GetFieldData(Keys::getKey(name.str(),"conducted_energy_source"), 
				   S_inter_->GetField(Keys::getKey(name.str(),"conducted_energy_source"))->owner()) = 
	     *S_next_->GetFieldData(Keys::getKey(name.str(),"conducted_energy_source"));
	  
	   *S_inter_->GetFieldData(Keys::getKey(name.str(),"water_source"), 
				   S_inter_->GetField(Keys::getKey(name.str(),"water_source"))->owner()) = 
	     *S_next_->GetFieldData(Keys::getKey(name.str(),"water_source"));

	   *S_inter_->GetFieldData(Keys::getKey(name_ss.str(),"water_source"), 
				   S_inter_->GetField(Keys::getKey(name_ss.str(),"water_source"))->owner()) = 
	     *S_next_->GetFieldData(Keys::getKey(name_ss.str(),"water_source"));

	   *S_inter_->GetFieldData(Keys::getKey(name.str(),"water_source_temperature"), 
				   S_inter_->GetField(Keys::getKey(name.str(),"water_source_temperature"))->owner()) = 
	    *S_next_->GetFieldData(Keys::getKey(name.str(),"water_source_temperature"));

	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"snow_temperature"), 
				   S_inter_->GetField(Keys::getKey(name.str(),"snow_temperature"))->owner()) = 
	    *S_next_->GetFieldData(Keys::getKey(name.str(),"snow_temperature"));
	  	  
	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"snow_age"), 
				   S_inter_->GetField(Keys::getKey(name.str(),"snow_age"))->owner()) = 
	    *S_next_->GetFieldData(Keys::getKey(name.str(),"snow_age"));
	    
	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"snow_density"), 
				   S_inter_->GetField(Keys::getKey(name.str(),"snow_density"))->owner()) = 
	    *S_next_->GetFieldData(Keys::getKey(name.str(),"snow_density"));

	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"stored_SWE"), 
				   S_inter_->GetField(Keys::getKey(name.str(),"stored_SWE"))->owner()) = 
	    *S_next_->GetFieldData(Keys::getKey(name.str(),"stored_SWE"));

	   *S_inter_->GetFieldData(Keys::getKey(name.str(),"surface_subsurface_flux"), 
				  S_inter_->GetField(Keys::getKey(name.str(),"surface_subsurface_flux"))->owner()) = 
	     *S_next_->GetFieldData(Keys::getKey(name.str(),"surface_subsurface_flux"));
	  
	   *S_inter_->GetFieldData(Keys::getKey(name.str(),"surface_subsurface_energy_flux"), 
				  S_inter_->GetField(Keys::getKey(name.str(),"surface_subsurface_energy_flux"))->owner()) = 
	     *S_next_->GetFieldData(Keys::getKey(name.str(),"surface_subsurface_energy_flux"));
	  
	   *S_inter_->GetFieldData(Keys::getKey(name_ss.str(),"saturation_liquid"), 
				   S_inter_->GetField(Keys::getKey(name_ss.str(),"saturation_liquid"))->owner()) = 
	     *S_next_->GetFieldData(Keys::getKey(name_ss.str(),"saturation_liquid"));
	   
	   *S_inter_->GetFieldData(Keys::getKey(name.str(),"unfrozen_fraction"), 
				   S_inter_->GetField(Keys::getKey(name.str(),"unfrozen_fraction"))->owner()) = 
	     *S_next_->GetFieldData(Keys::getKey(name.str(),"unfrozen_fraction"));
	   
	   *S_inter_->GetFieldData(Keys::getKey(name.str(),"porosity"), 
				   S_inter_->GetField(Keys::getKey(name.str(),"porosity"))->owner()) = 
	     *S_next_->GetFieldData(Keys::getKey(name.str(),"porosity"));

	   *S_inter_->GetFieldData(Keys::getKey(name.str(),"evaporative_flux"), 
				   S_inter_->GetField(Keys::getKey(name.str(),"evaporative_flux"))->owner()) = 
	     *S_next_->GetFieldData(Keys::getKey(name.str(),"evaporative_flux"));
	   

}

void
UpdateNextStateParameters(Teuchos::RCP<Amanzi::State>& S_next_, Teuchos::RCP<Amanzi::State>& S_inter_, int id){
  std::stringstream name, name_ss;
      
      name << "column_" << id <<"_surface";
      name_ss << "column_" << id;
      
  
	  *S_next_->GetFieldData(Keys::getKey(name_ss.str(),"pressure"), S_next_->GetField(Keys::getKey(name_ss.str(),"pressure"))->owner()) = 
	    *S_inter_->GetFieldData(Keys::getKey(name_ss.str(),"pressure"));
					     
	  
	  *S_next_->GetFieldData(Keys::getKey(name_ss.str(),"temperature"), S_next_->GetField(Keys::getKey(name_ss.str(),"temperature"))->owner()) = 
	    *S_inter_->GetFieldData(Keys::getKey(name_ss.str(),"temperature")); 
	  
	  *S_next_->GetFieldData(Keys::getKey(name.str(),"pressure"), S_next_->GetField(Keys::getKey(name.str(),"pressure"))->owner()) = 
	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"pressure"));
					     
	  
	  *S_next_->GetFieldData(Keys::getKey(name.str(),"temperature"), S_next_->GetField(Keys::getKey(name.str(),"temperature"))->owner()) = 
	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"temperature"));
	  
	  *S_next_->GetFieldData(Keys::getKey(name.str(),"snow_depth"), S_next_->GetField(Keys::getKey(name.str(),"snow_depth"))->owner()) = 
	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"snow_depth"));

	   *S_next_->GetFieldData(Keys::getKey(name.str(),"conducted_energy_source"), 
				   S_next_->GetField(Keys::getKey(name.str(),"conducted_energy_source"))->owner()) = 
	     *S_inter_->GetFieldData(Keys::getKey(name.str(),"conducted_energy_source"));
	  
	   *S_next_->GetFieldData(Keys::getKey(name.str(),"water_source"), 
				   S_next_->GetField(Keys::getKey(name.str(),"water_source"))->owner()) = 
	     *S_inter_->GetFieldData(Keys::getKey(name.str(),"water_source"));

	   *S_next_->GetFieldData(Keys::getKey(name_ss.str(),"water_source"), 
				   S_next_->GetField(Keys::getKey(name_ss.str(),"water_source"))->owner()) = 
	     *S_inter_->GetFieldData(Keys::getKey(name_ss.str(),"water_source"));

	   *S_next_->GetFieldData(Keys::getKey(name.str(),"water_source_temperature"), 
				   S_next_->GetField(Keys::getKey(name.str(),"water_source_temperature"))->owner()) = 
	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"water_source_temperature"));

	    *S_next_->GetFieldData(Keys::getKey(name.str(),"snow_temperature"), 
				   S_next_->GetField(Keys::getKey(name.str(),"snow_temperature"))->owner()) = 
	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"snow_temperature"));
	  	  
	    *S_next_->GetFieldData(Keys::getKey(name.str(),"snow_age"), 
				   S_next_->GetField(Keys::getKey(name.str(),"snow_age"))->owner()) = 
	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"snow_age"));
	    
	    *S_next_->GetFieldData(Keys::getKey(name.str(),"snow_density"), 
				   S_next_->GetField(Keys::getKey(name.str(),"snow_density"))->owner()) = 
	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"snow_density"));

	    *S_next_->GetFieldData(Keys::getKey(name.str(),"stored_SWE"), 
				   S_next_->GetField(Keys::getKey(name.str(),"stored_SWE"))->owner()) = 
	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"stored_SWE"));

	   *S_next_->GetFieldData(Keys::getKey(name.str(),"surface_subsurface_flux"), 
				  S_next_->GetField(Keys::getKey(name.str(),"surface_subsurface_flux"))->owner()) = 
	     *S_inter_->GetFieldData(Keys::getKey(name.str(),"surface_subsurface_flux"));
	  
	   *S_next_->GetFieldData(Keys::getKey(name.str(),"surface_subsurface_energy_flux"), 
				  S_next_->GetField(Keys::getKey(name.str(),"surface_subsurface_energy_flux"))->owner()) = 
	     *S_inter_->GetFieldData(Keys::getKey(name.str(),"surface_subsurface_energy_flux"));
	  
	  

	   *S_next_->GetFieldData(Keys::getKey(name_ss.str(),"saturation_liquid"), 
				  S_next_->GetField(Keys::getKey(name_ss.str(),"saturation_liquid"))->owner()) = 
	     *S_inter_->GetFieldData(Keys::getKey(name_ss.str(),"saturation_liquid"));

	   *S_next_->GetFieldData(Keys::getKey(name.str(),"unfrozen_fraction"), 
				   S_next_->GetField(Keys::getKey(name.str(),"unfrozen_fraction"))->owner()) = 
	     *S_inter_->GetFieldData(Keys::getKey(name.str(),"unfrozen_fraction"));
	   
	   *S_next_->GetFieldData(Keys::getKey(name.str(),"porosity"), 
				   S_next_->GetField(Keys::getKey(name.str(),"porosity"))->owner()) = 
	    *S_inter_->GetFieldData(Keys::getKey(name.str(),"porosity"));

	   *S_next_->GetFieldData(Keys::getKey(name.str(),"evaporative_flux"), 
				  S_next_->GetField(Keys::getKey(name.str(),"evaporative_flux"))->owner()) = 
	     *S_inter_->GetFieldData(Keys::getKey(name.str(),"evaporative_flux"));


}




}


