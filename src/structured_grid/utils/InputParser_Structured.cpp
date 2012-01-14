#include "InputParser_Structured.H"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

using Teuchos::Array;
using Teuchos::ParameterList;
using Teuchos::ParameterEntry;

namespace Amanzi {
  namespace AmanziInput {

    std::string underscore(const std::string& instring)
    {
      std::string s = instring;
      std::replace(s.begin(),s.end(),' ','_');
      return s;
    }

    //
    // convert parameterlist to format for structured code
    //
    ParameterList
    convert_to_structured(const ParameterList& parameter_list)
    {
      ParameterList struc_list = setup_structured();

      // determine spatial dimension of the problem
      ndim = parameter_list.sublist("Domain").get<int>("Spatial Dimension");
      
      // Note: Structured starts each run time 0 except for restart.
      //       There is only stop_time = End_Time-Start_Time.
      ParameterList eclist = parameter_list.sublist("Execution control");
      double Start_Time = eclist.get<double>("Start Time");
      double End_Time = eclist.get<double>("End Time");
      simulation_time = End_Time - Start_Time;
      struc_list.set("stop_time",End_Time-Start_Time);
      //
      // Mesh
      //
      convert_to_structured_mesh(parameter_list,struc_list);
      //
      // Execution control
      //
      convert_to_structured_control(parameter_list,struc_list);
      //
      // Regions
      //
      convert_to_structured_region(parameter_list, struc_list);
      //
      // Materials
      //
      convert_to_structured_material(parameter_list, struc_list);
      //
      // Components 
      //
      convert_to_structured_comp(parameter_list, struc_list);
      //
      // Tracers
      //
      convert_to_structured_tracer(parameter_list, struc_list);
      //
      // Output
      // 
      convert_to_structured_output(parameter_list, struc_list);

      return struc_list;
    }

    //
    // setup sublists of struc_list
    //
    ParameterList
    setup_structured()
    {
      ParameterList struc_list, empty_list;      
      struc_list.set("Mesh", empty_list);
      struc_list.set("amr" , empty_list);
      struc_list.set("cg"  , empty_list);
      struc_list.set("mg"  , empty_list);
      struc_list.set("mac" , empty_list);
      struc_list.set("comp", empty_list);
      struc_list.set("phase", empty_list);
      struc_list.set("press", empty_list);
      struc_list.set("prob" , empty_list); 
      struc_list.set("rock" , empty_list);
      struc_list.set("diffuse" , empty_list);
      struc_list.set("geometry", empty_list);
      struc_list.set("source"  , empty_list);
      struc_list.set("tracer"  , empty_list); 
      struc_list.set("observation", empty_list);

      return struc_list;

    }

    //
    // convert mesh to structured format
    //
    void
    convert_to_structured_mesh (const ParameterList& parameter_list, 
				ParameterList&       struc_list)
    {
      // Mesh
      struc_list.sublist("Mesh").set("Framework","Structured");
      const ParameterList& stlist = parameter_list.sublist("Mesh").sublist("Structured");
      for (ParameterList::ConstIterator j=stlist.begin(); j!=stlist.end(); ++j) 
	{
	  const std::string& rlabel = stlist.name(j);
	  const ParameterEntry& rentry = stlist.getEntry(rlabel);
	  if (rlabel == "Number of Cells")
	    struc_list.sublist("amr").setEntry("n_cell",rentry);
	  else if (rlabel == "Domain Low Corner")
	    domlo = stlist.get<Array<double> >(rlabel);
	  else if (rlabel == "Domain High Corner")
	    domhi = stlist.get<Array<double> >(rlabel);
	}

      const ParameterList& eclist = parameter_list.sublist("Execution control");
      int bfactor = 2;
      if (eclist.isSublist("amr"))
	if (eclist.sublist("amr").isParameter("blocking_factor"))
	  bfactor = eclist.sublist("amr").get<int>("blocking_factor");
     
      Array<int> n_cell = struc_list.sublist("amr").get<Array<int> >("n_cell");
      for (int i=0;i<ndim;i++) {
	if (n_cell[i]%bfactor > 0) {
	  std::cerr << "Number of Cells must be divisible by " << bfactor << std::endl;
	  throw std::exception();
	}
      }

      ParameterList& glist = struc_list.sublist("geometry");
      glist.set("prob_lo",domlo);
      glist.set("prob_hi",domhi);

      // these are by default
      Array<int> is_periodic(ndim,0);
      glist.set("is_periodic",is_periodic);
      glist.set("coord_sys","0");
    }

    //
    // convert execution control to structured format
    //
    void
    convert_to_structured_control(const ParameterList& parameter_list, 
				  ParameterList&       struc_list)
    {

      ParameterList& amr_list     = struc_list.sublist("amr");
      ParameterList& prob_list    = struc_list.sublist("prob");
      ParameterList& cg_list      = struc_list.sublist("cg");
      ParameterList& mg_list      = struc_list.sublist("mg");
      ParameterList& mac_list     = struc_list.sublist("mac");
      ParameterList& diffuse_list = struc_list.sublist("diffuse");
      
      const ParameterList& eclist =
	parameter_list.sublist("Execution control");

      if (eclist.isSublist("prob"))
	prob_list = eclist.sublist("prob");

      if (eclist.isSublist("amr")) 
	amr_list = eclist.sublist("amr");

      if (eclist.isSublist("mg"))
	mg_list = eclist.sublist("mg");

      if (eclist.isSublist("cg"))
	cg_list = eclist.sublist("cg");
    
      if (eclist.isSublist("mac"))
	mac_list = eclist.sublist("mac");

      if (eclist.isSublist("mac"))
	diffuse_list = eclist.sublist("diffuse");

      std::string flowmode = eclist.get<std::string>("Flow Mode");
      if (!flowmode.compare("transient single phase variably saturated flow") ||
	  !flowmode.compare("steady state single phase variably saturated flow"))
	{
	  prob_list.set("model_name","richard");
	  prob_list.set("have_capillary",1);
	  prob_list.set("visc_abs_tol",1.e-16);
	  prob_list.set("visc_tol",1.e-14);
	  prob_list.set("cfl",1);
	  prob_list.set("richard_max_dt",prob_list.get<double>("max_dt"));
	}
      else if (!flowmode.compare("transient single phase saturated flow"))
	{
	  prob_list.set("model_name","single-phase");
	  prob_list.set("cfl",0.95);
	  prob_list.set("do_simple",1);
	}
      else if (!flowmode.compare("steady state single phase saturated flow"))
	{
	  prob_list.set("model_name","single-phase");
	  prob_list.set("cfl",0.95);
	  prob_list.set("do_simple",2);
	}
      else if (!flowmode.compare("transient two phase flow"))
	{
	  prob_list.set("model_name","two-phase");
	  prob_list.set("cfl",0.75);
	}

      std::string chemmode = eclist.get<std::string>("Chemistry Mode");
      if (!chemmode.compare("none"))
	prob_list.set("do_chem",-1);

      //AMR
      const ParameterList& stlist = parameter_list.sublist("Mesh").sublist("Structured");
      amr_list.set<Array<int> >("n_cell",stlist.get<Array<int> >("Number of Cells"));
      int nlevel = amr_list.get<int>("max_level",0);
      if (nlevel == 0) nlevel = 1;
      Array<int> n_buf(nlevel,2);
      if (!amr_list.isParameter("ref_ratio"))
	amr_list.set("ref_ratio",n_buf);
    }

    //
    // convert region to structured format
    //
    void
    convert_to_structured_region(const ParameterList& parameter_list, 
				 ParameterList&       struc_list)
    {
      ParameterList& geom_list = struc_list.sublist("geometry");
      
      const ParameterList& rlist = parameter_list.sublist("Regions");
      Array<std::string> arrayregions;
      for (ParameterList::ConstIterator i=rlist.begin(); i!=rlist.end(); ++i) {
        
	std::string label = rlist.name(i);
	std::string _label = underscore(label);
	const ParameterEntry& entry = rlist.getEntry(label);
        
        if ( !entry.isList() ) {
          if (Teuchos::MPISession::getRank() == 0) {
            std::cerr << "Region section must define only regions  " 
		      << label << " is not a valid region definition." << std::endl;
          }
          throw std::exception();
        }
        
	ParameterList rsublist;
        const ParameterList& rslist = rlist.sublist(label);
        for (ParameterList::ConstIterator j=rslist.begin(); j!=rslist.end(); ++j) {
          
          const std::string& rlabel = rslist.name(j);
          const ParameterEntry& rentry = rslist.getEntry(rlabel);
          
          if (rentry.isList()) {
            if (rlabel=="Region: Color Function") {
	      const ParameterList& rsslist = rslist.sublist(rlabel);
	      rsublist.setEntry("color_file",rsslist.getEntry("File"));
	      rsublist.setEntry("color_value",rsslist.getEntry("Value"));
	      rsublist.set("purpose", "all");
	      rsublist.set("type", "color_function");
            }
	    else if (rlabel=="Region: Point"){	      
	      const ParameterList& rsslist = rslist.sublist(rlabel);
	      rsublist.setEntry("coordinate",rsslist.getEntry("Coordinate"));
	      rsublist.set("purpose", "all");
	      rsublist.set("type", "point");
            }
            else if (rlabel=="Region: Box") {
	      const ParameterList& rsslist = rslist.sublist(rlabel);
	      rsublist.setEntry("lo_coordinate",rsslist.getEntry("Low Coordinate"));
	      rsublist.setEntry("hi_coordinate",rsslist.getEntry("High Coordinate"));
	      rsublist.set("purpose", "all");
	      rsublist.set("type", "box");
            }
	    else if (rlabel=="Region: Plane") {
	      const ParameterList& rsslist = rslist.sublist(rlabel);
	      Array<double> low_coordinate(ndim),hi_coordinate(ndim);
	      Array<double> direction = rsslist.get<Array<double> >("Direction");
	      Array<double> coordinate = rsslist.get<Array<double> >("Coordinate");
	      for (int dir=0; dir<ndim; dir++) {
		if (abs(direction[dir] - 1.0) < 1.e-6 || abs(direction[dir] + 1) < 1.e-6){
		  low_coordinate[dir] = coordinate[dir];
		  hi_coordinate[dir]  = coordinate[dir];
		}
		else {
		  low_coordinate[dir] = domlo[dir];
		  hi_coordinate[dir]  = domhi[dir];
		}
		}
	      rsublist.set("low_coordinate",low_coordinate);
	      rsublist.set("hi_coordinate",hi_coordinate);
	      rsublist.set("purpose","all");
	      rsublist.set("type", "box");	      
	    }
	  }
	  else {
            std::cerr << rlabel << " is not a valid region type for structured.\n";
	    throw std::exception();
	  }
	  
	  geom_list.set(_label,rsublist);
	  // need to remove empty spaces
	  arrayregions.push_back(_label); 
	}
	geom_list.set("region",arrayregions);
      }
    }

    //
    // convert material to structured format
    //
    void
    convert_to_structured_material(const ParameterList& parameter_list, 
				   ParameterList&       struc_list)
    {
      ParameterList& rock_list = struc_list.sublist("rock");

      const ParameterList& rlist = parameter_list.sublist("Material Properties");
      Array<std::string> arrayrock;
      for (ParameterList::ConstIterator i=rlist.begin(); i!=rlist.end(); ++i) {
        
	std::string label = rlist.name(i);
	std::string _label = underscore(label);
	const ParameterEntry& entry = rlist.getEntry(label);
        
        if (entry.isList()) {
	  ParameterList rsublist;
	  const ParameterList& rslist = rlist.sublist(label);
	  for (ParameterList::ConstIterator j=rslist.begin(); j!=rslist.end(); ++j) {
	    
	    const std::string& rlabel = rslist.name(j);
	    const ParameterEntry& rentry = rslist.getEntry(rlabel);
	      
	    if (rentry.isList()) {
	      const ParameterList& rsslist = rslist.sublist(rlabel);
	      if (rlabel=="Porosity: Uniform"){
		rsublist.setEntry("porosity",rsslist.getEntry("Value"));
		rsublist.set("porosity_dist","uniform");
	      }
	      else if (rlabel=="Intrinsic Permeability: Anisotropic Uniform") {
		Array<double> array_p(2);
		array_p[0] = rsslist.get<double>("Horizontal");
		array_p[1] = rsslist.get<double>("Vertical");
		// convert from m^2 to mDa
		for (int k=0; k<2; k++) array_p[k] /= 9.869233e-16;
		rsublist.set("permeability",array_p);
		rsublist.set("permeability_dist","uniform");
	      }
	      else if (rlabel=="Capillary Pressure: van Genuchten") {
		Array<double> array_c(4);
		Array<double> array_k(3);
		// convert Pa^-1 to atm^-1
		double alpha = rsslist.get<double>("alpha");
	      array_c[1] = alpha*1.01325e5; 
	      array_c[0] = rsslist.get<double>("m");
	      array_c[2] = rsslist.get<double>("Sr");
	      array_c[3] = 0.0;
	      
	      array_k[0] = array_c[0];
	      array_k[1] = array_c[2];
	      array_k[2] = array_c[3];
	      int cpl_type = 3;
	      rsublist.set("cpl_type",cpl_type);
	      rsublist.set("kr_type",cpl_type);
	      rsublist.set("cpl_param",array_c);
	      rsublist.set("kr_param",array_k);
	      } 
	    }
	    else {
	      Array<std::string> tmp_region = rslist.get<Array<std::string> >(rlabel);
	      for (Array<std::string>::iterator it=tmp_region.begin();
		   it !=tmp_region.end(); it++)
		(*it) = underscore(*it);
	      rsublist.set("region",tmp_region);
	    }
	    rsublist.set("density",1e3);
	    rock_list.set(_label,rsublist);
	  }
	  // need to remove empty spaces
	  arrayrock.push_back(_label);
	}
      }
      rock_list.set("rock",arrayrock);
      std::string kp_file="kp";
      std::string pp_file="pp";

      if (rlist.isParameter("Permeability Output File"))
	kp_file = rlist.get<std::string>("Permeability Output File");
      if (rlist.isParameter("Porosity Output File"))
	pp_file = rlist.get<std::string>("Porosity Output File");
      rock_list.set("permeability_file",kp_file);
      rock_list.set("porosity_file",pp_file);
    }
 
    //
    // convert component to structured format
    //
    void
    convert_to_structured_comp(const ParameterList& parameter_list, 
			       ParameterList&       struc_list)
    {

      ParameterList& phase_list = struc_list.sublist("phase");
      ParameterList& comp_list  = struc_list.sublist("comp"); 
      ParameterList& press_list = struc_list.sublist("press");

      const ParameterList& rlist = parameter_list.sublist("Phase Definitions");
      Array<std::string> arrayphase;
      std::map<std::string,double> viscosity,density;

      // get Phase Properties
      for (ParameterList::ConstIterator i=rlist.begin(); i!=rlist.end(); ++i) {
        
	std::string label = rlist.name(i);
	std::string _label = underscore(label);
	const ParameterEntry& entry = rlist.getEntry(label);
        
        if ( !entry.isList() ) {
          if (Teuchos::MPISession::getRank() == 0) {
            std::cerr << "Phase section must define only phases  " 
		      << label << " is not a valid phase definition." << std::endl;
          }
          throw std::exception();
        }
	arrayphase.push_back(_label);	
        
        const ParameterList& rslist = rlist.sublist(label);
	if (rslist.isSublist("Phase Properties")) {
	  const ParameterList rsslist = rslist.sublist("Phase Properties");
	  for (ParameterList::ConstIterator k=rsslist.begin(); k!=rsslist.end(); ++k) {
	    const std::string& rrlabel = rsslist.name(k);
	    const ParameterEntry& rrentry = rsslist.getEntry(rrlabel);
	    
	    if (entry.isList()) {
	      const ParameterList& rssslist = rsslist.sublist(rrlabel);
	      if (rrlabel == "Viscosity: Uniform")
		{
		  viscosity[label] = rssslist.get<double>("Viscosity");
		  // convert to cP
		  viscosity[label] *= 1.e3;
		}
	      else if (rrlabel == "Density: Uniform")
		density[label] = rssslist.get<double>("Density");
	    }
	  }
	}
	else
	  std::cerr << "Phase properties are not defined\n";
      }
      phase_list.set("phase",arrayphase);

      // get Phase Components
      Array<std::string> array_comp;
      for (ParameterList::ConstIterator i=rlist.begin(); i!=rlist.end(); ++i) {        
	std::string label = rlist.name(i);
        const ParameterList& rslist = rlist.sublist(label);
	if (rslist.isSublist("Phase Components")) {
	  const ParameterList& rsslist = rslist.sublist("Phase Components");
	  for (ParameterList::ConstIterator j=rsslist.begin(); j!=rsslist.end(); ++j) {
	    ParameterList csublist;
	    std::string rlabel = rsslist.name(j);
	    std::string _rlabel = underscore(rlabel);
	    csublist.set("phase",label);
	    csublist.set("density",density[label]);
	    csublist.set("viscosity",viscosity[label]);
	    csublist.set("diffusivity",0.0);
	    comp_list.set(_rlabel,csublist);
	    array_comp.push_back(_rlabel);
	  }
	}
      }
      comp_list.set("comp",array_comp);

      // get initial conditions
      const ParameterList& ilist = parameter_list.sublist("Initial Conditions");
      Array<std::string> init_region;

      for (ParameterList::ConstIterator i=ilist.begin(); i!=ilist.end(); ++i) {        
	std::string label = ilist.name(i);
	std::string _label = underscore(label);
	init_region.push_back(_label);
	const ParameterEntry& entry = ilist.getEntry(label);
        
	ParameterList rsublist;
        const ParameterList& rslist = ilist.sublist(label);
        for (ParameterList::ConstIterator j=rslist.begin(); j!=rslist.end(); ++j) {
          const std::string& rlabel = rslist.name(j);	
	  Array<std::string > regions = rslist.get<Array<std::string> >("Assigned Regions");
	  Array<std::string> _regions;
	  for (Array<std::string>::iterator it=regions.begin(); it!=regions.end(); it++) {
	    _regions.push_back(underscore(*it));
	  }
	  rsublist.set("region",_regions);

          const ParameterEntry& rentry = rslist.getEntry(rlabel);
          if (rentry.isList()) {
	    const ParameterList& rsslist = rslist.sublist(rlabel);
	    if (rlabel=="IC: Linear Pressure"){
	      Array<double> ref_coor = rsslist.get<Array<double> >("Reference Coordinate");
	      rsublist.set("type","hydrostatic");
	      rsublist.set("water_table",ref_coor[ref_coor.size()-1]);
	    }
	    else if (rlabel=="IC: Uniform Saturation") {
	      rsublist.set("type","scalar");
	      double sat = rsslist.get<double>("Value");
	      rsublist.set(array_comp[0],density[array_comp[0]]);
 	    }
	    else if (rlabel=="IC: Flux") {
	      rsublist.set("type","zero_total_velocity");
	      rsublist.set("inflow",rsslist.get<double>("inflow velocity"));
	    }
	    else if (rlabel=="IC: Scalar") {
	      for (Array<std::string>::iterator k=array_comp.begin(); k!=array_comp.end(); ++k)
		rsublist.set(*k,rsslist.get<double>(*k));
	    }
	    else if (rlabel=="IC: Uniform Pressure") {
	      std::cerr << rlabel << " does not work with Amanzi-S yet.\n";
	      throw std::exception();
	    }
	      
	    comp_list.set(_label,rsublist);
	  }
	}
      }
      comp_list.set("init",init_region);

      // Boundary conditions
      // Require information related to the regions and materials.  
      ParameterList regionlist = parameter_list.sublist("Regions");

      Array<std::string> bc_name_array;
      Array<bool> hi_def(ndim,false), lo_def(ndim,false);
      Array<int> hi_bc(ndim,4), lo_bc(ndim,4);
      Array<int> p_hi_bc(ndim,1), p_lo_bc(ndim,1);
      Array<int> inflow_hi_bc(ndim,0), inflow_lo_bc(ndim,0);
      Array<double> inflow_hi_vel(ndim,0), inflow_lo_vel(ndim,0);
      Array<double> press_hi(ndim,0), press_lo(ndim,0);
      double water_table_hi=0;
      double water_table_lo=0;

      if (parameter_list.isSublist("Boundary Conditions"))
	{
	  const ParameterList& blist = parameter_list.sublist("Boundary Conditions");
	  for (ParameterList::ConstIterator i=blist.begin(); i!=blist.end(); ++i) {
        
	    std::string label = blist.name(i);
	    std::string _label = underscore(label);
	    const ParameterEntry& entry = blist.getEntry(label);
        
	    ParameterList rsublist;
	    const ParameterList& rslist = blist.sublist(label);
	    for (ParameterList::ConstIterator j=rslist.begin(); j!=rslist.end(); ++j) {
	      
	      const std::string& rlabel = rslist.name(j);
	      Array<std::string > regions = rslist.get<Array<std::string> >("Assigned Regions");
	      
	      Array<std::string > valid_region_def_lo, valid_region_def_hi;
	      valid_region_def_lo.push_back("XLOBC");
	      valid_region_def_hi.push_back("XHIBC");
	      valid_region_def_lo.push_back("YLOBC");
	      valid_region_def_hi.push_back("YHIBC");
	      if (ndim == 3) {
		valid_region_def_lo.push_back("ZLOBC");
		valid_region_def_hi.push_back("ZHIBC");
	      }
	      Array<std::string> _regions;
	      for (Array<std::string>::iterator it=regions.begin(); it!=regions.end(); it++) {
		_regions.push_back(underscore(*it));
		bool valid = false;
		for (int j = 0; j<ndim; j++) {
		  if (!(*it).compare(valid_region_def_lo[j]))
		    valid = true;
		  else if (!(*it).compare(valid_region_def_hi[j]))
		    valid = true;
		}
		if (!valid)
		  std::cerr << "Structured: boundary conditions can only be applied to "
			    << "regions labeled as XLOBC, XHIBC, YLOBC, YHIBC, ZLOBC and ZHIBC\n";
	      }
	      
	      rsublist.set("region",_regions);
	      const ParameterEntry& rentry = rslist.getEntry(rlabel);
	      if (rentry.isList()) {
		const ParameterList& rsslist = rslist.sublist(rlabel);
		if (rlabel=="BC: Flux"){
		  // don't do time-dependent flux
		  Array<double> flux;
		  if (rsslist.isParameter("Intensive Volumetric Flux"))
		    flux = rsslist.get<Array<double> >("Intensive Volumetric Flux");
		  else if (rsslist.isParameter("Intensive Mass Flux"))
		    {
		      flux = rsslist.get<Array<double> >("Intensive Mass Flux");
		      for (int i=0;i<flux.size();i++)
			flux[i] /= 1.e3;
		    }
		  if (!rsslist.isParameter("Material Type at Boundary"))
		    {
		      // this is temporary.  It will be removed once the structured code is updated.
		      const ParameterList& mat_list = parameter_list.sublist("Material Properties"); 
		      std::string first_material = mat_list.name(mat_list.begin());
		      rsublist.set("rock",first_material);
		    }
		  
		  rsublist.set("inflow",-flux[0]);
		  rsublist.set("type","zero_total_velocity");
		  comp_list.set(_label,rsublist);
		  bc_name_array.push_back(_label);
		  
		  for (Array<std::string>::iterator it=regions.begin(); 
		       it!=regions.end(); it++) {
		    for (int j=0; j<ndim; j++) {
		      if (!(*it).compare(valid_region_def_lo[j])) {
			lo_def[j] = true;
			lo_bc[j] = 1;
			p_lo_bc[j] = 1;
			inflow_lo_bc[j] = 1;
			inflow_lo_vel[j] = flux[0];
		      }
		      else if (!(*it).compare(valid_region_def_hi[j])) {
			hi_def[j] = true;
			hi_bc[j] = 1;
			p_hi_bc[j] = 1;
			inflow_hi_bc[j] = 1;
			inflow_hi_vel[j] = -flux[0];
		      }
		    }
		  }
		}
		else if (rlabel=="BC: Uniform Pressure") {
		  // set to saturated condition
		  rsublist.set("Water",density[arrayphase[0]]);
		  rsublist.set("type","scalar");
		  comp_list.set(_label,rsublist);
		  bc_name_array.push_back(_label);
		  
		  Array<double> values = rsslist.get<Array<double> >("Values");
		  // convert to atm with datum at atmospheric pressure
		  for (Array<double>::iterator it=values.begin();
		       it!=values.end(); it++) {
		    *it = *it/1.01325e5 - 1.e0;
		  }
		  for (Array<std::string>::iterator it=regions.begin(); 
		   it!=regions.end(); it++) {
		    for (int j=0; j<ndim; j++) {
		      if (!(*it).compare(valid_region_def_lo[j])) {
			lo_def[j] = true;
			lo_bc[j] = 1;
			p_lo_bc[j] = 2;
			press_lo[j] = values[0];
		      }
		      else if (!(*it).compare(valid_region_def_hi[j])) {
			hi_def[j] = true;
			hi_bc[j] = 1;
			p_hi_bc[j] = 2;
			press_hi[j] = values[0];
		      }
		    }
		  }
		} 
		else if (rlabel=="BC: Linear Pressure") {
		  // set to saturated condition
		  rsublist.set("Water",1000e0);
		  rsublist.set("type","scalar");
		  comp_list.set(_label,rsublist);
		  bc_name_array.push_back(_label);
		  
		  Array<double> ref_values = rsslist.get<Array<double> >("Reference Values");
		  Array<double> ref_coor   = rsslist.get<Array<double> >("Reference Coordinate");
		  // convert to atm with datum at atmospheric pressure
		  for (Array<double>::iterator it=ref_values.begin();
		       it!=ref_values.end(); it++) {
		    *it = *it/101325 - 1.e0;
		  }
		  for (Array<std::string>::iterator it=regions.begin(); 
		       it!=regions.end(); it++) {
		    for (int j=0; j<ndim; j++) {
		      if (!(*it).compare(valid_region_def_lo[j])) {
			lo_def[j] = true;
			lo_bc[j] = 1;
			p_lo_bc[j] = 2;
			press_lo[j] = ref_values[0];
			if (j == ndim-1) water_table_lo = ref_coor[j];
		      }
		      else if (!(*it).compare(valid_region_def_hi[j])) {
			hi_def[j] = true;
			hi_bc[j] = 1;
			p_hi_bc[j] = 2;
			press_hi[j] = ref_values[0];
			if (j == ndim-1) water_table_hi = ref_coor[j];
		      }
		    }
		  }
		}
		else if (rlabel=="BC: Zero Flow") {
		  for (Array<std::string>::iterator it=regions.begin(); 
		       it!=regions.end(); it++) {
		    for (int j=0; j<ndim; j++) {
		      if (!(*it).compare(valid_region_def_lo[j])) {
			lo_def[j] = true;
			lo_bc[j] = 4;
			p_lo_bc[j] = 4;
		      }
		      else if (!(*it).compare(valid_region_def_hi[j])) {
			hi_def[j] = true;
			hi_bc[j] = 4;
			p_hi_bc[j] = 4;
		      }
		    }
		  }
		}	      
	      }          
	    }
	  }
	}
      //default to no flow
      /*
	bool all_bc_are_defined = true;
	for (int j=0; j<ndim; j++) {
	if (!lo_def[j])
	all_bc_are_defined = false;
	}
	if (!all_bc_are_defined)
	std::cerr << "Structured: boundarys conditions must be defined to all the "
	<< "following regions: XLOBC, XHIBC, YLOBC, YHIBC, ZLOBC and ZHIBC\n";
      */
      comp_list.set("inflow",bc_name_array);      
      comp_list.set("lo_bc",lo_bc);
      comp_list.set("hi_bc",hi_bc);
      press_list.set("lo_bc",p_lo_bc);
      press_list.set("hi_bc",p_hi_bc);
      press_list.set("inflow_bc_lo",inflow_lo_bc);
      press_list.set("inflow_bc_hi",inflow_hi_bc);
      press_list.set("inflow_vel_lo",inflow_lo_vel);
      press_list.set("inflow_vel_hi",inflow_hi_vel);
      press_list.set("press_lo",press_lo);
      press_list.set("press_hi",press_hi);
      press_list.set("water_table_lo",water_table_lo);
      press_list.set("water_table_hi",water_table_hi);
    }

    //
    // convert tracer to structured format
    //
    void
    convert_to_structured_tracer(const ParameterList& parameter_list, 
				 ParameterList&       struc_list)
    {
      const ParameterList& eclist = parameter_list.sublist("Execution control");
      std::string tmode = "none";
      if (eclist.isParameter("Transport Mode"))
	tmode = eclist.get<std::string>("Transport Mode");
      bool do_tracer = true;
      if (!tmode.compare("none"))
	do_tracer = false;

      if (do_tracer) {
	const ParameterList& rlist = parameter_list.sublist("Phase Definitions");
	
	// get tracers
	Array<std::string> array_phase, array_comp, array_tracer;
	for (ParameterList::ConstIterator i=rlist.begin(); i!=rlist.end(); ++i) {        
	  std::string label = rlist.name(i);
	  array_phase.push_back(label);
	  const ParameterList& rslist = rlist.sublist(label);
	  if (rslist.isSublist("Phase Components")) {
	    const ParameterList& rsslist = rslist.sublist("Phase Components");
	    for (ParameterList::ConstIterator j=rsslist.begin(); j!=rsslist.end(); ++j) {
	      std::string rlabel = rsslist.name(j);
	      array_comp.push_back(rlabel);
	      const ParameterList& rssslist = rsslist.sublist(rlabel);
	      if (rssslist.isParameter("Component Solutes")) {
		Array<std::string> solutes = rssslist.get<Array<std::string> >("Component Solutes");
		for (Array<std::string>::iterator it=solutes.begin(); 
		     it!=solutes.end();it++) {
		  array_tracer.push_back(underscore(*it));
		}
	      }
	    }
	  }
	}

	if (array_tracer.size() > 0) {
	  ParameterList& tracer_list = struc_list.sublist("tracer");
	  tracer_list.set("tracer",array_tracer);
	  tracer_list.set("group","Total");
	  for (Array<std::string>::iterator it=array_tracer.begin(); 
	       it!=array_tracer.end(); it++) {
	    ParameterList sublist;
	    sublist.set("group","Total");
	    tracer_list.set(*it,sublist);
	  }
	  
	  // get initial conditions
	  const ParameterList& ilist = parameter_list.sublist("Initial Conditions");
	  Array<std::string> init_set;
	  
	  for (ParameterList::ConstIterator i=ilist.begin(); i!=ilist.end(); ++i) {        
	    std::string label = ilist.name(i);
	    std::string _label = underscore(label);
	    init_set.push_back(_label);
	    
	    ParameterList rsublist;
	    const ParameterList& rslist = ilist.sublist(label);
	    Array<std::string> regions = rslist.get<Array<std::string> >("Assigned Regions");
	    Array<std::string> _regions;
	    for (Array<std::string>::iterator it=regions.begin(); it!=regions.end(); it++) {
	      _regions.push_back(underscore(*it));
	    }
	    rsublist.set("region",_regions);
	    const ParameterList& rsslist = rslist.sublist("Solute IC");
	    for (Array<std::string>::iterator it=array_phase.begin(); 
		 it!=array_phase.end(); it++ ) {
	      for (Array<std::string>::iterator jt=array_comp.begin(); 
		   jt!=array_comp.end(); jt++ ) {
		for (Array<std::string>::iterator kt=array_tracer.begin(); 
		     kt!=array_tracer.end(); kt++ ) {
		  const ParameterList& rssslist = rsslist.sublist(*it).sublist(*jt).sublist(*kt);
		  for (ParameterList::ConstIterator j=rssslist.begin(); 
		       j!=rssslist.end(); ++j) {
		    std::string rlabel = rssslist.name(j);
		    const ParameterEntry& rentry = rssslist.getEntry(rlabel);
		    if (rentry.isList()) {
		      const ParameterList& rsssslist = rssslist.sublist(rlabel);
		      if (rlabel=="IC: Uniform") {
			if (!rsublist.isParameter("type"))
			  rsublist.set("type","scalar");
			else
			  if (rsublist.get<std::string>("scalar") != "scalar")
			    std::cerr << "Definition of IC for tracer is not consistent!\n";
			rsublist.set(*kt,rsssslist.get<double>("Value"));
		      }
		    }
		  }
		}     
	      }
	    }
	    tracer_list.set(_label,rsublist);
	  }
	  tracer_list.set("init",init_set);
	  
	  // Boundary conditions
	  // Require information related to the regions.  
	  ParameterList regionlist = parameter_list.sublist("Regions");
	  
	  if (parameter_list.isSublist("Boundary Conditions")) {
	    const ParameterList& blist = parameter_list.sublist("Boundary Conditions");
	    Array<std::string> bc_set;
	    for (ParameterList::ConstIterator i=blist.begin(); i!=blist.end(); ++i) {        
	      std::string label  = blist.name(i);
	      std::string _label = underscore(label);
	      ParameterList rsublist;
	      const ParameterList& rslist = blist.sublist(label);
	      Array<std::string> regions = rslist.get<Array<std::string> >("Assigned Regions");
	      Array<std::string> _regions;
	      Array<std::string> valid_region_def_lo, valid_region_def_hi;
	      valid_region_def_lo.push_back("XLOBC");
	      valid_region_def_hi.push_back("XHIBC");
	      valid_region_def_lo.push_back("YLOBC");
	      valid_region_def_hi.push_back("YHIBC");
	      if (ndim == 3) {
		valid_region_def_lo.push_back("ZLOBC");
		valid_region_def_hi.push_back("ZHIBC");
	      }
	      for (Array<std::string>::iterator it=regions.begin(); it!=regions.end(); it++) {
		_regions.push_back(underscore(*it));
		bool valid = false;
		for (int j = 0; j<ndim; j++) {
		  if (!(*it).compare(valid_region_def_lo[j]))
		    valid = true;
		  else if (!(*it).compare(valid_region_def_hi[j]))
		    valid = true;
		}
		if (!valid)
		  std::cerr << "Structured: boundary conditions can only be applied to "
			    << "regions labeled as XLOBC, XHIBC, YLOBC, YHIBC, ZLOBC and ZHIBC\n";
	      }
	      const ParameterList& rsslist = rslist.sublist("Solute BC");
	      for (Array<std::string>::iterator it=array_phase.begin(); 
		   it!=array_phase.end(); it++ ) {
		for (Array<std::string>::iterator jt=array_comp.begin(); 
		     jt!=array_comp.end(); jt++ ) {
		  for (Array<std::string>::iterator kt=array_tracer.begin(); 
		       kt!=array_tracer.end(); kt++ ) {
		    const ParameterList& rssslist = rsslist.sublist(*it).sublist(*jt).sublist(*kt);
		    for (ParameterList::ConstIterator j=rssslist.begin(); 
			 j!=rssslist.end(); ++j) {
		      std::string rlabel = rssslist.name(j);
		      const ParameterEntry& rentry = rssslist.getEntry(rlabel);
		      if (rentry.isList()) {
			const ParameterList& rsssslist = rssslist.sublist(rlabel);
			if (rlabel=="BC: Inflow") {
			  if (!rsublist.isParameter("type"))
			    rsublist.set("type","scalar");
			  else
			    if (rsublist.get<std::string>("scalar") != "scalar")
			      std::cerr << "Definition of BC for tracer is not consistent!\n";
			  rsublist.set(*kt,rsssslist.get<Array<double> >("Values"));
			}
		      }
		    }
		  }     
		}
	      }
	      if (rsublist.isParameter("type"))
		if (rsublist.get<std::string>("type") == "scalar") {
		  bc_set.push_back(_label);
		  rsublist.set("region",_regions);
		  tracer_list.set(_label,rsublist);
		}
	    }
	    tracer_list.set("inflow",bc_set);
	  }
	}
      }
    }

    //
    // convert output to structured format
    //
    void
    convert_to_structured_output(const ParameterList& parameter_list, 
				 ParameterList&       struc_list)
    {

      ParameterList& amr_list = struc_list.sublist("amr");
      ParameterList& obs_list = struc_list.sublist("observation");

      const ParameterList& rlist = parameter_list.sublist("Output");
      
      // cycle macros
      const ParameterList& clist = rlist.sublist("Cycle Macros");
      std::map<std::string,int> cycle_map;
      for (ParameterList::ConstIterator i=clist.begin(); i!=clist.end(); ++i) {
	std::string label = clist.name(i);
	const ParameterList& rslist = clist.sublist(label);
	Array<int> tmp = rslist.get<Array<int> >("Start_Period_Stop");
	cycle_map[label] = tmp[1];
      }

      // vis data
      const ParameterList& vlist = rlist.sublist("Visualization Data");
      amr_list.set("plot_file",vlist.get<std::string>("File Name Base"));
      amr_list.set("plot_int",cycle_map[vlist.get<std::string>("Cycle Macro")]);
      if (vlist.isParameter("Visualised Derives")) {
          Array<std::string> visNames = vlist.get<Array<std::string> >("Visualised Derives");
          for (int i=0; i<visNames.size(); ++i) {
              visNames[i] = underscore(visNames[i]);
          }
          amr_list.set<Array<std::string> >("derive_plot_vars",visNames);
      }
      amr_list.set<std::string>("plot_vars",""); // Shut off, per spec

      // check point
      const ParameterList& plist = rlist.sublist("Checkpoint Data");
      amr_list.set("check_file",plist.get<std::string>("File Name Base"));
      amr_list.set("check_int",cycle_map[plist.get<std::string>("Cycle Macro")]);

      // time macros
      const ParameterList& tlist = rlist.sublist("Time Macros");
      std::map<std::string,Array<double> > time_map;
      for (ParameterList::ConstIterator i=tlist.begin(); i!=tlist.end(); ++i) {
	std::string label = tlist.name(i);
	const ParameterList& rslist = tlist.sublist(label);
	if (rslist.isParameter("Values")) {
	    Array<double> tmp = rslist.get<Array<double> >("Values");
	    time_map[label] = tmp;
	}
	else if (rslist.isParameter("Start_Period_Stop")) {
	  Array<double> tmp = rslist.get<Array<double> >("Start_Period_Stop");
	  Array<double> timeseries;
	  timeseries.push_back(0.e0);
	  double timecount = 0.0;
	  while (timecount < simulation_time) {
	    timecount += tmp[1];
	    timeseries.push_back(timecount);
	  }
	  time_map[label] = timeseries;

	}
      }      

      // variable macros
      const ParameterList& varlist = rlist.sublist("Variable Macros");
      std::map<std::string,tridata> var_map;
      for (ParameterList::ConstIterator i=varlist.begin(); i!=varlist.end(); ++i) {
	std::string label = varlist.name(i);
	const ParameterList& rslist = varlist.sublist(label);
	tridata tri;
	tri.phase = rslist.get<std::string>("Phase");
	if (rslist.isParameter("Component")) {
	  tri.component = rslist.get<std::string>("Component");
	  if (rslist.isParameter("Solute"))
	    tri.solute = rslist.get<std::string>("Solute");
	}
	var_map[label] = tri;
      }
      
      // observation
      Array<std::string> arrayobs;
      const ParameterList& olist = rlist.sublist("Observation Data");
      ParameterList sublist;
      for (ParameterList::ConstIterator i=olist.begin(); i!=olist.end(); ++i) {
	std::string label = olist.name(i);
	std::string _label = underscore(label);
	const ParameterEntry& entry = olist.getEntry(label);
	if (entry.isList()) {
	  const ParameterList& rslist = olist.sublist(label);
	  std::string functional = rslist.get<std::string>("Functional");
	  std::string region_name = rslist.get<std::string>("Region");
	  if (functional == "Observation Data: Integral")
	    sublist.set("obs_type","integral");
	  else if (functional == "Observation Data: Point")
	    {
	      sublist.set("obs_type","point_sample");
	      const ParameterList& lregion = 
		parameter_list.sublist("Regions").sublist(region_name);
	      if (!lregion.isSublist("Region: Point"))
		{
		  std::cerr << label << " is a point observation and "
			    << region_name << " is not a point region.\n";
		  throw std::exception();
		}
	    }
	  sublist.set("region",region_name);
	  sublist.set("times",time_map[rslist.get<std::string>("Time Macro")]);
	  
	  Array<std::string > variables = rslist.get<Array<std::string> >("Variable Macro");
	  tridata tri = var_map[variables[0]];
	  if (tri.solute.empty()) {
	    sublist.set("var_type","comp");
	    sublist.set("var_id",tri.component);
	  }
	  else if (tri.component.empty()) {
	    sublist.set("var_type","phase");
	    sublist.set("var_id",tri.phase);
	  }
	  else {
	    sublist.set("var_type","tracer");
	    sublist.set("var_id",tri.solute);
	  }	  
	  obs_list.set(_label,sublist);
	  arrayobs.push_back(_label);
	}
      }
      obs_list.set("observation",arrayobs);
    }
  }
}


