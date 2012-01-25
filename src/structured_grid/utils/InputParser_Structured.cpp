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

    double atmToMKS = 101325;

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
	geom_list.set("regions",arrayregions);
      }
    }

      typedef std::map<std::string,bool> MTEST;
      std::vector<std::string> remaining_false(const MTEST& p) 
      {
          std::vector<std::string> ret;
          for (MTEST::const_iterator it=p.begin(); it!=p.end(); ++it) {
              if (!it->second) {
                  ret.push_back(it->first);
              }
          }
          return ret;
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
        
        for (ParameterList::ConstIterator i=rlist.begin(); i!=rlist.end(); ++i)
        {
            MTEST mtest;
            mtest.insert(MTEST::value_type("Porosity",false));
            mtest.insert(MTEST::value_type("Density",false));
            mtest.insert(MTEST::value_type("Intrinsic_Permeability",false));
            mtest.insert(MTEST::value_type("Capillary_Pressure",false));
            mtest.insert(MTEST::value_type("Relative_Permeability",false));
            mtest.insert(MTEST::value_type("Regions_Assigned",false));
        
            std::string label = rlist.name(i);
            std::string _label = underscore(label);
                    
            // Add this rock label to list of rocks
            arrayrock.push_back(_label);

            const ParameterEntry& entry = rlist.getEntry(label);
            
            if (entry.isList()) {
                ParameterList rsublist;
                const ParameterList& rslist = rlist.sublist(label);
                for (ParameterList::ConstIterator j=rslist.begin(); j!=rslist.end(); ++j) 
                {
                    const std::string& rlabel = rslist.name(j);
                    const ParameterEntry& rentry = rslist.getEntry(rlabel);
                    
                    if (rentry.isList())
                    {
                        const ParameterList& rsslist = rslist.sublist(rlabel);
                        if (rlabel=="Porosity: Uniform"){
                            rsublist.setEntry("porosity",rsslist.getEntry("Value"));
                            rsublist.set("porosity_dist","uniform");
                            mtest["Porosity"] = true;
                        }
                        else if (rlabel=="Intrinsic Permeability: Anisotropic Uniform") {
                            Array<double> array_p(2);
                            array_p[0] = rsslist.get<double>("Horizontal");
                            array_p[1] = rsslist.get<double>("Vertical");
                            // convert from m^2 to mDa
                            for (int k=0; k<2; k++) {
                                array_p[k] /= 9.869233e-16;
                            }
                            rsublist.set("permeability",array_p);
                            rsublist.set("permeability_dist","uniform");
                            mtest["Intrinsic_Permeability"] = true;
                        }
                        else if (rlabel=="Capillary Pressure: van Genuchten") {
                            int cpl_type = 3;
                            rsublist.set("cpl_type",cpl_type);
                            
                            double alpha = rsslist.get<double>("alpha");
                            
                            Array<double> array_c(4);
                            array_c[1] = alpha*1.01325e5; // convert Pa^-1 to atm^-1 
                            array_c[0] = rsslist.get<double>("m");
                            array_c[2] = rsslist.get<double>("Sr");
                            array_c[3] = 0.0;                  
                            rsublist.set("cpl_param",array_c);
                            mtest["Capillary_Pressure"] = true;

                            std::string krType = rsslist.get<std::string>("Relative Permeability");
                            if (krType=="Mualem") {
                                Array<double> array_k(3);
                                array_k[0] = array_c[0];
                                array_k[1] = array_c[2];
                                array_k[2] = array_c[3];
                                rsublist.set("kr_type",cpl_type);
                                rsublist.set("kr_param",array_k);
                                mtest["Relative_Permeability"] = true;
                            }
                            else {
                                std::cerr << "Unsupported Relative Permeability model: " << krType << std::endl;
                                throw std::exception();
                            }
                        } 
                    }
                    else if (rlabel=="Assigned Regions") {
                        Array<std::string> tmp_region = rslist.get<Array<std::string> >(rlabel);
                        for (Array<std::string>::iterator it=tmp_region.begin();
                             it !=tmp_region.end(); it++)
                            (*it) = underscore(*it);
                        rsublist.set("regions",tmp_region);
                        mtest["Regions_Assigned"] = true;
                    }
                    else if (rlabel=="Density") {
                        rsublist.set<double>("density",rslist.get<double>("Density"));
                        //rsublist.set("density",1e3);
                        mtest["Density"] = true;
                    }
                    else {
                        std::cerr << "Unrecognized rock parameter: " << rlabel << std::endl;
                        throw std::exception();
                    }
                }
                std::vector<std::string> region_check = remaining_false(mtest); 
                if (region_check.size()==0) {
                    rock_list.set(_label,rsublist);
                }
                else {
                    std::cerr << "Material not completely defined: " << label << std::endl;
                    std::cerr << "   unfilled: ";
                    for (int i=0; i<region_check.size(); ++i)
                        std::cerr << region_check[i] << " ";
                    std::cerr << '\n';
                    throw std::exception();
                }
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
      Array<std::string> arrayphase, phaseNames;
      std::map<std::string,double> viscosity,density;
      
      // get Phase Properties
      Array<std::string> array_comp;
      for (ParameterList::ConstIterator i=rlist.begin(); i!=rlist.end(); ++i) {
          
          const std::string& phaseLabel = rlist.name(i);
          std::string _phaseLabel = underscore(phaseLabel);
          const ParameterEntry& phaseEntry = rlist.getEntry(phaseLabel);
          
          if ( !phaseEntry.isList() ) {
              if (Teuchos::MPISession::getRank() == 0) {
                  std::cerr << "Phase section must define only phases  " 
                            << phaseLabel << " is not a valid phase definition." << std::endl;
              }
              throw std::exception();
          }
          arrayphase.push_back(_phaseLabel);	
          
          // FIXME: Set arbitrary ordering of phases for IC/BC value semantics
          phaseNames.push_back(phaseLabel); // remember Amanzi phase name for later

          ParameterList phasePLtr;
          double visc = -1;
          double rho = -1;
          double diff = 0;
          Array<std::string> compArray;
          const ParameterList& phasePL = rlist.sublist(phaseLabel);

          for (ParameterList::ConstIterator j=phasePL.begin(); j!=phasePL.end(); ++j) {
              const std::string& phasePLlabel = phasePL.name(j);
              const ParameterEntry& phasePLentry = phasePL.getEntry(phasePLlabel);

              if (phasePLentry.isList()) {

                  if ( phasePLlabel == "Phase Properties") {
                      const ParameterList phasePropPL = phasePL.sublist(phasePLlabel);
                      for (ParameterList::ConstIterator k=phasePropPL.begin(); k!=phasePropPL.end(); ++k) {
                          const std::string& pPropLabel = phasePropPL.name(k);
                          const ParameterEntry& pPropEntry = phasePropPL.getEntry(pPropLabel);
                          
                          if (pPropEntry.isList()) {
                              const ParameterList& pPropPLi = phasePropPL.sublist(pPropLabel);
                              if (pPropLabel == "Viscosity: Uniform") {
                                  visc = pPropPLi.get<double>("Viscosity");
                              }
                              else if (pPropLabel == "Density: Uniform") {
                                  rho = pPropPLi.get<double>("Density");
                              }
                              else if (pPropLabel == "Diffusivity: Uniform") {
                                  diff = pPropPLi.get<double>("Diffusivity");
                              }
                              else {
                                  std::cerr << "Invalid Phase Property entry: " << pPropLabel << std::endl;
                                  throw std::exception();
                              }
                          }
                      }
                  }
                  else if ( phasePLlabel=="Phase Components") {
                      const ParameterList compPL = phasePL.sublist(phasePLlabel);
                      for (ParameterList::ConstIterator L=compPL.begin(); L!=compPL.end(); ++L) {
                          const std::string& compLabel = compPL.name(L);
                          const ParameterEntry& compEntry = compPL.getEntry(compLabel);
                          if (compEntry.isList()) {
                              std::string _compLabel = underscore(compLabel);
                              compArray.push_back(_compLabel);
                              
                              array_comp.push_back(_compLabel); // Flattened array here
                          }
                          else {
                              std::cerr << "Phase components must be specified as lists";
                              throw std::exception();
                          }
                      }
                  }
                  else {
                      std::cerr << "Invalid entry in phase definition: " << phaseLabel << std::endl;
                      throw std::exception();
                  }
              }
              else {
                  std::cerr << "Invalid entry in phase definition: " << phaseLabel << std::endl;
                  throw std::exception();
              }
              
              if (compArray.size()==0) {
                  std::cerr << "No components found in phase " << phaseLabel << std::endl;
                  std::cerr << "  rho, visc: " << rho << ", " << visc << std::endl;
                  throw std::exception();
              }
          }
          phasePLtr.set("comps",compArray);
          
          if (rho<0) {
              std::cerr << "density not defined for phase=" << phaseLabel << std::endl;
              throw std::exception();
          }
          if (visc<0) {
              std::cerr << "viscosity not defined for phase=" << phaseLabel << std::endl;
              throw std::exception();
          }
          if (diff<0) {
              std::cerr << "diffusivity not defined for phase=" << phaseLabel << std::endl;
              throw std::exception();
          }
          visc *= rho;

          phasePLtr.set("density",rho);
          phasePLtr.set("viscosity",visc);
          phasePLtr.set("diffusivity",diff);
          phase_list.set(_phaseLabel,phasePLtr);
      }

      phase_list.set("phase",arrayphase);

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
          std::string _rlabel = underscore(rlabel); 
	  Array<std::string > regions = rslist.get<Array<std::string> >("Assigned Regions");
	  Array<std::string> _regions;
	  for (Array<std::string>::iterator it=regions.begin(); it!=regions.end(); it++) {
	    _regions.push_back(underscore(*it));
	  }
	  rsublist.set("regions",_regions);

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
          }  
          ParameterList thisIClist;
          thisIClist.set(_label,rsublist);
          comp_list.set("ics",thisIClist);
	}
      }
      comp_list.set("inits",init_region);

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
      //Array<std::string> inflow_hi_vel(ndim,0), inflow_lo_vel(ndim,0); // time function names
      //Array<std::string> press_hi(ndim,0), press_lo(ndim,0); // time function names
      double water_table_hi=0;
      double water_table_lo=0;

      if (parameter_list.isSublist("Boundary Conditions"))
      {
          ParameterList thisBClist;
	  const ParameterList& blist = parameter_list.sublist("Boundary Conditions");
	  for (ParameterList::ConstIterator i=blist.begin(); i!=blist.end(); ++i) {
        
              std::string label = blist.name(i);
              std::string _label = underscore(label);
              const ParameterEntry& entry = blist.getEntry(label);
              
              const ParameterList& rslist = blist.sublist(label);
              for (ParameterList::ConstIterator j=rslist.begin(); j!=rslist.end(); ++j) {
                  
                  const std::string& rlabel = rslist.name(j);
                  std::string _rlabel = underscore(rlabel);
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
	      
                  ParameterList rsublist;
                  rsublist.set("regions",_regions);

                  const ParameterEntry& rentry = rslist.getEntry(rlabel);
                  if (rentry.isList())
                  {
                      const ParameterList& rsslist = rslist.sublist(rlabel);
                      if (rlabel=="BC: Flux")
                      {
                          bool is_in_vol = rsslist.isParameter("Inward Volumetric Flux");
                          bool is_in_mass = rsslist.isParameter("Inward Mass Flux");
                          bool is_out_vol = rsslist.isParameter("Outward Volumetric Flux");
                          bool is_out_mass = rsslist.isParameter("Outward Mass Flux");
                          
                          bool is_mass = is_in_mass || is_out_mass;
                          bool is_out = is_out_vol || is_out_mass;
                          
                          Array<double> fluxvals, timevals;
                          Array<std::string> timeforms;
                          
                          if ( !(is_in_mass || is_out_mass || is_in_vol || is_out_vol) ) {
                              std::cout << "BC: Unsupported flux quantity" << std::endl;
                              throw std::exception();
                          }
                          
                          if (is_in_vol) {
                              fluxvals = rsslist.get<Array<double> >("Inward Volumetric Flux");
                          }
                          else if (is_out_vol) {
                              fluxvals = rsslist.get<Array<double> >("Outward Volumetric Flux");
                          }
                          else if (is_in_mass) {
                              fluxvals = rsslist.get<Array<double> >("Inward Mass Flux");
                          }
                          else if (is_out_mass) {
                              fluxvals = rsslist.get<Array<double> >("Outward Mass Flux");
                          }
                          else {
                              std::cerr << "Flux not specified in recognized form" << std::endl;
                              throw std::exception();
                          }
                          
                          if (is_mass) {
                              // FIXME: Assumes this is Aqueous
                              if (density.find("Aqueous") == density.end()) {
                                  std::cerr << "Aqueous phase not found, not scaling mass flux" << std::endl;
                                  throw std::exception();
                              }
                              else {
                                  for (int i=0; i<fluxvals.size(); ++i) {
                                      fluxvals[i] /= density["Aqueous"];
                                  }
                              }
                          }
                          
                          if (is_out) {
                              for (int i=0; i<fluxvals.size(); ++i) {
                                  fluxvals[i] = -fluxvals[i];
                              }
                          }
                          
                          if (!rsslist.isParameter("Material Type at Boundary"))
                          {
                              // this is temporary.  It will be removed once the structured code is updated.
                              const ParameterList& mat_list = parameter_list.sublist("Material Properties"); 
                              std::string first_material = mat_list.name(mat_list.begin());
                              rsublist.set("rock",first_material);
                          }
                          
                          
                          timevals = rsslist.get<Array<double> >("Times");
                          if (timevals.size() != fluxvals.size()) {
                              std::cout << "tv,fv: " << timevals.size() << ", " << fluxvals.size() << std::endl;
                              std::cerr << "Wrong number of Time Values specified for bc: \""
                                        << rlabel << "\"" << std::endl;
                              throw std::exception();                          
                          }
                          
                          if (timevals.size() > 0) {
                              timeforms = rsslist.get<Array<std::string> >("Time Functions");
                              if (timeforms.size() != fluxvals.size() - 1) {
                                  std::cerr << "Wrong number of Time Functions specified for bc: \"" 
                                            << rlabel << "\"" << std::endl;
                                  throw std::exception();                          
                              }
                          }
                          
                          rsublist.set<Array<double> >("aqueous_volumetric_flux",fluxvals);
                          rsublist.set<Array<double> >("inflowtimes",timevals);
                          rsublist.set<Array<std::string> >("inflowfncs",timeforms);
                          rsublist.set("type","zero_total_velocity");
                          bc_name_array.push_back(_label);
                          
                          for (Array<std::string>::iterator it=regions.begin(); 
                               it!=regions.end(); it++) {
                              for (int j=0; j<ndim; j++) {
                                  if (!(*it).compare(valid_region_def_lo[j])) {
                                      lo_def[j] = true;
                                      lo_bc[j] = 1;
                                      p_lo_bc[j] = 1;
                                      inflow_lo_bc[j] = 1;
                                      inflow_lo_vel[j] = fluxvals[0];
                                  }
                                  else if (!(*it).compare(valid_region_def_hi[j])) {
                                      hi_def[j] = true;
                                      hi_bc[j] = 1;
                                      p_hi_bc[j] = 1;
                                      inflow_hi_bc[j] = 1;
                                      inflow_hi_vel[j] = -fluxvals[0];
                                  }
                              }
                          }
                          thisBClist.set(_label,rsublist);
                      }                      
                      else if (rlabel=="BC: Uniform Pressure")
                      {
                          // set to saturated condition
                          rsublist.set("Water",density[arrayphase[0]]);
                          rsublist.set("type","scalar");
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
                          thisBClist.set(_label,rsublist);
                      } 
                      else if (rlabel=="BC: Linear Pressure") {
                          // set to saturated condition
                          rsublist.set("Water",density["Aqueous"]);
                          rsublist.set("type","scalar");
                          bc_name_array.push_back(_label);
                          
                          Array<double> ref_values = rsslist.get<Array<double> >("Reference Values");
                          Array<double> ref_coor   = rsslist.get<Array<double> >("Reference Coordinate");
                          // convert to atm with datum at atmospheric pressure
                          for (Array<double>::iterator it=ref_values.begin();
                               it!=ref_values.end(); it++) {
                              *it = *it/atmToMKS - 1.e0;
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
                          thisBClist.set(_label,rsublist);
                      }
                      else if (rlabel=="BC: No Flow") {
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
                          thisBClist.set(_label,rsublist);
                      }
                      else if (rlabel=="Solute BC") {
                          // Parse elsewhere
                      }
                      else {
                          std::cerr << "Unrecognized BC type: " << rlabel << std::endl;
                          throw std::exception();
                      }
                  }
              }
	  }
          comp_list.set("bcs",thisBClist);
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

      struct TRACER
      {
          TRACER() {};
          TRACER(const std::string& phase,
                 const std::string& comp,
                 const std::string& solute,
                 const Array<std::string>& regions,
                 const std::string& units,
                 double             val,
                 const std::string& AMR_type,
                 const std::string& Amanzi_type,
                 const std::string& label)
              : phase(phase), comp(comp), solute(solute), regions(regions), units(units), vals(1,val),
                AMR_type(AMR_type), Amanzi_type(Amanzi_type), label(label) {}
          
          TRACER(const std::string&   phase,
                 const std::string&   comp,
                 const std::string&   solute,
                 const Array<std::string>& regions,
                 const std::string&   units,
                 const Array<double>& vals,
                 const Array<double>& times,
                 const Array<std::string>& forms,
                 const std::string& AMR_type,
                 const std::string& Amanzi_type,
                 const std::string& label)
              : phase(phase), comp(comp), solute(solute), regions(regions), units(units), vals(vals),
                times(times), forms(forms), AMR_type(AMR_type), Amanzi_type(Amanzi_type), label(label) {}
          
          const std::string& Phase() const {return phase;}
          const std::string& Comp() const {return comp;}
          const std::string& Solute() const {return solute;}
          const std::string& AMR_Type() const {return AMR_type;}
          const std::string& Amanzi_Type() const {return Amanzi_type;}
          const std::string& Label() const {return label;}
          const Array<std::string>& Regions() const {return regions;}
          const Array<double>& Values() const {return vals;}
          const Array<double>& Times() const {return times;}
          const Array<std::string>& Forms() const {return forms;}
      protected:
          std::string AMR_type, Amanzi_type;
          std::string label;
          std::string phase;
          std::string comp;
          std::string solute;
          Array<std::string> regions;
          std::string units;
          Array<double> vals;
          Array<double> times;
          Array<std::string> forms;
      };
      
      
      
      struct COMP
      {
          std::map<std::string, TRACER> tracers; // tracerLabel->TRACER
      };
      
      struct PHASE
      {
          std::map<std::string, COMP> comps; // compLabel->COMP
          double density;
          double viscosity;
      };
      
      typedef std::map<std::string,TRACER> TracerMap;

      bool build_timeindep_function_data(double& val,
                                         const ParameterList& pl,
                                         const std::string& val_name)
      {
          bool retval = true;
          bool found_val = false;
          for (ParameterList::ConstIterator i=pl.begin(); i!=pl.end(); ++i) {
              const std::string& name = pl.name(i);
              if (name == val_name) {
                  found_val = true;
                  if (!pl.isSublist(name)) {
                      val = pl.get<double>(val_name);
                  }
                  else {
                      retval = false;
                      std::cerr << val_name << " cannot be in a ParameterList" << std::endl;
                      throw std::exception();
                  }
              } else {
                  retval = false;
                  std::cerr << "Unrecognized function parameter: " << name << std::endl;
                  throw std::exception();
              }
          }
          if (!found_val) {
              retval = false;
              std::cerr << "Function requires " << val_name << std::endl;
              throw std::exception();
          }
          return retval;
      }

      bool build_timedep_function_data(Array<double>& vals,
                                       Array<double>& times,
                                       Array<std::string>& forms,
                                       const ParameterList& pl,
                                       const std::string& vals_name)
      {
          bool retval = true;
          const std::string times_name="Times";
          const std::string forms_name="Time Functions";

          bool found_vals = false; vals.resize(0);
          bool found_times = false; times.resize(0);
          bool found_forms = false; forms.resize(0);
          for (ParameterList::ConstIterator i=pl.begin(); i!=pl.end(); ++i) {
              const std::string& name = pl.name(i);
              if (name == vals_name) {
                  found_vals = true;
                  if (!pl.isSublist(name)) {
                      vals = pl.get<Array<double> >(vals_name);
                  }
                  else {
                      retval = false;
                      std::cerr << "Time function parameter values (in " << vals_name << " cannot be in a ParameterList" << std::endl;
                      throw std::exception();
                  }
              } else if (name == times_name) {
                  found_times = true;
                  if (!pl.isSublist(times_name)) {
                      times = pl.get<Array<double> >(times_name);
                  }
                  else {
                      retval = false;
                      std::cerr << "Time function parameter values (in " << times_name << " cannot be in a ParameterList" << std::endl;
                      throw std::exception();
                  }
              } else if (name == forms_name) {
                  found_forms = true;
                  if (!pl.isSublist(times_name)) {
                      forms = pl.get<Array<std::string> >(forms_name);
                  }
                  else {
                      retval = false;
                      std::cerr << "Time function parameter values (in " << forms_name << " cannot be in a ParameterList" << std::endl;
                      throw std::exception();
                  }
              } else {
                  retval = false;
                  std::cerr << "Unrecognized Time function parameter: " << name << std::endl;
                  throw std::exception();
              }
          }

          if (!found_vals) {
              retval = false;
              std::cerr << "Time function requires " << vals_name << std::endl;
              throw std::exception();
          }
          if (!found_times) {
              retval = false;
              std::cerr << "Time function requires " << times_name << std::endl;
              throw std::exception();
          }
          if (!found_forms) {
              retval = false;
              std::cerr << "Time function requires " << forms_name << std::endl;
              throw std::exception();
          }
          return true;
      }

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
              
              TracerMap ic_tmap;
              const ParameterList& ilist = parameter_list.sublist("Initial Conditions");
              
              Array<std::string> ic_assigned_regions;
              std::string ic_phaseLabel,ic_compLabel,ic_soluteLabel, ic_units;
              std::string ic_AMR_type,ic_Amanzi_type,ic_label;
              double ic_value, ic_coord, ic_grad;
              
              for (ParameterList::ConstIterator i=ilist.begin(); i!=ilist.end(); ++i) {        
                  const std::string& IClabel = ilist.name(i);
                  if (ilist.isSublist(IClabel)) {
                      const ParameterList& ICList = ilist.sublist(IClabel);
                      for (ParameterList::ConstIterator j=ICList.begin(); j!=ICList.end(); ++j) {
                          const std::string& ICsolutePLLabel = ICList.name(j);
                          if (ICsolutePLLabel=="Assigned Regions") {
                              ic_assigned_regions = ICList.get<Array<std::string> >(ICsolutePLLabel);
                          }
                          else if (ICsolutePLLabel=="Solute IC" && ICList.isSublist(ICsolutePLLabel)) {
                              const ParameterList& ICsolutePL = ICList.sublist(ICsolutePLLabel);
                              for (ParameterList::ConstIterator k=ICsolutePL.begin(); k!=ICsolutePL.end(); ++k) {
                                  ic_phaseLabel = ICsolutePL.name(k);
                                  if (ICsolutePL.isSublist(ic_phaseLabel)) {
                                      const ParameterList& phasePL = ICsolutePL.sublist(ic_phaseLabel);
                                      for (ParameterList::ConstIterator L=phasePL.begin(); L!=phasePL.end(); ++L) {
                                          ic_compLabel = phasePL.name(L);
                                          if (phasePL.isSublist(ic_compLabel)) {
                                              const ParameterList& compPL = phasePL.sublist(ic_compLabel);
                                              for (ParameterList::ConstIterator M=compPL.begin(); M!=compPL.end(); ++M) {
                                                  ic_soluteLabel = compPL.name(M);
                                                  if (compPL.isSublist(ic_soluteLabel)) {
                                                      const ParameterList& soluteICPL = compPL.sublist(ic_soluteLabel);
                                                      for (ParameterList::ConstIterator N=soluteICPL.begin(); N!=soluteICPL.end(); ++N) {
                                                          ic_Amanzi_type = soluteICPL.name(N);

                                                          if ( (ic_Amanzi_type == "IC: Uniform Concentration") && soluteICPL.isSublist(ic_Amanzi_type) )
                                                          {
                                                              const ParameterList& funcPL = soluteICPL.sublist(ic_Amanzi_type); 
                                                              build_timeindep_function_data(ic_value,funcPL,"Value");
                                                          }
                                                          else if (ic_Amanzi_type == "Concentration Units") {
                                                              if (soluteICPL.isParameter(ic_Amanzi_type)) {
                                                                  ic_units = soluteICPL.get<std::string>("Concentration Units");
                                                              }
                                                              else {
                                                                  std::cerr << "Unsupported parameters to " << ic_Amanzi_type << std::endl;
                                                                  throw std::exception();
                                                              }
                                                          }
                                                          else if (ic_Amanzi_type == "IC: File Pressure"  || ic_Amanzi_type == "IC: File Saturation")
                                                          {
                                                              std::cerr << "Solute IC type " << ic_Amanzi_type << " not yet supported" << std::endl;
                                                              throw std::exception();
                                                          }
                                                          else {
                                                              std::cerr << "Unsupported solute IC: " << ic_Amanzi_type << std::endl;
                                                              std::cerr << " is list: " << soluteICPL.isSublist(ic_Amanzi_type) << std::endl;
                                                              throw std::exception();
                                                          }
                                                      }
                                                      
                                                  }
                                                  else {
                                                      std::cerr << "Solute IC must be a list" << std::endl;
                                                      throw std::exception();
                                                  }
                                              }
                                          }
                                          else {
                                              std::cerr << "Solute IC Comp must be a list" << std::endl;
                                              throw std::exception();
                                          }
                                      }
                                  }
                                  else {
                                      std::cerr << "Solute IC Phase must be a list" << std::endl;
                                      throw std::exception();
                                  }
                                  
                              }
                              
                          }
                      }

                      ic_AMR_type = "scalar";

                      ic_tmap[ic_Amanzi_type] = TRACER(ic_phaseLabel,ic_compLabel,ic_soluteLabel,ic_assigned_regions,ic_units,ic_value,
                                                       ic_AMR_type,underscore(ic_Amanzi_type),IClabel);

                  }
                  else {
                      std::cerr << "Solute IC must be a list" << std::endl;
                      throw std::exception();
                  }

              }
          
              TracerMap bc_tmap;
              const ParameterList& blist = parameter_list.sublist("Boundary Conditions");
              
              Array<std::string> bc_assigned_regions;
              std::string bc_phaseLabel,bc_compLabel,bc_soluteLabel, bc_units;
              Array<double> bc_values, bc_times, bc_ref_values, bc_coords, bc_gradient;
              Array<std::string> bc_funcs, bc_time_funcs;
              std::string bc_AMR_type, bc_Amanzi_type;

              for (ParameterList::ConstIterator i=blist.begin(); i!=blist.end(); ++i) {        
                  const std::string& BClabel = blist.name(i);
                  if (blist.isSublist(BClabel)) {
                      const ParameterList& BCList = blist.sublist(BClabel);
                      for (ParameterList::ConstIterator j=BCList.begin(); j!=BCList.end(); ++j) {
                          const std::string& BCsolutePLLabel = BCList.name(j);
                          if (BCsolutePLLabel=="Assigned Regions") {
                              bc_assigned_regions = BCList.get<Array<std::string> >(BCsolutePLLabel);
                          }
                          else if (BCsolutePLLabel=="Solute BC" && BCList.isSublist(BCsolutePLLabel)) {
                              const ParameterList& BCsolutePL = BCList.sublist(BCsolutePLLabel);
                              for (ParameterList::ConstIterator k=BCsolutePL.begin(); k!=BCsolutePL.end(); ++k) {
                                  bc_phaseLabel = BCsolutePL.name(k);
                                  if (BCsolutePL.isSublist(bc_phaseLabel)) {
                                      const ParameterList& phasePL = BCsolutePL.sublist(bc_phaseLabel);
                                      for (ParameterList::ConstIterator L=phasePL.begin(); L!=phasePL.end(); ++L) {
                                          bc_compLabel = phasePL.name(L);
                                          if (phasePL.isSublist(bc_compLabel)) {
                                              const ParameterList& compPL = phasePL.sublist(bc_compLabel);
                                              for (ParameterList::ConstIterator M=compPL.begin(); M!=compPL.end(); ++M) {
                                                  bc_soluteLabel = compPL.name(M);
                                                  if (compPL.isSublist(bc_soluteLabel)) {
                                                      const ParameterList& soluteBCPL = compPL.sublist(bc_soluteLabel);
                                                      for (ParameterList::ConstIterator N=soluteBCPL.begin(); N!=soluteBCPL.end(); ++N) {
                                                          bc_Amanzi_type = soluteBCPL.name(N);
                                                          if ( (bc_Amanzi_type == "BC: Uniform Concentration" || bc_Amanzi_type == "BC: Inflow")
                                                               && soluteBCPL.isSublist(bc_Amanzi_type) )
                                                          {
                                                              const ParameterList& funcPL = soluteBCPL.sublist(bc_Amanzi_type);
                                                              build_timedep_function_data(bc_values,bc_times,bc_time_funcs,funcPL,"Values");
                                                              bc_AMR_type = "scalar";
                                                          }
                                                          else if ( (bc_Amanzi_type == "BC: Outflow"   || bc_Amanzi_type == "BC: Zero Gradient")
                                                                    && soluteBCPL.isSublist(bc_Amanzi_type) ) { 
                                                              bc_values.resize(0);
                                                              bc_AMR_type = "outflow";
                                                          }
                                                          else if ( (bc_Amanzi_type == "BC: Zero Flow") && soluteBCPL.isSublist(bc_Amanzi_type) ) { 
                                                              bc_values.resize(0);
                                                              bc_AMR_type = "noflow";
                                                          }
                                                          else if (bc_Amanzi_type == "Concentration Units") {
                                                              if (soluteBCPL.isParameter(bc_Amanzi_type)) {
                                                                  bc_units = soluteBCPL.get<std::string>("Concentration Units");
                                                              }
                                                              else {
                                                                  std::cerr << "Unsupported parameters to " << bc_Amanzi_type << std::endl;
                                                                  throw std::exception();
                                                              }
                                                          }
                                                          else {
                                                              std::cerr << "Unsupported solute BC: " << bc_Amanzi_type << std::endl;
                                                              throw std::exception();
                                                          }
                                                      }
                                                      
                                                  }
                                                  else {
                                                      std::cerr << "Solute BC must be a list" << std::endl;
                                                      throw std::exception();
                                                  }
                                              }
                                          }
                                          else {
                                              std::cerr << "Solute BC Comp must be a list" << std::endl;
                                              throw std::exception();
                                          }
                                      }
                                  }
                                  else {
                                      std::cerr << "Solute BC Phase must be a list" << std::endl;
                                      throw std::exception();
                                  }
                                  
                              }
                              
                          }
                      }
                      bc_tmap[bc_Amanzi_type] = TRACER(bc_phaseLabel,bc_compLabel,bc_soluteLabel,bc_assigned_regions,bc_units,
                                                       bc_values,bc_times,bc_funcs,
                                                       bc_AMR_type,underscore(bc_Amanzi_type), BClabel);

                  }
                  else {
                      std::cerr << "Solute BC must be a list" << std::endl;
                      throw std::exception();
                  }

              }
          

              // Check that solutes consistent with phase definitions
              const ParameterList& rlist = parameter_list.sublist("Phase Definitions");

              std::map<std::string,PHASE> phases;
              for (ParameterList::ConstIterator i=rlist.begin(); i!=rlist.end(); ++i) {        
                  const std::string& phaseLabel = rlist.name(i);
                  if (rlist.isSublist(phaseLabel)) {
                      const ParameterList& phasePL = rlist.sublist(phaseLabel);
                      for (ParameterList::ConstIterator j=phasePL.begin(); j!=phasePL.end(); ++j) {
                          const std::string& compPLLabel = phasePL.name(j);
                          if (compPLLabel=="Phase Components") {
                              const ParameterList& compPLList = phasePL.sublist(compPLLabel);
                              for (ParameterList::ConstIterator k=compPLList.begin(); k!=compPLList.end(); ++k) {
                                  const std::string& compLabel = compPLList.name(k);
                                  if (compPLList.isSublist(compLabel)) {
                                      const ParameterList& solutePLList = compPLList.sublist(compLabel);
                                      for (ParameterList::ConstIterator M=solutePLList.begin(); M!=solutePLList.end(); ++M) {
                                          
                                          const std::string& solutePLLabel = solutePLList.name(M);
                                          if (solutePLLabel == "Component Solutes") {
                                              Array<std::string> solutes = solutePLList.get<Array<std::string> >("Component Solutes");
                                              for (int L=0; L<solutes.size(); ++L) {
                                                  const std::string& soluteInPhase = solutes[L];

                                                  bool tracer_found_in_ics = false;
                                                  for (TracerMap::const_iterator it=ic_tmap.begin(); it!=ic_tmap.end(); ++it) {
                                                      if (it->second.Solute() == soluteInPhase) {
                                                          tracer_found_in_ics = true;
                                                          if ( ! (it->second.Phase() == phaseLabel) && (it->second.Comp() == compLabel)) {
                                                              std::cerr << "Solute in IC definition for " << it->first << " found in Phase Definitions" << std::endl;
                                                              std::cerr << "but has inconsistent phase or comp value"  << std::endl;
                                                              throw std::exception();
                                                          }
                                                      }
                                                  }

                                                  if ( ! tracer_found_in_ics ) {
                                                      std::cerr << "no IC specified for solute: " << soluteInPhase << " declared in Phase definition" << std::endl;
                                                      throw std::exception();
                                                  }
                                                  

                                                  bool tracer_found_in_bcs = false;
                                                  for (TracerMap::const_iterator it=bc_tmap.begin(); it!=bc_tmap.end(); ++it) {
                                                      if (it->second.Solute() == soluteInPhase) {
                                                          tracer_found_in_bcs = true;
                                                          if ( !(it->second.Phase() == phaseLabel) && (it->second.Comp() == compLabel)) {
                                                              std::cerr << "Solute in BC definition for " << it->first << " found in Phase Definitions" << std::endl;
                                                              std::cerr << "but has inconsistent phase or comp value"  << std::endl;
                                                              throw std::exception();
                                                          }
                                                      }
                                                  }

                                                  if ( ! tracer_found_in_bcs ) {
                                                      std::cerr << "no BC specified for solute: " << soluteInPhase << " declared in Phase definition" << std::endl;
                                                      throw std::exception();
                                                  }

                                              }
                                          }
                                          else {
                                              std::cerr << "Can only define Component Solutes at this level of the Phase Definitions" << std::endl;
                                              throw std::exception();
                                          }
                                      }
                                  }
                                  else {
                                      std::cerr << "Can only define Components at this level of the Phase Definitions" << std::endl;
                                      throw std::exception();
                                  }
                              }
                          }
                      }
                  }
              }

              // Process Solute ICs
              Array<std::string> ta(ic_tmap.size()), _labels(ic_tmap.size());
              
              TracerMap::const_iterator it=ic_tmap.begin();
              int tcnt = 0;
              while (it != ic_tmap.end()) {
                  ta[tcnt] = it->first;
                  _labels[tcnt] = underscore(ic_tmap[ta[tcnt]].Label());
                  tcnt++; it++;
              }

              struc_list.sublist("tracer").set("inits",_labels);
              for (int i=0; i<ta.size(); ++i) {
                  const TRACER& t =ic_tmap[ta[i]]; 
                  ParameterList tPL;
                  tPL.set<Array<std::string> >("regions",t.Regions());
                  tPL.set<std::string>("type",t.AMR_Type());
                  tPL.set<Array<double> >("val",t.Values());
                  struc_list.sublist("tracer").set(underscore(t.Label()),tPL);
              }

              // Process Solute BCs
              Array<std::string> tba(bc_tmap.size()), _blabels(bc_tmap.size());
              
              TracerMap::const_iterator bit=bc_tmap.begin();
              int btcnt = 0;
              while (bit != bc_tmap.end()) {
                  tba[btcnt] = bit->first;
                  _blabels[btcnt] = underscore(bc_tmap[tba[btcnt]].Label());
                  btcnt++; bit++;
              }

              struc_list.sublist("tracer").set("bcs",_blabels);
              for (int i=0; i<tba.size(); ++i) {
                  const TRACER& t =bc_tmap[tba[i]]; 
                  ParameterList tPL;
                  tPL.set<Array<std::string> >("regions",t.Regions());
                  tPL.set<std::string>("type",t.AMR_Type());
                  if (t.Values().size()>0) {
                      tPL.set<Array<double> >("vals",t.Values());
                      if (t.Values().size()>1) {
                          tPL.set<Array<double> >("times",t.Times());
                          tPL.set<Array<std::string> >("forms",t.Forms());
                      }
                  }
                  struc_list.sublist("tracer").set(underscore(t.Label()),tPL);
              }


          } // if do_tracer
      } // end of function
              

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
      Array<std::string> visNames;
      if (vlist.isParameter("Variables")) {
          visNames = vlist.get<Array<std::string> >("Variables");
          for (int i=0; i<visNames.size(); ++i) {
              visNames[i] = underscore(visNames[i]);
          }
      }
      else {
          visNames.resize(1);
          visNames[0] = "ALL";
      }
      amr_list.set<Array<std::string> >("derive_plot_vars",visNames);
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


