#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <InputParser_Structured.H>
#include <PMAMR_Labels.H>

#include <BoxLib.H>

using Teuchos::Array;
using Teuchos::ParameterList;
using Teuchos::ParameterEntry;

namespace Amanzi {
    namespace AmanziInput {

        void MyAbort(const std::string& m) {
            std::cerr << m << std::endl;
            throw std::exception();
        }

        double atmToMKS = 101325;

        std::string underscore(const std::string& instring)
        {
            std::string s = instring;
            std::replace(s.begin(),s.end(),' ','_');
            AMR_to_Amanzi_label_map[s] = instring;
            return s;
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
            Array<std::string> reqP, reqL;

            // Mesh
            std::string mesh_str = "Mesh"; reqL.push_back(mesh_str);
            PLoptions Mopt(parameter_list,reqL,reqP,false,false); 

            std::string struc_str = "Structured"; reqL.clear(); reqL.push_back(struc_str);
            PLoptions MSopt(parameter_list.sublist(mesh_str),reqL,reqP,true,true); 

            struc_list.sublist("Mesh").set("Framework","Structured");
            const ParameterList& st_list = parameter_list.sublist(mesh_str).sublist(struc_str);

            reqL.clear(); reqP.clear();
            std::string NCells_str = "Number of Cells"; reqP.push_back(NCells_str);
            std::string ProbLo_str = "Domain Low Corner"; reqP.push_back(ProbLo_str);
            std::string ProbHi_str = "Domain High Corner"; reqP.push_back(ProbHi_str);
            PLoptions MSOopt(st_list,reqL,reqP,true,true); 

            Array<int> n_cell = st_list.get<Array<int> >(NCells_str);
            domlo = st_list.get<Array<double> >(ProbLo_str);
            domhi = st_list.get<Array<double> >(ProbHi_str);


            // FIXME: This multiplier should be input
            if (domlo.size()<ndim  || domhi.size()<ndim) {
                MyAbort("Domain size nonsensical");
            }
            double max_size = domhi[0]-domlo[0];
            for (int i=0; i<ndim; ++i) {
                max_size=std::max(max_size,domhi[i]-domlo[i]);
            }
            geometry_eps = 1.e-6*max_size;

            const ParameterList& eclist = parameter_list.sublist("Execution Control");
            int bfactor = 2;
            if (eclist.isSublist("amr"))
                if (eclist.sublist("amr").isParameter("blocking_factor"))
                    bfactor = eclist.sublist("amr").get<int>("blocking_factor");
     
            for (int i=0;i<ndim;i++) {
                if (n_cell[i]%bfactor > 0) {
                    MyAbort("Number of Cells must be divisible by " + bfactor);
                }
            }

            ParameterList& alist = struc_list.sublist("amr");
            alist.set("n_cell",n_cell);

            ParameterList& glist = struc_list.sublist("geometry");
            glist.set("prob_lo",domlo);
            glist.set("prob_hi",domhi);
            glist.set("geometry_eps",geometry_eps);

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
                                      ParameterList&       struc_out_list,
                                      bool&                do_tracer)
        {
            std::string ec_str = "Execution Control";
            std::string amr_str = "Adaptive Mesh Refinement Control";
            std::string prob_str = "Basic Algorithm Control";
            std::string cg_str = "Conjugate Gradient Algorithm";
            std::string mg_str = "Multigrid Algorithm";
            std::string mac_str = "Pressure Discretization Control";
            std::string diffuse_str = "Diffusion Discretization Control";

            const ParameterList& ec_list = parameter_list.sublist(ec_str);

            Array<std::string> reqL, reqP;
            std::string flow_str = "Flow Model"; reqP.push_back(flow_str);
            std::string trans_str = "Transport Model";reqP.push_back(trans_str);
            std::string chem_str = "Chemistry Model"; reqP.push_back(chem_str);
            std::string tim_str = "Time Integration Mode"; reqL.push_back(tim_str);
            std::string v_str = "Verbosity";
            PLoptions ECopt(ec_list,reqL,reqP,false,false); 

            ParameterList& amr_out_list     = struc_out_list.sublist("amr");
            ParameterList& prob_out_list    = struc_out_list.sublist("prob");
            ParameterList& cg_out_list      = struc_out_list.sublist("cg");
            ParameterList& mg_out_list      = struc_out_list.sublist("mg");
            ParameterList& mac_out_list     = struc_out_list.sublist("mac");
            ParameterList& diffuse_out_list = struc_out_list.sublist("diffuse");

            //
            // Set flow model
            //
            std::string model_name;
            std::string flow_mode = ec_list.get<std::string>(flow_str);
            if (flow_mode == "Off") {
                MyAbort("Flow Mode = \"" + flow_mode + "\" not supported");
                prob_out_list.set("do_simple",2);
            }
            else if (flow_mode == "Richards") {
                model_name = "richard";
                prob_out_list.set("have_capillary",1);
            }
            else if (flow_mode == "Single-phase") {
                model_name = "single-phase";
                prob_out_list.set("do_simple",1);
            }
            else if (flow_mode == "Multi-phase") {
                model_name = "two-phase";
                prob_out_list.set("cfl",0.75);
            }
            prob_out_list.set("model_name",model_name);

            //
            // Set transport model
            //
            std::string transport_mode = ec_list.get<std::string>(trans_str);
            do_tracer = (transport_mode == "Off"  ?  0  :  1);
            prob_out_list.set<int>("do_tracer",do_tracer);

            //
            // Set chemistry model
            //
            std::string chem_mode = ec_list.get<std::string>(chem_str);
            if (chem_mode == "Off") {
                prob_out_list.set("do_chem",-1);
            }
            else {
                MyAbort("Chemistry Mode must be \"Off\"");
            }

            //
            // Set time evolution mode
            //
            std::string steady_str = "Steady";
            std::string transient_str = "Transient";
            std::string init_to_steady_str = "Initialize To Steady";
            const ParameterList& t_list = ec_list.sublist(tim_str);
            if (t_list.isSublist(steady_str)) {
                MyAbort("\"" + steady_str + "\" not yet supported" );
                if (model_name == "single_phase") {
                    prob_out_list.set("do_simple",2);
                }
            }
            else if (t_list.isSublist(transient_str)) {
                const ParameterList& tran_list = t_list.sublist(transient_str);
               
                reqP.clear(); 
                std::string Start_str = "Start"; reqP.push_back(Start_str);
                std::string End_str = "End"; reqP.push_back(End_str);
                std::string Init_Time_Step_str = "Initial Time Step";
                std::string Time_Step_Mult_str = "Time Step Multiplier";
                PLoptions opt(tran_list,reqL,reqP,true,false); 
                const Array<std::string> optStr = opt.OptLists();

                struc_out_list.set<double>("start_time", tran_list.get<double>(Start_str));
                struc_out_list.set<double>("stop_time", tran_list.get<double>(End_str));
                for (int i=0; i<optStr.size(); ++i) {
                    if (optStr[i] == Init_Time_Step_str) {
                        struc_out_list.set<double>("dt_init", tran_list.get<double>(Init_Time_Step_str));
                    }
                    else {
                        MyAbort("Unrecognized option under \""+transient_str+"\"");
                    }
                }
            }
            else if (t_list.isSublist(init_to_steady_str)) {
                MyAbort("\"" + init_to_steady_str + "\" not yet supported" );
            }
            else {
                std::cout << t_list << std::endl;
                MyAbort("No recognizable value for \"" + tim_str + "\"");
            }


            // Deal with optional settings
            const Array<std::string> optL = ECopt.OptLists();
            for (int i=0; i<optL.size(); ++i) {
                if (optL[i] == amr_str) {
                    const ParameterList& amr_list = ec_list.sublist(amr_str);
                    std::string max_level_str = "Max AMR Level";
                    int max_level = 0;
                    if (amr_list.isParameter(max_level_str)) {
                        max_level = amr_list.get<int>(max_level_str);
                    }
                    amr_out_list.set<int>("max_level",max_level);
                    
                    std::string ref_ratio_str = "Refinement Ratio";
                    Array<int> ref_ratio(max_level,2);
                    if (amr_list.isParameter(ref_ratio_str)) {
                        ref_ratio = amr_list.get<Array<int> >(ref_ratio_str);
                    }
                    amr_out_list.set<Array<int> >("ref_ratio",ref_ratio);
                }
                else if (optL[i] == mg_str) {
                    const ParameterList& mg_list = ec_list.sublist(mg_str);
                    for (ParameterList::ConstIterator it=mg_list.begin(); it!=mg_list.end(); ++it) {
                        const std::string& name = mg_list.name(it);
                        mg_out_list.setEntry(name,mg_list.getEntry(name));
                    }
                }
                else if (optL[i] == cg_str) {
                    const ParameterList& cg_list = ec_list.sublist(cg_str);
                    for (ParameterList::ConstIterator it=cg_list.begin(); it!=cg_list.end(); ++it) {
                        const std::string& name = cg_list.name(it);
                        cg_out_list.setEntry(name,cg_list.getEntry(name));
                    }
                }
                else if (optL[i] == prob_str) {
                    const ParameterList& prob_list = ec_list.sublist(prob_str);
                    std::string Max_Time_Step_Change_str = "Maximum Time Step Change";
                    std::string Max_Time_Step_Size_str = "Maximum Time Step Size";
                    std::string Max_Step_str = "Maximum Cycle Number";
                    for (ParameterList::ConstIterator it=prob_list.begin(); it!=prob_list.end(); ++it) {
                        const std::string& name = prob_list.name(it);
                        std::cerr << "name: " << name << std::endl;
                        if (name == Max_Time_Step_Change_str) {
                            prob_out_list.set<double>("change_max", prob_list.get<double>(Max_Time_Step_Change_str));
                        }
                        else if (name == Max_Step_str) {
                            int max_step = prob_list.get<int>(Max_Step_str);
                            struc_out_list.set<double>("max_step", max_step);
                        }
                        else if (name == Max_Time_Step_Size_str) {
                            double max_dt = prob_list.get<double>(Max_Time_Step_Size_str);
                            struc_out_list.set<double>("max_dt", max_dt);
                        }
                    }
                }
                else {
                    MyAbort("Unrecognized optional parameter to \"" + ec_str + "\" list: \"" + optL[i] + "\"");
                }
            }

            const Array<std::string> optP = ECopt.OptParms();
            int prob_v, mg_v, cg_v, amr_v, diffuse_v;
            for (int i=0; i<optP.size(); ++i) {
                if (optP[i] == v_str) {
                    std::string v_val = "Medium";
                    if (ec_list.isParameter(v_str)) {
                        v_val = ec_list.get<std::string>(v_str);
                    }

                    if (v_val == "None") {
                        prob_v = 0; mg_v = 0; cg_v = 0; amr_v = 0; diffuse_v = 0;
                    }
                    else if (v_val == "Low") {
                        prob_v = 1; mg_v = 0; cg_v = 0; amr_v = 1;  diffuse_v = 0;
                    }
                    else if (v_val == "Medium") {
                        prob_v = 1; mg_v = 0; cg_v = 0; amr_v = 1;  diffuse_v = 1;
                    }
                    else if (v_val == "High") {
                        prob_v = 2; mg_v = 1; cg_v = 1; amr_v = 2;  diffuse_v = 1;
                    }
                    else if (v_val == "Extreme") {
                        prob_v = 3; mg_v = 3; cg_v = 3; amr_v = 3;  diffuse_v = 0;
                    }
                    prob_out_list.set("v",prob_v);
                    amr_out_list.set("v",amr_v);
                    mg_out_list.set("v",mg_v);
                    cg_out_list.set("v",cg_v);
                    diffuse_out_list.set("v",diffuse_v);
                }
                else {
                    MyAbort("Unrecognized optional parameter to \"" + ec_str + "\": \"" + optP[i] + "\"");
                }
            }

            // Set other basic algorithm defaults
            prob_out_list.set("visc_abs_tol",1.e-16);
            prob_out_list.set("visc_tol",1.e-14);
            prob_out_list.set("cfl",1);
            

        }

        static std::string dirStr[6] = {"-X", "-Y", "-Z", "+X", "+Y", "+Z"};
        std::pair<bool,int>orient(const std::string& dir)
        {

            bool is_lo = false;
            int coord = -1;
            for (int i=0 ; i<6 && coord<0; ++i) {
                if (dir == dirStr[i]) {
                    coord = i%3;
                    is_lo = i>2;
                }
            }
            if (coord<0  ||  (ndim<3 && coord>1)) {
                std::cerr << dir << " is not a valid Direction value.\n";
                throw std::exception();
            }
            return std::pair<bool,int>(is_lo,coord);
        }


        void convert_Region_ColorFunction(const ParameterList& rslist,
                                          const std::string&   rlabel,
                                          ParameterList&       rsublist)
        {
            const ParameterList& rsslist = rslist.sublist(rlabel);
            rsublist.setEntry("color_file",rsslist.getEntry("File"));
            rsublist.setEntry("color_value",rsslist.getEntry("Value"));
            rsublist.set("purpose", "all");
            rsublist.set("type", "color_function");
        }

        void convert_Region_Point(const ParameterList& rslist,
                                  const std::string&   rlabel,
                                  ParameterList&       rsublist)
        {
            const ParameterList& rsslist = rslist.sublist(rlabel);
            rsublist.setEntry("coordinate",rsslist.getEntry("Coordinate"));
            rsublist.set("purpose", "all");
            rsublist.set("type", "point");
        }

        void convert_Region_Box(const ParameterList& rslist,
                                const std::string&   rlabel,
                                ParameterList&       rsublist)
        {
            const ParameterList& rsslist = rslist.sublist(rlabel);
            Array<double> lo = rsslist.get<Array<double> >("Low Coordinate");
            Array<double> hi = rsslist.get<Array<double> >("High Coordinate");

            std::string purpose, type;
            for (int d=0; d<ndim; ++d) {
                if (std::abs(hi[d] - lo[d]) < geometry_eps) // This is a (ndim-1) dimensional region
                {
                    lo[d] = hi[d];
                    type = "surface";
                    purpose = "all";
                    
                    // Is this on the domain boundary?
                    if (lo[d] == domlo[d]) {
                        purpose = PMAMR::RpurposeDEF[d];
                    }
                    else if (lo[d] == domhi[d]) {
                        purpose = PMAMR::RpurposeDEF[d+3];
                    }
                }
                else {
                    type = "box";
                    purpose = "all";
                }
            }

            rsublist.set<Array<double> >("lo_coordinate",lo);
            rsublist.set<Array<double> >("hi_coordinate",hi);
            rsublist.set("purpose", purpose);
            rsublist.set("type", type);
        }

        void convert_Region_Plane(const ParameterList& rslist,
                                  const std::string&   rlabel,
                                  ParameterList&       rsublist)
        {
            const ParameterList& rsslist = rslist.sublist(rlabel);
            std::string dir = rsslist.get<std::string>("Direction");
            std::pair<bool,int> o = orient(dir); bool is_lo=o.first; int coord=o.second;

            // FIXME: Ignores sign of direction

            double loc = rsslist.get<double>("Location");
          
            Array<double> lo(domlo);
            Array<double> hi(domhi);
            lo[coord] = loc;
            hi[coord] = loc;
          
            rsublist.set("lo_coordinate",lo);
            rsublist.set("hi_coordinate",hi);
            rsublist.set("type", "surface");
            int iPurpose = 6;

            if (std::abs(lo[coord] - domlo[coord]) < geometry_eps) {
                iPurpose = (coord==0 ?  0 : (ndim==2 || coord==1 ? 1 : 2) );
            }
            else if (std::abs(hi[coord] - domhi[coord]) < geometry_eps) {
                iPurpose = (coord==0 ?  3 : (ndim==2 || coord==1 ? 4 : 5) );
            }
            std::string purpose = underscore(PMAMR::RpurposeDEF[iPurpose]);
            rsublist.set("purpose",purpose);
        }


        static std::string RlabelDEF[7] = {"XLOBC", "YLOBC", "ZLOBC", "XHIBC", "YHIBC", "ZHIBC", "ALL"};
        Array<std::string>
        generate_default_regions(ParameterList& rsublist)
        {
            Array<std::string> def_regionNames;
            ParameterList t1PL, t2PL, t3PL;
            t1PL.set<Array<double> >("Low Coordinate",domlo);
            t1PL.set<Array<double> >("High Coordinate",domhi);
            std::string blabel = "Region: Box";
            t2PL.set(blabel,t1PL);
            convert_Region_Box(t2PL,blabel,t3PL);
            rsublist.set(RlabelDEF[6],t3PL);
            def_regionNames.push_back(RlabelDEF[6]);

            std::string dir_name = "Direction";
            std::string loc_name = "Location";
            for (int i=0 ; i<6; ++i) {
                if (i%3<ndim) {
                    ParameterList rslist, rsslist, tmp;
                    rsslist.set(dir_name,dirStr[i]);
                    double loc = (i>2 ? domhi[i%3] : domlo[i]);
                    rsslist.set(loc_name,loc);
                    rslist.set("Region: Plane",rsslist);
                    convert_Region_Plane(rslist,"Region: Plane",tmp);
                    rsublist.set(RlabelDEF[i],tmp);
                    def_regionNames.push_back(RlabelDEF[i]);
                }
            }   
            return def_regionNames;
        }


        //
        // convert region to structured format
        //
        void
        convert_to_structured_region(const ParameterList& parameter_list, 
                                     ParameterList&       struc_list)
        {
            ParameterList& geom_list = struc_list.sublist("geometry");

            // Initialize with default region set
            Array<std::string> arrayregions = generate_default_regions(geom_list);

            const ParameterList& rlist = parameter_list.sublist("Regions");
            for (ParameterList::ConstIterator i=rlist.begin(); i!=rlist.end(); ++i) {
        
                std::string label = rlist.name(i);
                std::string _label = underscore(label);
                const ParameterEntry& entry = rlist.getEntry(label);
        
                if ( !entry.isList() ) {
                    if (Teuchos::MPISession::getRank() == 0) {
                        std::cerr << "Region section must define only regions. \"" 
                                  << label << "\" is not a valid region definition." << std::endl;
                    }
                    throw std::exception();
                }
        
                ParameterList rsublist;
                const ParameterList& rslist = rlist.sublist(label);

                // Add user-defined regions
                for (ParameterList::ConstIterator j=rslist.begin(); j!=rslist.end(); ++j) {
          
                    const std::string& rlabel = rslist.name(j);
                    const ParameterEntry& rentry = rslist.getEntry(rlabel);
          
                    if (rentry.isList()) {
                        if (rlabel=="Region: Color Function") {
                            convert_Region_ColorFunction(rslist,rlabel,rsublist);
                        }
                        else if (rlabel=="Region: Point"){	      
                            convert_Region_Point(rslist,rlabel,rsublist);
                        }
                        else if (rlabel=="Region: Box") {
                            convert_Region_Box(rslist,rlabel,rsublist);
                        }
                        else if (rlabel=="Region: Plane") {
                            convert_Region_Plane(rslist,rlabel,rsublist);
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
                //mtest.insert(MTEST::value_type("Density",false)); // This is not used anywhere
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
                            Array<std::string> tmp_regions = rslist.get<Array<std::string> >(rlabel);
                            for (int j=0; j<tmp_regions.size(); ++j) {
                                tmp_regions[j] = underscore(tmp_regions[j]);
                            }
                            rsublist.set("regions",tmp_regions);
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
      
        StateDef::StateDef(const ParameterList& parameter_list)
            : parameter_list(parameter_list)
        {
            build_state_def();
        }

        void
        StateDef::clear()
        {
            phase_comp_map.clear();
            state_ics.clear();
            state_bcs.clear();
        }

        void
        StateDef::build_state_def()
        {
            clear();

            Array<std::string> reqL, nullList;
            reqL.push_back("Phases");
            PLoptions opt(parameter_list,reqL,nullList,false,false);
            const ParameterList& plist = parameter_list.sublist(reqL[0]);

            PLoptions optP(plist,nullList,nullList,false,false); // each optional list is a phase 
            const Array<std::string>& phaseLabels = optP.OptLists();

            for (int i=0; i<phaseLabels.size(); ++i) {
                const std::string& phaseLabel = phaseLabels[i];
                const ParameterList& psublist = plist.sublist(phaseLabel);
                const ParameterList& pplist = plist.sublist(phaseLabel);

                Array<std::string> reqLp;
                reqLp.push_back("Phase Properties");
                reqLp.push_back("Phase Components");
                PLoptions optP1(pplist,reqLp,nullList,true,true);

                // Start with a phase def and then add components
                const ParameterList& ppsublist = pplist.sublist(reqLp[0]);

                PLoptions optPP(ppsublist,nullList,nullList,false,true);

                Array<std::string> optLpp = optPP.OptLists();
                
                double density=-1;
                double viscosity=-1; // FIXME: Assumes constant is only model so stores values
                double diffusivity=-1;
                for (int j=0; j<optLpp.size(); ++j) {
                    const std::string& propLabel = optLpp[j];

                    const ParameterList& ppolist = ppsublist.sublist(propLabel);

                    if (propLabel == "Phase Mass Density: Uniform") {
                        Array<std::string> reqP;
                        reqP.push_back("Value");
                        PLoptions optPD(ppolist,nullList,reqP,true,true); 
                        density = ppolist.get<double>(reqP[0]);
                    }
                    else if (propLabel == "Phase Viscosity: Uniform") {
                        Array<std::string> reqP;
                        reqP.push_back("Value");
                        PLoptions optPD(ppolist,nullList,reqP,true,true); 
                        viscosity= ppolist.get<double>(reqP[0]);
                    }
                    else {
                        std::cerr << "Unrecognized phase property parameter: " << propLabel << std::endl;
                        throw std::exception();
                    }
                    
                }
                if (density<0 || viscosity<0) {
                    std::cerr << "Must define density and viscosity for each phase present" << std::endl;
                    throw std::exception();
                }
                getPhases().insert(std::pair<std::string,PHASE>
                                   (phaseLabel,
                                    PHASE(density,viscosity,diffusivity)));
                

                const ParameterList& pclist = plist.sublist(phaseLabel).sublist(reqLp[1]);

                PLoptions optC(pclist,nullList,nullList,false,true);
                const Array<std::string>& cLabels = optC.OptLists(); // each optional list names a component
                for (int j=0; j<cLabels.size(); ++j) {
                    const std::string& compLabel = cLabels[j];

                    Array<std::string> reqLc, reqPc;
                    const ParameterList& slist = pclist.sublist(compLabel);

                    PLoptions optCC(slist,reqLc,reqPc,true,false); 
                    const Array<std::string>& sParams = optCC.OptParms();
                  
                    Array<std::string> sLabels;
                    for (int k=0; k<sParams.size(); ++k) {
                        if (sParams[k] == "Component Solutes") {
                            sLabels = pclist.sublist(compLabel).get<Array<std::string> >(sParams[k]);
                        }
                    }

                    COMP& c = (*this)[phaseLabel][compLabel];
                    for (int L=0; L<sLabels.size(); ++L) {
                        c.push_back(sLabels[L]);
                    }
                }
            }

            build_solute_funcs(state_bcs,"Boundary Conditions","Solute BC");
            build_solute_funcs(state_ics,"Initial Conditions","Solute IC");
          
        }


        void
        StateDef::build_solute_funcs(StateFuncMap& s,
                                     const std::string& top_level_label,
                                     const std::string& solute_section_label)
        {
            Array<std::string> nullList;
            Array<std::string> reqL, reqP;
            reqL.push_back(top_level_label);

            PLoptions optBC(parameter_list,reqL,nullList,false,false); 
            const ParameterList& plist = parameter_list.sublist(reqL[0]);

            PLoptions BCs(plist,nullList,nullList,false,true);  
            const Array<std::string>& bcLabels = BCs.OptLists();

            for (int i=0; i<bcLabels.size(); ++i)
            {
                //
                // Under BC:
                //  (1) Region assignment
                //  (2) Solute BC 
                //  (3) Phase/comp BC
                //
                const std::string& BClabel = bcLabels[i];
                const ParameterList& bc_plist = plist.sublist(BClabel);

                reqP.clear(); reqP.push_back("Assigned Regions");
                reqL.clear(); reqL.push_back(solute_section_label);
                PLoptions BCs(bc_plist,reqL,reqP,false,true);  

                // The optional list is for the phase/comp func
                const Array<std::string>& phaseBCfuncs = BCs.OptLists();
                if (phaseBCfuncs.size()!=1) {
                    std::cerr << "Phase/comp BCs and ICs required for BC/IC label: " << BClabel << std::endl;
                    throw std::exception();
                }
                const ParameterList& func_plist = bc_plist.sublist(phaseBCfuncs[0]);

                const Array<std::string>& assigned_regions = bc_plist.get<Array<std::string> >(reqP[0]);
                const std::string& Amanzi_type = phaseBCfuncs[0];              

                s[BClabel] = StateFunc(BClabel, Amanzi_type, func_plist, assigned_regions);

                // Now add solute data
                Array<std::string> nullList;
                const ParameterList& solute_plist = bc_plist.sublist(reqL[0]);
                PLoptions soluteOPT(solute_plist,nullList,nullList,false,false); // Expect only phase names here
                const Array<std::string>& phaseNames = soluteOPT.OptLists();

                for (int icp=0; icp<phaseNames.size(); ++icp) {
                    const ParameterList& phasePL = solute_plist.sublist(phaseNames[icp]);
                    PLoptions soluteOPTc(phasePL,nullList,nullList,false,true);  // Expect only comp names here
                    const Array<std::string>& compNames = soluteOPTc.OptLists();
                  
                    for (int icc=0; icc<compNames.size(); ++icc) {
                        const ParameterList& compPL = phasePL.sublist(compNames[icc]);
                        PLoptions soluteOPTs(compPL,nullList,nullList,false,true); // Expect only solute names here
                        const Array<std::string>& soluteNames = soluteOPTs.OptLists();

                        for (int ics=0; ics<soluteNames.size(); ++ics) {                          
                            const ParameterList& solutePL = compPL.sublist(soluteNames[ics]);
                          
                            Array<std::string> funcL, funcP;
                            funcP.push_back("Concentration Units");
                            PLoptions soluteOPTf(solutePL,nullList,funcP,false,true);
                          
                            // Get units
                            const std::string& units = solutePL.get<std::string>(funcP[0]);

                            // Get function name/list
                            const Array<std::string>& funcNames = soluteOPTf.OptLists();
                            if (funcNames.size()!=1) {
                                std::cout << "Each solute BC expects a single function" << std::endl;
                                throw std::exception();
                            }
                            const std::string& Amanzi_solute_type = funcNames[0];
                            ParameterList solute_func_plist = solutePL.sublist(Amanzi_solute_type);

                            s[BClabel][phaseNames[icp]][compNames[icc]][soluteNames[ics]]
                                = ICBCFunc(solute_func_plist,Amanzi_solute_type,units,BClabel);
                        }
                    }
                }
              

                // Confirm that embedded phase/comp and solute names are compatible with state
                StateFunc::PhaseFuncMap& pfm = s[BClabel].getPhaseFuncMap();
                for (StateFunc::PhaseFuncMap::iterator pit=pfm.begin(); pit!=pfm.end(); ++pit) {
                    const std::string& phaseName = pit->first;
                    PhaseFunc::CompFuncMap& cfm = pit->second.getCompFuncMap();
                    for (PhaseFunc::CompFuncMap::iterator cit=cfm.begin(); cit!=cfm.end(); ++cit) {
                        const std::string& compName = cit->first;
                        CompFunc::ICBCFuncMap& fm = cit->second.getICBCFuncMap();
                        for (CompFunc::ICBCFuncMap::iterator fit=fm.begin(); fit!=fm.end(); ++fit) {
                            const std::string& soluteName = fit->first;
                            const Array<std::string>& ds=(*this)[phaseName][compName].getTracerArray();
                            bool found = false;
                            for (int it=0; it<ds.size() && !found; ++it) {
                                found = ds[it] == soluteName;
                            }
                            if (!found) {
                                std::cerr << "function: phase/comp/solute not in Phase Definition: "
                                          << fit->second.Amanzi_Type() << std::endl;
                                throw std::exception();
                            }                          
                        }
                    }
                }
            }

        }



        void convert_ICSaturation(const ParameterList& fPLin,
                                  const std::string&   Amanzi_type,
                                  ParameterList&       fPLout)
        {
            fPLout.set<std::string>("type","saturation");
            Array<std::string> nullList, reqP;
            const std::string val_name="Value"; reqP.push_back(val_name);
            PLoptions opt(fPLin,nullList,reqP,true,true); 
            fPLout.set<double>("val",fPLin.get<double>(val_name));
        }


        void convert_ICPressure(const ParameterList& fPLin,
                                const std::string&   Amanzi_type,
                                ParameterList&       fPLout)
        {
            const std::string val_name="Reference Value";
            const std::string grad_name="Gradient Value";
            const std::string ref_name="Reference Coordinate";

            Array<std::string> reqP, nullList;
            reqP.push_back(val_name);
            if (Amanzi_type == "IC: Linear Pressure") {
                reqP.push_back(grad_name);
                reqP.push_back(ref_name);
            }
            PLoptions opt(fPLin,nullList,reqP,true,true);  
    
            fPLout.set<std::string>("type","hydrostatic");
            fPLout.set<double>("val",fPLin.get<double>(val_name));
    
            if (Amanzi_type == "IC: Linear Pressure") {
                fPLout.set<Array<double> >("grad",fPLin.get<Array<double> >(grad_name));
                const Array<double>& water_table = fPLin.get<Array<double> >(ref_name);
                int coord = water_table.size()-1;
                fPLout.set<double>("water_table_height",water_table[coord]);                      
            }
        }

        void convert_ICHydrostatic(const ParameterList& fPLin,
                                   const std::string&   Amanzi_type,
                                   ParameterList&       fPLout)
        {
            Array<std::string> reqP, nullList;
            const std::string ref_name="Water Table Height"; reqP.push_back(ref_name);
            PLoptions opt(fPLin,nullList,reqP,true,true);  
            
            fPLout.set<double>("water_table_height",fPLin.get<double>(ref_name));
            fPLout.set<std::string>("type","hydrostatic");
        }

        void convert_solute_ICConcentration(const ICBCFunc& solute_ic,
                                            ParameterList&  fPLout)
        {
            const ParameterList& fPLin = solute_ic.PList();
            const std::string& solute_ic_Amanzi_type = solute_ic.Amanzi_Type();
            const std::string& solute_ic_label = solute_ic.Label();
            const std::string& solute_ic_units = solute_ic.Units();

            Array<std::string> reqP, nullList;
            const std::string val_name="Value"; reqP.push_back(val_name);
            PLoptions opt(fPLin,nullList,reqP,true,true);  
            fPLout.set<double>("val",fPLin.get<double>(val_name));
    
            // Adjust dimensions of data
            if (solute_ic_units=="Molar Concentration" 
                || solute_ic_units=="Molal Concentration")
            {
                std::cerr << "IC label \"" << solute_ic_label
                          << "\" function: \"" << solute_ic_Amanzi_type
                          << "\" requests unsupported units: \"" << solute_ic_units
                          << "\"" << std::endl;
                throw std::exception();
            }
            else if (solute_ic_units=="Specific Concentration") {
                // This is the units expected by the structured code
            }
            else {
                std::cerr << "Unsupported Solute IC function: \"" << solute_ic_units << "\"" << std::endl;
                throw std::exception();
            }
    
            fPLout.set<std::string>("type","concentration");
        }

        void
        convert_ics(const ParameterList& parameter_list, 
                    ParameterList&       struc_list,
                    StateDef&            stateDef)
        {
            ParameterList& comp_list  = struc_list.sublist("comp"); 

            Array<std::string> ic_label_list;
            Array<std::string> reqL, reqP, nullList;
            StateFuncMap& state_ics = stateDef.IC();


            ParameterList icPLout_master;
            for (StateFuncMap::iterator ic_it = state_ics.begin(); ic_it!=state_ics.end(); ++ic_it)
            {
                const std::string& IClabel = ic_it->first;
                ic_label_list.push_back(underscore(IClabel));
                StateFunc& state_ic = ic_it->second;
        
                const Array<std::string>& regions = state_ic.Regions();
                Array<std::string> _regions;
                for (int i=0; i<regions.size(); ++i) {
                    _regions.push_back(regions[i]);
                }
                ParameterList fPLout;
                fPLout.set("regions",_regions);
        
                // Phase/comp ICs
                const ParameterList& fPLin= state_ic.FuncPList();
                //const std::string& units = state_ic.units();
                const std::string& Amanzi_type = state_ic.Amanzi_Type();
        
                if (Amanzi_type == "IC: Uniform Saturation"
                    || Amanzi_type == "IC: Linear Saturation")
                {
                    convert_ICSaturation(fPLin,Amanzi_type,fPLout);
                }
                else if ( (Amanzi_type == "IC: Uniform Pressure"
                           || Amanzi_type == "IC: Linear Pressure") )
                {
                    convert_ICPressure(fPLin,Amanzi_type,fPLout);
                }
                else if ( Amanzi_type == "IC: Hydrostatic" )
                {
                    convert_ICHydrostatic(fPLin,Amanzi_type,fPLout);
                }
                else {
                    std::cerr << "Unsupported IC: " << Amanzi_type << std::endl;
                    throw std::exception();
                }           

                icPLout_master.set(underscore(IClabel),fPLout);

            }

            comp_list.set<Array<std::string> >("ic_labels",ic_label_list);
            comp_list.set("ics",icPLout_master);

        }

        void convert_BCNoFlow(const ParameterList& fPLin,
                              const std::string&   Amanzi_type,
                              ParameterList&       fPLout)
        {
            Array<std::string> nullList;
            PLoptions opt(fPLin,nullList,nullList,true,true); 
            fPLout.set<std::string>("type","noflow");
        }

        void convert_BCSaturation(const ParameterList& fPLin,
                                  const std::string&   Amanzi_type,
                                  ParameterList&       fPLout)
        {
            // FIXME: Assumes single Aqueous saturation value specified
            fPLout.set<std::string>("type","saturation");
            Array<std::string> nullList, reqP;
            const std::string val_name="Values"; reqP.push_back(val_name);
            const std::string time_name="Times"; reqP.push_back(time_name);
            const std::string form_name="Time Functions"; reqP.push_back(form_name);
            PLoptions opt(fPLin,nullList,reqP,true,true); 
            const Array<double>& vals = fPLin.get<Array<double> >(val_name);
            fPLout.set<Array<double> >("vals",vals);
            if (vals.size()>1) {
                fPLout.set<Array<double> >("times",fPLin.get<Array<double> >(time_name));
                fPLout.set<Array<double> >("forms",fPLin.get<Array<double> >(form_name));
            }
        }


        void convert_BCPressure(const ParameterList& fPLin,
                                const std::string&   Amanzi_type,
                                ParameterList&       fPLout)
        {
            Array<std::string> nullList, reqP;
            fPLout.set<std::string>("type","pressure");

            if (Amanzi_type == "BC: Linear Pressure")
            {
                const std::string val_name="Reference Values"; reqP.push_back(val_name);
                const std::string grad_name="Gradient Value";reqP.push_back(grad_name);
                const std::string ref_name="Reference Coordinate"; reqP.push_back(ref_name);
                PLoptions opt(fPLin,nullList,reqP,true,true);

                fPLout.set<Array<double> >("vals",fPLin.get<Array<double> >(val_name));
                fPLout.set<Array<double> >("grad",fPLin.get<Array<double> >(grad_name));
                const Array<double>& water_table = fPLin.get<Array<double> >(ref_name);
                int coord = water_table.size()-1;
                fPLout.set<double>("water_table",water_table[coord]);
            }
            else if (Amanzi_type == "BC: Uniform Pressure")
            {
                const std::string val_name="Values"; reqP.push_back(val_name);
                const std::string time_name="Times"; reqP.push_back(time_name);
                const std::string form_name="Time Functions"; reqP.push_back(form_name);
                PLoptions opt(fPLin,nullList,reqP,true,true);  

                Array<double> vals = fPLin.get<Array<double> >(val_name);
                fPLout.set<Array<double> >("vals",vals);
                if (vals.size()>1) {
                    fPLout.set<Array<double> >("times",fPLin.get<Array<double> >(time_name));
                    fPLout.set<Array<double> >("forms",fPLin.get<Array<double> >(form_name));
                }
            }
        }

        void convert_BCFlux(const ParameterList& fPLin,
                            const std::string&   Amanzi_type,
                            ParameterList&       fPLout)
        {
            // FIXME: Assumes that the flux specified is that of the Aqueous phase
            bool is_in_vol = fPLin.isParameter("Inward Volumetric Flux");
            bool is_in_mass = fPLin.isParameter("Inward Mass Flux");
            bool is_out_vol = fPLin.isParameter("Outward Volumetric Flux");
            bool is_out_mass = fPLin.isParameter("Outward Mass Flux");
    
            bool is_mass = is_in_mass || is_out_mass;
            bool is_out = is_out_vol || is_out_mass;
    
            Array<std::string> reqL, reqP;
            std::string val_name;
            if (is_in_vol) {
                val_name = "Inward Volumetric Flux";
            }
            else if (is_out_vol) {
                val_name = "Outward Volumetric Flux";
            }
            else if (is_in_mass) {
                val_name = "Inward Mass Flux";
            }
            else if (is_out_mass) {
                val_name = "Outward Mass Flux";
            }
            else {
                std::cerr << "Flux not specified in recognized form" << std::endl;
                throw std::exception();
            }
    
            reqP.clear(); reqP.push_back(val_name);
            const std::string time_name="Times"; reqP.push_back(time_name);
            const std::string form_name="Time Functions"; reqP.push_back(form_name);
            const std::string rock_name = "Material Type at Boundary";
            bool require_rock = true;
            if (require_rock) {                
                reqP.push_back(rock_name);
            }
    
            PLoptions opt(fPLin,reqL,reqP,false,true);
    
    
            Array<double> times, fluxvals = fPLin.get<Array<double> >(val_name);
            Array<std::string> forms;
            if (fluxvals.size()>1) {
                times = fPLin.get<Array<double> >(time_name);
                forms = fPLin.get<Array<std::string> >(form_name);
            }
    
            // Convert mass flux to volumetric flux
            if (is_mass) {
                std::cerr << "Mass fluxes not yet supported" << std::endl;
                throw std::exception();
            }
    
            // Convert to inward flux
            if (is_out) {
                for (int i=0; i<fluxvals.size(); ++i) {
                    fluxvals[i] = -fluxvals[i];
                }
            }
    
            std::string rock_label;
            if (require_rock) {
                rock_label = fPLin.get<std::string>(rock_name);
            }
            else{ 
                const Array<std::string>& optional_lists = opt.OptLists();
                if (optional_lists.size() > 0) {
            
                    if (optional_lists.size()!=1 || optional_lists[0] != rock_name) {
                        std::cerr << "BC: Flux - invalid optional arg(s): ";
                        for (int i=0; i<optional_lists.size(); ++i) {
                            std::cerr << "\"" << optional_lists[i] << "\" ";
                        }
                        std::cerr << std::endl;
                    }
                    else {
                        rock_label = fPLin.get<std::string>(rock_name);
                    }
                }
            }
    
            fPLout.set<Array<double> >("aqueous_vol_flux",fluxvals);
            if (fluxvals.size() > 1) {
                fPLout.set<Array<double> >("inflowtimes",fPLin.get<Array<double> >(time_name));
                fPLout.set<Array<std::string> >("inflowfncs",fPLin.get<Array<std::string> >(form_name));
            }
            if (require_rock) {
                fPLout.set<std::string>("rock",rock_label);
            }
            fPLout.set<std::string>("type","zero_total_velocity");
        }

        void convert_solute_BCConcentration(const ICBCFunc&    solute_bc,
                                            ParameterList&     fPLout)
        {
            const ParameterList& fPLin = solute_bc.PList();
            const std::string& solute_bc_Amanzi_type = solute_bc.Amanzi_Type();
            const std::string& solute_bc_label = solute_bc.Label();
            const std::string& solute_bc_units = solute_bc.Units();

            Array<std::string> reqP, nullList;
            const std::string val_name="Values"; reqP.push_back(val_name);
            const std::string time_name="Times"; reqP.push_back(time_name);
            const std::string form_name="Time Functions"; reqP.push_back(form_name);
    
            PLoptions opt(fPLin,nullList,reqP,true,true);  
            Array<double> vals = fPLin.get<Array<double> >(val_name);
            fPLout.set<Array<double> >("vals",vals);
            if (vals.size() > 1) {
                fPLout.set<Array<double> >("times",fPLin.get<Array<double> >(time_name));
                fPLout.set<Array<std::string> >("forms",fPLin.get<Array<std::string> >(form_name));
            }
    
            // Adjust dimensions of data
            if (solute_bc_units=="Molar Concentration" 
                || solute_bc_units=="Molal Concentration")
            {
                std::cerr << "BC label \"" << solute_bc_label
                          << "\" function: \"" << solute_bc_Amanzi_type
                          << "\" requests unsupported units: \"" << solute_bc_units
                          << "\"" << std::endl;
                throw std::exception();
            }
            else if (solute_bc_units=="Specific Concentration") {
                // This is the units expected by the structured code
            } 
            else {
                std::cerr << "Solute BC - invalid units: \"" << solute_bc_units << "\"" << std::endl;
                throw std::exception();
            }
    
            fPLout.set<std::string>("type","concentration");
        }


        void convert_solute_BCOutflow(const ICBCFunc& solute_bc,
                                      ParameterList&  fPLout)
        {
            Array<std::string> nullList;
            const ParameterList& fPLin = solute_bc.PList();
            PLoptions opt(fPLin,nullList,nullList,true,true);  
            fPLout.set<std::string>("type","outflow");
        }

        void convert_solute_BCNoflow(const ICBCFunc& solute_bc,
                                     ParameterList&  fPLout)
        {
            Array<std::string> nullList;
            const ParameterList& fPLin = solute_bc.PList();
            PLoptions opt(fPLin,nullList,nullList,true,true);  
            fPLout.set<std::string>("type","noflow");
        }

        void
        convert_bcs(const ParameterList& parameter_list, 
                    ParameterList&       struc_list,
                    StateDef&            stateDef)
        {
            ParameterList& comp_list  = struc_list.sublist("comp"); 
            ParameterList& press_list  = struc_list.sublist("press"); 

            Array<std::string> bc_label_list;
            Array<std::string> reqL, reqP, nullList;
            StateFuncMap& state_bcs = stateDef.BC();

            ParameterList bcPLout_master;
            for (StateFuncMap::iterator bc_it = state_bcs.begin(); bc_it!=state_bcs.end(); ++bc_it)
            {
                const std::string& BClabel = bc_it->first;
                bc_label_list.push_back(underscore(BClabel));
                StateFunc& state_bc = bc_it->second;
        
                const Array<std::string>& regions = state_bc.Regions();
                Array<std::string> _regions;
                for (int i=0; i<regions.size(); ++i) {
                    _regions.push_back(underscore(regions[i]));
                }
                ParameterList fPLout;
                fPLout.set("regions",_regions);
        
                // Phase/comp BCs
                const ParameterList& fPLin= state_bc.FuncPList();
                //const std::string& units = state_bc.units();
                const std::string& Amanzi_type = state_bc.Amanzi_Type();

                if (Amanzi_type == "BC: Uniform Saturation"
                    || Amanzi_type == "BC: Linear Saturation")
                {
                    convert_BCSaturation(fPLin,Amanzi_type,fPLout);
                }
                else if ( (Amanzi_type == "BC: Uniform Pressure"
                           || Amanzi_type == "BC: Linear Pressure") )
                {
                    convert_BCPressure(fPLin,Amanzi_type,fPLout);
                }
                else if (Amanzi_type == "BC: Flux")
                {
                    convert_BCFlux(fPLin,Amanzi_type,fPLout);
                }
                else if (Amanzi_type == "BC: No Flow")
                {
                    convert_BCNoFlow(fPLin,Amanzi_type,fPLout);
                }
                else {
                    std::cerr << "Unsupported BC: \"" << Amanzi_type << "\"" << std::endl;
                    throw std::exception();
                }           

                bcPLout_master.set(underscore(BClabel),fPLout);

            }
            comp_list.set<Array<std::string> >("bc_labels",bc_label_list);
            comp_list.set("bcs",bcPLout_master);

            // After all bcs set, scan through to check that each orientation has only a single "type" of bc
            // and then set this type into the translated PL

    
            ParameterList& geom_list = struc_list.sublist("geometry");
            typedef std::pair<std::string,std::string> Spair;
            std::map<std::string,Array<Spair> > orient_RT_map;
            for (int i=0; i<bc_label_list.size(); ++i) {
                const std::string& bc_label = bc_label_list[i];
                const ParameterList& bc_sublist = comp_list.sublist("bcs").sublist(bc_label);

                const std::string& bc_type = bc_sublist.get<std::string>("type");

                const Array<std::string>& bc_regions = bc_sublist.get<Array<std::string> >("regions");
                for (int j=0; j<bc_regions.size(); ++j) {
                    const std::string& regionName = bc_regions[j];
                    if (geom_list.isSublist(regionName)) {
                        const std::string& purpose = geom_list.sublist(regionName).get<std::string>("purpose");
                        orient_RT_map[purpose].push_back(Spair(regionName,bc_type));
                    }
                    else {
                        std::cerr << "BC: \"" << AMR_to_Amanzi_label_map[bc_label] 
                                  << "\" refers to undefined region \""
                                  << AMR_to_Amanzi_label_map[regionName] << "\"" << std::endl;
                        throw std::exception();                        
                    }
                }
            }

            Array<int> inflow_hi_bc(ndim), inflow_lo_bc(ndim);
            Array<int> hi_bc(ndim), lo_bc(ndim);
            Array<int> phi_bc(ndim), plo_bc(ndim);
            Array<double> press_lo(ndim), press_hi(ndim);
            Array<double> inflow_lo_vel(ndim), invlow_hi_vel(ndim);

            for (int i=0; i<6; ++i)
            {
                int k = i%3;
                if (k < ndim) {
                    const std::string& Amanzi_purpose = AMR_to_Amanzi_label_map[PMAMR::RpurposeDEF[i]];
                    const Array<Spair>& orient_RTs = orient_RT_map[Amanzi_purpose];
                    if (orient_RTs.size()==0) {
                        std::cerr << "No boundary conditions found for " << Amanzi_purpose << std::endl;
                        throw std::exception();
                    }

                    const std::string& orient_type = orient_RTs[0].second;
                    for (int j=1; j<orient_RTs.size(); ++j) {
                        if (orient_type != orient_RTs[j].second) {
                            std::cerr << "Structured grid requires that all BCs on "
                                      << Amanzi_purpose << " be of the same type" << std::endl;
                            throw std::exception();
                        }
                    }

                    int& sat_bc      = (i<3  ?        lo_bc[k] :  hi_bc[k]);
                    int& pressure_bc = (i<3  ?       plo_bc[k] : phi_bc[k]);
                    int& inflow_bc   = (i<3  ? inflow_lo_bc[k] : inflow_hi_bc[k]);

                    if (orient_type == "noflow") {
                        sat_bc      = 4; // No flow for saturation
                        pressure_bc = 4; // No flow for p
                        inflow_bc   = 0; // Automatically set to zero
                    }
                    else if (orient_type == "pressure" || orient_type == "hydrostatic") {
                        // Must set components by name, and phase press set by press_XX(scalar) or hydro 
                        sat_bc      = 1; // Dirichlet for saturation,
                        pressure_bc = 2; // Dirichlet for p
                        inflow_bc   = 0; // Unused
                    }
                    else if (orient_type == "saturation") {
                        // Must set components by name, and phase press set by press_XX(scalar) or hydro 
                        sat_bc      = 1; // Dirichlet for saturation,
                        pressure_bc = 2; // Dirichlet for p
                        inflow_bc   = 0; // Unused
                    }
                    else if (orient_type == "zero_total_velocity") {
                        sat_bc      = 1; // Dirichlet for saturation (but will compute values based on p)
                        pressure_bc = 1; // Dirichlet for p
                        inflow_bc   = 1; // Requires inflow_XX_vel velocity values for nphase-1 phases
                    }
                    else {
                        std::cerr << "Structured grid inputs translator generated unrecognized BCs "
                                  << orient_type << std::endl;
                        throw std::exception();
                    }
                }
            }

            comp_list.set("lo_bc",lo_bc);
            comp_list.set("hi_bc",hi_bc);
            press_list.set("lo_bc",plo_bc);
            press_list.set("hi_bc",phi_bc);
            press_list.set("inflow_bc_lo",inflow_lo_bc);
            press_list.set("inflow_bc_hi",inflow_hi_bc);

        }

        typedef std::multimap<std::string,ParameterList> SolutePLMMap;

        SolutePLMMap
        convert_solute_bcs(StateDef& stateDef)
        {
            SolutePLMMap solute_to_BClabel;
            StateFuncMap& state_bcs = stateDef.BC();    
            for (StateFuncMap::iterator bc_it = state_bcs.begin(); bc_it!=state_bcs.end(); ++bc_it)
            {
                const std::string& BClabel = bc_it->first;
                StateFunc& state_bc = bc_it->second;
                const Array<std::string>& regions = state_bc.Regions();
                Array<std::string> _regions;
                for (int i=0; i<regions.size(); ++i) {
                    _regions.push_back(underscore(regions[i]));
                }
                const std::string& Amanzi_type = state_bc.Amanzi_Type();
        
                //
                // Scan through all phases, comps to find solute BCs organized by solute
                //
                StateFunc::PhaseFuncMap& pfm = state_bcs[BClabel].getPhaseFuncMap();
                for (StateFunc::PhaseFuncMap::iterator pit=pfm.begin(); pit!=pfm.end(); ++pit) {
                    const std::string& phaseName = pit->first;
                    PhaseFunc::CompFuncMap& cfm = pit->second.getCompFuncMap();
                    for (PhaseFunc::CompFuncMap::iterator cit=cfm.begin(); cit!=cfm.end(); ++cit) {
                        const std::string& compName = cit->first;
                        CompFunc::ICBCFuncMap& fm = cit->second.getICBCFuncMap();
                        for (CompFunc::ICBCFuncMap::iterator fit=fm.begin(); fit!=fm.end(); ++fit) {
                            const std::string& soluteName = fit->first;
                    
                            ICBCFunc& solute_bc = state_bc[phaseName][compName][soluteName];
                    
                            const ParameterList& fPLin = solute_bc.PList();
                            const std::string& solute_bc_Amanzi_type = solute_bc.Amanzi_Type();
                            const std::string& solute_bc_label = solute_bc.Label();
                            const std::string& solute_bc_units = solute_bc.Units();
                    
                            ParameterList fPL;
                            if (solute_bc_Amanzi_type == "BC: Uniform Concentration")
                            {
                                convert_solute_BCConcentration(solute_bc,fPL);
                            }
                            else if (solute_bc_Amanzi_type == "BC: Zero Gradient"
                                     || solute_bc_Amanzi_type == "BC: Outflow") 
                            {
                                convert_solute_BCOutflow(solute_bc,fPL);
                            }
                            else if (solute_bc_Amanzi_type == "BC: No Flow") 
                            {
                                convert_solute_BCNoflow(solute_bc,fPL);
                            }
                            else
                            {
                                std::cerr << "Unsupported Solute BC function: \"" 
                                          << solute_bc_Amanzi_type << "\"" << std::endl;
                                throw std::exception();
                            }


                            ParameterList typePL;
                            fPL.set("regions",_regions);
                            fPL.setName(underscore(solute_bc_label));
                            solute_to_BClabel.insert(
                                std::pair<std::string,ParameterList>(underscore(soluteName),fPL));
                        }
                    }
                }
            }
            return solute_to_BClabel;
        }


        SolutePLMMap
        convert_solute_ics(StateDef& stateDef)
        {
            SolutePLMMap solute_to_IClabel;
            StateFuncMap& state_ics = stateDef.IC();    
            for (StateFuncMap::iterator ic_it = state_ics.begin(); ic_it!=state_ics.end(); ++ic_it)
            {
                const std::string& IClabel = ic_it->first;
                StateFunc& state_ic = ic_it->second;
                const Array<std::string>& regions = state_ic.Regions();
                Array<std::string> _regions;
                for (int i=0; i<regions.size(); ++i) {
                    _regions.push_back(regions[i]);
                }
                const std::string& Amanzi_type = state_ic.Amanzi_Type();
                //
                // Scan through all phases, comps to find solute ICs organized by solute
                //
                StateFunc::PhaseFuncMap& pfm = state_ics[IClabel].getPhaseFuncMap();
                for (StateFunc::PhaseFuncMap::iterator pit=pfm.begin(); pit!=pfm.end(); ++pit) {
                    const std::string& phaseName = pit->first;
                    PhaseFunc::CompFuncMap& cfm = pit->second.getCompFuncMap();
                    for (PhaseFunc::CompFuncMap::iterator cit=cfm.begin(); cit!=cfm.end(); ++cit) {
                        const std::string& compName = cit->first;
                        CompFunc::ICBCFuncMap& fm = cit->second.getICBCFuncMap();
                        for (CompFunc::ICBCFuncMap::iterator fit=fm.begin(); fit!=fm.end(); ++fit) {
                            const std::string& soluteName = fit->first;
                    
                            ICBCFunc& solute_ic = state_ic[phaseName][compName][soluteName];
                            const std::string& solute_ic_Amanzi_type = solute_ic.Amanzi_Type();
                            const std::string& solute_ic_label= solute_ic.Label();
                    
                            ParameterList fPL;
                            if (solute_ic_Amanzi_type == "IC: Uniform Concentration") 
                            {
                                convert_solute_ICConcentration(solute_ic,fPL);
                            }
                            else
                            {
                                std::cerr << "Unsuppoprted Solute IC function: \"" 
                                          << solute_ic_Amanzi_type << "\"" << std::endl;
                                throw std::exception();
                            }
                    
                            fPL.set("regions",_regions);
                            ParameterList icPL;
                            icPL.set(underscore(solute_ic_label),fPL);
                            icPL.setName(underscore(solute_ic_label));
                            solute_to_IClabel.insert(
                                std::pair<std::string,ParameterList>(underscore(soluteName),icPL));
                        }
                    }
                }
            }
            return solute_to_IClabel;
        }

        void
        convert_to_structured_state(const ParameterList& parameter_list, 
                                    ParameterList&       struc_list,
                                    StateDef&            stateDef,
                                    bool                 do_tracer)
        {
            ParameterList& phase_list  = struc_list.sublist("phase");
            ParameterList& comp_list   = struc_list.sublist("comp"); 
            ParameterList& solute_list = struc_list.sublist("tracer"); 
            ParameterList& press_list  = struc_list.sublist("press");
    
            typedef StateDef::PhaseCompMap PhaseCompMap;
            typedef StateDef::CompMap  CompMap;
    
            // FIXME: Flattens the hierarchy, as expected for PMAMR
            Array<std::string> arrayphase;
            Array<std::string> arraysolute;  
            Array<double> arraydensity;  
            Array<double> arrayviscosity;  
            Array<double> arraydiffusivity;  // FIXME: No in current spec
    
            const PhaseCompMap phase_map = stateDef.getPhaseCompMap();
            if (phase_map.size() != 1) {
                std::cerr << "Only single-phase currently supported" << std::endl;
                throw std::exception();
            }
    
            StateDef::Phases& phases = stateDef.getPhases();
            for (StateDef::Phases::const_iterator pit = phases.begin(); pit!=phases.end(); ++pit) 
            {
                const std::string& phaseLabel = pit->first;                
                std::string _phaseLabel = underscore(phaseLabel);
                PHASE& phase = stateDef.getPhases()[phaseLabel];

                arrayphase.push_back(phaseLabel);
                arraydensity.push_back(phase.Density());
                arrayviscosity.push_back(phase.Viscosity());
                arraydiffusivity.push_back(phase.Diffusivity());
        
                Array<std::string> arraycomp;  
                const CompMap comp_map = stateDef[phaseLabel];
                for (CompMap::const_iterator cit = comp_map.begin(); cit!=comp_map.end(); ++cit) 
                {
                    const std::string& compLabel = cit->first;
                    arraycomp.push_back(compLabel);

                    const Array<std::string>& soluteNames = cit->second.getTracerArray();
                    for (int i=0; i<soluteNames.size(); ++i)
                    {
                        arraysolute.push_back(soluteNames[i]);
                    }
                }

                ParameterList phasePLtr;
                phasePLtr.set("comps",arraycomp);
                phasePLtr.set("density",arraydensity);
                phasePLtr.set("viscosity",arrayviscosity);
                phasePLtr.set("diffusivity",arraydiffusivity);
                phase_list.set(_phaseLabel,phasePLtr);
            }
    
            phase_list.set("phases",arrayphase);
            solute_list.set("tracers",arraysolute);

            // Convert ICs, BCs for phase/comp
            convert_ics(parameter_list,struc_list,stateDef);    
            convert_bcs(parameter_list,struc_list,stateDef);    

            if (do_tracer) 
            {
                SolutePLMMap solute_to_ictype = convert_solute_ics(stateDef);
                SolutePLMMap solute_to_bctype = convert_solute_bcs(stateDef);

                typedef SolutePLMMap::const_iterator SPLit;
                SPLit it;
                std::pair<SPLit,SPLit> retIC, retBC;

                for (int i=0; i<arraysolute.size(); ++i) {
                    const std::string& soluteName = arraysolute[i];

                    ParameterList tmp;
                    Array<std::string> icLabels, bcLabels;
                    retIC = solute_to_ictype.equal_range(soluteName);
                    for (it=retIC.first; it!=retIC.second; ++it) {
                        const ParameterList& pl=it->second;
                        for (ParameterList::ConstIterator pit=pl.begin(); pit!=pl.end(); ++pit) {
                            const std::string& name = pl.name(pit);
                            tmp.setEntry(name,pl.getEntry(name));
                        }
                        icLabels.push_back(it->second.name());
                    }
     
                    retBC = solute_to_bctype.equal_range(soluteName);
                    for (it=retBC.first; it!=retBC.second; ++it) {
                        tmp.set(it->second.name(),it->second);
                        bcLabels.push_back(it->second.name());
                    }
                    
                    Array<std::string> regions;
                    tmp.set<Array<std::string> >("regions",regions);
                    tmp.set<Array<std::string> >("tinits",icLabels);
                    tmp.set<Array<std::string> >("tbcs",bcLabels);

                    // 
                    // FIXME: Solute groups not yet in spec, default set here
                    std::string group_name = "Total";
                    tmp.set<std::string>("group",group_name);

                    solute_list.set(soluteName,tmp);
            
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

            // Create list of available field quantities.  All these must be recognized
            //  by name inside AmrLevel::derive
            user_derive_list.push_back(underscore("Material ID"));
            user_derive_list.push_back(underscore("Capillary Pressure"));
            user_derive_list.push_back(underscore("Volumetric Water Content"));
            user_derive_list.push_back(underscore("Porosity"));
            user_derive_list.push_back(underscore("Aqueous Saturation"));
            user_derive_list.push_back(underscore("Aqueous Pressure"));
            amr_list.set<Array<std::string> >("user_derive_list",user_derive_list);

            const ParameterList& rlist = parameter_list.sublist("Output");
      
            // time macros
            std::set<std::string> time_macros;
            const ParameterList& tlist = rlist.sublist("Time Macros");
            ParameterList tmPL;
            for (ParameterList::ConstIterator i=tlist.begin(); i!=tlist.end(); ++i) {
                std::string label = underscore(tlist.name(i));
                const ParameterList& rslist = tlist.sublist(label);
                Array<double> times, vals;
                
                ParameterList tPL;
                for (ParameterList::ConstIterator ii=rslist.begin(); ii!=rslist.end(); ++ii) {
                    const std::string& name = rslist.name(ii);
                    
                    if (name == "Times") {
                        times = rslist.get<Array<double> >("Times");
                    }
                    else if (name == "Start_Period_Stop") {
                        vals = rslist.get<Array<double> >("Start_Period_Stop");
                    }
                    else {
                        std::cerr << "Invalid parameter in time macro  \"" << name
                                  << "\""  << std::endl;
                        throw std::exception();
                    }
                    
                    std::string type;
                    if (times.size()) {
                        if (vals.size()) {
                            std::cerr << "Cannot specify both time values and periods for time macro  \"" << name
                                      << "\""  << std::endl;
                            throw std::exception();
                        }
                        type = "times";
                        tPL.set<Array<double> >("times",times);
                    }
                    else {
                        if (vals.size()!=3) {
                            std::cerr << "Incorrect number of values in time macro:  \"" << name
                                      << "\""  << std::endl;
                            throw std::exception();
                        }
                        type = "period";
                        tPL.set<double>("start",vals[0]);
                        tPL.set<double>("period",vals[1]);
                        tPL.set<double>("stop",vals[2]);
                    }
                    tPL.set<std::string>("type",type);
                    tmPL.set(label,tPL);
                    time_macros.insert(label);
                }
            }
            amr_list.set("time_macro",tmPL);

            Array<std::string> tma(time_macros.size());
            int cnt = 0;
            for (std::set<std::string>::const_iterator it=time_macros.begin(); it!=time_macros.end(); ++it) {
                tma[cnt++] = *it;
            }
            amr_list.set<Array<std::string> >("time_macros",tma);
            


            // cycle macros
            std::set<std::string> cycle_macros;
            const ParameterList& clist = rlist.sublist("Cycle Macros");
            ParameterList cmPL;
            for (ParameterList::ConstIterator i=clist.begin(); i!=clist.end(); ++i) {
                std::string label = underscore(clist.name(i));
                const ParameterList& rslist = clist.sublist(clist.name(i));
                Array<int> cycles, vals;
                
                ParameterList cPL;
                for (ParameterList::ConstIterator ii=rslist.begin(); ii!=rslist.end(); ++ii) {
                    const std::string& name = rslist.name(ii);
                    
                    if (name == "Cycles") {
                        cycles = rslist.get<Array<int> >(name);
                    }
                    else if (name == "Start_Period_Stop") {
                        vals = rslist.get<Array<int> >(name);
                    }
                    else {
                        std::cerr << "Invalid parameter in cycle macro  \"" << name
                                  << "\""  << std::endl;
                        throw std::exception();
                    }
                    
                    if (cycles.size()) {
                        if (vals.size()) {
                            std::cerr << "Cannot specify both cycle values and periods for cycle macro  \"" << name
                                      << "\""  << std::endl;
                            throw std::exception();
                        }
                        cPL.set<std::string>("type","cycles");
                        cPL.set<Array<int> >("cycles",cycles);
                    }
                    else {
                        if (vals.size()!=3) {
                            std::cerr << "Incorrect number of values in cycle macro:  \"" << name
                                      << "\""  << std::endl;
                            throw std::exception();
                        }
                        cPL.set<std::string>("type","period");
                        cPL.set<int>("start",vals[0]);
                        cPL.set<int>("period",vals[1]);
                        cPL.set<int>("stop",vals[2]);
                    }
                }
                cmPL.set(label,cPL);
                cycle_macros.insert(label);
            }
            amr_list.set("cycle_macro",cmPL);

            Array<std::string> cma(cycle_macros.size());
            int ccnt = 0;
            for (std::set<std::string>::const_iterator it=cycle_macros.begin(); it!=cycle_macros.end(); ++it) {
                cma[ccnt++] = *it;
            }
            amr_list.set<Array<std::string> >("cycle_macros",cma);
            


            // vis data
            const ParameterList& vlist = rlist.sublist("Visualization Data");
            amr_list.set("plot_file",vlist.get<std::string>("File Name Base"));
            Array<std::string> visNames;

            if (vlist.isParameter("Variables")) {
                visNames = vlist.get<Array<std::string> >("Variables");
                for (int i=0; i<visNames.size(); ++i) {
                    std::string _visName = underscore(visNames[i]);
                    bool found = false;
                    for (int j=0; j<user_derive_list.size(); ++j) {
                        if (_visName == user_derive_list[j]) {
                            found = true;
                        }
                    }
                    if (!found) {
                        std::cerr << "Invalid variable (\"" << visNames[i]
                                  << "\") in \"Visualization Data\" -> \"Variables\"" << std::endl;
                        throw std::exception();
                    }
                    visNames[i] = underscore(visNames[i]);
                }
            }
            else {
                std::cerr << "Must select variables to put into visualization files from the list: (\"";
                for (int j=0; j<user_derive_list.size(); ++j) {
                    std::cout << "\""<< AMR_to_Amanzi_label_map[user_derive_list[j]] << "\" ";
                }
                std::cout << std::endl;
                throw std::exception();
            }
            amr_list.set<Array<std::string> >("derive_plot_vars",visNames);

            Array<std::string> vis_cMacroNames, vis_tMacroNames;
            if (vlist.isParameter("Cycle Macros")) {
                const Array<std::string>& vcMacros = vlist.get<Array<std::string> >("Cycle Macros");
                vis_cMacroNames.resize(vcMacros.size());
                for (int j=0; j<vcMacros.size(); ++j) {
                    std::string label = underscore(vcMacros[j]);
                    if (cycle_macros.find(label) != cycle_macros.end()) {
                        vis_cMacroNames[j] = label;
                    }
                    else {
                        std::cerr << "Unrecognized cycle macro in Visualization Data: \""
                                  << vcMacros[j] << "\"" << std::endl;

                        for (std::set<std::string>::const_iterator it=cycle_macros.begin(); it!=cycle_macros.end(); ++it) {
                            std::cout << *it << " " << std::endl;
                        }

                        throw std::exception();
                    }
                }
            }
            if (vlist.isParameter("Time Macros")) {
                const Array<std::string>& vtMacros = vlist.get<Array<std::string> >("Time Macros");
                vis_tMacroNames.resize(vtMacros.size());
                for (int j=0; j<vtMacros.size(); ++j) {
                    std::string label = underscore(vtMacros[j]);
                    if (time_macros.find(label) != time_macros.end()) {
                        vis_tMacroNames[j] = label;
                    }
                    else {
                        std::cerr << "Unrecognized time macro in Visualization Data: \""
                                  << vtMacros[j] << "\"" << std::endl;

                        for (std::set<std::string>::const_iterator it=time_macros.begin(); it!=time_macros.end(); ++it) {
                            std::cout << *it << " " << std::endl;
                        }

                        throw std::exception();
                    }
                }
            }

            amr_list.set<Array<std::string> >("vis_cycle_macros",vis_cMacroNames);
            amr_list.set<Array<std::string> >("vis_time_macros",vis_tMacroNames);

            //amr_list.set<std::string>("plot_vars",""); // Shut off, per spec



            // check point data
            const ParameterList& chlist = rlist.sublist("Checkpoint Data");
            amr_list.set("check_file",chlist.get<std::string>("File Name Base"));
            Array<std::string> chkNames;

            Array<std::string> chk_cMacroNames;
            if (chlist.isParameter("Cycle Macros")) {
                const Array<std::string>& ccMacros = chlist.get<Array<std::string> >("Cycle Macros");
                chk_cMacroNames.resize(ccMacros.size());
                for (int j=0; j<ccMacros.size(); ++j) {
                    std::string label = underscore(ccMacros[j]);
                    if (cycle_macros.find(label) != cycle_macros.end()) {
                        chk_cMacroNames[j] = label;
                    }
                    else {
                        std::cerr << "Unrecognized cycle macro in Checkpoint Data: \""
                                  << ccMacros[j] << "\"" << std::endl;
                        throw std::exception();
                    }
                }
            }
            else {
                std::cerr << "Must provide \"Cycle Macros\" in Checkpoint Data" << std::endl;
                throw std::exception();
            }
            amr_list.set<Array<std::string> >("chk_cycle_macros",chk_cMacroNames);


        
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
                    if (functional == "Observation Data: Integral") {
                        sublist.set("obs_type","integral");
                    }
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

                    const std::string& timeMacro = rslist.get<std::string>("Time Macro");
                    if (time_macros.find(timeMacro) == time_macros.end()) {
                        std::cerr << "Unrecognized time macro: \"" << timeMacro
                                  << "\" for observation data: \"" << label << "\"" << std::endl;
                        throw std::exception();                        
                    }
                    else {
                        sublist.set("time_macro",timeMacro);
                    }
	  
                    Array<std::string> arrayvariables;
                    Array<std::string> variables = rslist.get<Array<std::string> >("Variables");
                    if (variables.size()!=1) {
                        std::cerr << "Currently must provide a single Variable per observation" << std::endl;
                            throw std::exception();                        
                    }
                    for (int j=0; j<variables.size(); ++j) {
                        std::string _variable = underscore(variables[j]);
                        bool found = false;
                        for (int k=0; k<user_derive_list.size() && !found; ++k) {
                            if (_variable == user_derive_list[k]) {
                                found = true;
                            }
                        }
                        if (found) {
                            arrayvariables.push_back(_variable);
                        } else 
                        {
                            std::cerr << variables[j] 
                                      << " is not a valid derive variable name. Must be on of: ";
                            for (int k=0; k<user_derive_list.size() && !found; ++k) {
                                std::cerr << "\"" << user_derive_list[k] << "\" ";
                            }
                            std::cerr << std::endl;
                            throw std::exception();
                        }

                    }

                    sublist.set<std::string>("field",arrayvariables[0]);

                    obs_list.set(_label,sublist);
                    arrayobs.push_back(_label);
                }
            }
            obs_list.set("observation",arrayobs);
        }

        //
        // convert parameterlist to format for structured code
        //
        ParameterList
        convert_to_structured(const ParameterList& parameter_list)
        {
            ParameterList struc_list = setup_structured();
            bool do_tracer;
            //
            // determine spatial dimension of the problem
            // 
            ndim = parameter_list.sublist("Domain").get<int>("Spatial Dimension");
            //
            // Mesh
            //
            convert_to_structured_mesh(parameter_list,struc_list);
            //
            // Execution control
            //
            convert_to_structured_control(parameter_list,struc_list, do_tracer);
            //
            // Regions
            //
            convert_to_structured_region(parameter_list, struc_list);
            //
            // Materials
            //
            convert_to_structured_material(parameter_list, struc_list);
            //
            // State
            //
            StateDef stateDef(parameter_list);
            convert_to_structured_state(parameter_list, struc_list, stateDef, do_tracer);
            //
            // Output
            // 
            convert_to_structured_output(parameter_list, struc_list);
            return struc_list;
        }

    }
}


