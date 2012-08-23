#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include <InputParser_Structured.H>
#include <PMAMR_Labels.H>

#include <BoxLib.H>

using Teuchos::Array;
using Teuchos::ParameterList;
using Teuchos::ParameterEntry;

namespace Amanzi {
    namespace AmanziInput {

        std::map<std::string,SolidChem::SorptionIsothermData> SolidChem::sorption_isotherms; // One for all materials, indexed on solute name

        void MyAbort(const std::string& m) {
            if (Teuchos::MPISession::getRank() == 0) {
                std::cerr << m << std::endl;
                throw std::exception();
            }
        }

        void MyAbort(const Array<std::string>& m) {
            if (Teuchos::MPISession::getRank() == 0) {
                for (int i=0; i<m.size(); ++i) {
                    std::cerr << m[i] << " ";
                }
                std::cerr << std::endl;
                throw std::exception();
            }
        }

        static std::map<std::string,std::string> AMR_to_Amanzi_label_map;

        double atmToMKS = 101325;

        std::string underscore(const std::string& instring)
        {
            std::string s = instring;
            std::replace(s.begin(),s.end(),' ','_');
            AMR_to_Amanzi_label_map[s] = instring;
            return s;
        }

        Array<std::string> underscore(const Array<std::string>& instrings)
        {
            Array<std::string> ss(instrings.size());
            for (int i=0; i<instrings.size(); ++i) {
                ss[i] = instrings[i];
                std::replace(ss[i].begin(),ss[i].end(),' ','_');
                AMR_to_Amanzi_label_map[ss[i]] = instrings[i];
            }
            return ss;
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

            if (domlo.size()<ndim  || domhi.size()<ndim) {
                MyAbort("Domain size nonsensical");
            }
            double max_size = domhi[0]-domlo[0];
            for (int i=0; i<ndim; ++i) {
                max_size=std::max(max_size,domhi[i]-domlo[i]);
            }
            // FIXME: This multiplier should be input
            geometry_eps = 1.e-6*max_size;

            const ParameterList& eclist = parameter_list.sublist("Execution Control");

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

        std::pair<std::string,bool> one_picked(const std::map<std::string,bool>& b)
        {
            std::pair<std::string,bool> res("",false);
            if (b.size()) {
                for (std::map<std::string,bool>::const_iterator it=b.begin(); it!=b.end(); ++it)
                {
                    if (res.second) {
                        if (it->second) return std::pair<std::string,bool>("",false);
                    }
                    else {
                        res = *it;
                    }
                }
            }
            return res;
        }


        void process_expert_options(const ParameterList& pl,
                                    ParameterList&       out_pl)
        {
	  const std::string expert_str = "Expert Settings";
	  if (pl.isSublist(expert_str)) {
	    const ParameterList& expert_list = pl.sublist(expert_str);
	    for (ParameterList::ConstIterator it=expert_list.begin(); it!=expert_list.end(); ++it) {
	      const std::string& name = expert_list.name(it);
	      out_pl.setEntry(name,expert_list.getEntry(name));
	    }
	  }
        }

        //
        // convert execution control to structured format
        //
        void
        convert_to_structured_control(const ParameterList& parameter_list, 
                                      ParameterList&       struc_out_list,
                                      bool&                do_tracer_transport,
                                      bool&                do_chem)
        {
            std::string ec_str = "Execution Control";
            std::string amr_str = "Adaptive Mesh Refinement Control";
            std::string prob_str = "Basic Algorithm Control";
            std::string it_str = "Iterative Linear Solver Control";
            std::string cg_str = "Conjugate Gradient Algorithm";
            std::string mg_str = "Multigrid Algorithm";
            std::string mac_str = "Pressure Discretization Control";
            std::string diffuse_str = "Diffusion Discretization Control";

            std::string chem_str = "Chemistry";
            const ParameterList& ec_list = parameter_list.sublist(ec_str);

            Array<std::string> reqL, reqP;
            std::string flow_str = "Flow Model"; reqP.push_back(flow_str);
            std::string trans_str = "Transport Model";reqP.push_back(trans_str);
            std::string chem_mod_str = "Chemistry Model"; reqP.push_back(chem_mod_str);
            std::string tim_str = "Time Integration Mode"; reqL.push_back(tim_str);
            std::string v_str = "Verbosity";
            std::string num_str = "Numerical Control Parameters";
            PLoptions ECopt(ec_list,reqL,reqP,false,false); 

            ParameterList& amr_out_list     = struc_out_list.sublist("amr");
            ParameterList& prob_out_list    = struc_out_list.sublist("prob");
            ParameterList& cg_out_list      = struc_out_list.sublist("cg");
            ParameterList& mg_out_list      = struc_out_list.sublist("mg");
            ParameterList& mac_out_list     = struc_out_list.sublist("mac");
            ParameterList& diffuse_out_list = struc_out_list.sublist("diffuse");

            ParameterList& chem_out_list    = prob_out_list.sublist("amanzi");

            bool echo_inputs = false;
            std::string echo_str = "Echo Inputs";
            if (parameter_list.isParameter(echo_str)) {
                echo_inputs = parameter_list.get<bool>(echo_str);
            }
            struc_out_list.set<bool>("echo_inputs",echo_inputs);

            //
            // Set flow model
            //
            std::string model_name;
            std::string flow_mode = ec_list.get<std::string>(flow_str);
            if (flow_mode == "Off") {
                MyAbort("\"" + flow_str + "\" = \"" + flow_mode + "\" not supported");
                prob_out_list.set("do_simple",2);
            }
            else if (flow_mode == "Richards") {
                model_name = "richard";
                prob_out_list.set("have_capillary",1);
                prob_out_list.set("cfl",-1);
            }
            else if (flow_mode == "Steady State Saturated") {
                model_name = "steady-saturated";
                prob_out_list.set("have_capillary",0);
                prob_out_list.set("cfl",-1);
            }
            else if (flow_mode == "Single-phase") {
                model_name = "single-phase";
                prob_out_list.set("do_simple",1);
            }
            else if (flow_mode == "Multi-phase") {
                model_name = "two-phase";
                prob_out_list.set("cfl",0.75);
            }
            else {
                MyAbort("\"" + flow_str + "\" = \"" + flow_mode + "\" not supported");
            }
            prob_out_list.set("model_name",model_name);

            //
            // Set transport model
            //
            std::string transport_mode = ec_list.get<std::string>(trans_str);
            do_tracer_transport = (transport_mode == "Off"  ?  0  :  1);
            prob_out_list.set<int>("do_tracer_transport",do_tracer_transport);

            //
            // Set chemistry model
            //
            std::string chem_mode = ec_list.get<std::string>(chem_mod_str);
            if (chem_mode == "Off") {
                prob_out_list.set("do_chem",0);
                do_chem = false;
            }
            else if (chem_mode == "On") {
                prob_out_list.set("do_chem",1);
                do_chem = true;
                const ParameterList& chem_list = parameter_list.sublist(chem_str);
                reqP.clear(); reqL.clear();
                PLoptions CHopt(chem_list,reqL,reqP,true,false); 
                const Array<std::string>& CHoptP = CHopt.OptParms();
                for (int i=0; i<CHoptP.size(); ++i) {
                    const std::string& name = CHoptP[i];
                    std::string _name = underscore(name);
                    if (name=="Thermodynamic Database Format") {
                        chem_out_list.set(_name,chem_list.get<std::string>(name));
                    }
                    else if (name=="Thermodynamic Database File") {
                        chem_out_list.set("chem_database_file",chem_list.get<std::string>(name));
                    }
                    else if (name=="Verbosity") {
                        chem_out_list.set("verbose_chemistry_init",chem_list.get<std::string>(name));
                    }
                    else if (name=="Activity Model") {
                        chem_out_list.set(_name,chem_list.get<std::string>(name));
                    }
                    else if (name=="Tolerance") {
                        chem_out_list.set(_name,chem_list.get<double>(name));
                    }
                    else if (name=="Maximum Newton Iterations") {
                        chem_out_list.set(_name,chem_list.get<int>(name));
                    }
                    else if (name=="Output File Name") {
                        chem_out_list.set(_name,chem_list.get<std::string>(name));
                    }
                    else if (name=="Use Standard Out") {
                        chem_out_list.set(_name,chem_list.get<bool>(name));
                    }
                    else if (name=="Auxiliary Data") {
                        chem_out_list.set(_name,chem_list.get<Array<std::string> >(name));
                    }
                }
            }
            else {
                MyAbort("Chemistry Model \"" + chem_mode + "\" not yet supported" );
            }

            //
            // Set time evolution mode
            //
            std::string steady_str = "Steady";
            std::string transient_str = "Transient";
            std::string init_to_steady_str = "Initialize To Steady";
            const ParameterList& t_list = ec_list.sublist(tim_str);
	    if (model_name == "single_phase") {
	      prob_out_list.set("do_simple",2);
            }
            
	    if (t_list.isSublist(transient_str))
            {
                const ParameterList& tran_list = t_list.sublist(transient_str);
                reqP.clear(); reqL.clear();
                std::string Start_str = "Start"; reqP.push_back(Start_str);
                std::string End_str = "End"; reqP.push_back(End_str);
                std::string Init_Time_Step_str = "Initial Time Step";
                std::string Time_Step_Grow_Max_str = "Maximum Time Step Grow";
                std::string Time_Step_Shrink_Max_str = "Maximum Time Step Shrink";
                std::string Init_Time_Step_Mult_str = "Initial Time Step Multiplier";
                std::string Max_Time_Step_Size_str = "Maximum Time Step Size";
                std::string Max_Step_str = "Maximum Cycle Number";

                PLoptions Topt(tran_list,reqL,reqP,true,false); 
                const Array<std::string> ToptP = Topt.OptParms();
                struc_out_list.set<double>("strt_time", tran_list.get<double>(Start_str));
                struc_out_list.set<double>("stop_time", tran_list.get<double>(End_str));
                double dt_init = -1;
                double dt_init_mult = -1;
                double dt_grow_max = -1;
                double dt_shrink_max = -1;
                int step_max = -1;
                double dt_max = -1;

		struc_out_list.set<std::string>("execution_mode", "transient");

                for (int i=0; i<ToptP.size(); ++i) {
                    if (ToptP[i] == Init_Time_Step_str) {
                        dt_init = tran_list.get<double>(Init_Time_Step_str);
                    }
                    else if (ToptP[i] == Time_Step_Shrink_Max_str) {
                        dt_shrink_max = tran_list.get<double>(Time_Step_Shrink_Max_str);
                    }
                    else if (ToptP[i] == Time_Step_Grow_Max_str) {
                        dt_grow_max = tran_list.get<double>(Time_Step_Grow_Max_str);
                    }
                    else if (ToptP[i] == Init_Time_Step_Mult_str) {
                        dt_init_mult = tran_list.get<double>(Init_Time_Step_Mult_str);
                    }
                    else if (ToptP[i] == Max_Time_Step_Size_str) {
                        dt_max = tran_list.get<double>(Max_Time_Step_Size_str);
                    }
                    else if (ToptP[i] == Max_Step_str) {
                        step_max = tran_list.get<int>(Max_Step_str);
                    }
                    else {
                        MyAbort("Unrecognized option under \""+transient_str+"\": \""+ToptP[i]+"\"" );
                    }
                }

                if (dt_init > 0) {
                    prob_out_list.set<double>("dt_init", dt_init);
                }

                if (dt_init_mult > 0) {
                    prob_out_list.set<double>("init_shrink", dt_init_mult);
                } 

                if (dt_shrink_max > 0) {                        
                    prob_out_list.set<double>("dt_shrink_max", dt_shrink_max);
                }

                if (dt_grow_max > 0) {                        
                    prob_out_list.set<double>("dt_grow_max", dt_grow_max);
                }

                if (step_max>=0) {
                    struc_out_list.set<int>("max_step", step_max);
                }
                
                if (dt_max > 0) {
		  prob_out_list.set<double>("max_dt", dt_max);
		}
            }
            else if (t_list.isSublist(init_to_steady_str))
            {
                const ParameterList& tran_list = t_list.sublist(init_to_steady_str);

                reqP.clear(); reqL.clear();
                std::map<std::string,double> optPd;
                std::map<std::string,int> optPi;
                std::map<std::string,bool> optPb;
                std::string Start_str = "Start"; reqP.push_back(Start_str);
                std::string Switch_str = "Switch"; reqP.push_back(Switch_str);
                std::string End_str = "End"; reqP.push_back(End_str);
                std::string Steady_Init_Time_Step_str = "Steady Initial Time Step"; 
		reqP.push_back(Steady_Init_Time_Step_str);
                std::string Transient_Init_Time_Step_str = "Transient Initial Time Step"; 
		reqP.push_back(Transient_Init_Time_Step_str);

                // Set defaults for optional parameters
                std::string Init_Time_Step_Mult_str = "Initial Time Step Multiplier"; optPd[Init_Time_Step_Mult_str] = -1;
                std::string Max_Step_str = "Maximum Cycle Number"; optPi[Max_Step_str] = 100000;
		std::string Use_Picard_str = "Use Picard"; optPb[Use_Picard_str] = true; // This will be ignored
                std::string Steady_Max_Time_Step_Size_str = "Steady Maximum Time Step Size";
                optPd[Steady_Max_Time_Step_Size_str] = -1; // <0 means inactive
                std::string Transient_Max_Time_Step_Size_str = "Transient Maximum Time Step Size";
                optPd[Transient_Max_Time_Step_Size_str] = -1; // <0 means inactive
                std::string Time_Step_Grow_Max_str = "Maximum Time Step Grow";
                optPd[Time_Step_Grow_Max_str] = -1; // <0 means inactive
                std::string Time_Step_Shrink_Max_str = "Maximum Time Step Shrink";
                optPd[Time_Step_Shrink_Max_str] = -1; // <0 means inactive


                struc_out_list.set<double>("strt_time", tran_list.get<double>(Start_str));
                struc_out_list.set<double>("stop_time", tran_list.get<double>(End_str));
                struc_out_list.set<double>("switch_time", tran_list.get<double>(Switch_str));
		prob_out_list.set<double>(underscore("steady_max_pseudo_time"),tran_list.get<double>(Switch_str));
		prob_out_list.set<double>(underscore("steady_init_time_step"),tran_list.get<double>(Steady_Init_Time_Step_str));
		prob_out_list.set<double>(underscore("dt_init"),tran_list.get<double>(Transient_Init_Time_Step_str));

		struc_out_list.set<std::string>("execution_mode", "init_to_steady");

                // Extract optional parameters
                PLoptions Topt(tran_list,reqL,reqP,true,false); 
                const Array<std::string> ToptP = Topt.OptParms();
                for (int i=0; i<ToptP.size(); ++i) {
                    if (optPd.find(ToptP[i]) != optPd.end()) {
                        optPd[ToptP[i]] = tran_list.get<double>(ToptP[i]); // replace default value
                    }
                    else if (optPi.find(ToptP[i]) != optPi.end()) {
                        optPi[ToptP[i]] = tran_list.get<int>(ToptP[i]); // replace default value
                    }
                    else if (optPb.find(ToptP[i]) != optPb.end()) {
                        optPb[ToptP[i]] = tran_list.get<bool>(ToptP[i]); // replace default value
                    }
                    else {
                        MyAbort("Unrecognized option under \""+init_to_steady_str+"\": \""+ToptP[i]+"\"" );
                    }
                }

                for (std::map<std::string,double>::const_iterator it=optPd.begin(); it!=optPd.end(); ++it) {
                    if (it->first==Time_Step_Grow_Max_str) {
                        prob_out_list.set<double>("dt_grow_max", it->second);
                    } else if (it->first==Time_Step_Shrink_Max_str) {
                        prob_out_list.set<double>("dt_shrink_max", it->second);
                    } else if (it->first==Steady_Max_Time_Step_Size_str) {
                        prob_out_list.set<double>("steady_richard_max_dt", it->second);
                    } else if (it->first==Transient_Max_Time_Step_Size_str) {
                        prob_out_list.set<double>("transient_richard_max_dt", it->second);
                    } else {
                        prob_out_list.set<double>(underscore(it->first), it->second);
                    }
                }
                for (std::map<std::string,int>::const_iterator it=optPi.begin(); it!=optPi.end(); ++it) {
		  if (it->first==Max_Step_str) {
		    struc_out_list.set<int>("max_step", it->second);
		  } else {
                    prob_out_list.set<int>(underscore(it->first), it->second);
		  }
                }
                for (std::map<std::string,bool>::const_iterator it=optPb.begin(); it!=optPb.end(); ++it) {
		  prob_out_list.set<bool>(underscore(it->first), it->second);
                }
            }
	    else if (t_list.isSublist(steady_str)) {
                const ParameterList& tran_list = t_list.sublist(transient_str);
                reqP.clear(); reqL.clear();
                PLoptions Topt(tran_list,reqL,reqP,true,true); 

		struc_out_list.set<std::string>("execution_mode", "steady");
		prob_out_list.set<double>("stop_time",0);
		prob_out_list.set<double>("start_time",0);
	    }
            else {
                std::cout << t_list << std::endl;
                MyAbort("No recognizable value for \"" + tim_str + "\"");
            }

            // Deal with optional settings
            const Array<std::string> optL = ECopt.OptLists();
            const Array<std::string> optP = ECopt.OptParms();

            //
            // Basic Algorithm Defaults
            //
            double visc_abs_tol = 1.e-12; prob_out_list.set<double>("visc_abs_tol",visc_abs_tol);
            double visc_tol = 1.e-6; prob_out_list.set<double>("visc_tol",visc_tol);

            //
            // Verbosity Default
            //
            std::string v_val = "Medium";

            //
            // AMR gridding control defaults
            //
            int num_levels = 1;
            int max_level = num_levels-1;
            bool do_amr_subcycling = false;                            
            int ref_ratio_DEF = 2;
            Array<int> ref_ratio(max_level,ref_ratio_DEF);
            int regrid_int_DEF = 2;
            Array<int> regrid_int(num_levels,regrid_int_DEF);
            int blocking_factor_DEF = 8;
            Array<int> blocking_factor(num_levels,blocking_factor_DEF);
            int n_err_buf_DEF = 1;
            Array<int> n_err_buf(max_level,n_err_buf_DEF);
            int max_grid_DEF = (ndim==2 ? 128  :  32);
            Array<int> max_grid(num_levels,max_grid_DEF);

            // 
            // Optional parameters
            //
            std::string restart_file_str = "Restart from Checkpoint Data File";

            for (int i=0; i<optP.size(); ++i) {
                if (optP[i] == v_str) {
                    //
                    // Verbosity level
                    //
                    v_val = ec_list.get<std::string>(v_str);
                }
                else if (optP[i] == restart_file_str) {
                    //
                    // Restart from checkpoint?
                    //
                    amr_out_list.set<std::string>("restart",ec_list.get<std::string>(restart_file_str));
                }
                else {
                    MyAbort("Unrecognized optional parameter to \"" + ec_str + "\": \"" + optP[i] + "\"");
                }
            }

            //
            // Verbosity implementation
            //
            int prob_v, mg_v, cg_v, amr_v, diffuse_v;
            if (v_val == "None") {
                prob_v = 0; mg_v = 0; cg_v = 0; amr_v = 0; diffuse_v = 0;
            }
            else if (v_val == "Low") {
                prob_v = 1; mg_v = 0; cg_v = 0; amr_v = 1;  diffuse_v = 0;
            }
            else if (v_val == "Medium") {
                prob_v = 1; mg_v = 0; cg_v = 0; amr_v = 2;  diffuse_v = 1;
            }
            else if (v_val == "High") {
                prob_v = 2; mg_v = 1; cg_v = 1; amr_v = 3;  diffuse_v = 1;
            }
            else if (v_val == "Extreme") {
                prob_v = 3; mg_v = 2; cg_v = 2; amr_v = 3;  diffuse_v = 1;
            }

            // 
            // Optional lists
            //
            for (int i=0; i<optL.size(); ++i)
            {
                if (optL[i] == num_str)
                {
                    const ParameterList& num_list = ec_list.sublist(num_str);
                    Array<std::string> nL, nP;
                    PLoptions NUMopt(num_list,nL,nP,false,true); 
                    const Array<std::string> NUMoptL = NUMopt.OptLists();

                    for (int j=0; j<NUMoptL.size(); ++j)
                    {
                        if (NUMoptL[j] == amr_str)
                        {
                            //
                            // AMR Options
                            //
                            const ParameterList& amr_list = num_list.sublist(amr_str);
                            std::string num_level_str = "Number Of AMR Levels";
                            if (amr_list.isParameter(num_level_str)) {
                                num_levels = amr_list.get<int>(num_level_str);
                            }
                            if (num_levels < 1) {
                                MyAbort("Must have at least 1 AMR level");
                            }
                            max_level = num_levels - 1;

                            std::string ref_ratio_str = "Refinement Ratio";
                            ref_ratio.resize(max_level,2);
                            if (amr_list.isParameter(ref_ratio_str)) {
                                ref_ratio = amr_list.get<Array<int> >(ref_ratio_str);
                            }
                            if (ref_ratio.size() < max_level) {
                                MyAbort("Must provide a refinement ratio for each refined level");
                            }

                            for (int k=0; k<max_level; ++k) {
                                if (ref_ratio[k] != 2 && ref_ratio[k]!=4) {
                                    MyAbort("\"Refinement Ratio\" values must be 2 or 4");
                                }
                            }

                            std::string regrid_int_str = "Regrid Interval";
                            regrid_int.resize(max_level,2);
                            if (amr_list.isParameter(regrid_int_str)) {
                                regrid_int = amr_list.get<Array<int> >(regrid_int_str);
                            }
                            if (max_level > 0) {
                                if (do_amr_subcycling) {
                                    if (regrid_int.size() < max_level) {
                                        MyAbort("Must provide a regridding interval for each refined level");
                                    }
                                }
                                else {
                                    if (regrid_int.size() != 1) {
                                        MyAbort("Subcycling is disabled, only a single regridding interval is supported");
                                    }                                
                                }

                                for (int k=0; k<regrid_int.size(); ++k) {
                                    if (regrid_int[k] <= 0) {
                                        MyAbort("Each value in \"" + regrid_int_str
                                                + "\" must be values must be a postive integer");
                                    }
                                }
                            }

                            
                            std::string blocking_factor_str = "Blocking Factor";
                            blocking_factor.resize(max_level+1,blocking_factor_DEF);
                            if (amr_list.isParameter(blocking_factor_str)) {
                                blocking_factor = amr_list.get<Array<int> >(blocking_factor_str);
                            }
                            if (blocking_factor.size() < max_level+1) {
                                MyAbort("If provided, value of \"" + blocking_factor_str
                                        + "\" required for each level");
                            }
                            for (int k=0; k<blocking_factor.size(); ++k) {
                                double twoPower = std::log(blocking_factor[k])/std::log(2);
                                if (twoPower != (int)(twoPower)) {
                                    MyAbort("\"" + blocking_factor_str + "\" must be a power of two");
                                }
                            }

                            std::string n_err_buf_str = "Number Error Buffer Cells";
                            int n_err_buf_DEF = 1;
                            n_err_buf.resize(max_level+1,n_err_buf_DEF);
                            if (amr_list.isParameter(n_err_buf_str)) {
                                n_err_buf = amr_list.get<Array<int> >(n_err_buf_str);
                            }
                            if (n_err_buf.size() < max_level) {
                                MyAbort("If provided, value of \"" + n_err_buf_str
                                        + "\" required for each refined level");
                            }
                            for (int k=0; k<n_err_buf.size(); ++k) {
                                if (n_err_buf[k] < 0) {
                                    MyAbort("\"" + n_err_buf_str + "\" must be > 0");
                                }
                            }

                            
                            std::string max_grid_str = "Maximum Grid Size";
                            max_grid.resize(max_level+1,max_grid_DEF);
                            if (amr_list.isParameter(max_grid_str)) {
                                max_grid = amr_list.get<Array<int> >(max_grid_str);
                            }
                            if (max_grid.size() < max_level+1) {
                                MyAbort("If provided, value of \"" + max_grid_str
                                        + "\" required for each level");
                            }
                            for (int k=0; k<max_grid.size(); ++k) {
                                if (max_grid[k] < blocking_factor[k]) {
                                    MyAbort("\"" + max_grid_str + "\" must be > \"" + blocking_factor_str + "\"");
                                }
                            }


                            std::string refineNames_str = "Refinement Indicators";
                            if (amr_list.isParameter(refineNames_str)) {
                                const Array<std::string>& refineNames = 
                                    amr_list.get<Array<std::string> >(refineNames_str);
                                Array<std::string> names(refineNames.size());
                                for (int k=0; k<refineNames.size(); ++k) {
                                    names[k] = underscore(refineNames[k]);
                                    const ParameterList& ref_list = amr_list.sublist(refineNames[k]);

                                    std::string val_greater_str = "Value Greater";
                                    std::string val_less_str = "Value Less";
                                    std::string diff_greater_str = "Adjacent Difference Greater";
                                    std::string in_region_str = "Inside Region";

                                    bool do_greater = false;
                                    bool do_less    = false;
                                    bool do_diff    = false;
                                    bool do_region  = false;
                                    std::map<std::string,bool> pick_one;
                                    if (ref_list.isParameter(val_greater_str)) {
                                        pick_one["do_greater"] = true;
                                    }
                                    if (ref_list.isParameter(val_less_str)) {
                                        pick_one["do_less"] = true;
                                    }
                                    if (ref_list.isParameter(diff_greater_str)) {
                                        pick_one["do_diff"] = true;
                                    }
                                    if (ref_list.isParameter(in_region_str)) {
                                        pick_one["do_region"] = true;
                                    }

                                    std::pair<std::string,bool> ch = one_picked(pick_one);
                                    if (! ch.second)
                                    {
                                        MyAbort("Refinement indicator \"" + refineNames[k]
                                                + "\" must specify one condition from the list: "
                                                + "\"" + val_greater_str + "\", "
                                                + "\"" + val_less_str + "\", "
                                                + "\"" + diff_greater_str + "\", "
                                                + "\"" + in_region_str + "\"" );
                                    }
                                    
                                    std::string fieldName_str = "Field Name";
                                    std::string regName_str = "Regions";

                                    std::string ref_type = ch.first;
                                    ParameterList ref_out_list;
                                    if (ref_type != "do_region") {
                                        ref_out_list.set<std::string>
                                            ("field",underscore(ref_list.get<std::string>(fieldName_str)));

                                        if (ref_type == "do_greater") {
                                            ref_out_list.set<double>(
                                                "val_greater_than",ref_list.get<double>(val_greater_str));
                                        }
                                        else if (ref_type == "do_less") {
                                            ref_out_list.set<double>(
                                                "val_less_than",ref_list.get<double>(val_less_str));
                                        }
                                        else if (ref_type == "do_diff") {
                                            ref_out_list.set<double>(
                                                "diff_greater_than",ref_list.get<double>(diff_greater_str));
                                        }
                                    }
                                    else {
                                        ref_out_list.set<bool>("in_region","TRUE");
                                    }

                                    std::string maxLev_str = "Maximum Refinement Level";
                                    int this_max_level = max_level;
                                    if (ref_list.isParameter(maxLev_str)) {
                                        this_max_level = ref_list.get<int>(maxLev_str);
                                    }
                                    ref_out_list.set<int>("max_level",this_max_level);

                                    ref_out_list.set<Array<std::string> >(
                                        "regions",underscore(ref_list.get<Array<std::string> >(regName_str)));
                                    
                                    std::string start_str = "Start Time";
                                    std::string end_str = "End Time";
                                    if (ref_list.isParameter(start_str)) {
                                        ref_out_list.set<double>("start_time",ref_list.get<double>(start_str));
                                    }
                                    if (ref_list.isParameter(end_str)) {
                                        ref_out_list.set<double>("end_time",ref_list.get<double>(end_str));
                                    }

                                    amr_out_list.set(names[k],ref_out_list);
                                }

                                amr_out_list.set<Array<std::string> >("refinement_indicators",names);
                            }
                        }
                    }
                }
            }
                    
            Array<int> n_cell = amr_out_list.get<Array<int> >("n_cell");
            for (int i=0;i<blocking_factor.size();i++) {
                
                for (int n=0;n<ndim;n++) {
                    if (n_cell[n]%blocking_factor[i] > 0) {
                        std::stringstream buf;
                        for (int j=0;j<blocking_factor.size();j++) {
                            buf << blocking_factor[j] << " ";
                        }
                        std::string m; buf >> m;
                        MyAbort("Number of Cells must be divisible by blocking_factor = " + m);
                    }
                }
            }
                            
            amr_out_list.set<int>("max_level",max_level);
            amr_out_list.set<Array<int> >("ref_ratio",ref_ratio);
            amr_out_list.set<Array<int> >("regrid_int",regrid_int);
            amr_out_list.set<Array<int> >("blocking_factor",blocking_factor);
            amr_out_list.set<Array<int> >("n_error_buf",n_err_buf);
            amr_out_list.set<Array<int> >("max_grid_size",max_grid);
            int amr_nosub = ( do_amr_subcycling ? 0 : 1);
            amr_out_list.set<int>("nosub", amr_nosub);

            amr_out_list.set("v",amr_v);
            mg_out_list.set("v",mg_v);
            cg_out_list.set("v",cg_v);
            prob_out_list.set("v",prob_v);
            diffuse_out_list.set("v",diffuse_v);
            
            for (int i=0; i<optL.size(); ++i)
            {
                if (optL[i] == num_str)
                {
                    const ParameterList& num_list = ec_list.sublist(num_str);
                    Array<std::string> nL, nP;
                    PLoptions NUMopt(num_list,nL,nP,false,true); 
                    const Array<std::string> NUMoptL = NUMopt.OptLists();

                    //
                    // Now process expert lists to overwrite any settings directly
                    //
                    for (int j=0; j<NUMoptL.size(); ++j)
                    {
                        if (NUMoptL[j] == amr_str)
                        {
                            process_expert_options(num_list.sublist(amr_str),amr_out_list);
                        }
                        else if (NUMoptL[j] == it_str)
                        {
                            const ParameterList& it_list = num_list.sublist(it_str);
                            for (ParameterList::ConstIterator k=it_list.begin(); k!=it_list.end(); ++k) 
                            {
                                const std::string& itname = it_list.name(k);
                                if (itname == mg_str)
                                {
                                    process_expert_options(it_list.sublist(mg_str),mg_out_list);
                                }
                                else if (itname == cg_str)
                                {
                                    process_expert_options(it_list.sublist(cg_str),cg_out_list);
                                }
                                else
                                {
                                    MyAbort("Unrecognized optional parameter to \"" + it_str + "\" list: \"" + itname + "\"");
                                }
                            }
                        }
                        else if (NUMoptL[j] == prob_str)
                        {
                            process_expert_options(num_list.sublist(prob_str),prob_out_list);
                        }
                        else if (NUMoptL[j] == diffuse_str)
                        {
                            process_expert_options(num_list.sublist(diffuse_str),diffuse_out_list);
                        }
                        else
                        {
                            MyAbort("Unrecognized optional parameter to \"" + num_str + "\" list: \"" + NUMoptL[j] + "\"");
                        }
                    }
                }
                else {
                    MyAbort("Unrecognized optional parameter to \"" + ec_str + "\" list: \"" + optL[i] + "\"");
                }
            }
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
                                       ParameterList&       struc_list,
                                       StateDef&            state)
        {
            ParameterList& rock_list = struc_list.sublist("rock");
        
            const ParameterList& rlist = parameter_list.sublist("Material Properties");

            std::string mineralogy_str = "Mineralogy";
            std::string complexation_str = "Surface Complexation Sites";
            std::string isotherm_str = "Sorption Isotherms";
            std::string cation_exchange_str = "Cation Exchange Capacity";
            state.getSolid().has_cation_exchange = false; // until we encounter the keyword for a material

            Array<std::string> arrayrock;

            bool add_chemistry_properties = false;

            std::map<std::string,SolidChem> solid_chem;
            std::map<std::string,double> cation_exchange_capacity;
            std::map<std::string,ParameterList> rsublist_mat;
            std::string kp_file="kp";
            std::string pp_file="pp";

            Array<std::string> reqL, reqP, nullList;
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
                const ParameterEntry& entry = rlist.getEntry(label);
            
                std::string _label = underscore(label);
                if (entry.isList()) {
                    
                    // Add this rock label to list of rocks
                    arrayrock.push_back(_label);

                    ParameterList& rsublist = rsublist_mat[label];

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
                            else if (rlabel=="Intrinsic Permeability: Anisotropic Uniform"  || 
				     rlabel=="Intrinsic Permeability: Uniform") {
                                Array<double> array_p(2);
				if (rlabel=="Intrinsic Permeability: Anisotropic Uniform") {
				  array_p[1] = rsslist.get<double>("Vertical");
				  array_p[0] = rsslist.get<double>("Horizontal");
				}
				else {
				  array_p[0] = rsslist.get<double>("Value");
				  array_p[1] = array_p[0];
				}
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
                            else if (rlabel==mineralogy_str) {

                                add_chemistry_properties = true;

                                PLoptions minP(rsslist,nullList,nullList,false,false); // each optional list is a mineral
                                const Array<std::string>& minLabels = minP.OptLists();
                                for (int k=0; k<minLabels.size(); ++k) {
                                    const std::string& minLabel = minLabels[k];
                                    if (state.getSolid().IsAMineral(minLabel)) {
                                        const ParameterList& minSL = rsslist.sublist(minLabel);
                                        PLoptions minP1(minSL,nullList,nullList,true,false);
                                        const Array<std::string>& minLabels1 = minP1.OptParms();
                                        for (int L=0; L<minLabels1.size(); ++L) {
                                            const std::string& minLabel1 = minLabels1[L];
                                            if (minLabel1 == "Volume Fraction" ) {
                                                solid_chem[label].Mineral(minLabel).volume_frac = minSL.get<double>(minLabel1);
                                            }
                                            else if (minLabel1 == "Specific Surface Area" ) {
                                                solid_chem[label].Mineral(minLabel).specific_surface_area = minSL.get<double>(minLabel1);
                                            }
                                            else {
                                                std::cerr << "Unsupported Mineralogy condition for "
                                                          << minLabel << ": " << minLabel1 << std::endl;
                                                throw std::exception();
                                            }
                                        }
                                    }
                                    else {
                                        std::cerr << "Unknown mineral in " << rlabel << ": " << minLabel << std::endl;
                                        throw std::exception();                                                
                                    }
                                }
                            }
                            else if (rlabel==complexation_str) {

                                add_chemistry_properties = true;

                                PLoptions scsP(rsslist,nullList,nullList,false,true); // each optional list is a mineral
                                const Array<std::string>& scsLabels = scsP.OptLists();
                                for (int k=0; k<scsLabels.size(); ++k) {
                                    const std::string& scsLabel = scsLabels[k];
                                    if (state.getSolid().IsASorptionSite(scsLabel)) {
                                        const ParameterList& scsSL = rsslist.sublist(scsLabel);
                                        PLoptions scsP1(scsSL,nullList,nullList,true,false);
                                        const Array<std::string>& scsLabels1 = scsP1.OptParms();
                                        for (int L=0; L<scsLabels1.size(); ++L) {
                                            const std::string& scsLabel1 = scsLabels1[L];
                                            if (scsLabel1 == "Site Density" ) {
                                                solid_chem[label].SorptionSite(scsLabel).site_density = scsSL.get<double>(scsLabel1);
                                            }
                                            else {
                                                std::cerr << "Unsupported Surface Complexation condition for "
                                                          << scsLabel << ": " << scsLabel1 << std::endl;
                                                throw std::exception();
                                            }
                                        }
                                    }
                                    else {
                                        std::cerr << "Unknown Sorption Site in " << rlabel << ": " << scsLabel << std::endl;
                                        throw std::exception();                                                
                                    }
                                }
                            }
                            else if (rlabel==isotherm_str) {

                                add_chemistry_properties = true;

                                PLoptions sipP(rsslist,nullList,nullList,false,true); // each optional list is a phase
                                const Array<std::string>& sipLabels = sipP.OptLists();
                                for (int k=0; k<sipLabels.size(); ++k) {
                                    const std::string& sipLabel = sipLabels[k];
                                    std::string _sipLabel = underscore(sipLabel);

                                    StateDef::PhaseCompMap& pc_map = state.getPhaseCompMap();
                                    if (pc_map.find(sipLabel)==pc_map.end()) {
                                        std::cerr << "Unknown phase " << sipLabel << " in " << rlabel << " for " << label << std::endl;
                                        throw std::exception();                                                
                                    }

                                    const ParameterList& sipSL = rsslist.sublist(sipLabel);
                                    PLoptions sipcP(sipSL,nullList,nullList,false,true); // each optional list is a component
                                    const Array<std::string>& sipcLabels = sipcP.OptLists();
                                    for (int L=0; L<sipcLabels.size(); ++L) {
                                        const std::string& sipcLabel = sipcLabels[L];
                                        std::string _sipcLabel = underscore(sipcLabel);
                                        
                                        PHASE::CompMap& c_map = pc_map[sipLabel];
                                        if (c_map.find(sipcLabel)==c_map.end()) {
                                            std::cerr << "Unknown component " << sipcLabel
                                                      << " in phase " << sipLabel << " for " <<
                                                rlabel << " in " << label << std::endl;
                                            throw std::exception();                                                
                                        }
                                        
                                        const ParameterList& sipcSL = sipSL.sublist(sipcLabel);
                                        PLoptions sipcsP(sipcSL,nullList,nullList,false,true); // each optional list is a solute
                                        const Array<std::string>& sipcsLabels = sipcsP.OptLists();
                                        for (int M=0; M<sipcsLabels.size(); ++M) {
                                            const std::string& sipcsLabel = sipcsLabels[M];
                                            std::string _sipcsLabel = underscore(sipcsLabel);
                                            const ParameterList& sipcsSL = sipcSL.sublist(sipcsLabel);
                                            
                                            if ( !(c_map[sipcLabel].HasTracer(sipcsLabel)) ) {
                                                std::cerr << "Unknown solute " << sipcsLabel << " in component "
                                                          << sipcLabel << " in phase " << sipLabel << " for " <<
                                                    rlabel << " in " << label << std::endl;
                                                throw std::exception();
                                            }
                                            

                                            SolidChem::SorptionIsothermData iso;

                                            reqP.clear();
                                            std::string Kd_str("Kd");
                                            std::string Lb_str("Langmuir b");
                                            std::string Fn_str("Freundlich n");
                                            PLoptions siP(sipcsSL,nullList,reqP,true,false);
                                            const Array<std::string>& siLabels = siP.OptParms();
					    bool n_found = false;
					    bool b_found = false;
                                            for (int N=0; N<siLabels.size(); ++N) {
                                                if (siLabels[N] == Kd_str) {
                                                    iso.Kd = sipcsSL.get<double>(siLabels[N]);
                                                }
                                                else if (siLabels[N] == Lb_str) {
                                                    iso.Langmuir_b = sipcsSL.get<double>(siLabels[N]);
                                                    iso.Freundlich_n = -1.0;
						    n_found = true;
                                                }
                                                else if (siLabels[N] == Fn_str) {
                                                    iso.Freundlich_n = sipcsSL.get<double>(siLabels[N]);
                                                    iso.Langmuir_b = -1.0;
						    b_found = true;
                                                }
                                                else {
                                                    std::cerr << "Unknown parameter for \"" << rlabel << "\": \"" << siLabels[N] 
                                                              << "\" in solute \"" << sipcsLabel << "\" of component \""
                                                              << sipcLabel << "\" in phase \"" << sipLabel << "\"." << std::endl;
                                                    throw std::exception();
                                                }
					    }
					    if (b_found && n_found) {
					      std::cerr << "Only one \"" << Lb_str << "\" and \"" << Fn_str 
							<< "\" may be specified for each solute.  Both given for \"" 
							<< sipcsLabel << "\" of component \""
							<< sipcLabel << "\" in phase \"" << sipLabel << "\"." << std::endl;
					      throw std::exception();
                                            }

                                            bool successfully_inserted = SolidChem::InsertSorptionIsotherm(sipcsLabel,iso);
                                            if (!successfully_inserted) {
                                                bool compatible = iso.IsFreundlich() ^ SolidChem::SorptionIsotherm(sipcsLabel).IsFreundlich();
                                                if (!compatible) {
                                                    std::cerr << "Only one \"" << Lb_str << "\" and \"" << Fn_str 
                                                              << "\" may be specified for each solute.  Both given for \"" 
                                                              << sipcsLabel << "\" in different materials." << std::endl;
                                                    throw std::exception();
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            else {
                                std::cerr << "Unrecognized Material Property: " << rlabel << std::endl;
                                throw std::exception();
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
                            mtest["Density"] = true;
                        }
                        else if (rlabel==cation_exchange_str) {
                            add_chemistry_properties = true;                                
                            cation_exchange_capacity[label] = rslist.get<double>(rlabel);
                            state.getSolid().has_cation_exchange = true;
                        }
                        else {

                            std::cerr << "Unrecognized Material Property: " << rlabel << std::endl;
                            throw std::exception();
                        }
                    }

                    // check for complete
                    std::vector<std::string> region_check = remaining_false(mtest); 
                    if (region_check.size()) {
		      for (int i=0; i<region_check.size(); ++i) {
			if (region_check[i]=="Capillary_Pressure") {
			  mtest[region_check[i]] = true;
			  rsublist.set("cpl_type",0);
			}
			else if (region_check[i]=="Relative_Permeability") {
			  mtest[region_check[i]] = true;
			  rsublist.set("kr_type",0);
			}
		      }
		      region_check = remaining_false(mtest); 
		      if (region_check.size()) {
                        std::cerr << "Material not completely defined: " << label << std::endl;
                        std::cerr << "   unfilled: ";
                        for (int i=0; i<region_check.size(); ++i)
			  std::cerr << region_check[i] << " ";
                        std::cerr << '\n';
                        throw std::exception();
		      }
		    }
		}
                else {

		    std::string sat_str = "Saturation Threshold For vg Kr";
                    if (rlist.isParameter(sat_str)) {
                        double saturation_threshold_for_vg_Kr;
                        saturation_threshold_for_vg_Kr = rlist.get<double>(sat_str);
                        rock_list.set(underscore(sat_str),saturation_threshold_for_vg_Kr);
                    }
                    
		    std::string shift_str = "Use Shifted Kr Eval";
                    if (rlist.isParameter(shift_str)) {
                        bool use_shifted_kr_eval = rlist.get<bool>(shift_str);
                        rock_list.set(underscore(shift_str),(int)use_shifted_kr_eval);
                    }
                    
                    if (rlist.isParameter("Permeability Output File"))
                        kp_file = rlist.get<std::string>("Permeability Output File");
                    if (rlist.isParameter("Porosity Output File"))
                        pp_file = rlist.get<std::string>("Porosity Output File");
                }

            }

            if (add_chemistry_properties)
            {
                
                const Array<std::string>& minerals = state.getSolid().mineral_names;
                const Array<std::string>& sorption_sites = state.getSolid().sorption_site_names;

                for (int k=0; k<arrayrock.size(); ++k)
                {
                    const std::string& _label = arrayrock[k];
                    const std::string& label = AMR_to_Amanzi_label_map[_label];

                    ParameterList& rsublist = rsublist_mat[label];

                    // Add chemistry data, if necessary
                    if (minerals.size()>0) {
                        ParameterList mPL;
                        for (int i=0; i<minerals.size(); ++i) {
                            const std::string& name = minerals[i];
                            SolidChem::MineralData md = solid_chem[label].Mineral(name); // Will call defctr if not set by inputs
                            ParameterList minPL = md.BuildPL();
                            mPL.set(underscore(name),minPL);
                        }
                        rsublist.set(underscore(mineralogy_str),mPL);
                    }

                    if (sorption_sites.size()>0) {
                        ParameterList sPL;
                        for (int i=0; i<sorption_sites.size(); ++i) {
                            const std::string& name = sorption_sites[i];
                            SolidChem::SorptionSiteData ssd = solid_chem[label].SorptionSite(name); 
                            ParameterList ssdPL = ssd.BuildPL();
                            sPL.set(underscore(name),ssdPL);
                        }
                        rsublist.set(underscore(complexation_str),sPL);
                    }
                        
		    if (cation_exchange_capacity.count(label)) {
		      rsublist.set(underscore(cation_exchange_str),cation_exchange_capacity[label]);
		    }

                    // if ntracers>0...
                    //
                    // Must "flatten" hierarchy of phases/comps to be compatible with current Amanzi-S 
                    // in fact, tracers are listed flat requiring that there be no name clashes across phase/comp
                    ParameterList siPL;
                    StateDef::Phases& phases = state.getPhases();
                    std::set<std::string> solutes_with_isotherms;
                    for (StateDef::Phases::iterator pit=phases.begin(); pit!=phases.end(); ++pit) {
                        const std::string& p=pit->first;
                        StateDef::CompMap& comps = state[p];
                        for (StateDef::CompMap::iterator cit=comps.begin(); cit!=comps.end(); ++cit) {
                            const std::string& c=cit->first;
                            const Array<std::string>& solutes = cit->second.getTracerArray();
                            for (int i=0; i<solutes.size(); ++i) {
                                const std::string& s=solutes[i];
                                if (SolidChem::HasSorptionIsotherm(s)) {
                                    SolidChem::SorptionIsothermData sid = SolidChem::SorptionIsotherm(s);
                                    ParameterList sitPL = sid.BuildPL();
                                    siPL.set(underscore(s),sitPL);
                                    solutes_with_isotherms.insert(s);
                                }
                            }
                        }
                    }
                    if (solutes_with_isotherms.size()>0) {
                        rsublist.set(underscore(isotherm_str),siPL);
                        state.getSolid().sorption_isotherm_names.resize(solutes_with_isotherms.size());
                        for (std::set<std::string>::const_iterator it=solutes_with_isotherms.begin(), End=solutes_with_isotherms.end(); it!=End; ++it) {
                            state.getSolid().sorption_isotherm_names.push_back(*it);
                        }
                    }
                }
            }

            // Now we can add the sublist the result
            for (int k=0; k<arrayrock.size(); ++k)
            {
                const std::string& _label = arrayrock[k];
                const std::string& label = AMR_to_Amanzi_label_map[_label];
                rock_list.set(_label,rsublist_mat[label]);
            }
            
            rock_list.set("rock",arrayrock);
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
            reqL.push_back("Phase Definitions");
            PLoptions opt(parameter_list,reqL,nullList,false,false);
            const ParameterList& plist = parameter_list.sublist(reqL[0]);

            PLoptions optP(plist,nullList,nullList,false,false); // each optional list is a phase 
            const Array<std::string>& phaseLabels = optP.OptLists();

            for (int i=0; i<phaseLabels.size(); ++i) {
                const std::string& phaseLabel = phaseLabels[i];
                const ParameterList& psublist = plist.sublist(phaseLabel);
                const ParameterList& pplist = plist.sublist(phaseLabel);

                if (phaseLabel=="Solid") {
                    PLoptions optP1(pplist,nullList,nullList,true,false);
                    Array<std::string> optPpp = optP1.OptParms();
                    for (int j=0; j<optPpp.size(); ++j) {
                        const std::string& name = optPpp[j];
                        if (name == "Minerals") {
                            getSolid().mineral_names = pplist.get<Array<std::string> >(name);
                        }
                        else if (name == "Sorption Sites") {
                            getSolid().sorption_site_names = pplist.get<Array<std::string> >(name);
                        }
                        else {
                            std::cerr << "Unrecognized Solid phase parameter: " << name << std::endl;
                            throw std::exception();
                        }
                    }
                }
                else {
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
                        
                        if (propLabel == "Density: Uniform") {
                            Array<std::string> reqP;
                            reqP.push_back("Density");
                            PLoptions optPD(ppolist,nullList,reqP,true,true); 
                            density = ppolist.get<double>(reqP[0]);
                        }
                        else if (propLabel == "Viscosity: Uniform") {
                            Array<std::string> reqP;
                            reqP.push_back("Viscosity");
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
            Array<std::string> nullList, reqP;
            const std::string val_name="Value"; reqP.push_back(val_name);
            PLoptions opt(fPLin,nullList,reqP,true,true); 
            // FIXME: Assumes Water exists, and that this is what was intended....
            std::string _name = underscore("Water");
            fPLout.set<double>(_name,fPLin.get<double>(val_name));
            fPLout.set<std::string>("type","saturation");
        }


        void convert_IC_ConstPressure(const ParameterList& fPLin,
				      const std::string&   Amanzi_type,
				      ParameterList&       fPLout)
        {
            Array<std::string> nullList, reqP;
	    const std::string phase_name="Phase";reqP.push_back(phase_name);
            const std::string val_name="Value";reqP.push_back(val_name);
            PLoptions opt(fPLin,nullList,reqP,true,true);      
	    fPLout.set<std::string>("type","pressure");
	    fPLout.set<double>("val",fPLin.get<double>(val_name));
	    fPLout.set<std::string>("phase",fPLin.get<std::string>(phase_name));
        }

        void convert_IC_LinPressure(const ParameterList& fPLin,
				    const std::string&   Amanzi_type,
				    ParameterList&       fPLout)
        {
            Array<std::string> reqP, nullList;
	    const std::string phase_name="Phase";reqP.push_back(phase_name);
            const std::string rval_name="Reference Value";reqP.push_back(rval_name);
            const std::string grad_name="Gradient Value";reqP.push_back(grad_name);
            const std::string ref_name="Reference Coordinate";reqP.push_back(ref_name);
            PLoptions opt(fPLin,nullList,reqP,true,true);  
    
	    fPLout.set<std::string>("type","hydrostatic");
	    fPLout.set<double>("val",fPLin.get<double>(rval_name));
	    const Array<double>& grad = fPLin.get<Array<double> >(grad_name);
	    const Array<double>& water_table = fPLin.get<Array<double> >(ref_name);
	    int coord = water_table.size()-1;
	    fPLout.set<double>("water_table_height",water_table[coord]); 
	    fPLout.set<double>("grad",grad[coord]);
        }

        void convert_ICFlow(const ParameterList& fPLin,
			    const std::string&   Amanzi_type,
			    ParameterList&       fPLout)
        {
	    const std::string phase_name="Phase";
            const std::string val_name="Reference Value";
            const std::string grad_name="Gradient Value";
            const std::string ref_name="Reference Coordinate";
            const std::string vel_name="Aqueous Volumetric Flux";

            Array<std::string> reqP, nullList;
            reqP.push_back(val_name);
            if (Amanzi_type == "IC: Flow") {
	        reqP.push_back(phase_name);
                reqP.push_back(grad_name);
                reqP.push_back(ref_name);
		reqP.push_back(vel_name);
            }
            PLoptions opt(fPLin,nullList,reqP,true,true);  
    
            fPLout.set<std::string>("type","zero_total_velocity");
            fPLout.set<double>("val",fPLin.get<double>(val_name));
            if (Amanzi_type == "IC: Flow") {
                const Array<double>& grad = fPLin.get<Array<double> >(grad_name);
                const Array<double>& water_table = fPLin.get<Array<double> >(ref_name);
		double AqVolFlux = fPLin.get<double>(vel_name);
                int coord = water_table.size()-1;
                fPLout.set<double>("water_table_height",water_table[coord]); 
                fPLout.set<double>("grad",grad[coord]);
                fPLout.set<double>("aqueous_vol_flux",AqVolFlux); 
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
                                            ParameterList&  fPLout,
                                            int             do_chem)
        {
            const ParameterList& fPLin = solute_ic.PList();
            const std::string& solute_ic_Amanzi_type = solute_ic.Amanzi_Type();
            const std::string& solute_ic_label = solute_ic.Label();
            const std::string& solute_ic_units = solute_ic.Units();

            Array<std::string> reqP, nullList;
            const std::string val_name="Value"; reqP.push_back(val_name);
            const std::string ion_name="Free Ion Guess";
            if (do_chem) {
                reqP.push_back(ion_name);
            }
            PLoptions opt(fPLin,nullList,reqP,true,true);  
            fPLout.set<double>("val",fPLin.get<double>(val_name));
            if (do_chem) {
	        fPLout.set<double>(underscore(ion_name),fPLin.get<double>(ion_name));
            }
            // Adjust dimensions of data
            if (solute_ic_units=="Molar Concentration" 
                || solute_ic_units=="Molal Concentration")
            {
	      //std::cerr << "IC label \"" << solute_ic_label
              //            << "\" function: \"" << solute_ic_Amanzi_type
              //            << "\" requests unsupported units: \"" << solute_ic_units
              //            << "\"" << std::endl;
              //  throw std::exception();
            }
            else if (solute_ic_units=="Specific Concentration") {
                // This is the units expected by the structured code
            }
            else {
	        //std::cerr << "Unsupported Solute IC function: \"" << solute_ic_units << "\"" << std::endl;
	        //throw std::exception();
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
                else if (Amanzi_type == "IC: Uniform Pressure")
                {
                    convert_IC_ConstPressure(fPLin,Amanzi_type,fPLout);
                }
                else if ( Amanzi_type == "IC: Linear Pressure" )
                {
                    convert_IC_LinPressure(fPLin,Amanzi_type,fPLout);
                }
                else if ( Amanzi_type == "IC: Hydrostatic" )
                {
                    convert_ICHydrostatic(fPLin,Amanzi_type,fPLout);
                }
                else if ( Amanzi_type == "IC: Flow" )
                {
		     convert_ICFlow(fPLin,Amanzi_type,fPLout);
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
                fPLout.set<Array<std::string> >("forms",fPLin.get<Array<std::string> >(form_name));
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
                    fPLout.set<Array<std::string> >("forms",fPLin.get<Array<std::string> >(form_name));
                }
            }
        }

        void convert_BCFlux(const ParameterList& fPLin,
                            const std::string&   Amanzi_type,
                            ParameterList&       fPLout,
			    StateDef&            stateDef)
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
            PLoptions opt(fPLin,reqL,reqP,false,true);
    
            Array<double> times, fluxvals = fPLin.get<Array<double> >(val_name);
            Array<std::string> forms;
            if (fluxvals.size()>1) {
                times = fPLin.get<Array<double> >(time_name);
                forms = fPLin.get<Array<std::string> >(form_name);
            }
    
            // Convert mass flux to volumetric flux
            if (is_mass) {
  	        double density = stateDef.getPhases()["Aqueous"].Density();
		for (int i=0; i<fluxvals.size(); ++i) {
		  fluxvals[i] *= 1/density;
		}
            }

            // Convert to inward flux
            if (is_out) {
                for (int i=0; i<fluxvals.size(); ++i) {
                    fluxvals[i] = -fluxvals[i];
                }
            }

	    const Array<std::string>& optional_lists = opt.OptLists();
	    if (optional_lists.size() > 0) {
	      
	      if (optional_lists.size()!=1) {
		std::cerr << "BC: Flux - invalid optional arg(s): ";
		for (int i=0; i<optional_lists.size(); ++i) {
		  std::cerr << "\"" << optional_lists[i] << "\" ";
		}
		std::cerr << std::endl;
	      }
	    }
    
            fPLout.set<Array<double> >("aqueous_vol_flux",fluxvals);
            if (fluxvals.size() > 1) {
                fPLout.set<Array<double> >("inflowtimes",fPLin.get<Array<double> >(time_name));
                fPLout.set<Array<std::string> >("inflowfncs",fPLin.get<Array<std::string> >(form_name));
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
	      //std::cerr << "BC label \"" << solute_bc_label
	      //           << "\" function: \"" << solute_bc_Amanzi_type
	      //          << "\" requests unsupported units: \"" << solute_bc_units
	      //          << "\"" << std::endl;
	      //throw std::exception();
            }
            else if (solute_bc_units=="Specific Concentration") {
                // This is the units expected by the structured code
            } 
            else {
	      //std::cerr << "Solute BC - invalid units: \"" << solute_bc_units << "\"" << std::endl;
	      //throw std::exception();
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
		  convert_BCFlux(fPLin,Amanzi_type,fPLout,stateDef);
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

		    int& sat_bc      = (i<3  ?        lo_bc[k] :  hi_bc[k]);
		    int& pressure_bc = (i<3  ?       plo_bc[k] : phi_bc[k]);
		    int& inflow_bc   = (i<3  ? inflow_lo_bc[k] : inflow_hi_bc[k]);

                    if (orient_RTs.size()!=0) {
  		        const std::string& orient_type = orient_RTs[0].second;
			for (int j=1; j<orient_RTs.size(); ++j) {
			    if (orient_type != orient_RTs[j].second) {
			      std::cerr << "Structured grid requires that all BCs on "
					<< Amanzi_purpose << " be of the same type" << std::endl;
			      throw std::exception();
			    }
			}


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
			    pressure_bc = 1; // Neumann for p
			    inflow_bc   = 1; // Requires inflow_XX_vel velocity values for nphase-1 phases
			}
			else {
			  std::cerr << "Structured grid inputs translator generated unrecognized BCs "
				    << orient_type << std::endl;
			  throw std::exception();
			}
		    }
		    else {
  		        sat_bc = 3;
			pressure_bc = 4;
			inflow_bc = 4;
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
        convert_solute_ics(StateDef& stateDef, int do_chem)
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
                                convert_solute_ICConcentration(solute_ic,fPL,do_chem);
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
                                    bool                 do_tracer_transport,
                                    bool                 do_chem)
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
            Array<double> arraydiffusivity;  // FIXME: Not in current spec
    
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

                arrayphase.push_back(_phaseLabel);
                arraydensity.push_back(phase.Density());
                arrayviscosity.push_back(phase.Viscosity());
                arraydiffusivity.push_back(phase.Diffusivity());
        
                Array<std::string> arraycomp;  
                const CompMap comp_map = stateDef[phaseLabel];
                for (CompMap::const_iterator cit = comp_map.begin(); cit!=comp_map.end(); ++cit) 
                {
                    const std::string& compLabel = cit->first;
                    std::string _compLabel = underscore(compLabel);
                    arraycomp.push_back(_compLabel);

                    const Array<std::string>& soluteNames = cit->second.getTracerArray();
                    for (int i=0; i<soluteNames.size(); ++i)
                    {
                        std::string _soluteLabel = underscore(soluteNames[i]);
                        arraysolute.push_back(_soluteLabel);
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

            std::map<std::string,ParameterList> solutePLs;

            typedef SolutePLMMap::const_iterator SPLit;
            SPLit it;
            SolutePLMMap solute_to_ictype = convert_solute_ics(stateDef,do_chem);

            for (int i=0; i<arraysolute.size(); ++i) {
                const std::string& soluteName = arraysolute[i];
                
                Array<std::string> icLabels;
                std::pair<SPLit,SPLit> retIC = solute_to_ictype.equal_range(soluteName);
                for (it=retIC.first; it!=retIC.second; ++it) {
                    const ParameterList& pl=it->second;
                    for (ParameterList::ConstIterator pit=pl.begin(); pit!=pl.end(); ++pit) {
                        const std::string& name = pl.name(pit);
                        solutePLs[soluteName].setEntry(name,pl.getEntry(name));
                    }
                    icLabels.push_back(it->second.name());
                }
     
                Array<std::string> regions;
                solutePLs[soluteName].set<Array<std::string> >("regions",regions);
                solutePLs[soluteName].set<Array<std::string> >("tinits",icLabels);
            }

            // Only do solute BCs if do_tracer_transport
            if (do_tracer_transport)
            {
                SolutePLMMap solute_to_bctype = convert_solute_bcs(stateDef);

                for (int i=0; i<arraysolute.size(); ++i) {
                    const std::string& soluteName = arraysolute[i];
                    Array<std::string> bcLabels;
                    std::pair<SPLit,SPLit> retBC = solute_to_bctype.equal_range(soluteName);
                    for (it=retBC.first; it!=retBC.second; ++it) {
                        solutePLs[soluteName].set(it->second.name(),it->second);
                        bcLabels.push_back(it->second.name());
                    }
                    solutePLs[soluteName].set<Array<std::string> >("tbcs",bcLabels);
                    
                }
            }

            // Now add solute info
            for (int i=0; i<arraysolute.size(); ++i) {
                const std::string& soluteName = arraysolute[i];
                solute_list.set(soluteName,solutePLs[soluteName]);
            }


            if (do_chem) 
            {
                if (stateDef.HasSolidChem()) {
                    const Array<std::string>& mineral_names = stateDef.getSolid().mineral_names;
                    if (mineral_names.size() > 0) {
                        struc_list.sublist("mineral").set("minerals",underscore(mineral_names));
                    }
                    
                    const Array<std::string>& sorption_site_names = stateDef.getSolid().sorption_site_names;
                    if (sorption_site_names.size() > 0) {
                        struc_list.sublist("sorption_site").set("sorption_sites",underscore(sorption_site_names));
                    }
                }
            }
        }

        //
        // convert output to structured format
        //
        void
        convert_to_structured_output(const ParameterList& parameter_list, 
                                     ParameterList&       struc_list,
                                     StateDef&            state,
                                     bool                 do_chem)
        {
            ParameterList& amr_list = struc_list.sublist("amr");
            ParameterList& obs_list = struc_list.sublist("observation");

            // Create list of available field quantities.  All these must be recognized
            //  by name inside AmrLevel::derive
            user_derive_list.push_back(underscore("Material ID"));
            user_derive_list.push_back(underscore("Grid ID"));
            user_derive_list.push_back(underscore("Core ID"));
            user_derive_list.push_back(underscore("Cell ID"));
            user_derive_list.push_back(underscore("Capillary Pressure"));
            user_derive_list.push_back(underscore("Volumetric Water Content"));
            user_derive_list.push_back(underscore("Porosity"));
            user_derive_list.push_back(underscore("Aqueous Saturation"));
            user_derive_list.push_back(underscore("Aqueous Pressure"));
            user_derive_list.push_back(underscore("Aqueous Volumetric Flux X"));
            user_derive_list.push_back(underscore("Aqueous Volumetric Flux Y"));
#if BL_SPACEDIM==3
            user_derive_list.push_back(underscore("Aqueous Volumetric Flux Z"));
#endif

            if (struc_list.isSublist("tracer")) {
                const Array<std::string>& solute_names = struc_list.sublist("tracer").get<Array<std::string> >("tracers");
                for (int i=0; i<solute_names.size(); ++i) {
                    const std::string& name = solute_names[i];
                    user_derive_list.push_back(underscore("Aqueous "+name+" Concentration"));
                }
            }
            
            if (do_chem && state.HasSolidChem()) {
                if (struc_list.isSublist("tracer")) {
                    const Array<std::string>& solute_names = struc_list.sublist("tracer").get<Array<std::string> >("tracers");
                    for (int i=0; i<solute_names.size(); ++i) {
                        const std::string& name = solute_names[i];
                        if (state.getSolid().UsingSorption()) {
                            user_derive_list.push_back(underscore("Total Sorbed "+name));
                        }
                        if (SolidChem::HasSorptionIsotherm(name)) {
                            user_derive_list.push_back(underscore("Kd "+name));
                            if (SolidChem::SorptionIsotherm(name).IsFreundlich()) {
                                user_derive_list.push_back(underscore("Freundlich n "+name));
                            }
                            else if (SolidChem::SorptionIsotherm(name).IsLangmuir()) 
                            {
                                user_derive_list.push_back(underscore("Langmuir b "+name));
                            }
                        }
                        user_derive_list.push_back(underscore("Free Ion Guess "+name));
                        user_derive_list.push_back(underscore("Activity Coefficient "+name));
                    }
                }

                if (state.getSolid().has_cation_exchange) {
                    user_derive_list.push_back(underscore("Cation Exchange Capacity"));
                }

                const Array<std::string>& mineral_names = state.getSolid().mineral_names;
                for (int i=0; i<mineral_names.size(); ++i) {
                    const std::string& name = mineral_names[i];
                    user_derive_list.push_back(underscore("Volume Fraction "+name));
                    user_derive_list.push_back(underscore("Specific Surface Area "+name));
                }

                const Array<std::string>& sorption_site_names = state.getSolid().sorption_site_names;
                for (int i=0; i<sorption_site_names.size(); ++i) {
                    const std::string& name = sorption_site_names[i];
                    user_derive_list.push_back(underscore("Site Density "+name));
                }
            }

            amr_list.set<Array<std::string> >("user_derive_list",user_derive_list);

            const ParameterList& rlist = parameter_list.sublist("Output");
      
            // time macros
            std::set<std::string> time_macros;
            const ParameterList& tlist = rlist.sublist("Time Macros");
            ParameterList tmPL;
            for (ParameterList::ConstIterator i=tlist.begin(); i!=tlist.end(); ++i) {
                std::string label = tlist.name(i);
                std::string _label = underscore(label);
                const ParameterList& rslist = tlist.sublist(label);
                Array<double> times, vals;
                
                ParameterList tPL;
                for (ParameterList::ConstIterator ii=rslist.begin(); ii!=rslist.end(); ++ii) {
                    const std::string& name = underscore(rslist.name(ii));
                    
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
                    tmPL.set(_label,tPL);
                    time_macros.insert(_label);
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
	    Array<std::string> cma;
	    if (rlist.isSublist("Cycle Macros")) {
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
	      
	      for (std::set<std::string>::const_iterator it=cycle_macros.begin(); it!=cycle_macros.end(); ++it) {
                cma.push_back(*it);
	      }
            }
	    amr_list.set<Array<std::string> >("cycle_macros",cma);


            // vis data
            const std::string vis_data_str = "Visualization Data";
            const std::string vis_file_str = "File Name Base";
            const std::string vis_var_str = "Variables";
            const std::string vis_cycle_str = "Cycle Macros";
            const std::string vis_time_str = "Time Macros";
            bool vis_vars_set = false;
            Array<std::string> visNames, vis_cMacroNames, vis_tMacroNames;
            int vis_digits = 5;
            std::string vis_file = "plt";
            if (rlist.isSublist(vis_data_str)) {                
                const ParameterList& vlist = rlist.sublist(vis_data_str);
                for (ParameterList::ConstIterator i=vlist.begin(); i!=vlist.end(); ++i)
                {
                    const std::string& name = vlist.name(i);

                    if (name == vis_file_str) {
                        vis_file = vlist.get<std::string>(vis_file_str);
                    }
                    else if (name == vis_var_str)
                    {                        
                        visNames = vlist.get<Array<std::string> >(vis_var_str);
                        for (int i=0; i<visNames.size(); ++i) {
                            std::string _visName = underscore(visNames[i]);
                            bool found = false;
                            for (int j=0; j<user_derive_list.size(); ++j) {
                                if (_visName == user_derive_list[j]) {
                                    found = true;
                                }
                            }
                            if (!found) {
                                MyAbort("Invalid variable (\""+visNames[i]+
                                        "\") in \"Visualization Data\" -> \"Variables\"");
                            }
                            visNames[i] = underscore(visNames[i]);
                        }
                        vis_vars_set = true;
                    }
                    else if (name == vis_cycle_str)
                    {
                        const Array<std::string>& vcMacros = vlist.get<Array<std::string> >(vis_cycle_str);
                        for (int i=0; i<vcMacros.size(); ++i) {
                            std::string label = underscore(vcMacros[i]);
                            if (cycle_macros.find(label) != cycle_macros.end()) {
                                vis_cMacroNames.push_back(label);
                            }
                            else {
                                if (Teuchos::MPISession::getRank() == 0) {
                                    std::cerr << "Unrecognized cycle macro in \""+vis_data_str+"\": \""
                                              << vcMacros[i] << "\"" << std::endl;
                                    
                                    for (std::set<std::string>::const_iterator it=cycle_macros.begin();
                                         it!=cycle_macros.end(); ++it) {
                                        std::cout << *it << " " << std::endl;
                                    }
                                    throw std::exception();
                                }
                            }
                        }
                    }
                    else if (name == vis_time_str)
                    {
                        const Array<std::string>& vtMacros = vlist.get<Array<std::string> >(vis_time_str);
                        for (int i=0; i<vtMacros.size(); ++i) {
                            std::string label = underscore(vtMacros[i]);
                            if (time_macros.find(label) != time_macros.end()) {
                                vis_tMacroNames.push_back(label);
                            }
                            else {
                                if (Teuchos::MPISession::getRank() == 0) {
                                    std::cerr << "Unrecognized time macro in \""+vis_data_str+"\": \""
                                              << vtMacros[i] << "\"" << std::endl;
                                    
                                    for (std::set<std::string>::const_iterator it=time_macros.begin();
                                         it!=time_macros.end(); ++it) {
                                        std::cout << *it << " " << std::endl;
                                    }
                                    throw std::exception();
                                }
                            }
                        }
                    }
                    else {
                        MyAbort("Unrecognized entry in \""+vis_data_str+"\" parameter list: \""+name+"\"");
                    }
                }
            }

            //
            // Set default to dump all known fields
            if (!vis_vars_set) {
                visNames.resize(user_derive_list.size());
                for (int j=0; j<user_derive_list.size(); ++j) {
                    visNames.push_back(underscore(user_derive_list[j])); 
                }
            }


            amr_list.set<Array<std::string> >("derive_plot_vars",visNames);
            amr_list.set<Array<std::string> >("vis_cycle_macros",vis_cMacroNames);
            amr_list.set<Array<std::string> >("vis_time_macros",vis_tMacroNames);
            amr_list.set<std::string>("plot_file",vis_file);


            // chk data
            const std::string chk_data_str = "Checkpoint Data";
            const std::string chk_file_str = "File Name Base";
            const std::string chk_cycle_str = "Cycle Macros";
            const std::string chk_digits_str = "File Name Digits";
            bool cycle_macro_set = false;
            std::string chk_file = "chk";
            Array<std::string> chk_cMacroNames;
            if (rlist.isSublist(chk_data_str)) {                
                const ParameterList& clist = rlist.sublist(chk_data_str);
                for (ParameterList::ConstIterator i=clist.begin(); i!=clist.end(); ++i)
                {
                    const std::string& name = clist.name(i);
                    if (name == chk_file_str) {
                        chk_file = clist.get<std::string>(chk_file_str);
                    }
                    else if (name == chk_cycle_str) {
                        const Array<std::string>& ccMacros = clist.get<Array<std::string> >(chk_cycle_str);
                        for (int i=0; i<ccMacros.size(); ++i) {
                            std::string label = underscore(ccMacros[i]);
                            if (cycle_macros.find(label) != cycle_macros.end()) {
                                chk_cMacroNames.push_back(label);
                            }
                            else {
                                if (Teuchos::MPISession::getRank() == 0) {
                                    std::cerr << "Unrecognized cycle macro in \""+chk_data_str+"\": \""
                                              << ccMacros[i] << "\"" << std::endl;
                                    
                                    for (std::set<std::string>::const_iterator it=cycle_macros.begin();
                                         it!=cycle_macros.end(); ++it) {
                                        std::cout << *it << " " << std::endl;
                                    }
                                    throw std::exception();
                                }
                            }
                        }
                    }
                    else {
                        MyAbort("Unrecognized entry in \""+chk_data_str+"\" parameter list: \""+name+"\"");
                    }
                }
            }
            amr_list.set<Array<std::string> >("chk_cycle_macros",chk_cMacroNames);
            amr_list.set<std::string>("check_file",chk_file);

            const std::string file_name_digits_str = "File Name Digits";
            int file_name_digits = 5;
            if (rlist.isParameter(file_name_digits_str)) {
                file_name_digits = rlist.get<int>(file_name_digits_str);
                if (file_name_digits<=0) {
                    MyAbort("Output -> \""+file_name_digits_str+"\" must be > 0");
                }
            }
            amr_list.set<int>("file_name_digits",file_name_digits);


        
            // observation
            Array<std::string> arrayobs;
            std::string obs_str = "Observation Data";
	    if (rlist.isSublist(obs_str)) {
	      const ParameterList& olist = rlist.sublist(obs_str);
	      ParameterList sublist;
	      std::string obs_file_str = "Observation Output Filename";
	      std::string obs_file="observation.out";
	      for (ParameterList::ConstIterator i=olist.begin(); i!=olist.end(); ++i) {                
                std::string label = olist.name(i);
                const ParameterEntry& entry = olist.getEntry(label);

                if (entry.isList()) {
		  std::string _label = underscore(label);
		  const ParameterList& rslist = olist.sublist(label);
		  std::string functional = rslist.get<std::string>("Functional");
		  const std::string& region_name = rslist.get<std::string>("Region");
		  std::string _region_name = underscore(region_name);
		  const std::string& variable = rslist.get<std::string>("Variable");
		  std::string _variable = underscore(variable);

		  if (functional == "Observation Data: Integral") {
		    sublist.set("obs_type","integral");
		  }
		  else if (functional == "Observation Data: Point") {
                      sublist.set("obs_type","point_sample");

                      // Check our list of regions to ensure it exists and is the correct type
		      const ParameterList& lregion = 
                          parameter_list.sublist("Regions").sublist(region_name);

		      if (!lregion.isSublist("Region: Point"))
                      {
			  std::cerr << label << " is a point observation and "
				    << region_name << " is not a point region.\n";
			  throw std::exception();
                        }
                    }
		  sublist.set("region",_region_name);

                  std::string Time_Macro_str = "Time Macro";
                  std::string Cycle_Macro_str = "Cycle Macro";
                  if (rslist.isParameter(Time_Macro_str)) {
                      const std::string& _timeMacro = underscore(rslist.get<std::string>(Time_Macro_str));
                      if (time_macros.find(_timeMacro) == time_macros.end()) {
                          std::cerr << "Unrecognized time macro: \"" << AMR_to_Amanzi_label_map[_timeMacro]
                                    << "\" for observation data: \"" << label << "\"" << std::endl;
                          throw std::exception();                        
                      }
                      sublist.set("time_macro",_timeMacro);
                  }
		  else if (rslist.isParameter(Cycle_Macro_str)) {
                      const std::string& _cycleMacro = underscore(rslist.get<std::string>(Cycle_Macro_str));
                      if (cycle_macros.find(_cycleMacro) == cycle_macros.end()) {
                          std::cerr << "Unrecognized cycle macro: \"" << AMR_to_Amanzi_label_map[_cycleMacro]
                                    << "\" for observation data: \"" << label << "\"" << std::endl;
                          throw std::exception();                        
                      }
                      sublist.set("cycle_macro",_cycleMacro);
		  }
                  else {
                      std::cerr << "Must specify either time or cycle macro forobservation data: \"" 
                                << label << "\"" << std::endl;
                      throw std::exception();                 
                  }
	  
		  bool found = false;
		  for (int k=0; k<user_derive_list.size() && !found; ++k) {
		    if (_variable == user_derive_list[k]) {
		      found = true;
		    }
		  }
		  if (!found) {
		    MyAbort(variable + " is not a valid derive variable name");
		  }

		  sublist.set<std::string>("field",_variable);

		  obs_list.set(_label,sublist);
		  arrayobs.push_back(_label);
                }
                else {
		  if (label == obs_file_str) {
		    obs_file = underscore(olist.get<std::string>(obs_file_str));
		  }
		  else {
		    MyAbort("Unrecognized option under \""+obs_str+"\": \""+label+"\"" );
		  }
                }
	      }
	      obs_list.set("output_file",obs_file);
	      obs_list.set("observation",arrayobs);
	    }
	}

        //
        // convert parameterlist to format for structured code
        //
        ParameterList
        convert_to_structured(const ParameterList& parameter_list)
        {
            ParameterList struc_list = setup_structured();
            bool do_tracer_transport, do_chem;
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
            convert_to_structured_control(parameter_list,struc_list, do_tracer_transport, do_chem);
            //
            // Regions
            //
            convert_to_structured_region(parameter_list, struc_list);
            //
            // State 
            //
            StateDef stateDef(parameter_list);
            convert_to_structured_state(parameter_list, struc_list, stateDef, do_tracer_transport, do_chem);
            //
            // Materials
            //
            convert_to_structured_material(parameter_list, struc_list, stateDef);
            //
            // Output
            // 
            convert_to_structured_output(parameter_list, struc_list,stateDef,do_chem);

            std::string dump_str = "Dump ParmParse Table";
            if (parameter_list.isParameter(dump_str)) {
                struc_list.set<std::string>("dump_parmparse_table",parameter_list.get<std::string>(dump_str));
            }
            return struc_list;
        }

    }
}


/*
Mesh
  Structured
    Number of Cells
    Domain Low Corner
    Domain High Corner
     ... sets geometry_eps = 1.e-6*max_size;

Domain
  Spatial Dimension

Execution Control

  Flow Model:
     Off (do_simple=2)
     Richards (... sets model_name=richard, have_capillary=1, cfl=-1)
     Single-phase (... sets model_name=single-phase, do_simple=1)
     Multi-phase (... sets model_name=two-phase, cfl=0.75)

  Transport Model
     Off (... sets do_tracer_transport=0)
     On  (... sets do_tracer_transport=1)

  Chemistry Model
     Off (... sets do_chem=-1)

  Time Integration Mode
     Steady
     Transient
       Start (start_time)
       End (stop_time)
       Initial Time Step (dt_init)
       Maximum Time Step Change (change_max)
       Initial Time Step Multiplier (init_shrink)
       Maximum Time Step Size (dt_max)
       Maximum Cycle Number (max_step)
     Initialize To Steady

  Verbosity
    None     ... sets prob_v = 0; mg_v = 0; cg_v = 0; amr_v = 0; diffuse_v = 0
    Low      ... sets prob_v = 1; mg_v = 0; cg_v = 0; amr_v = 1; diffuse_v = 0
    Medium   ... sets prob_v = 1; mg_v = 0; cg_v = 0; amr_v = 2; diffuse_v = 1
    High     ... sets prob_v = 2; mg_v = 1; cg_v = 1; amr_v = 3; diffuse_v = 1
    Extreme  ... sets prob_v = 3; mg_v = 2; cg_v = 2; amr_v = 3; diffuse_v = 1

  Restart from Checkpoint Data File (amr.restart)
  Echo Inputs (TRUE-> echo_inputs=true)
  Dump ParmParse Table

  Numerical Control Parameters

    Adaptive Mesh Refinement Control
      Number Of AMR Levels (1)
      Refinement Ratio (2)
      Regrid Interval (2)
      Blocking Factor (8)
      Number Error Buffer Cells (1)
      Maximum Grid Size (128, 32)
      Refinement Indicators
        Value Greater, Value Less, Adjacent Difference Greater, Inside Region
        Field Name
        Regions
        Maximum Refinement Level
        Start Time
        End Time
      Expert Settings

      ....set amr.nosub, is_periodic=0

    Basic Algorithm Control
      Expert Settings
      ...set visc_abs_tol=1.e-16, visc_tol=1.e-14

    Iterative Linear Solver Control
      Conjugate Gradient Algorithm
        Expert Settings
      Multigrid Algorithm
        Expert Settings

    Pressure Discretization Control
      Expert Settings

    Diffusion Discretization Control
      Expert Settings

Regions
   NAME
     REGION-FUNC     

Material Properties
  NAME
    Porosity: Uniform
      Value
       ... set porosity_dist=uniform
    Intrinsic Permeability: Anisotropic Uniform
      Horizontal
      Vertical
       ... set permeability_dist=uniform
    Capillary_Pressure: van Genuchten
      alpha
      m
      Sr
      ... set cpl_type=3, vals=[m,alpha*Pa/atm,Sr,0]
    Relative_Permeability
      ... assume Mualem, set kr_type=3, set vals=[cpl[0,2,3]]
    Regions_Assigned
    Density 
      ... set kp (def="Permeability Output File"), pp (def="Porosity Output File")

Phase Definitions
  Phase Properties
    Density: Uniform
      Density
    Viscosity: Uniform
      Viscosity
  Phase Components
    NAME
      Component Solutes
       ...sets group=Total

Boundary Conditions
  NAME
    Assigned Regions
    BC-FUNC
    Solute BC
      NAME
        SOLUTE-BC-FUNC
        Concentration Units

Initial Conditions
  NAME
    IC: Uniform Saturation, IC: Linear Saturation, IC: Uniform Pressure,IC: Linear Pressure, IC: Hydrostatic
    Assigned Regions
    IC-FUNC
    Solute IC
      NAME
        SOLUTE-IC-FUNC
        Concentration Units


Output
  Time Macros
    NAME
      Times, Start_Period_Stop

  Cycle Macros
    NAME
      Cycles
      Start_Period_Stop

  Visualization Data
    File Name Base
    Variables
    Cycle Macros, Time Macros, File Name Digits

  Checkpoint Data
    File Name Base
    Variables
    Cycle Macros, File Name Digits

  Observation Data
    Functional
    Region
    Observation Data: Integral, Observation Data: Point
    Time Macro
    Variables



REGION-FUNC:
------------

Region: Color Function
  File
  Value
   ... sets purpose=all, type=color_function

Region: Point
  Coordinate
   ... sets purpose=all, type=point

Region: Box
  Low Coordinate
  High Coordinate
   if (length(d)<geometry_eps) {
     ... sets type=surface, purpose=all
   }
   else {
     ... sets type=box, purpose=all
   }

Region: Plane
  Direction
  Location
   FIXME: ignores sign of direction
  ... sets type=surface, purpose=0-5 for orientation




BC-FUNC:
---------

BC: No Flow
  .. sets type=noflow

BC: Uniform Saturation
  ... sets type = saturation
  Values
  Times
  Time Functions

BC: Linear Saturation
  ... sets type = saturation
  Values
  Times
  Time Functions
  FIXME: no gradient info processed

BC: Linear Pressure
  ... sets type = pressure
  Reference Values
  Gradient Value
  Reference Coordinate

BC: Uniform Pressure
  Values
  Times
  Time Functions

BC: Flux
  Inward Volumetric Flux, Inward Mass Flux, Outward Volumetric Flux, Outward Mass Flux
  Times
  Time Functions
    ...sets aqueous_vol_flux, inflowtimes, inflowfncs, and type=zero_total_velocity

BC: Uniform Concentration
  Values
  Times
  Time Functions
  Molar Concentration, Molal Concentration, Specific Concentration
    ..sets type=concentration

BC: Zero Gradient or BC: Outflow
  ... sets type=outflow

BC: No Flow
  ... sets type=noflow


Additionally:

  Hard_set bcs in PMAMR (for sat, pressure, inflow)
    noflow:                  (4,4,0)
    pressure or hydrostatic: (1,2,0)
    saturation:              (1,2,0)
    zero_total_velocity:     (1,1,1)
    else:                    (4,4,4)




IC-FUNC:
--------

IC: Uniform Saturation
  Value (sets Water saturation value, type=saturation)

IC: Linear Saturation
  Value (sets Water saturation value, type=saturation)
  FIXME: no gradient info processed

IC: Uniform Pressure
  Phase
  Value
  ... sets type=hydrostatic, val, water_table_height

IC: Linear Pressure
  Phase
  Reference Value
  Gradient Value
  Reference Coordinate
  ... sets type=hydrostatic, val, grad, water_table_height

IC: Uniform Concentration
  Value
  Molar Concentration, Molal Concentration, Specific Concentration, 
   ...sets type=concentration


user_derives: "Material ID", "Grid ID", "Core ID", "Cell ID", "Capillary Pressure", "Volumetric Water Content",
              "Porosity", "Aqueous Saturation", "Aqueous Pressure", "Aqueous Volumetric Flux X",
              "Aqueous Volumetric Flux Y", "Aqueous Volumetric Flux Z", "Aqueous TRACER Concentration"


*/
