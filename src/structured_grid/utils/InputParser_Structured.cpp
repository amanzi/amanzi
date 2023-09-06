/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include <InputParser_Structured.H>
#include <PMAMR_Labels.H>

#include <BoxLib.H>

#include <algorithm>

using Teuchos::Array;
using Teuchos::ParameterList;
using Teuchos::ParameterEntry;

namespace Amanzi {
  namespace AmanziInput {

    std::map<std::string,std::string> GlobalData::AMR_to_Amanzi_label_map;
    std::map<std::string,std::string>& AMR_to_Amanzi_label_map = AMRToAmanziLabelMap();
    //std::map<std::string,SolidChem::SorptionIsothermData> SolidChem::sorption_isotherms; // One for each material, indexed on solute name


    static int require_static_velocity = -1; // <0 means it has not yet been set

    void MyAbort(const std::string& m) {
      if (Teuchos::GlobalMPISession::getRank() == 0) {
        std::cerr << m << std::endl;
        throw std::exception();
      }
    }

    void MyWarning(const std::string& m) {

      if (Teuchos::GlobalMPISession::getRank() == 0) {
        std::cerr << "WARNING::" << m << "!!!" << std::endl;
      }
    }

    void MyAbort(const Array<std::string>& m) {
      if (Teuchos::GlobalMPISession::getRank() == 0) {
        for (int i=0; i<m.size(); ++i) {
          std::cerr << m[i] << " ";
        }
        std::cerr << std::endl;
        throw std::exception();
      }
    }

    static double atmToMKS = 101325;
    static double gravity_mag_DEF = 9.807;
    static int gravity_dir_DEF = BL_SPACEDIM - 1;
#if BL_SPACEDIM == 2
    static double z_location_DEF = 0;
#endif

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

    std::string lc(const std::string& in)
    {
      std::string out;
      std::transform(in.begin(), in.end(), std::back_inserter(out), ::tolower);
      return out;
    }

    std::string uc(const std::string& in)
    {
      std::string out;
      std::transform(in.begin(), in.end(), std::back_inserter(out), ::toupper);
      return out;
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
      std::string ProbLo_str = "Domain Low Coordinate"; reqP.push_back(ProbLo_str);
      std::string ProbHi_str = "Domain High Coordinate"; reqP.push_back(ProbHi_str);
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
                                  bool&                do_tracer_advection,
                                  bool&                do_tracer_diffusion,
                                  bool&                do_chem)
    {
      std::string ec_str = "Execution Control";
      std::string tpc_str = "Time Period Control";
      std::string amr_str = "Adaptive Mesh Refinement Control";
      //EIB: this is the actual name of the list, is there a reason for a difference????
      //std::string amr_str = "Adaptive Mesh Refinement";
      //EIB: change due to 1.2.2 update
      //std::string prob_str = "Basic Algorithm Control";
      std::string str_str = "Structured Algorithm";
      std::string exp_str = "Expert Settings";
      std::string io_str = "IO Control";
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
      ParameterList& io_out_list      = struc_out_list.sublist("vismf");
      ParameterList& fabarr_out_list  = struc_out_list.sublist("fabarray");
      ParameterList& fab_out_list     = struc_out_list.sublist("fab");

      bool echo_inputs = false;
      std::string echo_str = "Echo Inputs";
      if (parameter_list.isParameter(echo_str)) {
        echo_inputs = parameter_list.get<bool>(echo_str);
      }
      struc_out_list.set<bool>("echo_inputs",echo_inputs);

      int prob_v, mg_v, cg_v, amr_v, diffuse_v, io_v, fab_v;
      //
      // Set flow model
      //
      std::string model_name;
      std::string flow_mode = ec_list.get<std::string>(flow_str);
      if (flow_mode == "Off") {
        model_name = "steady-saturated";
        prob_out_list.set("have_capillary",0);
        prob_out_list.set("cfl",1);
      }
      else if (flow_mode == "Richards") {
        model_name = "richards";
        prob_out_list.set("have_capillary",0);
        prob_out_list.set("cfl",-1);
      }
      else if (flow_mode == "Single Phase") {
        model_name = "saturated"; // Note: May later be overwritten as steady-saturated, depending on time integration mode
        prob_out_list.set("have_capillary",0);
        prob_out_list.set("cfl",-1);
      }
      else {
        MyAbort("\"" + flow_str + "\" = \"" + flow_mode + "\" not supported");
      }
      // prob_out_list.set("model_name",model_name); // NOTE: This now set AFTER discovering time integration mode

      //
      // Set transport model
      //
      std::string transport_mode = ec_list.get<std::string>(trans_str);
      if (transport_mode == "Off") {
        do_tracer_advection = 0;
      }
      else if (transport_mode == "On") {
        do_tracer_advection = 1;
      }
      else {
        MyAbort("\"" + trans_str + "\" = \"" + transport_mode + "\" not supported");
      }
      do_tracer_diffusion = false;

      //
      // Set chemistry model
      //
      std::string chem_mode = ec_list.get<std::string>(chem_mod_str);
      ParameterList chem_out_list;
      prob_out_list.set("chemistry_model",chem_mode);
      if (chem_mode == "Off") {
	// FIXME: Do we need chemistry info even if no reactions?
        do_chem = false;
      }
      else if (chem_mode == "Alquimia" || chem_mode == "Amanzi") {
        do_chem = true;
        const ParameterList& chem_list = parameter_list.sublist(chem_str);

        const std::string ThermoDB_str = "Thermodynamic Database";
        const std::string ThermoDB_Fmt_str = "Format";
        const std::string ThermoDB_File_str = "File";
        const std::string Chemistry_Engine_str = "Engine";
        const std::string Chemistry_Engine_Input_str = "Engine Input File";
        const std::string Chemistry_Verbosity_str = "Verbosity";
        const std::string Chemistry_Activity_Model_str = "Activity Model";
        const std::string Chemistry_Tol_str = "Tolerance";
        const std::string Chemistry_Newton_str = "Maximum Newton Iterations";
        const std::string Chemistry_Aux_str = "Auxiliary Data";
        const std::string Chemistry_Max_Step_str = "Max Time Step (s)";
        const std::string Chemistry_Max_Step_str_tr = "Max Time Step";

        reqP.clear(); reqL.clear();
        bool thermoDB_reqd = chem_mode == "Amanzi";
        if (thermoDB_reqd) {
          reqL.push_back(ThermoDB_str);
        }
        if (chem_mode == "Alquimia") {
          reqP.push_back(Chemistry_Engine_str);
          reqP.push_back(Chemistry_Engine_Input_str);
        }
        reqP.push_back(Chemistry_Tol_str);
        reqP.push_back(Chemistry_Newton_str);
        reqP.push_back(Chemistry_Max_Step_str);

        PLoptions CHopt(chem_list,reqL,reqP,true,false);
        const Array<std::string>& CHoptP = CHopt.OptParms();

        if (chem_list.isSublist(ThermoDB_str)) {
          const Teuchos::ParameterList& thermoPL = chem_list.sublist(ThermoDB_str);
          Array<std::string> treqL, treqP;
          treqP.push_back(ThermoDB_Fmt_str);
          treqP.push_back(ThermoDB_File_str);
          PLoptions ThDBopt(thermoPL,treqL,treqP,true,true);
          chem_out_list.set<std::string>(underscore(ThermoDB_str)+"_"+underscore(ThermoDB_Fmt_str), underscore(thermoPL.get<std::string>(ThermoDB_Fmt_str)));
          chem_out_list.set<std::string>(underscore(ThermoDB_str)+"_"+underscore(ThermoDB_File_str), thermoPL.get<std::string>(ThermoDB_File_str));
        }

        if (chem_mode=="Alquimia") {
          chem_out_list.set<std::string>(underscore(Chemistry_Engine_str), underscore(chem_list.get<std::string>(Chemistry_Engine_str)));
          chem_out_list.set<std::string>(underscore(Chemistry_Engine_Input_str), underscore(chem_list.get<std::string>(Chemistry_Engine_Input_str)));
        }
        for (int i=0; i<CHoptP.size(); ++i) {
          const std::string& name = CHoptP[i];
          std::string _name = underscore(name);
          if (name==Chemistry_Aux_str) {
            chem_out_list.set(_name,underscore(chem_list.get<Array<std::string> >(name)));
          }
          else if (name==Chemistry_Verbosity_str) {
            chem_out_list.set<std::string>(underscore(Chemistry_Verbosity_str), underscore(chem_list.get<std::string>(Chemistry_Verbosity_str)));
          }
          else if (name==Chemistry_Activity_Model_str) {
            chem_out_list.set<std::string>(underscore(Chemistry_Activity_Model_str), underscore(chem_list.get<std::string>(Chemistry_Activity_Model_str)));
          }
          else if (name==Chemistry_Tol_str) {
            chem_out_list.set<double>(underscore(Chemistry_Tol_str), chem_list.get<double>(Chemistry_Tol_str));
          }
          else if (name==Chemistry_Newton_str) {
            chem_out_list.set<int>(underscore(Chemistry_Newton_str), chem_list.get<int>(Chemistry_Newton_str));
          }
          else if (name==Chemistry_Max_Step_str) {
            chem_out_list.set<double>(underscore(Chemistry_Max_Step_str_tr), chem_list.get<double>(Chemistry_Max_Step_str));
          }
          else {
            MyAbort("Uncrecognized option in \"Chemistry\" parameter list: \"" + name + "\"");
          }
        }
      }
      else {
        MyAbort("\"" + chem_mod_str + "\" = \"" + chem_mode + "\" not supported");
      }
      struc_out_list.set(underscore(chem_str),chem_out_list);

      //
      // Set time evolution mode
      //
      std::string steady_str = "Steady";
      std::string transient_str = "Transient";
      std::string transient_with_static_flow_str = "Transient with Static Flow";
      std::string init_to_steady_str = "Initialize To Steady";
      const ParameterList& t_list = ec_list.sublist(tim_str);
      if (model_name == "single_phase") {
        prob_out_list.set("do_simple",2);
      }

      require_static_velocity = 0;
      if (t_list.isSublist(transient_str) || t_list.isSublist(transient_with_static_flow_str))
      {
        const std::string& sublist_name = (t_list.isSublist(transient_str) ? transient_str : transient_with_static_flow_str);
        if (sublist_name == transient_with_static_flow_str) {
          require_static_velocity = 1;
          if (model_name != "saturated" && model_name != "steady-saturated") {
            MyAbort("\"" + tim_str + "\" = \"" +  transient_with_static_flow_str + "\" only applicable to saturated flow" );
          }
          model_name = "steady-saturated";
        }
        const ParameterList& tran_list = t_list.sublist(sublist_name);
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
          prob_out_list.set<double>("transient_max_dt", dt_max);
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
        std::string Steady_Init_Time_Step_str = "Steady Initial Time Step";
        optPd[Steady_Init_Time_Step_str] = -1; // <0 means inactive
        std::string Transient_Init_Time_Step_str = "Transient Initial Time Step";
        optPd[Transient_Init_Time_Step_str] = -1; // <0 means inactive

        double strt_time = tran_list.get<double>(Start_str);
        double stop_time = tran_list.get<double>(End_str);
        double switch_time = tran_list.get<double>(Switch_str);
        struc_out_list.set<double>("strt_time", strt_time);
        struc_out_list.set<double>("stop_time", stop_time);
        struc_out_list.set<double>("switch_time", switch_time);
        prob_out_list.set<double>(underscore("steady_max_pseudo_time"),tran_list.get<double>(Switch_str));
        if (switch_time < stop_time) {
          prob_out_list.set<double>(underscore("dt_init"),tran_list.get<double>(Transient_Init_Time_Step_str));
        }
        prob_out_list.set<double>(underscore("steady_init_time_step"),tran_list.get<double>(Steady_Init_Time_Step_str));

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
            prob_out_list.set<double>("steady_max_dt", it->second);
          } else if (it->first==Transient_Max_Time_Step_Size_str) {
            prob_out_list.set<double>("transient_max_dt", it->second);
          } else if (it->first==Transient_Init_Time_Step_str) {
            prob_out_list.set<double>("dt_init", it->second);
          }
          else {
            prob_out_list.set<double>(underscore(it->first), it->second);
          }
        }
        for (std::map<std::string,int>::const_iterator it=optPi.begin(); it!=optPi.end(); ++it) {
          if (it->first==Max_Step_str) {
            int max_step = it->second;
            if (max_step<0) {
              std::cerr << "Negative value specified for \""+Max_Step_str+"\": " << max_step << std::endl;
              MyAbort("");
            }
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
        const ParameterList& steady_list = t_list.sublist(steady_str);
        reqP.clear(); reqL.clear();
        std::string Start_str = "Start"; reqP.push_back(Start_str);
        std::string End_str = "End"; reqP.push_back(End_str);
        std::string Init_Time_Step_str = "Initial Time Step"; reqP.push_back(Init_Time_Step_str);

        PLoptions Sopt(steady_list,reqL,reqP,true,false);
        struc_out_list.set<double>("strt_time", steady_list.get<double>(Start_str));
        struc_out_list.set<double>("stop_time", steady_list.get<double>(End_str));
        prob_out_list.set<double>(underscore("steady_init_time_step"),steady_list.get<double>(Init_Time_Step_str));

        struc_out_list.set<std::string>("execution_mode", "steady");


        // Set defaults for optional parameters
        std::map<std::string,double> optPd;
        std::map<std::string,int> optPi;
        std::string Max_Time_Step_Size_str = "Maximum Time Step Size";     optPd[Max_Time_Step_Size_str] = -1; // <0 means inactive
        std::string Time_Step_Grow_Max_str = "Maximum Time Step Grow";     optPd[Time_Step_Grow_Max_str] = -1; // <0 means inactive
        std::string Time_Step_Shrink_Max_str = "Maximum Time Step Shrink"; optPd[Time_Step_Shrink_Max_str] = -1; // <0 means inactive
	optPd[Init_Time_Step_str] = -1; // <0 means inactive

        // Extract optional parameters
        const Array<std::string> SoptP = Sopt.OptParms();
        for (int i=0; i<SoptP.size(); ++i) {
          if (optPd.find(SoptP[i]) != optPd.end()) {
            optPd[SoptP[i]] = steady_list.get<double>(SoptP[i]); // replace default value
          }
          else if (optPi.find(SoptP[i]) != optPi.end()) {
            optPi[SoptP[i]] = steady_list.get<int>(SoptP[i]); // replace default value
          }
          else {
            MyAbort("Unrecognized option under \""+steady_str+"\": \""+SoptP[i]+"\"" );
          }
        }

        for (std::map<std::string,double>::const_iterator it=optPd.begin(); it!=optPd.end(); ++it) {
          if (it->first==Time_Step_Grow_Max_str) {
            prob_out_list.set<double>("dt_grow_max", it->second);
          } else if (it->first==Time_Step_Shrink_Max_str) {
            prob_out_list.set<double>("dt_shrink_max", it->second);
          } else if (it->first==Max_Time_Step_Size_str) {
            prob_out_list.set<double>("steady_max_dt", it->second);
          } else if (it->first==Init_Time_Step_str) {
            prob_out_list.set<double>("dt_init", it->second);
          }
          else {
            prob_out_list.set<double>(underscore(it->first), it->second);
          }
        }

        struc_out_list.set<int>("max_step", 1); // For steady flow, we take 1 step
        for (std::map<std::string,int>::const_iterator it=optPi.begin(); it!=optPi.end(); ++it) {
          prob_out_list.set<int>(underscore(it->first), it->second);
        }
      }
      else {
        std::cerr << t_list << std::endl;
        MyAbort("No recognizable value for \"" + tim_str + "\"");
      }

      // NOTE: Can now set this AFTER discovering time integration mode
      prob_out_list.set("model_name",model_name);

      // Generate a useful warning

      double start_time = struc_out_list.get<double>("strt_time");
      double stop_time = struc_out_list.get<double>("stop_time");
      if (stop_time <= start_time) {
        MyWarning("End <= Start, code will halt immediately after initializing data");
      }

      // Deal with optional settings
      const Array<std::string> optL = ECopt.OptLists();
      const Array<std::string> optP = ECopt.OptParms();

      //
      // Basic Algorithm Defaults
      //
      double visc_abs_tol = 1.e-12; prob_out_list.set<double>("visc_abs_tol",visc_abs_tol);
      double visc_tol = 1.e-12; prob_out_list.set<double>("visc_tol",visc_tol);

      //
      // Verbosity Default
      //
      std::string v_val = "Medium";

      //
      // AMR gridding control defaults
      //
      int num_levels = 1;
      int max_level = num_levels-1;
      bool do_amr_subcycling = true;
      int ref_ratio_DEF = 2;
      Array<int> ref_ratio(max_level,ref_ratio_DEF);
      int regrid_int_DEF = 2;
      Array<int> regrid_int(num_levels,regrid_int_DEF);
      int blocking_factor_DEF = 1;
      Array<int> blocking_factor(num_levels,blocking_factor_DEF);
      int n_err_buf_DEF = 1;
      Array<int> n_err_buf(max_level,n_err_buf_DEF);
      int max_grid_DEF = (ndim==2 ? 128  :  32);
      Array<int> max_grid(num_levels,max_grid_DEF);

      //
      // Optional parameters
      //
      for (int i=0; i<optP.size(); ++i) {
        if (optP[i] == v_str) {
          //
          // Verbosity level
          //
          v_val = ec_list.get<std::string>(v_str);
        }
        else {
          MyAbort("Unrecognized optional parameter to \"" + ec_str + "\": \"" + optP[i] + "\"");
        }
      }

      //
      // Verbosity implementation
      //
      if (lc(v_val) == "none") {
        prob_v = 0; mg_v = 0; cg_v = 0; amr_v = 0; diffuse_v = 0; io_v = 0; fab_v = 0;
      }
      else if (lc(v_val) == "low") {
        prob_v = 1; mg_v = 0; cg_v = 0; amr_v = 1;  diffuse_v = 0; io_v = 0; fab_v = 0;
      }
      else if (lc(v_val) == "medium") {
        prob_v = 1; mg_v = 0; cg_v = 0; amr_v = 2;  diffuse_v = 0; io_v = 0; fab_v = 0;
      }
      else if (lc(v_val) == "high") {
        prob_v = 2; mg_v = 0; cg_v = 0; amr_v = 3;  diffuse_v = 0; io_v = 0; fab_v = 0;
      }
      else if (lc(v_val) == "extreme") {
        prob_v = 3; mg_v = 2; cg_v = 2; amr_v = 3;  diffuse_v = 1; io_v = 1; fab_v = 1;
      }
      //
      // Optional lists
      //
      const std::string restart_str = "Restart";

      for (int i=0; i<optL.size(); ++i)
      {

        if (optL[i] == restart_str) {

          const std::string restart_file_str = "File Name";
          const ParameterList& restart_list = ec_list.sublist(restart_str);
          Array<std::string> nL, nP;
          nP.push_back(restart_file_str);
          PLoptions Resopt(restart_list,nL,nP,true,true);
          amr_out_list.set<std::string>("restart",restart_list.get<std::string>(restart_file_str));

        }
        else if (optL[i] == tpc_str)
        {
          const std::string tpc_start_times_str = "Start Times";
          const std::string tpc_initial_time_steps_str = "Initial Time Step";
          const std::string tpc_initial_time_step_multipliers_str = "Initial Time Step Multiplier";
          const std::string tpc_maximum_time_steps_str = "Maximum Time Step";
          const ParameterList& tpc_list = ec_list.sublist(tpc_str);
          Array<std::string> nL, nP;
          PLoptions TPCopt(tpc_list,nL,nP,true,false);
          const Array<std::string> TPCoptP = TPCopt.OptParms();

          Array<double> tpc_start_times, tpc_initial_time_steps, tpc_initial_time_step_multipliers, tpc_maximum_time_steps;
          for (int j=0; j<TPCoptP.size(); ++j)
          {
            const std::string& name = TPCoptP[j];
            if (name == tpc_start_times_str) {
              tpc_start_times = tpc_list.get<Array<double> >(tpc_start_times_str);
            }
            else if (name == tpc_initial_time_steps_str) {
              tpc_initial_time_steps = tpc_list.get<Array<double> >(tpc_initial_time_steps_str);
            }
            else if (name == tpc_initial_time_step_multipliers_str) {
              tpc_initial_time_step_multipliers = tpc_list.get<Array<double> >(tpc_initial_time_step_multipliers_str);
            }
            else if (name == tpc_maximum_time_steps_str) {
              tpc_maximum_time_steps = tpc_list.get<Array<double> >(tpc_maximum_time_steps_str);
            }
          }

          int num_periods = tpc_start_times.size();
          if (tpc_initial_time_steps.size() != num_periods) {
            if (tpc_initial_time_steps.size() != 0) {
              MyAbort("If specified, number of \""+tpc_initial_time_steps_str+"\" entries must equal number of \""+tpc_start_times_str+"\" entries");
            }
            else {
              tpc_initial_time_steps.resize(num_periods,-1);
            }
          }
          if (tpc_initial_time_step_multipliers.size() != num_periods) {
            if (tpc_initial_time_step_multipliers.size() != 0) {
              MyAbort("If specified, number of \""+tpc_initial_time_step_multipliers_str+"\" entries must equal number of \""+tpc_start_times_str+"\" entries");
            }
            else {
              tpc_initial_time_step_multipliers.resize(num_periods,1);
            }
          }
          if (tpc_maximum_time_steps.size() != num_periods) {
            if (tpc_maximum_time_steps.size() != 0) {
              MyAbort("If specified, number of \""+tpc_maximum_time_steps_str+"\" entries must equal number of \""+tpc_start_times_str+"\" entries");
            }
            else {
              tpc_maximum_time_steps.resize(num_periods,-1);
            }
          }

          if (num_periods>0) {
            prob_out_list.set<Array<double> >(underscore("TPC "+tpc_start_times_str), tpc_start_times);
            prob_out_list.set<Array<double> >(underscore("TPC "+tpc_initial_time_steps_str), tpc_initial_time_steps);
            prob_out_list.set<Array<double> >(underscore("TPC "+tpc_initial_time_step_multipliers_str), tpc_initial_time_step_multipliers);
            prob_out_list.set<Array<double> >(underscore("TPC "+tpc_maximum_time_steps_str), tpc_maximum_time_steps);
          }
        }
        else if (optL[i] == num_str)
        {
	  const ParameterList& num_list = ec_list.sublist(num_str);
          Array<std::string> nL, nP;
          PLoptions NUMopt(num_list,nL,nP,false,true);
          const Array<std::string> NUMoptL = NUMopt.OptLists();

	  std::set<std::string> ncp_allowed_lists;
	  ncp_allowed_lists.insert(str_str);
          for (int j=0; j<NUMoptL.size(); ++j)
          {
	    if (ncp_allowed_lists.find(NUMoptL[j]) == ncp_allowed_lists.end()) {
	      std::string str = "Unrecognized parameter list to \"" + num_str + "\" list: \"" + NUMoptL[j] + "\".\n";
	      str += "Acceptable lists:\n";
	      for (std::set<std::string>::const_iterator it=ncp_allowed_lists.begin(); it!=ncp_allowed_lists.end(); ++it) {
		str += "\t\"" + (*it) + "\"\n";
	      }
	      MyAbort(str);
	    }

	    if (NUMoptL[j] == str_str) {

	      Array<std::string> sL, sP;
	      const ParameterList& str_list = num_list.sublist(str_str);
	      PLoptions STRopt(str_list,sL,sP,false,true);

	      std::set<std::string> str_allowed_lists;
	      str_allowed_lists.insert(amr_str);
	      str_allowed_lists.insert(it_str);
	      str_allowed_lists.insert(diffuse_str);
	      str_allowed_lists.insert(exp_str);

	      const Array<std::string> STRoptL = STRopt.OptLists();

	      for (int k=0; k<STRoptL.size(); ++k)
	      {
		if (str_allowed_lists.find(STRoptL[k]) == str_allowed_lists.end()) {
		  std::string str = "Unrecognized parameter list to \"" + str_str + "\" list: \"" + STRoptL[k] + "\".\n";
		  str += "Acceptable lists:\n";
		  for (std::set<std::string>::const_iterator it=str_allowed_lists.begin(); it!=str_allowed_lists.end(); ++it) {
		    str += "\t\"" + (*it) + "\"\n";
		  }
		  MyAbort(str);
		}
		//
		// AMR Options
		//
		if (STRoptL[k] == amr_str)
		{
		  const ParameterList& amr_list = str_list.sublist(amr_str);
		  std::string amr_subcycling_str = "Do AMR Subcycling";
		  if (amr_list.isParameter(amr_subcycling_str)) {
		    do_amr_subcycling = amr_list.get<bool>(amr_subcycling_str);
		  }

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

		  for (int kk=0; kk<max_level; ++kk) {
		    if (ref_ratio[kk] != 2 && ref_ratio[kk]!=4) {
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
		      if (regrid_int.size() != max_level && regrid_int.size() != 1) {
			MyAbort("Must provide either a single regridding interval for all refinement levels, or one for each");
		      }
		    }
		    else {
		      if (regrid_int.size() != 1) {
			MyAbort("Subcycling is disabled, only a single regridding interval is supported");
		      }
		    }

		    for (int kk=0; kk<regrid_int.size(); ++kk) {
		      if (regrid_int[kk] <= 0) {
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
		  for (int kk=0; kk<blocking_factor.size(); ++kk) {
		    double twoPower = std::log(blocking_factor[kk])/std::log(2);
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
		  for (int kk=0; kk<n_err_buf.size(); ++kk) {
		    if (n_err_buf[kk] < 0) {
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
		  for (int kk=0; kk<max_grid.size(); ++kk) {
		    if (max_grid[kk] < blocking_factor[kk]) {
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
		} // End of AMR option processing
	      } // End of Structured Options optional lists processing
	    } // End of Structured Options processing
	  } // End of Numerical Control Parameters optional lists processing
	} // End of Numerical Control Parameters processing
      } // End of Execution Controls optional lists processing

      Array<int> n_cell = amr_out_list.get<Array<int> >("n_cell");

      if (blocking_factor.size() > 0) {
        for (int n=0;n<ndim;n++) {
          if (n_cell[n]%blocking_factor[0] > 0) {
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
      std::string amr_subcycling_mode = ( do_amr_subcycling ? "Auto" : "None");
      amr_out_list.set<std::string>("subcycling_mode", amr_subcycling_mode);

      amr_out_list.set("v",amr_v);
      mg_out_list.set("v",mg_v);
      cg_out_list.set("v",cg_v);
      prob_out_list.set("v",prob_v);
      io_out_list.set("v",io_v);
      fabarr_out_list.set("verbose",fab_v);

      for (int i=0; i<optL.size(); ++i)
      {
        if (optL[i] == num_str)
        {
          const ParameterList& num_list = ec_list.sublist(num_str);
          Array<std::string> nL, nP;
          PLoptions NUMopt(num_list,nL,nP,false,true);
          const Array<std::string> NUMoptL = NUMopt.OptLists();

	  std::set<std::string> ncp_allowed_lists;
	  ncp_allowed_lists.insert(str_str);
          //
          // Now process expert lists to overwrite any settings directly
          //
          for (int j=0; j<NUMoptL.size(); ++j)
          {
	    if (ncp_allowed_lists.find(NUMoptL[j]) == ncp_allowed_lists.end()) {
	      std::string str = "Unrecognized parameter list to \"" + num_str + "\" list: \"" + NUMoptL[j] + "\".\n";
	      str += "Acceptable lists:\n";
	      for (std::set<std::string>::const_iterator it=ncp_allowed_lists.begin(); it!=ncp_allowed_lists.end(); ++it) {
		str += "\t\"" + (*it) + "\"\n";
	      }
	      MyAbort(str);
	    }

	    if (NUMoptL[j] == str_str) {

	      Array<std::string> sL, sP;
	      const ParameterList& str_list = num_list.sublist(str_str);
	      PLoptions STRopt(str_list,sL,sP,false,true);

	      std::set<std::string> str_allowed_lists;
	      str_allowed_lists.insert(amr_str);
	      str_allowed_lists.insert(it_str);
	      str_allowed_lists.insert(diffuse_str);
	      str_allowed_lists.insert(exp_str);

	      const Array<std::string> STRoptL = STRopt.OptLists();

	      for (int k=0; k<STRoptL.size(); ++k)
	      {
		if (str_allowed_lists.find(STRoptL[k]) == str_allowed_lists.end()) {
		  std::string str = "Unrecognized parameter list to \"" + str_str + "\" list: \"" + STRoptL[k] + "\".\n";
		  str += "Acceptable lists:\n";
		  for (std::set<std::string>::const_iterator it=str_allowed_lists.begin(); it!=str_allowed_lists.end(); ++it) {
		    str += "\t\"" + (*it) + "\"\n";
		  }
		  MyAbort(str);
		}

		if (STRoptL[k] == exp_str)
		{
		  ParameterList p1_out_list;
		  process_expert_options(str_list,p1_out_list);
		  prob_out_list.setParameters(p1_out_list);
		}
		else if (STRoptL[k] == amr_str)
		{
		  process_expert_options(str_list.sublist(amr_str),amr_out_list);
		}
		else if (STRoptL[k] == it_str)
		{
		  const ParameterList& it_list = str_list.sublist(it_str);
		  Array<std::string> iL, iP;
		  PLoptions ITopt(it_list,iL,iP,false,true);

		  std::set<std::string> it_allowed_lists;
		  it_allowed_lists.insert(mg_str);
		  it_allowed_lists.insert(cg_str);

		  const Array<std::string> IToptL = ITopt.OptLists();

		  for (int L=0; L<IToptL.size(); ++L)
		  {
		    if (it_allowed_lists.find(IToptL[L]) == it_allowed_lists.end()) {
		      std::string str = "Unrecognized parameter list to \"" + it_str + "\" list: \"" + IToptL[L] + "\".\n";
		      str += "Acceptable lists:\n";
		      for (std::set<std::string>::const_iterator it=it_allowed_lists.begin(); it!=it_allowed_lists.end(); ++it) {
			str += "\t\"" + (*it) + "\"\n";
		      }
		      MyAbort(str);
		    }
		    else {
		      if (IToptL[L] == mg_str)
		      {
			process_expert_options(it_list.sublist(mg_str),mg_out_list);
		      }
		      else if (IToptL[L] == cg_str)
		      {
			process_expert_options(it_list.sublist(cg_str),cg_out_list);
		      }
		    }
		  } // End processing IT optional lists
		} // End processing IT parameter list
	      } // End processing Structured Algorithm optional lists
	    } // End processing Structured Algorithm parameter list
	  } // End processing Numerical Control Parameters optional lists
	} // End processing Numerical Control Parameters parameter list
      } // End of Execution Controls optional lists processing
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

      std::string purpose="all", type="box";
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

    void convert_Region_Logical(const ParameterList& rslist,
                                const std::string&   rlabel,
                                ParameterList&       rsublist)
    {
      const std::string Op_str = "Operation";
      const std::string Comp_str = "Complement";
      const std::string Union_str = "Union";
      const std::string Sub_str = "Subtraction";
      const std::string Isect_str = "Intersection";
      const std::string CompReg_str = "Region";
      const std::string Regs_str = "Regions";

      const ParameterList& rsslist = rslist.sublist(rlabel);

      if (!rsslist.isParameter(Op_str)) {
        MyAbort("Keyword \""+Op_str+"\" required for \""+rlabel + "\"");
      }
      const std::string& op = rsslist.get<std::string>(Op_str);
      if (op==Comp_str) {
        if (!rsslist.isParameter(CompReg_str) && !rsslist.isParameter(Regs_str)) {
          MyAbort("\""+CompReg_str+"\" or \""+Regs_str+"\"required with "+"\""+Op_str+"\" = \"" + Comp_str + "\" "
                  "for \"" + rlabel + "\"");
        }
	if (rsslist.isParameter(CompReg_str)) {
	  const std::string& comp_reg = rsslist.get<std::string>(CompReg_str);
	  rsublist.set("operation","complement");
	  rsublist.set("region",underscore(comp_reg));
	} else {
	  const std::string& comp_reg = rsslist.get<Array<std::string> >(Regs_str)[0];
	  rsublist.set("operation","complement");
	  rsublist.set("region",underscore(comp_reg));
	}
      }
      else if (op==Union_str) {
        if (!rsslist.isParameter(Regs_str)) {
          MyAbort("\""+Regs_str+"\" required with "+"\""+Op_str+"\" = \"" + Union_str + "\" "
                  "for \"" + rlabel + "\"");
        }
        const Array<std::string>& u_reg = rsslist.get<Array<std::string> >(Regs_str);
        rsublist.set("operation","union");
        rsublist.set("regions",underscore(u_reg));
      }
      else if (op==Sub_str) {
        if (!rsslist.isParameter(Regs_str)) {
          MyAbort("\""+Regs_str+"\" required with "+"\""+Op_str+"\" = \"" + Sub_str + "\" "
                  "for \"" + rlabel + "\"");
        }
        const Array<std::string>& s_reg = rsslist.get<Array<std::string> >(Regs_str);
        rsublist.set("operation","subtraction");
        rsublist.set("regions",underscore(s_reg));
      }
      else if (op==Isect_str) {
        if (!rsslist.isParameter(Regs_str)) {
          MyAbort("\""+Regs_str+"\" required with "+"\""+Op_str+"\" = \"" + Isect_str + "\" "
                  "for \"" + rlabel + "\"");
        }
        const Array<std::string>& s_reg = rsslist.get<Array<std::string> >(Regs_str);
        rsublist.set("operation","intersection");
        rsublist.set("regions",underscore(s_reg));
      }
      else {
        MyAbort("Unknown \""+Op_str+"\" for \"" + rlabel + "\" (\""+op+"\")");
      }

      int iPurpose = 6;
      std::string purpose = underscore(PMAMR::RpurposeDEF[iPurpose]);
      rsublist.set("purpose",purpose);
      rsublist.set("type", "logical");
    }

#if BL_SPACEDIM == 2
    void convert_Region_Polygon(const ParameterList& rslist,
                                const std::string&   rlabel,
                                ParameterList&       rsublist)
    {
      const std::string V1_str = "VerticesV1";
      const std::string V2_str = "VerticesV2";
      const std::string Extent_str = "Extent";
      const std::string Plane_str = "Plane";
      const std::string Reg_str = "Region: Polygon";

      const ParameterList& rsslist = rslist.sublist(rlabel);
      const Array<double>& verticesV1 = rsslist.get<Array<double> >(V1_str);
      const Array<double>& verticesV2 = rsslist.get<Array<double> >(V2_str);
      if (verticesV1.size() != verticesV2.size()) {
        MyAbort("\""+V1_str+"\", and \""+V2_str+"\" must have same number of points"+
                "for \""+Reg_str+"\" region \"" + rlabel);
      }

      rsublist.set<Array<double> >("v1",verticesV1);
      rsublist.set<Array<double> >("v2",verticesV2);
      rsublist.set("type", "polygon");
      int iPurpose = 6;
      std::string purpose = underscore(PMAMR::RpurposeDEF[iPurpose]);
      rsublist.set("purpose",purpose);
    }
    void convert_Region_Ellipse(const ParameterList& rslist,
                                const std::string&   rlabel,
                               ParameterList&       rsublist)
    {
      const std::string C_str = "Center";
      const std::string R_str = "Radius";
      const std::string Reg_str = "Region: Ellipse";

      const ParameterList& rsslist = rslist.sublist(rlabel);
      const Array<double>& c = rsslist.get<Array<double> >(C_str);
      if (c.size() < BL_SPACEDIM) {
        MyAbort("\""+C_str+"\" length must match problem dimensions "+
                "for \""+Reg_str+"\" region \"" + rlabel);
      }

      rsublist.set<Array<double> >("center",c);

      const Array<double>& r = rsslist.get<Array<double> >(R_str);
      if (r.size() < BL_SPACEDIM) {
        MyAbort("\""+R_str+"\" length must match problem dimensions "+
                "for \""+Reg_str+"\" region \"" + rlabel);
      }

      rsublist.set<Array<double> >("radius",r);
      rsublist.set("type", "ellipse");
      int iPurpose = 6;
      std::string purpose = underscore(PMAMR::RpurposeDEF[iPurpose]);
      rsublist.set("purpose",purpose);
    }
#else
    void convert_Region_SweptPolygon(const ParameterList& rslist,
                                     const std::string&   rlabel,
                                     ParameterList&       rsublist)
    {
      const std::string V1_str = "VerticesV1";
      const std::string V2_str = "VerticesV1";
      const std::string Extent_str = "Extent";
      const std::string Plane_str = "Plane";
      const std::string Reg_str = "Region: Swept Polygon";

      const ParameterList& rsslist = rslist.sublist(rlabel);
      const Array<double>& verticesV1 = rsslist.get<Array<double> >(V1_str);
      const Array<double>& verticesV2 = rsslist.get<Array<double> >(V2_str);
      if (verticesV1.size() != verticesV2.size()) {
        MyAbort("\""+V1_str+"\", and \""+V2_str+"\" must have same number of points"+
                "for \""+Reg_str+"\" region \"" + rlabel);
      }
      const std::string plane = rsslist.get<std::string>(Plane_str);
      if (plane != "XY" && plane != "YZ" && plane != "XZ") {
        MyAbort("\""+Plane_str+"\" must be \"XY\", \"YZ\" or \"XZ\" "+
                "for \""+Reg_str+"\" region \"" + rlabel);
      }
      const Array<double>& extent = rsslist.get<Array<double> >(Extent_str);
      if (extent.size() != 2) {
        MyAbort("\""+Extent_str+"\" must be 2-element Array of doubles"+
                "for \""+Reg_str+"\" region \"" + rlabel);
      }

      rsublist.set<Array<double> >("v1",verticesV1);
      rsublist.set<Array<double> >("v2",verticesV2);
      rsublist.set<std::string>("plane",plane);
      rsublist.set("type", "swept_polygon");
      int iPurpose = 6;
      std::string purpose = underscore(PMAMR::RpurposeDEF[iPurpose]);
      rsublist.set("purpose",purpose);
    }

    void convert_Region_RotatedPolygon(const ParameterList& rslist,
                                       const std::string&   rlabel,
                                       ParameterList&       rsublist)
    {
      const std::string V1_str = "VerticesV1";
      const std::string V2_str = "VerticesV2";
      const std::string Axis_str = "Axis";
      const std::string Plane_str = "Plane";
      const std::string Ref_str = "Reference Point";
      const std::string Reg_str = "Region: Rotated Polygon";

      const ParameterList& rsslist = rslist.sublist(rlabel);
      const Array<double>& verticesV1 = rsslist.get<Array<double> >(V1_str);
      const Array<double>& verticesV2 = rsslist.get<Array<double> >(V2_str);
      if (verticesV1.size() != verticesV2.size()) {
        MyAbort("\""+V1_str+"\", and \""+V2_str+"\" must have same number of points"+
                "for \""+Reg_str+"\" region \"" + rlabel);
      }
      const std::string plane = rsslist.get<std::string>(Plane_str);
      if (plane != "XY" && plane != "YZ" && plane != "XZ") {
        MyAbort("\""+Plane_str+"\" must be \"XY\", \"YZ\" or \"XZ\" "+
                "for \""+Reg_str+"\" region \"" + rlabel);
      }

      const std::string axis = rsslist.get<std::string>(Axis_str);
      if (axis != "X" && axis != "Y" && axis != "Z") {
        MyAbort("\""+Axis_str+"\" must be \"X\", \"Y\" or \"Z\" "+
                "for \""+Reg_str+"\" region \"" + rlabel);
      }

      const Array<double>& refPt = rsslist.get<Array<double> >(Ref_str);
      if (refPt.size() != BL_SPACEDIM) {
        const std::string d = BL_SPACEDIM==2 ? "2" : "3";
        MyAbort("\""+Ref_str+"\" must be a "+d+"-element Array of doubles"+
                "for \""+Reg_str+"\" region \"" + rlabel);
      }

      rsublist.set<Array<double> >("v1",verticesV1);
      rsublist.set<Array<double> >("v2",verticesV2);
      rsublist.set<std::string>("plane",plane);
      rsublist.set<std::string>("axis",axis);
      rsublist.set("type", "rotated_polygon");
      rsublist.set<Array<double> >("reference_pt",refPt);
      int iPurpose = 6;
      std::string purpose = underscore(PMAMR::RpurposeDEF[iPurpose]);
      rsublist.set("purpose",purpose);
    }
#endif

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
      rsublist.set(PMAMR::RlabelDEF[6],t3PL);
      def_regionNames.push_back(PMAMR::RlabelDEF[6]);

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
          rsublist.set(PMAMR::RlabelDEF[i],tmp);
          def_regionNames.push_back(PMAMR::RlabelDEF[i]);
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
          if (Teuchos::GlobalMPISession::getRank() == 0) {
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

          bool fail = false;
          if (rentry.isList()) {
            if (rlabel=="Region: Color Function") {
              convert_Region_ColorFunction(rslist,rlabel,rsublist);
            }
            else if (rlabel=="Region: Point") {
              convert_Region_Point(rslist,rlabel,rsublist);
            }
            else if (rlabel=="Region: Box") {
              convert_Region_Box(rslist,rlabel,rsublist);
            }
            else if (rlabel=="Region: Plane") {
              convert_Region_Plane(rslist,rlabel,rsublist);
            }
            else if (rlabel=="Region: Logical") {
              convert_Region_Logical(rslist,rlabel,rsublist);
            }
#if BL_SPACEDIM==2
            else if (rlabel=="Region: Polygon") {
              convert_Region_Polygon(rslist,rlabel,rsublist);
            }
            else if (rlabel=="Region: Ellipse") {
              convert_Region_Ellipse(rslist,rlabel,rsublist);
            }
#else
            else if (rlabel=="Region: Swept Polygon") {
              convert_Region_SweptPolygon(rslist,rlabel,rsublist);
            }
            else if (rlabel=="Region: Rotated Polygon") {
              convert_Region_RotatedPolygon(rslist,rlabel,rsublist);
            }
#endif
            else {
              fail = true;
            }
          }
          else {
            fail = true;
          }

          if (fail) {
            std::cerr << rlabel << " is not a valid region type for structured.\n";
            throw std::exception();
          }

          geom_list.set(_label,rsublist);
          // need to remove empty spaces
          arrayregions.push_back(_label);
        }
      }
      geom_list.set("regions",arrayregions);
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

    void convert_PorosityUniform(const ParameterList& fPLin,
                                 ParameterList&       fPLout)
    {
      fPLout.set<std::string>("distribution_type","uniform");
      Array<std::string> nullList, reqP;
      if (fPLin.isParameter("Values")) {
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
      else if (fPLin.isParameter("Value")) {
        const std::string val_name="Value"; reqP.push_back(val_name);
        PLoptions opt(fPLin,nullList,reqP,true,true);
        double val = fPLin.get<double>(val_name);
        fPLout.set<double>("vals",val);
      }
      else {
        std::string str = "Unrecognized porosity function parameters";
        std::cerr << fPLin << std::endl;
        BoxLib::Abort(str.c_str());
      }
    }

    void convert_PorosityGSLib(const ParameterList& fPLin,
                               ParameterList&       fPLout)
    {
      fPLout.set<std::string>("distribution_type","gslib");
      Array<std::string> nullList, reqP;
      const std::string param_file_name="GSLib Parameter File";
      const std::string data_file_name="GSLib Data File"; reqP.push_back(data_file_name);
      const std::string value_name="Value"; reqP.push_back(value_name);
      PLoptions opt(fPLin,nullList,reqP,false,false);

      fPLout.set<std::string>("gslib_data_file",fPLin.get<std::string>(data_file_name));
      fPLout.set<double>("val",fPLin.get<double>(value_name));
      const Array<std::string>& optP = opt.OptParms();
      for (int i=0; i<optP.size(); ++i) {
	if (optP[i] == param_file_name) {
	  fPLout.set<std::string>("gslib_param_file",fPLin.get<std::string>(param_file_name));
	}
	else {
          std::string str = "Unrecognized \"Porosity: GSLib\" option";
          std::cerr << fPLin << std::endl;
          BoxLib::Abort(str.c_str());
	}
      }
    }

    void convert_PermeabilityGSLib(const ParameterList& fPLin,
				   ParameterList&       fPLout)
    {
      fPLout.set<std::string>("distribution_type","gslib");
      Array<std::string> nullList, reqP;
      const std::string param_file_name="GSLib Parameter File";
      const std::string data_file_name="GSLib Data File"; reqP.push_back(data_file_name);
      const std::string value_name="Value"; reqP.push_back(value_name);
      PLoptions opt(fPLin,nullList,reqP,false,false);

      fPLout.set<std::string>("gslib_data_file",fPLin.get<std::string>(data_file_name));

      fPLout.set<double>("val",fPLin.get<double>(value_name) * 1.01325e15); // convert from m^2 to mDa
      const Array<std::string>& optP = opt.OptParms();
      for (int i=0; i<optP.size(); ++i) {
	if (optP[i] == param_file_name) {
	  fPLout.set<std::string>("gslib_param_file",fPLin.get<std::string>(param_file_name));
	}
	else {
          std::string str = "Unrecognized \"Permeability: GSLib\" option";
          std::cerr << fPLin << std::endl;
          BoxLib::Abort(str.c_str());
	}
      }
    }

    void convert_PermeabilityAnisotropic(const ParameterList& fPLin,
                                         ParameterList&       fPLout,
                                         double               scale)
    {
      /* Handle isotropic and anisotropic permeabilities here */
      fPLout.set<std::string>("distribution_type","uniform");
      Array<std::string> nullList, reqP;
      const std::string vertical_str = (BL_SPACEDIM==3 ? "z" : "y");
      const std::string horizontal_str = "x";
#if BL_SPACEDIM==3
      const std::string horizontal1_str = "y";
#endif
      const std::string uniform_value_str = "Value";
      const std::string values_str = "Values";
      const std::string times_str = "Times";
      const std::string forms_str = "Time Functions";

      if (fPLin.isParameter(uniform_value_str)) {
        /* Isotropic */
        Array<double> local_val(1); local_val[0] = fPLin.get<double>(uniform_value_str);
        local_val[0] *= 1.01325e15 * scale; // convert from m^2 to mDa
        ParameterList vlist, hlist, h1list;
        vlist.set<Array<double> >("vals",local_val);
        fPLout.set("vertical",vlist);
        hlist.set<Array<double> >("vals",local_val);
        fPLout.set("horizontal",hlist);
        h1list.set<Array<double> >("vals",local_val);
        fPLout.set("horizontal1",h1list);
      } else {
        /* Ansisotropic */
        Array<double> local_vvals(1);
        ParameterList vlist;
        if (fPLin.isSublist(vertical_str)) {
          const ParameterList& vPLin = fPLin.sublist(vertical_str);
          reqP.push_back(values_str);
          reqP.push_back(times_str);
          reqP.push_back(forms_str);
          PLoptions opt(vPLin,nullList,reqP,true,true);
          local_vvals = vPLin.get<Array<double> >(values_str);
          if (local_vvals.size()>1) {
            vlist.set<Array<double> >("times",vPLin.get<Array<double> >(times_str));
            vlist.set<Array<std::string> >("forms",vPLin.get<Array<std::string> >(forms_str));
          }
        }
        else if (fPLin.isParameter(vertical_str)) {
          local_vvals[0] = fPLin.get<double>(vertical_str);
        } else {
          std::string str = "Unrecognized vertical permeability function parameters";
          std::cerr << fPLin << std::endl;
          BoxLib::Abort(str.c_str());
        }

        for (int k=0; k<local_vvals.size(); k++) {
          local_vvals[k] *= 1.01325e15 * scale; // convert from m^2 to mDa
        }
        vlist.set<Array<double> >("vals",local_vvals);
        fPLout.set("vertical",vlist);

        Array<double> local_hvals(1); // local copy because we need to convert from m^2 to mD here
        ParameterList hlist;
        if (fPLin.isSublist(horizontal_str)) {
          const ParameterList& hPLin = fPLin.sublist(horizontal_str);
          reqP.push_back(values_str);
          reqP.push_back(times_str);
          reqP.push_back(forms_str);
          PLoptions opt(hPLin,nullList,reqP,true,true);
          local_hvals = hPLin.get<Array<double> >(values_str);
          if (local_hvals.size()>1) {
            hlist.set<Array<double> >("times",hPLin.get<Array<double> >(times_str));
            hlist.set<Array<std::string> >("forms",hPLin.get<Array<std::string> >(forms_str));
          }
        }
        else if (fPLin.isParameter(horizontal_str)) {
          local_hvals[0] = fPLin.get<double>(horizontal_str);
        }
        else {
          std::string str = "Unrecognized horizontal permeability function parameters";
          std::cerr << fPLin << std::endl;
          BoxLib::Abort(str.c_str());
        }

        for (int k=0; k<local_hvals.size(); k++) {
          local_hvals[k] *= 1.01325e15 * scale; // convert from m^2 to mDa
        }
        hlist.set<Array<double> >("vals",local_hvals);
        fPLout.set("horizontal",hlist);

#if BL_SPACEDIM==3
        Array<double> local_h1vals(1); // local copy because we need to convert from m^2 to mD here
        ParameterList h1list;
        if (fPLin.isSublist(horizontal1_str)) {
          const ParameterList& h1PLin = fPLin.sublist(horizontal1_str);
          reqP.push_back(values_str);
          reqP.push_back(times_str);
          reqP.push_back(forms_str);
          PLoptions opt(h1PLin,nullList,reqP,true,true);
          local_h1vals = h1PLin.get<Array<double> >(values_str);
          if (local_h1vals.size()>1) {
            h1list.set<Array<double> >("times",h1PLin.get<Array<double> >(times_str));
            h1list.set<Array<std::string> >("forms",h1PLin.get<Array<std::string> >(forms_str));
          }
        }
        else if (fPLin.isParameter(horizontal1_str)) {
          local_h1vals[0] = fPLin.get<double>(horizontal1_str);
        }
        else {
          std::string str = "Unrecognized horizontal1 permeability function parameters";
          std::cerr << fPLin << std::endl;
          BoxLib::Abort(str.c_str());
        }

        for (int k=0; k<local_h1vals.size(); k++) {
          local_h1vals[k] *= 1.01325e15 * scale; // convert from m^2 to mDa
        }
        h1list.set<Array<double> >("vals",local_h1vals);
        fPLout.set("horizontal1",h1list);
#endif
      }
    }

    void convert_TortuosityUniform(const ParameterList& fPLin,
                                   ParameterList&       fPLout)
    {
      Array<std::string> nullList, reqP;
      if (fPLin.isParameter("Value")) {
        const std::string val_name="Value"; reqP.push_back(val_name);
        PLoptions opt(fPLin,nullList,reqP,true,true);
        double val = fPLin.get<double>(val_name);
        fPLout.set<double>("val",val);
      }
      else {
        std::string str = "Unrecognized tortuosity parameters";
        std::cerr << fPLin << std::endl;
        BoxLib::Abort(str.c_str());
      }
    }

    void convert_SpecificStorageUniform(const ParameterList& fPLin,
                                        ParameterList&       fPLout)
    {
      Array<std::string> nullList, reqP;
      if (fPLin.isParameter("Value")) {
        const std::string val_name="Value"; reqP.push_back(val_name);
        PLoptions opt(fPLin,nullList,reqP,true,true);
        double val = fPLin.get<double>(val_name);
        fPLout.set<double>("val",val);
      }
      else {
        std::string str = "Unrecognized specific_storage parameters";
        std::cerr << fPLin << std::endl;
        BoxLib::Abort(str.c_str());
      }
    }

    void convert_SpecificYieldUniform(const ParameterList& fPLin,
                                      ParameterList&       fPLout)
    {
      Array<std::string> nullList, reqP;
      if (fPLin.isParameter("Value")) {
        const std::string val_name="Value"; reqP.push_back(val_name);
        PLoptions opt(fPLin,nullList,reqP,true,true);
        double val = fPLin.get<double>(val_name);
        fPLout.set<double>("val",val);
      }
      else {
        std::string str = "Unrecognized specific_yield parameters";
        std::cerr << fPLin << std::endl;
        BoxLib::Abort(str.c_str());
      }
    }

    void convert_ParticleDensityUniform(const ParameterList& fPLin,
                                        ParameterList&       fPLout)
    {
      Array<std::string> nullList, reqP;
      if (fPLin.isParameter("Value")) {
        const std::string val_name="Value"; reqP.push_back(val_name);
        PLoptions opt(fPLin,nullList,reqP,true,true);
        double val = fPLin.get<double>(val_name);
        fPLout.set<double>("val",val);
      }
      else {
        std::string str = "Unrecognized particle_density parameters";
        std::cerr << fPLin << std::endl;
        BoxLib::Abort(str.c_str());
      }
    }

    bool convert_DispersionTensorUniform(const ParameterList& fPLin,
                                         ParameterList&       fPLout)
    {
      const std::string alphaL_str = "alphaL";
      const std::string alphaT_str = "alphaT";
      Array<std::string> nullList, reqP;
      reqP.push_back(alphaL_str);
      reqP.push_back(alphaT_str);
      bool is_nonzero = true;
      PLoptions opt(fPLin,nullList,reqP,true,true);
      for (int i=0; i<reqP.size(); ++i) {
        double val = fPLin.get<double>(reqP[i]);
        is_nonzero &= (val != 0);
        fPLout.set<double>(reqP[i],val);
      }
      return is_nonzero;
    }

    static double gravity_magnitude(const ParameterList& parameter_list)
    {
      double gravity_mag = gravity_mag_DEF;
      if (parameter_list.isSublist("Execution Control")) {
        const ParameterList& ec_list = parameter_list.sublist("Execution Control");
        if (ec_list.isSublist("Numerical Control Parameters")) {
          const ParameterList& ncp_list = ec_list.sublist("Numerical Control Parameters");
          //if (ncp_list.isSublist("Basic Algorithm Control")) {
            //const ParameterList& bac_list = ncp_list.sublist("Basic Algorithm Control");
          if (ncp_list.isSublist("Structured Algorithm")) {
            const ParameterList& bac_list = ncp_list.sublist("Structured Algorithm");
            if (bac_list.isSublist("Expert Settings")) {
              const ParameterList& es_list = bac_list.sublist("Expert Settings");
              if (es_list.isParameter("gravity")) {
                gravity_mag = es_list.get<double>("gravity");
              }
            }
          }
        }
      }
      return gravity_mag;
    }

    static bool gravity_is_nonzero(const ParameterList& parameter_list)
    {
      return gravity_magnitude(parameter_list) != 0;
    }

    static int gravity_dir(const ParameterList& parameter_list)
    {
      int gravity_dir = gravity_dir_DEF;
      if (parameter_list.isSublist("Execution Control")) {
        const ParameterList& ec_list = parameter_list.sublist("Execution Control");
        if (ec_list.isSublist("Numerical Control Parameters")) {
          const ParameterList& ncp_list = ec_list.sublist("Numerical Control Parameters");
          //if (ncp_list.isSublist("Basic Algorithm Control")) {
            //const ParameterList& bac_list = ncp_list.sublist("Basic Algorithm Control");
          if (ncp_list.isSublist("Structured Algorithm")) {
            const ParameterList& bac_list = ncp_list.sublist("Structured Algorithm");
            if (bac_list.isSublist("Expert Settings")) {
              const ParameterList& es_list = bac_list.sublist("Expert Settings");
              if (es_list.isParameter("gravity_dir")) {
                gravity_dir = es_list.get<int>("gravity_dir");
              }
            }
          }
        }
      }
      return gravity_dir;
    }

#if BL_SPACEDIM == 2
    static double z_location(const ParameterList& parameter_list)
    {
      double z_loc = z_location_DEF;
      if (parameter_list.isSublist("Execution Control")) {
        const ParameterList& ec_list = parameter_list.sublist("Execution Control");
        if (ec_list.isSublist("Numerical Control Parameters")) {
          const ParameterList& ncp_list = ec_list.sublist("Numerical Control Parameters");
          //if (ncp_list.isSublist("Basic Algorithm Control")) {
            //const ParameterList& bac_list = ncp_list.sublist("Basic Algorithm Control");
          if (ncp_list.isSublist("Structured Algorithm")) {
            const ParameterList& bac_list = ncp_list.sublist("Structured Algorithm");
            if (bac_list.isSublist("Expert Settings")) {
              const ParameterList& es_list = bac_list.sublist("Expert Settings");
              if (es_list.isParameter("z_location")) {
                z_loc = es_list.get<double>("z_location");
              }
            }
          }
        }
      }
      return z_loc;
    }
#endif

    //
    // convert material to structured format
    //
    void
    convert_to_structured_material(const ParameterList& parameter_list,
                                   ParameterList&       struc_list,
                                   StateDef&            state,
                                   bool&                do_tracer_advection,
                                   bool&                do_tracer_diffusion)
    {
      ParameterList& rock_list = struc_list.sublist("rock");

      const ParameterList& rlist = parameter_list.sublist("Material Properties");

      const std::string mineralogy_str = "Mineralogy";
      const std::string complexation_str = "Surface Complexation Sites";
      const std::string isotherm_str = "Sorption Isotherms";
      const std::string cation_exchange_str = "Cation Exchange Capacity";
      state.getSolid().has_cation_exchange = false; // until we encounter the keyword for a material

      Array<std::string> arrayrock;

      bool add_chemistry_properties = false;

      std::map<std::string,SolidChem> solid_chem;
      std::map<std::string,double> cation_exchange_capacity;
      std::map<std::string,ParameterList> rsublist_mat;

      const std::string porosity_uniform_str = "Porosity: Uniform";
      const std::string porosity_gslib_str = "Porosity: GSLib";
      const std::string hydraulic_conductivity_uniform_str = "Hydraulic Conductivity: Uniform";
      const std::string perm_file_str = "Intrinsic Permeability: Input File";
      const std::string perm_uniform_str = "Intrinsic Permeability: Uniform";
      const std::string perm_anisotropic_uniform_str = "Intrinsic Permeability: Anisotropic Uniform";
      const std::string perm_gslib_str = "Intrinsic Permeability: GSLib";
      const std::string tortuosity_str = "Tortuosity Aqueous: Uniform";
      const std::string dispersivity_str = "Dispersion Tensor: Uniform Isotropic";
      const std::string specific_storage_uniform_str = "Specific Storage: Uniform";
      const std::string specific_yield_uniform_str = "Specific Yield: Uniform";
      const std::string particle_density_uniform_str = "Particle Density: Uniform";

      std::string kp_file_in, kp_file_out, pp_file_in, pp_file_out;
      std::string porosity_plotfile_in, porosity_plotfile_out;
      std::string permeability_plotfile_in, permeability_plotfile_out;
      bool kp_file_in_set, kp_file_out_set, pp_file_in_set, pp_file_out_set;
      bool kp_plotfile_in_set, kp_plotfile_out_set, pp_plotfile_in_set, pp_plotfile_out_set;
      kp_file_in_set = kp_file_out_set = pp_file_in_set = pp_file_out_set = false;
      kp_plotfile_in_set = kp_plotfile_out_set = pp_plotfile_in_set = pp_plotfile_out_set = false;

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
        if (entry.isList())
        {
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
              if (rlabel==porosity_gslib_str || rlabel==porosity_uniform_str) {
                if (mtest["Porosity"]) {
                  std::string str = "More than one of: (\""+porosity_gslib_str+"\", \""+porosity_uniform_str+
                    "\") specified for material \""+label+"\"";
                  BoxLib::Abort(str.c_str());
                }
                ParameterList psublist;
                if (rlabel==porosity_uniform_str) {
                  convert_PorosityUniform(rsslist,psublist);
                } else if (rlabel==porosity_gslib_str) {
                  convert_PorosityGSLib(rsslist,psublist);
                }
                rsublist.set("porosity",psublist);
                mtest["Porosity"] = true;
              }
              else if (rlabel==perm_uniform_str
		       || rlabel==perm_anisotropic_uniform_str
		       || rlabel==perm_gslib_str) {
                if (mtest["Intrinsic_Permeability"]) {
                  std::string str = "More than one of: (\""+perm_uniform_str
                    +"\", \""+perm_anisotropic_uniform_str
                    +"\", \""+perm_gslib_str
                    +"\", \""+hydraulic_conductivity_uniform_str
                    +"\") specified for material \""+label+"\"";
                  BoxLib::Abort(str.c_str());
                }
                ParameterList psublist;
		if (rlabel == perm_gslib_str) {
                  convert_PermeabilityGSLib(rsslist,psublist);
		} else {
		  convert_PermeabilityAnisotropic(rsslist,psublist,1.0);
		}
		rsublist.set("permeability",psublist);
                mtest["Intrinsic_Permeability"] = true;
              }
              else if (rlabel==hydraulic_conductivity_uniform_str) {
                if (mtest["Intrinsic_Permeability"]) {
                  std::string str = "More than one of: (\""+perm_uniform_str
                    +"\", \""+perm_anisotropic_uniform_str
                    +"\", \""+hydraulic_conductivity_uniform_str
                    +"\") specified for material \""+label+"\"";
                  BoxLib::Abort(str.c_str());
                }
                const std::string aq = "Aqueous";
                if (state.getPhases().count(aq) == 0) {
                  std::string str = "Hydraulic conductivity specified for material \""+label
                    +"\" but phase \""+aq+"\" does not exist";
                  BoxLib::Abort(str.c_str());
                }
                double density = state.getPhases()[aq].Density();
                if (density == 0) {
                  std::string str = "Hydraulic conductivity specified for material \""+label
                    +"\" but density of \""+aq+"\"  = 0";
                  BoxLib::Abort(str.c_str());
                }
                double viscosity = state.getPhases()[aq].Viscosity();
                if (viscosity == 0) {
                  std::string str = "Hydraulic conductivity specified for material \""+label
                    +"\" but viscosity of \""+aq+"\"  = 0";
                  BoxLib::Abort(str.c_str());
                }
                double gravity_mag = gravity_magnitude(parameter_list);
                if (gravity_mag == 0) {
                  std::string str = "Hydraulic conductivity specified for material \""+label
                    +"\" but gravity magnitude = 0";
                  BoxLib::Abort(str.c_str());
                }
                double factor = viscosity / (density * gravity_mag);
                ParameterList psublist;
                convert_PermeabilityAnisotropic(rsslist,psublist,factor);
                rsublist.set("permeability",psublist);
                rsublist.set("permeability_dist","uniform");
                mtest["Intrinsic_Permeability"] = true;
              }
              else if (rlabel==tortuosity_str) {
                ParameterList dsublist;
                convert_TortuosityUniform(rsslist,dsublist);
                rsublist.set("tortuosity",dsublist);
              }
              else if (rlabel==dispersivity_str) {
                ParameterList dsublist;
                bool is_nonzero = convert_DispersionTensorUniform(rsslist,dsublist);
                if (is_nonzero) {
                  rsublist.set("dispersivity",dsublist);
                  do_tracer_diffusion = do_tracer_advection;
                }
              }
              else if (rlabel==specific_storage_uniform_str) {
                ParameterList ssublist;
                convert_SpecificStorageUniform(rsslist,ssublist);
                rsublist.set("specific_storage",ssublist);
              }
              else if (rlabel==specific_yield_uniform_str) {
                ParameterList ssublist;
                convert_SpecificYieldUniform(rsslist,ssublist);
                rsublist.set("specific_yield",ssublist);
              }
              else if (rlabel==particle_density_uniform_str) {
                ParameterList ssublist;
                convert_ParticleDensityUniform(rsslist,ssublist);
                rsublist.set("particle_density",ssublist);
              }
              else if (rlabel=="Capillary Pressure: van Genuchten" || rlabel=="Capillary Pressure: Brooks Corey") {
                double alpha = rsslist.get<double>("alpha");
                ParameterList cpl_pl;
                double ell;
                if (rlabel=="Capillary Pressure: van Genuchten") {
                  cpl_pl.set("type","VanGenuchten");
                  cpl_pl.set("m",rsslist.get<double>("m"));
                  ell = 0.5;
                } else {
                  cpl_pl.set("type","BrooksCorey");
                  cpl_pl.set("lambda",rsslist.get<double>("lambda"));
                  ell = 2;
                }
                cpl_pl.set("Sr",rsslist.get<double>("Sr"));
                cpl_pl.set("alpha",alpha*1.01325e5); // convert Pa^-1 to atm^-1
                if (rsslist.isParameter("ell")) {
                  ell = rsslist.get<double>("ell");
                }
                cpl_pl.set("Kr_ell",ell);
                rsublist.set("cpl",cpl_pl);
                double Kr_smoothing_max_pcap = -1;
                if (rsslist.isParameter("krel smoothing interval")) {
                  Kr_smoothing_max_pcap = rsslist.get<double>("krel smoothing interval");
                }
                rsublist.set("Kr_smoothing_max_pcap",Kr_smoothing_max_pcap);

                const std::string WRM_plot_file_str = "WRM Plot File";
                if (rsslist.isParameter(WRM_plot_file_str)) {
                  std::string WRM_plot_file = rsslist.get<std::string>(WRM_plot_file_str);
                  rsublist.set("WRM_plot_file",WRM_plot_file);

                  const std::string WRM_plot_file_num_pts_str = "WRM Plot File Number Of Points";
                  int WRM_plot_file_numPts = 500; // Default
                  if (rsslist.isParameter(WRM_plot_file_num_pts_str)) {
                    WRM_plot_file_numPts = rsslist.get<int>(WRM_plot_file_num_pts_str);
                  }
                  rsublist.set("WRM_plot_file_num_pts",WRM_plot_file_numPts);
                }
                mtest["Capillary_Pressure"] = true;

                std::string krType = rsslist.get<std::string>("Relative Permeability");
                if (krType=="Mualem") {
                  rsublist.set("Kr_model",krType);
                  double Kr_ell = 0.5; // Default
                  if (rsslist.isParameter("ell")) {
                    Kr_ell = rsslist.get<double>("ell");
                  }
                  rsublist.set("Kr_ell",Kr_ell);
                  mtest["Relative_Permeability"] = true;
                }
                else if (krType=="Burdine") {
                  rsublist.set("Kr_model",krType);
                  double Kr_ell = 2.0; // Default
                  if (rsslist.isParameter("ell")) {
                    Kr_ell = rsslist.get<double>("ell");
                  }
                  rsublist.set("Kr_ell",Kr_ell);
                  mtest["Relative_Permeability"] = true;
                }
                else {
                  std::cerr << "Unsupported Relative Permeability model: " << krType << std::endl;
                  throw std::exception();
                }
              }
              else if (rlabel=="Capillary Pressure: None") {
                ParameterList cpl_pl;
		cpl_pl.set("type","None");
                rsublist.set("cpl",cpl_pl);
                mtest["Capillary_Pressure"] = true;
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

		// FIXME: Short-circuit all this stuff, since the powers that
		//        be decided to special-case all this for Aqueous Water
                //PLoptions sipP(rsslist,nullList,nullList,false,true); // each optional list is a phase
                //const Array<std::string>& sipLabels = sipP.OptLists();
		Array<std::string> sipLabels; sipLabels.push_back("Aqueous");
                for (int k=0; k<sipLabels.size(); ++k) {
                  const std::string& sipLabel = sipLabels[k];
                  std::string _sipLabel = underscore(sipLabel);

                  StateDef::PhaseCompMap& pc_map = state.getPhaseCompMap();
                  if (pc_map.find(sipLabel)==pc_map.end()) {
                    std::cerr << "Unknown phase " << sipLabel << " in " << rlabel << " for " << label << std::endl;
                    throw std::exception();
                  }
                  //const ParameterList& sipSL = rsslist.sublist(sipLabel);

		  // FIXME: Short-circuit all this stuff, since the powers that
		  //        be decided to special-case all this for Aqueous Water
                  //PLoptions sipcP(sipSL,nullList,nullList,false,true); // each optional list is a component
                  //const Array<std::string>& sipcLabels = sipcP.OptLists();
		  Array<std::string> sipcLabels; sipcLabels.push_back("Water");
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

		    // FIXME: Short-circuit all this stuff, since the powers that
		    //        be decided to special-case all this for Aqueous Water
                    //PLoptions sipcsP(sipcSL,nullList,nullList,false,true); // each optional list is a solute
                    //const ParameterList& sipcSL = sipSL.sublist(sipcLabel);

		    PLoptions sipcsP(rsslist,nullList,nullList,false,true); // each optional list is a solute

                    const Array<std::string>& sipcsLabels = sipcsP.OptLists();
                    for (int M=0; M<sipcsLabels.size(); ++M) {
                      const std::string& sipcsLabel = sipcsLabels[M];
                      std::string _sipcsLabel = underscore(sipcsLabel);
                      //const ParameterList& sipcsSL = sipcSL.sublist(sipcsLabel);
                      const ParameterList& sipcsSL = rsslist.sublist(sipcsLabel);

                      if ( !(c_map[sipcLabel].HasTracer(sipcsLabel)) ) {
                        std::cerr << "Unknown solute " << sipcsLabel << " in component "
                                  << sipcLabel << " in phase " << sipLabel << " for " <<
                          rlabel << " in " << label << std::endl;
                        throw std::exception();
                      }


                      //SolidChem::SorptionIsothermData& iso = solid_chem[label].sorption_isotherms(sipcsLabel);

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
                          solid_chem[label].SorptionIsotherm(sipcsLabel).Kd = sipcsSL.get<double>(siLabels[N]);
                        }
                        else if (siLabels[N] == Lb_str) {
                          solid_chem[label].SorptionIsotherm(sipcsLabel).Langmuir_b = sipcsSL.get<double>(siLabels[N]);
                          solid_chem[label].SorptionIsotherm(sipcsLabel).Freundlich_n = -1.0;
                          n_found = true;
                        }
                        else if (siLabels[N] == Fn_str) {
                          solid_chem[label].SorptionIsotherm(sipcsLabel).Freundlich_n = sipcsSL.get<double>(siLabels[N]);
                          solid_chem[label].SorptionIsotherm(sipcsLabel).Langmuir_b = -1.0;
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
                      //else
                      //solid_chem[label].sorption_isotherms(sipcsLabel) = iso;
                      //bool successfully_inserted = SolidChem::InsertSorptionIsotherm(sipcsLabel,iso);
                      //if (!successfully_inserted) {
                      //    bool compatible = iso.IsFreundlich() ^ SolidChem::SorptionIsotherm(sipcsLabel).IsFreundlich();
                      //    if (!compatible) {
                      //        std::cerr << "Only one \"" << Lb_str << "\" and \"" << Fn_str
                      //                  << "\" may be specified for each solute.  Both given for \""
                      //                  << sipcsLabel << "\" in different materials." << std::endl;
                      //        throw std::exception();
                      //    }
                      //}
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
          std::cerr << "Unsupported material property: " << label << std::endl;
          throw std::exception();
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
              const Array<TRACER>& solutes = cit->second.getTracerArray();
              for (int i=0; i<solutes.size(); ++i) {
                const std::string& s=solutes[i].name;
                if (solid_chem[label].HasSorptionIsotherm(s)) {
                  SolidChem::SorptionIsothermData sid = solid_chem[label].SorptionIsotherm(s);
                  ParameterList sitPL = sid.BuildPL();
                  siPL.set(underscore(s),sitPL);
                  solutes_with_isotherms.insert(s);
                }
              }
            }
          }
          if (solutes_with_isotherms.size()>0) {
            rsublist.set("sorption_isotherms",siPL);
            state.getSolid().sorption_isotherm_names.resize(solutes_with_isotherms.size());
            for (std::set<std::string>::const_iterator it=solutes_with_isotherms.begin(),
                   End=solutes_with_isotherms.end(); it!=End; ++it) {
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

      rock_list.set("rocks",arrayrock);
      if (kp_file_out_set) {
        rock_list.set("permeability_output_file",kp_file_out);
      }
      if (kp_file_in_set) {
        rock_list.set("permeability_input_file",kp_file_in);
      }
      if (pp_file_out_set) {
        rock_list.set("porosity_output_file",pp_file_out);
      }
      if (pp_file_in_set) {
        rock_list.set("porosity_input_file",pp_file_in);
      }
      if (pp_plotfile_in_set) {
        rock_list.set("porosity_plotfile_in",porosity_plotfile_in);
      }
      if (pp_plotfile_out_set) {
        rock_list.set("porosity_plotfile_out",porosity_plotfile_out);
      }
      if (kp_plotfile_in_set) {
        rock_list.set("permeability_plotfile_in",permeability_plotfile_in);
      }
      if (kp_plotfile_out_set) {
        rock_list.set("permeability_plotfile_out",permeability_plotfile_out);
      }
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
        const ParameterList& pplist = plist.sublist(phaseLabel);
        const ParameterList& psublist = plist.sublist(phaseLabel);

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
            PLoptions optCC(slist,reqLc,reqPc,false,false);
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

            const Array<std::string>& sLists = optCC.OptLists();
            for (int k=0; k<sLists.size(); ++k) {
              int iSolute = -1;
              for (int L=0; L<sLabels.size(); ++L) {
                if (sLabels[L] == sLists[k]) {
                  iSolute = L;
                }
              }
              if (iSolute < 0) {
                std::cerr << "Solute ParameterList has name not in \"Component Solutes\" list: " << sLists[k] << std::endl;
                throw std::exception();
              }

              PLoptions optS(slist.sublist(sLists[k]),nullList,nullList,true,false);
              const Array<std::string>& slParams = optS.OptParms();
              const std::string molec_diff_str = "Molecular Diffusivity";
              const std::string fodecay_constant_str = "First Order Decay Constant";

              TRACER& s = c.getTracerArray()[iSolute];
              for (int M=0; M<slParams.size(); ++M) {
                if (slParams[M] == molec_diff_str) {
                  s.molecularDiffusivity = pclist.sublist(compLabel).sublist(sLabels[iSolute]).get<double>(molec_diff_str);
                }
                else if (slParams[M] == fodecay_constant_str) {
                  s.firstOrderDecayConstant = pclist.sublist(compLabel).sublist(sLabels[iSolute]).get<double>(fodecay_constant_str);
                }
                else {
                  std::cerr << "Solute ParameterList contains unrecognized parameter: " << slParams[M] << std::endl;
                  throw std::exception();
                }
              }
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

        Array<std::string> reqPbc, reqLbc;
        reqPbc.push_back("Assigned Regions");
        PLoptions BCs(bc_plist,reqLbc,reqPbc,false,false);

        // The optional list is for the phase/comp func
        Array<std::string> phaseBCfuncs;
        const Array<std::string>& phaseBCentries = BCs.OptLists();
        bool has_solute_bcs = false;
        for (int j=0; j<phaseBCentries.size(); ++j) {
          if (phaseBCentries[j] == solute_section_label) {
            has_solute_bcs = true;
          } else {
            phaseBCfuncs.append(phaseBCentries[j]);
          }
        }

        if (phaseBCfuncs.size()!=1) {
          std::cerr << "StateDef::build_solute_funcs: Only one Phase/comp BC and IC allowed for BC/IC label: " << BClabel << std::endl;
          throw std::exception();
        }
        const ParameterList& func_plist = bc_plist.sublist(phaseBCfuncs[0]);
        const Array<std::string>& _assigned_regions = bc_plist.get<Array<std::string> >(reqPbc[0]);

        Array<std::string> assigned_regions(_assigned_regions.size());
        for (int i=0; i<_assigned_regions.size(); ++i) {
          assigned_regions[i] = underscore(_assigned_regions[i]);
        }


        const std::string& Amanzi_type = phaseBCfuncs[0];

        s[BClabel] = StateFunc(BClabel, Amanzi_type, func_plist, assigned_regions);

        // Now add solute data
        if (has_solute_bcs) {
          Array<std::string> nullList;
          const ParameterList& solute_plist = bc_plist.sublist(solute_section_label);
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
                PLoptions soluteOPTf(solutePL,nullList,funcP,false,true);

                // Get function name/list
                const Array<std::string>& funcNames = soluteOPTf.OptLists();
                if (funcNames.size()!=1) {
                  std::cerr << "Each solute BC expects a single function" << std::endl;
                  throw std::exception();
                }
                const std::string& Amanzi_solute_type = funcNames[0];
                ParameterList solute_func_plist = solutePL.sublist(Amanzi_solute_type);
                s[BClabel][phaseNames[icp]][compNames[icc]][soluteNames[ics]]
                  = ICBCFunc(solute_func_plist,Amanzi_solute_type,BClabel);
              }
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
              const Array<TRACER>& ds=(*this)[phaseName][compName].getTracerArray();
              bool found = false;
              for (int it=0; it<ds.size() && !found; ++it) {
                found = ds[it].name == soluteName;
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
      {
        std::cerr << "IC: Saturation functions no longer supported " << std::endl;
        throw std::exception();
      }
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
      //const std::string phase_name="Phase";reqP.push_back(phase_name);
      const std::string val_name="Value";reqP.push_back(val_name);
      PLoptions opt(fPLin,nullList,reqP,true,true);
      fPLout.set<std::string>("type","pressure");
      fPLout.set<double>("val",fPLin.get<double>(val_name));
      //fPLout.set<std::string>("phase",fPLin.get<std::string>(phase_name));
    }

    void convert_IC_LinPressure(const ParameterList& fPLin,
                                const std::string&   Amanzi_type,
                                ParameterList&       fPLout)
    {
      Array<std::string> reqP, nullList;
      //const std::string phase_name="Phase";reqP.push_back(phase_name);
      const std::string rval_name="Reference Value";reqP.push_back(rval_name);
      const std::string grad_name="Gradient Value";reqP.push_back(grad_name);
      const std::string ref_name="Reference Point";reqP.push_back(ref_name);
      PLoptions opt(fPLin,nullList,reqP,true,true);

      fPLout.set<std::string>("type","linear_pressure");
      fPLout.set<double>("val",fPLin.get<double>(rval_name));
      fPLout.set<Array<double> >("grad",fPLin.get<Array<double> >(grad_name));
      fPLout.set<Array<double> >("ref_coord",fPLin.get<Array<double> >(ref_name));
    }

    void convert_ICFlow(const ParameterList& fPLin,
                        const std::string&   Amanzi_type,
                        ParameterList&       fPLout)
    {
      //const std::string phase_name="Phase";
      const std::string val_name="Reference Value";
      const std::string grad_name="Gradient Value";
      const std::string ref_name="Reference Point";
      const std::string vel_name="Aqueous Volumetric Flux";

      Array<std::string> reqP, nullList;
      reqP.push_back(val_name);
      if (Amanzi_type == "IC: Flow") {
        //reqP.push_back(phase_name);
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
                               ParameterList&       fPLout,
			       const ParameterList& entire_PL,
			       StateDef&            stateDef)
    {
      int gdir = gravity_dir(entire_PL);

      Array<std::string> reqP, nullList;
      const std::string wth_name="Water Table Height";reqP.push_back(wth_name);
      PLoptions opt(fPLin,nullList,reqP,true,true);
      double wth = fPLin.get<double>(wth_name);

      Array<double> pgrad(BL_SPACEDIM,0);

      if (gdir < 0 || gdir > 2) {
	MyAbort("gravity_dir cannot be < 0 or > 2");
      }

      double gval = gravity_magnitude(entire_PL);
      double density = stateDef.getPhases()["Aqueous"].Density();

      Array<double> ref_coord(BL_SPACEDIM,0);
      if (gdir < BL_SPACEDIM) {
	pgrad[gdir] = density * gval;
	ref_coord[gdir] = wth;
      }

      double val = atmToMKS; // Value at water table = 1atm
      if (gdir >= BL_SPACEDIM) {
	double z_loc = z_location(entire_PL);
	val += (wth - z_loc) * density * gval;
      }

      fPLout.set<std::string>("type","linear_pressure");
      fPLout.set<double>("val",val);
      fPLout.set<Array<double> >("grad",pgrad);
      fPLout.set<Array<double> >("ref_coord",ref_coord);
    }

    void convert_ICVel(const ParameterList& fPLin,
                       const std::string&   Amanzi_type,
                       ParameterList&       fPLout)
    {
      Array<std::string> reqP, nullList;
      const std::string vel_name="Velocity Vector"; reqP.push_back(vel_name);
      PLoptions opt(fPLin,nullList,reqP,true,true);

      fPLout.set<Array<double> >(underscore(vel_name),fPLin.get<Array<double> >(vel_name));
      fPLout.set<std::string>("type","constant_velocity");
    }

    void convert_solute_ICConcentration(const ICBCFunc& solute_ic,
                                        ParameterList&  fPLout,
                                        int             do_chem)
    {
      const ParameterList& fPLin = solute_ic.PList();
      const std::string& solute_ic_Amanzi_type = solute_ic.Amanzi_Type();
      const std::string& solute_ic_label = solute_ic.Label();

      const std::string geo_name="Geochemical Condition";
      const std::string val_name="Value";
      const std::string time_name="Times";
      const std::string form_name="Time Functions";
      const std::string ion_name="Free Ion Guess";

      if (fPLin.isParameter(geo_name)) {
        fPLout.set<std::string>("geochemical_condition",fPLin.get<std::string>(geo_name));
      } else if (fPLin.isParameter(val_name)) {
        fPLout.set<double>("val",fPLin.get<double>(val_name));
      } else {
        const std::string str = "Solute Concentration for IC \""+solute_ic_label+
          "\" must be specified using \""+geo_name+"\" or \""+val_name+"\"";
        MyAbort(str);
      }

      fPLout.set<std::string>("type","concentration");

      if (do_chem && fPLin.isParameter(ion_name)) {
        fPLout.set<double>(underscore(ion_name),fPLin.get<double>(ion_name));
      }
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

        if (require_static_velocity < 0) {
          std::cerr << "Initialization out of order internally, must determine "
                    << "if fixed velocity required prior to reading phase/component BCs" << std::endl;
          throw std::exception();
        }

        if (require_static_velocity == 1) {
          if ( Amanzi_type != "IC: Uniform Velocity" ) {
            std::cerr << "For this time integration mode, a velocity field must be specified "
                      << "rather than any other IC type" << std::endl;
            throw std::exception();
          }

          convert_ICVel(fPLin,Amanzi_type,fPLout);
        }
        else {
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
            convert_ICHydrostatic(fPLin,Amanzi_type,fPLout,parameter_list,stateDef);
          }
          else if ( Amanzi_type == "IC: Flow" )
          {
            convert_ICFlow(fPLin,Amanzi_type,fPLout);
          }
          else {
            std::cerr << "Unsupported IC: " << Amanzi_type << std::endl;
            throw std::exception();
          }
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

    void convert_length_units(const ParameterList& fPLin,
                              Array<double>&       vals)
    {
      const std::string units_name = "Units";
      std::map<std::string,double> conversion;
      conversion["m"] = 1;
      conversion["ft"] = 12*2.54*.01;

      Array<std::string> nullList;
      PLoptions opt(fPLin,nullList,nullList,false,false);
      const Array<std::string>& optional_params = opt.OptParms();
      if (optional_params.size()>0) {
        for (int i=0; i<optional_params.size(); ++i) {
          if (optional_params[i] == units_name) {
            const std::string& units = fPLin.get<std::string>(units_name);
            std::map<std::string,double>::const_iterator it=conversion.find(units);
            if (it != conversion.end()) {
              double factor = it->second;
              for (int j=0; j<vals.size(); ++j) {
                vals[j] *= factor;
              }
            }
            else {
              std::cerr << " Unrecognized units (\"" << units << "\") in " << fPLin << std::endl;
              throw std::exception();
            }
          }
        }
      }
    }


    void convert_BCPressure(const ParameterList& fPLin,
                            const std::string&   Amanzi_type,
                            ParameterList&       fPLout)
    {
      Array<std::string> nullList, reqP;

      if (Amanzi_type == "BC: Linear Pressure")
      {
        const std::string val_name="Reference Value"; reqP.push_back(val_name);
        const std::string grad_name="Gradient Value";reqP.push_back(grad_name);
        const std::string ref_name="Reference Point"; reqP.push_back(ref_name);
        PLoptions opt(fPLin,nullList,reqP,true,true);

        fPLout.set<double>("val",fPLin.get<double>(val_name));
        fPLout.set<Array<double> >("grad",fPLin.get<Array<double> >(grad_name));
        fPLout.set<Array<double> >("loc",fPLin.get<Array<double> >(ref_name));
        fPLout.set<std::string>("type","linear_pressure");
      }
      else if ((Amanzi_type == "BC: Uniform Hydraulic Head")
               || (Amanzi_type == "BC: Hydrostatic"))
      {
        const std::string val_name
          = (Amanzi_type == "BC: Uniform Hydraulic Head" ?  "Values" : "Water Table Height");
        reqP.push_back(val_name);
        const std::string time_name="Times"; reqP.push_back(time_name);
        const std::string form_name="Time Functions"; reqP.push_back(form_name);
        PLoptions opt(fPLin,nullList,reqP,true,false);
        Array<double> vals = fPLin.get<Array<double> >(val_name);
        convert_length_units(fPLin,vals);
        const std::string Coordinate_System_str = "Coordinate_System";
        const std::string Coordinate_System_Absolute_str = "Absolute";
        const std::string Coordinate_System_Relative_str = "Relative";
        std::string absolute_or_relative = Coordinate_System_Absolute_str; // Default value
        if (Amanzi_type == "Water Table Height") {
          const Array<std::string>& optParams = opt.OptParms();
          for (int i=0; i<optParams.size(); ++i) {
            const std::string& optParam = optParams[i];
            if (optParam == Coordinate_System_str) {
              absolute_or_relative = fPLin.get<std::string>(Coordinate_System_str);
              if ( ! (absolute_or_relative==Coordinate_System_Absolute_str
                      || absolute_or_relative==Coordinate_System_Relative_str)) {
                std::string j = "Value for \"BC: Hydrostatic\" parameter: \""+Coordinate_System_str+"\" (if provided) must be"
                  + " either \""+Coordinate_System_Absolute_str+"\" or \""+Coordinate_System_Absolute_str+"\" (\""
                  + absolute_or_relative +"\" given)";
              }
            }
          }
        }

        fPLout.set<Array<double> >("vals",vals);
        if (vals.size()>1) {
          fPLout.set<Array<double> >("times",fPLin.get<Array<double> >(time_name));
          fPLout.set<Array<std::string> >("forms",fPLin.get<Array<std::string> >(form_name));
        }
        fPLout.set<std::string>("type","hydraulic_head");
        fPLout.set<std::string>("normalization",absolute_or_relative);
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
        fPLout.set<std::string>("type","pressure");
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

      const std::string geo_name="Geochemical Condition";
      const std::string geos_name="Geochemical Conditions";
      const std::string val_name="Value";
      const std::string vals_name="Values";
      const std::string time_name="Times";
      const std::string form_name="Time Functions";

      if (fPLin.isParameter(geo_name)) {
        fPLout.set<std::string>("geochemical_conditions",fPLin.get<std::string>(geo_name));
      } else if (fPLin.isParameter(geos_name)) {
        const Array<std::string>& geo_names = fPLin.get<Array<std::string> >(geos_name);
        if (geo_names.size() > 1) {
          fPLout.set<Array<double> >("times",fPLin.get<Array<double> >(time_name));
          fPLout.set<Array<std::string> >("forms",fPLin.get<Array<std::string> >(form_name));
        }
        fPLout.set<Array<std::string> >("geochemical_conditions",geo_names);
      } else if (fPLin.isParameter(val_name)) {
        fPLout.set<double>("vals",fPLin.get<double>(val_name));
      } else if (fPLin.isParameter(vals_name)) {
        Array<double> vals = fPLin.get<Array<double> >(vals_name);
        fPLout.set<Array<double> >("vals",vals);
        if (vals.size() > 1) {
          fPLout.set<Array<double> >("times",fPLin.get<Array<double> >(time_name));
          fPLout.set<Array<std::string> >("forms",fPLin.get<Array<std::string> >(form_name));
        }
      } else {
        const std::string str = "Solute Concentration for BC \""+solute_bc_label
          +"\" must be specified using \""+geo_name+"\", \""+val_name+"\" or \""+vals_name;
        MyAbort(str);
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
        else if ( (Amanzi_type == "BC: Uniform Pressure")
                  || (Amanzi_type == "BC: Uniform Hydraulic Head")
                  || (Amanzi_type == "BC: Hydrostatic")
                  || (Amanzi_type == "BC: Linear Pressure") )
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
            std::string orient_type = orient_RTs[0].second;
            for (int j=1; j<orient_RTs.size(); ++j) {
              if (orient_type != orient_RTs[j].second) {
                std::cerr << "Structured grid requires that all BCs on "
                          << Amanzi_purpose << " be of the same type" << std::endl;
                throw std::exception();
              }
            }
            if (orient_type == "pressure"
                     || orient_type == "hydraulic_head"
                     || orient_type == "linear_pressure"
                     || orient_type == "hydrostatic") {
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
            else if (orient_type == "zero_total_velocity" || orient_type == "noflow") {
              sat_bc      = 4; // Neumann for saturation
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
            sat_bc = 4;
            pressure_bc = 1;
            inflow_bc = 1;
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
    convert_solute_bcs(ParameterList& struc_list,
                       StateDef&      stateDef)
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
              //const std::string& solute_bc_units = solute_bc.Units();

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

      ParameterList& geom_list = struc_list.sublist("geometry");
      std::map<std::string,std::string> orient_type_map;
      for (SolutePLMMap::const_iterator it = solute_to_BClabel.begin(),
             End=solute_to_BClabel.end(); it!=End; ++it) {

        const std::string& trac_name = it->first;
        const ParameterList& trac_bc_pl = it->second;
        const std::string& trac_bc_type = trac_bc_pl.get<std::string>("type");

        const Array<std::string>& trac_bc_regions = trac_bc_pl.get<Array<std::string> >("regions");
        for (int j=0; j<trac_bc_regions.size(); ++j) {
          const std::string& regionName = trac_bc_regions[j];
          if (geom_list.isSublist(regionName)) {
            const std::string& purpose = geom_list.sublist(regionName).get<std::string>("purpose");
            std::map<std::string,std::string>::const_iterator it = orient_type_map.find(purpose);
            if (it!=orient_type_map.end()) {
              if (trac_bc_type != it->second) {
                std::cerr << "All solutes must have the same boundary condition type a side\n";
                throw std::exception();
              }
            }
            else {
              orient_type_map[purpose] = trac_bc_type;
            }
          }
          else {
            std::cerr << "Unknown region " << regionName << " in boundary condition for " << trac_name << " \n";
            throw std::exception();
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

    void convert_Sources(const ParameterList& fPLin,
			 const std::string&   Amanzi_type,
			 ParameterList&       fPLout)
    {
      const std::string Source_Uniform_str = "Source: Uniform";
      const std::string Source_Volume_Weighted_str = "Source: Volume Weighted";
      const std::string Source_Permeability_Weighted_str = "Source: Permeability Weighted";
      const std::string Source_Point_str = "Source: Point";

      Array<std::string> nullList, reqP;

      if (Amanzi_type == Source_Uniform_str
	  || Amanzi_type == Source_Volume_Weighted_str
	  || Amanzi_type == Source_Permeability_Weighted_str
	  || Amanzi_type == Source_Point_str) {
        const std::string val_name="Values"; reqP.push_back(val_name);
        const std::string time_name="Times"; reqP.push_back(time_name);
        const std::string form_name="Time Functions"; reqP.push_back(form_name);
        PLoptions opt(fPLin,nullList,reqP,true,false);
        Array<double> vals = fPLin.get<Array<double> >(val_name);
        fPLout.set<Array<double> >("vals",vals);
        if (vals.size()>1) {
          fPLout.set<Array<double> >("times",fPLin.get<Array<double> >(time_name));
          fPLout.set<Array<std::string> >("forms",fPLin.get<Array<std::string> >(form_name));
        }
	if (Amanzi_type == Source_Uniform_str) {
	  fPLout.set<std::string>("type","uniform");
	} else if (Amanzi_type == Source_Volume_Weighted_str) {
	  fPLout.set<std::string>("type","volume_weighted");
	} else if (Amanzi_type == Source_Point_str) {
	  fPLout.set<std::string>("type","point");
	} else {
	  fPLout.set<std::string>("type","permeability_weighted");
	}
      }
      else {
        MyAbort("Unrecognized component source function: "+Amanzi_type);
      }
    }

    void convert_Solute_Sources(const ParameterList& fPLin,
				const std::string&   Amanzi_type,
				ParameterList&       fPLout)
    {
      const std::string Source_Uniform_str = "Source: Uniform Concentration";
      const std::string Source_Flow_Weighted_str = "Source: Flow Weighted Concentration";
      const std::string Source_DiffDomRel_str = "Source: Diffusion Dominated Release Model";
      const std::string Total_Inventory_str = "Total Inventory";
      const std::string Mixing_Length_str = "Mixing Length";
      const std::string Eff_Diff_Coeff_str = "Effective Diffusion Coefficient";
      const std::string Source_Point_str = "Source: Point";
      const std::string Values_str="Values";
      const std::string Times_str="Times";
      const std::string Forms_str="Time Functions";

      Array<std::string> nullList, reqP;
      if (Amanzi_type == Source_Uniform_str
	  || Amanzi_type == Source_Flow_Weighted_str
	  || Amanzi_type == Source_Point_str) {
        reqP.push_back(Values_str);
        reqP.push_back(Times_str);
        reqP.push_back(Forms_str);
        PLoptions opt(fPLin,nullList,reqP,true,false);
        Array<double> vals = fPLin.get<Array<double> >(Values_str);
        fPLout.set<Array<double> >("vals",vals);
        if (vals.size()>1) {
          fPLout.set<Array<double> >("times",fPLin.get<Array<double> >(Times_str));
          fPLout.set<Array<std::string> >("forms",fPLin.get<Array<std::string> >(Forms_str));
        }
	if (Amanzi_type == Source_Uniform_str) {
	  fPLout.set<std::string>("type","uniform");
	} else if (Amanzi_type == Source_Flow_Weighted_str) {
	  fPLout.set<std::string>("type","flow_weighted");
        } else if (Amanzi_type == Source_Point_str) {
	  fPLout.set<std::string>("type","point");
	}
      }
      else if (Amanzi_type == Source_DiffDomRel_str) {
        reqP.push_back(Total_Inventory_str);
        reqP.push_back(Mixing_Length_str);
        reqP.push_back(Eff_Diff_Coeff_str);
        reqP.push_back(Times_str);
        PLoptions opt(fPLin,nullList,reqP,true,true);
        fPLout.set<std::string>("type","diffusion_dominated_release_model");
        fPLout.set<double>("total_inventory",fPLin.get<double>(Total_Inventory_str));
        fPLout.set<double>("mixing_length",fPLin.get<double>(Mixing_Length_str));
        fPLout.set<double>("effective_diffusion_coef",fPLin.get<double>(Eff_Diff_Coeff_str));
        fPLout.set<double>("time_scale",100);
        const Array<double>& times = fPLin.get<Array<double> >(Times_str);
        fPLout.set<double>("start_time",times[0]);
        fPLout.set<double>("end_time",times[1]);
      }
      else {
        MyAbort("Unrecognized solute source function: "+Amanzi_type);
      }
    }

    void
    convert_to_structured_state(const ParameterList& parameter_list,
                                ParameterList&       struc_list,
                                StateDef&            stateDef,
                                bool&                do_tracer_advection,
                                bool&                do_tracer_diffusion,
                                bool&                do_chem)
    {
      ParameterList& phase_list  = struc_list.sublist("phase");
      ParameterList& comp_list   = struc_list.sublist("comp");
      ParameterList& solute_list = struc_list.sublist("tracer");

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

          const Array<TRACER>& solutes = cit->second.getTracerArray();
          for (int i=0; i<solutes.size(); ++i)
          {
            std::string _soluteLabel = underscore(solutes[i].name);
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

        // Find solute in phase defs in order to extract solute-specific properties
        bool found_solute = false;
        for (StateDef::Phases::const_iterator pit = phases.begin(); pit!=phases.end() && !found_solute; ++pit)
        {
          const std::string& phaseLabel = pit->first;
          PHASE& phase = stateDef.getPhases()[phaseLabel];
          const CompMap comp_map = stateDef[phaseLabel];
          for (CompMap::const_iterator cit = comp_map.begin(); cit!=comp_map.end() && !found_solute; ++cit)
          {
            const Array<TRACER>& solutes = cit->second.getTracerArray();
            for (int i=0; i<solutes.size(); ++i)
            {
              std::string _soluteLabel = underscore(solutes[i].name);
              if (_soluteLabel == soluteName) {
                found_solute = true;

                double D = solutes[i].molecularDiffusivity;
                if (D != 0) {
                  solutePLs[soluteName].set<double>("molecularDiffusivity",D);
                  do_tracer_diffusion = do_tracer_advection;
                }

                double lambda = solutes[i].firstOrderDecayConstant;
                if (lambda != 0) {
                  solutePLs[soluteName].set<double>("firstOrderDecayConstant",lambda);
                }
              }
            }
          }
        }
      }

      // Only do solute BCs if do_tracer_transport
      if (do_tracer_advection)
      {
        SolutePLMMap solute_to_bctype = convert_solute_bcs(struc_list,stateDef);

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

      /*
	Note: Looking to translate xml block that looks something like:


	<ParameterList name="Sources">
  	  <ParameterList name="SOURCE 1">
  	    <Parameter name="Assigned Regions" type="Array(string)" value="{All}"/>
            <ParameterList name="Water">
  	      <ParameterList name="Source: Uniform">
  	        <Parameter name="Values" type="Array(double)" value="{19}"/>
  	      </ParameterList>
  	      <ParameterList name="Solute SOURCE">
  	        <ParameterList name="Tc-99">
  	          <ParameterList name="Source: Flow Weighted Concentration">
  	    	      <Parameter name="Values" type="Array(double)" value="{19}"/>
  	          </ParameterList>
  	        </ParameterList>
  	        <ParameterList name="U 237">
  	          <ParameterList name="Source: Uniform Concentration">
  	    	      <Parameter name="Values" type="Array(double)" value="{21}"/>
  	          </ParameterList>
  	        </ParameterList>
  	      </ParameterList>
  	    </ParameterList>
          </ParameterList>

	FIXME: Although current version implemented to skip the component section entirely
	(so it is hardwired to "Water", respectively)

       */

      // Add source info
      const std::string Sources_str = "Sources";
      const std::string Assigned_Regions_str = "Assigned Regions";
      const std::string Solute_Source_str = "Solute SOURCE";

      // FIXME: NEED TO VERIFY PHASE/COMP/SOLUTE IS VALID
      bool do_source_term = false;
      if (parameter_list.isSublist(Sources_str)) {
        const ParameterList& src_list = parameter_list.sublist(Sources_str);
        ParameterList struct_src_list;
        Array<std::string> nullList;
        PLoptions s_opt(src_list,nullList,nullList,false,true);
        const Array<std::string> src_labels = s_opt.OptLists(); // src labels
	struct_src_list.set<Array<std::string> >("sources",underscore(src_labels));
        for (int i=0; i< src_labels.size(); ++i) {
          const std::string& src_label = src_labels[i];
          ParameterList struct_src_label_list;
          const ParameterList& src_label_pl = src_list.sublist(src_label);
          Array<std::string> src_label_reqd_params; src_label_reqd_params.push_back(Assigned_Regions_str);
          PLoptions sl_opt(src_label_pl,nullList,src_label_reqd_params,false,true);
          const Array<std::string>& regions = src_label_pl.get<Array<std::string> >(Assigned_Regions_str);
          struct_src_label_list.set("regions",underscore(regions));
          const Array<std::string>& src_f_or_s_labels = sl_opt.OptLists(); // must be a known src func or "Solute SOURCE"
          for (int j=0; j<src_f_or_s_labels.size(); ++j) {
            const std::string& src_f_or_s_label = src_f_or_s_labels[j];

            if (src_f_or_s_label == Solute_Source_str) {
              // Do solute sources

              ParameterList struct_src_solute_phase_list;
              const ParameterList& src_solute_phase_pl = src_label_pl.sublist(Solute_Source_str);
              PLoptions sl_solute_phase_opt(src_solute_phase_pl,nullList,nullList,false,true);
              const Array<std::string>& src_solute_phase_labels = sl_solute_phase_opt.OptLists(); // must be known phase labels
              for (int L=0; L<src_solute_phase_labels.size(); ++L) {
                const std::string& src_solute_phase_label = src_solute_phase_labels[L];

                ParameterList struct_src_solute_phase_list;


                const ParameterList& src_solute_phase_comp_pl = src_solute_phase_pl.sublist(src_solute_phase_label);
                PLoptions src_solute_phase_comp_opt(src_solute_phase_comp_pl,nullList,nullList,false,true);
                const Array<std::string>& src_solute_phase_clabels = src_solute_phase_comp_opt.OptLists();
                for (int M=0; M<src_solute_phase_clabels.size(); ++M) {
                  const std::string& src_solute_phase_clabel = src_solute_phase_clabels[M];
                  const ParameterList& src_solute_pl = src_solute_phase_comp_pl.sublist(src_solute_phase_clabel);
                  PLoptions src_solute_opt(src_solute_pl,nullList,nullList,false,true);
                  const Array<std::string>& src_solute_labels = src_solute_opt.OptLists();

		  ParameterList struct_src_solute_phase_comp_list;


                  for (int N=0; N<src_solute_labels.size(); ++N) {
                    const std::string& src_solute_label = src_solute_labels[N];
                    const ParameterList& src_solute_func_pl = src_solute_pl.sublist(src_solute_label);

                    PLoptions src_solute_func_opt(src_solute_func_pl,nullList,nullList,false,false);
                    const Array<std::string>& src_solute_opt_params = src_solute_func_opt.OptParms();

		    ParameterList struct_src_solute_phase_comp_solute_list;
		    for (int LL=0; LL<src_solute_opt_params.size(); ++LL) {
		      const std::string& src_solute_opt_param = src_solute_opt_params[LL];
		      if (src_solute_opt_param == "Concentration Units") {
			const std::string& param = underscore(src_solute_func_pl.get<std::string>(src_solute_opt_param));
			struct_src_solute_phase_comp_solute_list.set<std::string>(underscore(src_solute_opt_param),param);
		      }
		      else {
			MyAbort("Unrecognized option for Solute source function: \""+src_solute_opt_param+"\"");
		      }
		    }
                    const Array<std::string>& src_solute_func_labels = src_solute_func_opt.OptLists();
                    if (src_solute_func_labels.size() == 1) {
                      const std::string& Amanzi_type = src_solute_func_labels[0];
                      const ParameterList& fPLin = src_solute_func_pl.sublist(Amanzi_type);
		      convert_Solute_Sources(fPLin,Amanzi_type,struct_src_solute_phase_comp_solute_list);
		    } else {
		      MyAbort("Exactly one source function allowed for "+Sources_str+"->"+src_label+"->"
			      +Solute_Source_str+"->"+src_solute_phase_clabel);
		    }
		    struct_src_solute_phase_comp_list.set(underscore(src_solute_label),struct_src_solute_phase_comp_solute_list);
		  }
		  struct_src_solute_phase_comp_list.set("tracers_with_sources",underscore(src_solute_labels));
		  struct_src_solute_phase_list.set(src_solute_phase_clabel,struct_src_solute_phase_comp_list);
                }

		struct_src_label_list.set(src_solute_phase_label,struct_src_solute_phase_list);
              }
            }
            else {
              // Must be a src function
              const std::string& Amanzi_type = src_f_or_s_label;
              const ParameterList& src_func_pl = src_label_pl.sublist(Amanzi_type);
              convert_Sources(src_func_pl,Amanzi_type,struct_src_label_list);
            }
          }
          struct_src_list.set(underscore(src_label),struct_src_label_list);
        }
	struc_list.set("source",struct_src_list);
	do_source_term = true;
      }
      struc_list.sublist("prob").set<bool>("do_source_term",do_source_term);

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
      user_derive_list.push_back(underscore("Volumetric Water Content"));
      user_derive_list.push_back(underscore("Porosity"));
      user_derive_list.push_back(underscore("Aqueous Saturation"));
      user_derive_list.push_back(underscore("Aqueous Pressure"));

      if (gravity_is_nonzero(parameter_list)) {
        user_derive_list.push_back(underscore("Hydraulic Head"));
      }
      user_derive_list.push_back(underscore("Aqueous Volumetric Flux X"));
      user_derive_list.push_back(underscore("Aqueous Volumetric Flux Y"));
#if BL_SPACEDIM==3
      user_derive_list.push_back(underscore("Aqueous Volumetric Flux Z"));
#endif

      user_derive_list.push_back(underscore("Tortuosity X"));
      user_derive_list.push_back(underscore("Tortuosity Y"));
#if BL_SPACEDIM==3
      user_derive_list.push_back(underscore("Tortuosity Z"));
#endif

      user_derive_list.push_back(underscore("Specific Storage"));
      user_derive_list.push_back(underscore("Specific Yield"));
      user_derive_list.push_back(underscore("Particle Density"));

      user_derive_list.push_back(underscore("Intrinsic Permeability X"));
      user_derive_list.push_back(underscore("Intrinsic Permeability Y"));
#if BL_SPACEDIM==3
      user_derive_list.push_back(underscore("Intrinsic Permeability Z"));
#endif

      if (struc_list.isSublist("tracer")) {
        const Array<std::string>& solute_names = struc_list.sublist("tracer").get<Array<std::string> >("tracers");
        for (int i=0; i<solute_names.size(); ++i) {
          const std::string& name = solute_names[i];
          user_derive_list.push_back(underscore(name+" Aqueous Concentration"));
          user_derive_list.push_back(underscore("Volumetric_" + name + "_Content"));
        }
      }

      if (1 && do_chem && state.HasSolidChem()) { // FIXME: All this data currently managed by Alquimia, work interface later...
        if (struc_list.isSublist("tracer")) {
          const Array<std::string>& solute_names = struc_list.sublist("tracer").get<Array<std::string> >("tracers");
          for (int i=0; i<solute_names.size(); ++i) {
            const std::string& name = solute_names[i];

            if (state.getSolid().UsingSorption()) {
              user_derive_list.push_back(underscore(name + " Sorbed Concentration"));
            }

            /*
            if (state.getSolid().HasSorptionIsotherm(name)) {
              user_derive_list.push_back(underscore(name + " Isotherm Kd"));
              if (state.getSolid().SorptionIsotherm(name).IsFreundlich()) {
                user_derive_list.push_back(underscore(name + "Isotherm Freundlich n "));
              }
              else if (state.getSolid().SorptionIsotherm(name).IsLangmuir())
              {
                user_derive_list.push_back(underscore(name + " Isotherm Langmuir b "));
              }
            }
            */
            user_derive_list.push_back(underscore(name+" Free Ion Guess"));
            user_derive_list.push_back(underscore(name+" Activity Coefficient"));
          }
        }

        if (state.getSolid().has_cation_exchange) {
          user_derive_list.push_back(underscore("Cation Exchange Capacity"));
        }

        const Array<std::string>& mineral_names = state.getSolid().mineral_names;
        for (int i=0; i<mineral_names.size(); ++i) {
          const std::string& name = mineral_names[i];
          user_derive_list.push_back(underscore(name + " Volume Fraction"));
          user_derive_list.push_back(underscore(name + " Specific Surface Area"));
        }

        const Array<std::string>& sorption_site_names = state.getSolid().sorption_site_names;
        for (int i=0; i<sorption_site_names.size(); ++i) {
          const std::string& name = sorption_site_names[i];
          user_derive_list.push_back(underscore(name + " Surface Site Density"));
        }
      }

      std::string output_str = "Output";
      if (!parameter_list.isSublist(output_str)) {
        MyAbort("Must have an \""+output_str+"\" section in the XML input file");
      }
      const ParameterList& rlist = parameter_list.sublist("Output");

      // time macros
      std::set<std::string> time_macros;
      std::string time_macros_str = "Time Macros";
      if (rlist.isSublist(time_macros_str)) {
        const ParameterList& tlist = rlist.sublist(time_macros_str);
        ParameterList tmPL;
        for (ParameterList::ConstIterator i=tlist.begin(); i!=tlist.end(); ++i) {
          std::string label = tlist.name(i);
          std::string _label = underscore(label);
          if (!tlist.isSublist(label)) {
            MyAbort("Time macro \""+label+"\" must be a ParameterList");
          }
          const ParameterList& rslist = tlist.sublist(label);
          Array<double> times, vals;

          ParameterList tPL;
          for (ParameterList::ConstIterator ii=rslist.begin(); ii!=rslist.end(); ++ii) {
            const std::string& name = underscore(rslist.name(ii));

            if (name == "Values") {
              times = rslist.get<Array<double> >("Values");
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
      }

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

            if (name == "Values") {
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
      const std::string vis_cycles_str = "Cycle Macros";
      const std::string vis_cycle_str = "Cycle Macro";
      const std::string vis_times_str = "Time Macros";
      const std::string vis_time_str = "Time Macro";
      const std::string vis_digits_str = "File Name Digits";
      const std::string vis_write_regions_str = "Write Regions";
      bool vis_vars_set = false;
      Array<std::string> visNames, vis_cMacroNames, vis_tMacroNames, wrNames;
      std::string vis_file = "plt";
      int vis_digits = 5;
      ParameterList wrPL;
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
          else if (name == vis_cycles_str || name == vis_cycle_str)
          {
            Array<std::string> vcMacros;
            if (name == vis_cycles_str) {
              const Array<std::string> tmp = vlist.get<Array<std::string> >(vis_cycles_str);
              vcMacros.resize(tmp.size());
              for (int i=0; i<tmp.size(); ++i) {
                vcMacros[i] = tmp[i];
              }
            }
            else {
              vcMacros.resize(1);
              vcMacros[0] = vlist.get<std::string>(vis_cycle_str);
            }

            for (int i=0; i<vcMacros.size(); ++i) {
              std::string label = underscore(vcMacros[i]);
              if (cycle_macros.find(label) != cycle_macros.end()) {
                vis_cMacroNames.push_back(label);
              }
              else {
                if (Teuchos::GlobalMPISession::getRank() == 0) {
                  std::cerr << "Unrecognized cycle macro in \""+vis_data_str+"\": \""
                            << vcMacros[i] << "\"" << std::endl;

                  for (std::set<std::string>::const_iterator it=cycle_macros.begin();
                       it!=cycle_macros.end(); ++it) {
                    std::cerr << *it << " " << std::endl;
                  }
                  throw std::exception();
                }
              }
            }
          }
          else if (name == vis_times_str || name == vis_time_str)
          {
            Array<std::string> vtMacros;
            if (name == vis_times_str) {
              const Array<std::string> tmp = vlist.get<Array<std::string> >(vis_times_str);
              vtMacros.resize(tmp.size());
              for (int i=0; i<tmp.size(); ++i) {
                vtMacros[i] = tmp[i];
              }
            }
            else {
              vtMacros.resize(1);
              vtMacros[0] = vlist.get<std::string>(vis_time_str);
            }
            for (int i=0; i<vtMacros.size(); ++i) {
              std::string label = underscore(vtMacros[i]);
              if (time_macros.find(label) != time_macros.end()) {
                vis_tMacroNames.push_back(label);
              }
              else {
                if (Teuchos::GlobalMPISession::getRank() == 0) {
                  std::cerr << "Unrecognized time macro in \""+vis_data_str+"\": \""
                            << vtMacros[i] << "\"" << std::endl;

                  std::cerr << "Known macros: ";
                  for (std::set<std::string>::const_iterator it=time_macros.begin();
                       it!=time_macros.end(); ++it) {
                    std::cerr << *it << " ";
                  }
                  std::cerr << std::endl;
                }
                throw std::exception();
              }
            }
          }
          else if (name == vis_digits_str) {
            vis_digits = vlist.get<int>(vis_digits_str);
            if (vis_digits<=0) {
              MyAbort("Output -> \""+vis_digits_str+"\" must be > 0");
            }
          }
	  else if (name == vis_write_regions_str) {
	    const ParameterList& wrlist = vlist.sublist(vis_write_regions_str);
	    for (ParameterList::ConstIterator i=wrlist.begin(); i!=wrlist.end(); ++i) {
	      const std::string& label = wrlist.name(i);
	      if (!wrlist.isSublist(label)) {
		Array<std::string> _regions = underscore(wrlist.get<Array<std::string> >(label));
		std::string _label = underscore(label);
		wrPL.set<Array<std::string> >(_label, _regions);
		wrNames.push_back(_label);
	      }
	      else {
		MyAbort("The \"Write Regions\" parameter list cannot take a sublist");
	      }
	    }
	  }
          else {
            MyAbort("Unrecognized entry in \""+vis_data_str+"\" parameter list: \""+name+"\"");
          }
        }
      }

      //
      amr_list.set<Array<std::string> >("vis_cycle_macros",vis_cMacroNames);
      amr_list.set<Array<std::string> >("vis_time_macros",vis_tMacroNames);
      amr_list.set<std::string>("plot_file",vis_file);
      amr_list.set<int>("plot_file_digits",vis_digits);
      if (wrNames.size() > 0) {
	amr_list.set<Array<std::string> >("write_regions",wrNames);
	amr_list.sublist("write_region") = wrPL;
	user_derive_list.insert(user_derive_list.end(),wrNames.begin(),wrNames.end());
      }

      amr_list.set<Array<std::string> >("user_derive_list",user_derive_list);
      if (!vis_vars_set) {
        for (int j=0; j<user_derive_list.size(); ++j) {
          visNames.push_back(underscore(user_derive_list[j]));
        }
      }
      amr_list.set<Array<std::string> >("derive_plot_vars",visNames);

      // chk data
      const std::string chk_data_str = "Checkpoint Data";
      const std::string chk_file_str = "File Name Base";
      const std::string chk_cycles_str = "Cycle Macros";
      const std::string chk_cycle_str = "Cycle Macro";
      const std::string chk_digits_str = "File Name Digits";
      std::string chk_file = "chk";
      int check_digits = 5;
      Array<std::string> chk_cMacroNames;
      if (rlist.isSublist(chk_data_str)) {
        const ParameterList& clist = rlist.sublist(chk_data_str);
        for (ParameterList::ConstIterator i=clist.begin(); i!=clist.end(); ++i)
        {
          const std::string& name = clist.name(i);
          if (name == chk_file_str) {
            chk_file = clist.get<std::string>(chk_file_str);
          }
          else if (name == chk_cycles_str || name == chk_cycle_str)
          {
            Array<std::string> ccMacros;
            if (name == chk_cycles_str) {
              const Array<std::string> tmp = clist.get<Array<std::string> >(chk_cycles_str);
              ccMacros.resize(tmp.size());
              for (int i=0; i<tmp.size(); ++i) {
                ccMacros[i] = tmp[i];
              }
            }
            else {
              ccMacros.resize(1);
              ccMacros[0] = clist.get<std::string>(chk_cycle_str);
            }
            for (int i=0; i<ccMacros.size(); ++i) {
              std::string label = underscore(ccMacros[i]);
              if (cycle_macros.find(label) != cycle_macros.end()) {
                chk_cMacroNames.push_back(label);
              }
              else {
                if (Teuchos::GlobalMPISession::getRank() == 0) {
                  std::cerr << "Unrecognized cycle macro in \""+chk_data_str+"\": \""
                            << ccMacros[i] << "\"" << std::endl;

                  for (std::set<std::string>::const_iterator it=cycle_macros.begin();
                       it!=cycle_macros.end(); ++it) {
                    std::cerr << *it << " " << std::endl;
                  }
                  throw std::exception();
                }
              }
            }
          }
          else if (name == chk_digits_str) {
            check_digits = clist.get<int>(chk_digits_str);
            if (check_digits<=0) {
              MyAbort("Output -> \""+chk_digits_str+"\" must be > 0");
            }
          }
          else {
            MyAbort("Unrecognized entry in \""+chk_data_str+"\" parameter list: \""+name+"\"");
          }
        }
      }
      amr_list.set<Array<std::string> >("chk_cycle_macros",chk_cMacroNames);
      amr_list.set<std::string>("check_file",chk_file);
      amr_list.set<int>("chk_file_digits",check_digits);


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
            else if (functional == "Observation Data: Mean") {
              sublist.set("obs_type","average");
            }
            else if (functional == "Observation Data: Point") {
              sublist.set("obs_type","point_sample");

              // Check our list of regions to ensure it exists and is the correct type
              const ParameterList& lregion =
                parameter_list.sublist("Regions").sublist(region_name);

              if (!lregion.isSublist("Region: Point")) {
                std::cerr << label << " is a point observation and "
                          << region_name << " is not a point region.\n";
                throw std::exception();
              }
            }
            else if (functional == "Observation Data: Peak Value") {
              sublist.set("obs_type","peak_value");
            }
            else {
              MyAbort("Unsupported functional for observation \""+label+"\": functional = \""+functional+"\"");
            }
            sublist.set("region",_region_name);

            std::string Time_Macros_str = "Time Macros";
            std::string Cycle_Macros_str = "Cycle Macros";
            if (rslist.isParameter(Time_Macros_str)) {
              const Array<std::string>& _timeMacros = underscore(rslist.get<Array<std::string> >(Time_Macros_str));
	      for (int k=0; k<_timeMacros.size(); ++k) {
		const std::string& _timeMacro = _timeMacros[k];
		if (time_macros.find(_timeMacro) == time_macros.end()) {
		  std::cerr << "Unrecognized time macro: \"" << AMR_to_Amanzi_label_map[_timeMacro]
			    << "\" for observation data: \"" << label << "\"" << std::endl;
		  throw std::exception();
		}
              }
              sublist.set("time_macros",_timeMacros);
            }
            else if (rslist.isParameter(Cycle_Macros_str)) {
              const Array<std::string>& _cycleMacros = underscore(rslist.get<Array<std::string> >(Cycle_Macros_str));
	      for (int k=0; k<_cycleMacros.size(); ++k) {
		const std::string& _cycleMacro = _cycleMacros[k];
		if (cycle_macros.find(_cycleMacro) == cycle_macros.end()) {
		  std::cerr << "Unrecognized cycle macro: \"" << AMR_to_Amanzi_label_map[_cycleMacro]
			    << "\" for observation data: \"" << label << "\"" << std::endl;
		  throw std::exception();
		}
              }
              sublist.set("cycle_macros",_cycleMacros);
            }
            else {
              std::cerr << "Must specify either time or cycle macro for observation data: \""
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
              std::cerr << "Available derived variables: ";
              for (int k=0; k<user_derive_list.size(); ++k) {
                std::cerr << "\"" << AMR_to_Amanzi_label_map[user_derive_list[k]] << "\"";
                if (k<=user_derive_list.size()-2) {
                  std::cerr << ", ";
                }
              }
              std::cerr << std::endl;
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
      bool do_tracer_advection, do_tracer_diffusion, do_chem;
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
      convert_to_structured_control(parameter_list,struc_list, do_tracer_advection, do_tracer_diffusion, do_chem);
      //
      // Regions
      //
      convert_to_structured_region(parameter_list, struc_list);
      //
      // State
      //
      StateDef stateDef(parameter_list);
      convert_to_structured_state(parameter_list, struc_list, stateDef,
                                  do_tracer_advection, do_tracer_diffusion, do_chem);
      //
      // Materials
      //
      convert_to_structured_material(parameter_list, struc_list, stateDef,
                                     do_tracer_advection, do_tracer_diffusion);
      //
      // Output
      //
      convert_to_structured_output(parameter_list,struc_list,stateDef,do_chem);

      ParameterList& prob_out_list = struc_list.sublist("prob");
      prob_out_list.set("do_tracer_advection",do_tracer_advection);
      prob_out_list.set("do_tracer_diffusion",do_tracer_diffusion);
      prob_out_list.set("gravity",gravity_magnitude(parameter_list));
      prob_out_list.set("gravity_dir",gravity_dir(parameter_list));
      prob_out_list.set("z_location",z_location(parameter_list));

      std::string dump_str = "Structured Native Input File";
      if (parameter_list.isParameter(dump_str)) {
        struc_list.set<std::string>("dump_parmparse_table",parameter_list.get<std::string>(dump_str));
      }

      return struc_list;
    }

  }
}
