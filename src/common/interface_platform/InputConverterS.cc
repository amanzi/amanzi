/*
  Input Converter

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson (jnjohnson@lbl.gov)
*/

#ifdef ENABLE_Structured

#include <float.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "BoxLib.H"
#include "Utility.H"

#include "errors.hh"
#include "exceptions.hh"
#include "StringExt.hh"
#include "InputConverterS.hh"

namespace Amanzi {
namespace AmanziInput {

// For saving a bit of typing.
using namespace std;
using namespace xercesc;

// Internally-used stuff.
namespace {

// This converts a time (contained in a string) to its value in seconds.
string ConvertTimeToSeconds(const string& time_string)
{
  vector<string> tokens = split(time_string, " ,;");
  if (tokens.size() == 1) {
    return tokens[0];
  } else {
    double value = atof(tokens[0].c_str());
    double factor;
    switch(tokens[1][0])
    {
      case 'h': factor = 3600.0; break;
      case 'd': factor = 3600.0 * 24.0; break;
      case 'y': factor = 3600.0 * 24.0 * 365.25; break;
      default: factor = 1.0;
    }
    value *= factor;
    stringstream s;
    s.precision(15);
    s << value;
    return s.str();
  }
}

typedef list<ParmParse::PP_entry>::iterator list_iterator;
typedef list<ParmParse::PP_entry>::const_iterator const_list_iterator;

// Trims strings and replaces spaces with underscores in a string.
string MangleString(const string& s)
{
  string s1 = s;
  trim(s1);
  for (int i = 0; i < s.length(); ++i)
  {
    if (s1[i] == ' ')
      s1[i] = '_';
  }
  return s1;
}

// Construct a ParmParse prefix name from sets of strings.
string MakePPPrefix(const string& s1)
{
  return MangleString(s1);
}

string MakePPPrefix(const string& s1, const string& s2)
{
  return MangleString(s1) + string(".") + MangleString(s2);
}

string MakePPPrefix(const string& s1, const string& s2, const string& s3)
{
  return MangleString(s1) + string(".") + MangleString(s2) + string(".") + MangleString(s3);
}

string MakePPPrefix(const string& s1, const string& s2, const string& s3, const string& s4)
{
  return MangleString(s1) + string(".") + MangleString(s2) + string(".") + MangleString(s3) + string(".") + MangleString(s4);
}

string MakePPPrefix(const string& s1, const string& s2, const string& s3, const string& s4, const string& s5)
{
  return MangleString(s1) + string(".") + MangleString(s2) + string(".") + MangleString(s3) + string(".") + MangleString(s4) + string(".") + MangleString(s5);
}

// Construct a ParmParse entry from sets of strings/values.
list<string> MakePPEntry(const string& s1)
{
  list<string> pp;
  pp.push_back(MangleString(s1));
  return pp;
}

list<string> MakePPEntry(double d)
{
  list<string> pp;
  stringstream s;
  s << setprecision(15);
  s << d;
  pp.push_back(s.str());
  return pp;
}

list<string> MakePPEntry(int i)
{
  list<string> pp;
  stringstream s;
  s << i;
  pp.push_back(s.str());
  return pp;
}

list<string> MakePPEntry(const vector<string>& ss)
{
  list<string> pp;
  for (size_t i = 0; i < ss.size(); ++i)
    pp.push_back(MangleString(ss[i]));
  return pp;
}

list<string> MakePPEntry(const vector<double>& ds)
{
  list<string> pp;
  for (size_t i = 0; i < ds.size(); ++i)
  {
    stringstream s;
    s << setprecision(15);
    s << ds[i];
    pp.push_back(s.str());
  }
  return pp;
}

list<string> MakePPEntry(const vector<int>& is)
{
  list<string> pp;
  for (size_t i = 0; i < is.size(); ++i)
  {
    stringstream s;
    s << is[i];
    pp.push_back(s.str());
  }
  return pp;
}


// Shortcut for putting entries into tables.
void AddToTable(list<ParmParse::PP_entry>& table,
                const string& entry_name,
                list<string> entry)
{
  table.push_back(ParmParse::PP_entry(entry_name, entry));
}

} // anonymous namespace


// Helper for parsing a mechanical property.
bool InputConverterS::ParseMechProperty_(DOMElement* mech_prop_node, 
                                         const string& material_name, 
                                         const string& property_name,
                                         list<ParmParse::PP_entry>& table,
                                         bool required)
{
  bool found;
  bool non_zero_diffusion = false;
  DOMElement* property = GetChildByName_(mech_prop_node, property_name, found, required);
  if (found)
  {
    if (property_name == "dispersion_tensor") // Weirdo!
    {
      string type = GetAttributeValueS_(property, "type");
      if (type == "uniform_isotropic")
      {
        string alpha_l = GetAttributeValueS_(property, "alpha_l");
        string alpha_t = GetAttributeValueS_(property, "alpha_t");
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "alpha_l"),
                   MakePPEntry(alpha_l));
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "alpha_t"),
                   MakePPEntry(alpha_t));

	non_zero_diffusion = (atof(alpha_l.c_str())!=0)
	  || (atof(alpha_t.c_str())!=0);
      }
      else if (type == "burnett_frind")
      {
        string alpha_l = GetAttributeValueS_(property, "alpha_l");
        string alpha_th = GetAttributeValueS_(property, "alpha_th");
        string alpha_tv = GetAttributeValueS_(property, "alpha_tv");
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "alpha_l"),
                   MakePPEntry(alpha_l));
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "alpha_th"),
                   MakePPEntry(alpha_th));
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "alpha_tv"),
                   MakePPEntry(alpha_tv));
	non_zero_diffusion = (atof(alpha_l.c_str())!=0)
	  || (atof(alpha_th.c_str())!=0)
	  || (atof(alpha_tv.c_str()));
      }
      else if (type == "lichtner_kelkar_robinson")
      {
        string alpha_lh = GetAttributeValueS_(property, "alpha_lh");
        string alpha_lv = GetAttributeValueS_(property, "alpha_lv");
        string alpha_th = GetAttributeValueS_(property, "alpha_th");
        string alpha_tv = GetAttributeValueS_(property, "alpha_tv");
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "alpha_lh"),
                   MakePPEntry(alpha_lh));
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "alpha_lv"),
                   MakePPEntry(alpha_lv));
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "alpha_th"),
                   MakePPEntry(alpha_th));
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "alpha_tv"),
                   MakePPEntry(alpha_tv));
	non_zero_diffusion = (atof(alpha_lh.c_str())!=0)
	  || (atof(alpha_lv.c_str())!=0)
	  || (atof(alpha_th.c_str())!=0)
	  || (atof(alpha_tv.c_str()));
      }
      else if (type == "file")
      {
        string filename = GetAttributeValueS_(property, "filename");
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "filename"),
                   MakePPEntry(filename));
      }
      else
        ThrowErrorIllformed_("materials->mechanical_properties", "type", property_name);
      AddToTable(table, MakePPPrefix("rock", material_name, property_name, "type"),
                 MakePPEntry(type));
    }
    else
    {
      string value = GetAttributeValueS_(property, "value", TYPE_NUMERICAL, false);
      string type = GetAttributeValueS_(property, "type", TYPE_NUMERICAL, false);
      if (!value.empty() && ((property_name != "porosity") || (type != "gslib")))
      {
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "vals"),
                   MakePPEntry(value));
        AddToTable(table, MakePPPrefix("rock", material_name, property_name, "distribution_type"),
                   MakePPEntry("uniform"));
      }
      else
      {
        if (!type.empty())
        {
          if (type == "file")
          {
            string filename = GetAttributeValueS_(property, "filename");
            AddToTable(table, MakePPPrefix("rock", material_name, property_name, "filename"),
                       MakePPEntry(filename));
            AddToTable(table, MakePPPrefix("rock", material_name, property_name, "type"),
                       MakePPEntry("file"));
          }
          else if ((property_name == "porosity") && (type == "gslib")) // special considerations!
          {
            string parameter_file = GetAttributeValueS_(property, "parameter_file");
            AddToTable(table, MakePPPrefix("rock", material_name, property_name, "parameter_file"),
                       MakePPEntry(parameter_file));
            AddToTable(table, MakePPPrefix("rock", material_name, property_name, "type"),
                       MakePPEntry("gslib"));

            string data_file = GetAttributeValueS_(property, "data_file", TYPE_NONE, false);
            if (!data_file.empty())
            {
              AddToTable(table, MakePPPrefix("rock", material_name, property_name, "data_file"),
                         MakePPEntry(data_file));
            }
          }
        }
      }
    }
  }
  return non_zero_diffusion;
}

InputConverterS::InputConverterS(const string& input_filename):
  InputConverter(input_filename)
{
}

InputConverterS::InputConverterS(const string& input_filename, xercesc::DOMDocument* input_doc):
  InputConverter(input_filename, input_doc)
{
}

InputConverterS::~InputConverterS()
{
}

void InputConverterS::ParseUnits_()
{
  // Units are supposedly not supported yet. I guess we'll find out!
}

void InputConverterS::ParseDefinitions_()
{
  list<ParmParse::PP_entry> table;
  bool found;
  DOMNode* constants_block = GetUniqueElementByTagsString_("definitions, constants", found);
  if (found)
  {
    vector<DOMNode*> constants = GetChildren_(constants_block, "constant", found);
    if (found)
    {
      for (size_t i = 0; i < constants.size(); ++i)
      {
        DOMElement* constant = static_cast<DOMElement*>(constants[i]);
        string name = GetAttributeValueS_(constant, "name");
        string type = GetAttributeValueS_(constant, "type");
        string value = GetAttributeValueS_(constant, "value");
        if (type == "time")
          labeled_times_[name] = value;
        else if (type == "numerical")
          labeled_numbers_[name] = value;
        else if (type == "area_mass_flux")
          labeled_area_mass_fluxes_[name] = value;
      }
    }
  }

  DOMNode* macros = GetUniqueElementByTagsString_("definitions, macros", found);
  if (found)
  {
    vector<DOMNode*> time_macros = GetChildren_(macros, "time_macro", found);
    vector<string> time_macro_names;
    if (found)
    {
      for (size_t i = 0; i < time_macros.size(); ++i)
      {
        DOMElement* time_macro = static_cast<DOMElement*>(time_macros[i]);
        string macro_name = GetAttributeValueS_(time_macro, "name");
        time_macro_names.push_back(macro_name);

        // Before we look for specific times, check for other stuff.
        string start = GetChildValueS_(time_macro, "start", found);
        if (found)
        {
          // We've got one of the interval-based time macros.
          // Get the other required elements.
          string timestep_interval = GetChildValueS_(time_macro, "timestep_interval", found, true);
          string stop = GetChildValueS_(time_macro, "stop", found, true);
          
          char* endptr = const_cast<char*>(start.c_str());
          strtod(start.c_str(), &endptr);
          if (endptr == start.c_str()) // Not parsed!
            ThrowErrorIllformed_("definitions->macros", "start", "time_macro");
          endptr = const_cast<char*>(stop.c_str());
          strtod(stop.c_str(), &endptr);
          if (endptr == stop.c_str()) // Not parsed!
            ThrowErrorIllformed_("definitions->macros", "stop", "time_macro");
          endptr = const_cast<char*>(timestep_interval.c_str());
          strtod(timestep_interval.c_str(), &endptr);
          if (endptr == timestep_interval.c_str()) // Not parsed!
            ThrowErrorIllformed_("definitions->macros", "timestep_interval", "time_macro");

          // Shove this macro into our table.
          AddToTable(table, MakePPPrefix("amr", "time_macro", macro_name, "type"),
                                         MakePPEntry("period"));
          AddToTable(table, MakePPPrefix("amr", "time_macro", macro_name, "start"),
                                         MakePPEntry(start));
          AddToTable(table, MakePPPrefix("amr", "time_macro", macro_name, "stop"),
                                         MakePPEntry(stop));
          AddToTable(table, MakePPPrefix("amr", "time_macro", macro_name, "period"),
                                         MakePPEntry(timestep_interval));
        }
        else
        {
          // We're just looking for times.
          vector<DOMNode*> time_nodes = GetChildren_(time_macro, "time", found, true);
          vector<string> times;
          for (size_t j = 0; j < time_nodes.size(); ++j)
            times.push_back(XMLString::transcode(time_nodes[j]->getTextContent()));

          AddToTable(table, MakePPPrefix("amr", "time_macro", macro_name, "type"),
                                         MakePPEntry("times"));
          AddToTable(table, MakePPPrefix("amr", "time_macro", macro_name, "times"),
                                         MakePPEntry(times));
        }
      }
    }
    // List of all time macro names.
    AddToTable(table, MakePPPrefix("amr", "time_macros"), 
                                   MakePPEntry(time_macro_names));

    vector<DOMNode*> cycle_macros = GetChildren_(macros, "cycle_macro", found);
    vector<string> cycle_macro_names;
    if (found)
    {
      for (size_t i = 0; i < cycle_macros.size(); ++i)
      {
        DOMElement* cycle_macro = static_cast<DOMElement*>(cycle_macros[i]);
        string macro_name = GetAttributeValueS_(cycle_macro, "name");
        cycle_macro_names.push_back(macro_name);
        vector<string> cycles;

        string start = GetChildValueS_(cycle_macro, "start", found, true);
        string timestep_interval = GetChildValueS_(cycle_macro, "timestep_interval", found, true);
        string stop = GetChildValueS_(cycle_macro, "stop", found, true);
         
        // Verify the entries.
        char* endptr = const_cast<char*>(start.c_str());
        strtol(start.c_str(), &endptr, 10);
        if (endptr == start.c_str()) // Not parsed!
          ThrowErrorIllformed_("definitions->macros", "start", "cycle_macro");
        endptr = const_cast<char*>(stop.c_str());
        strtol(stop.c_str(), &endptr, 10);
        if (endptr == stop.c_str()) // Not parsed!
          ThrowErrorIllformed_("definitions->macros", "stop", "cycle_macro");
        endptr = const_cast<char*>(timestep_interval.c_str());
        strtol(timestep_interval.c_str(), &endptr, 10);
        if (endptr == timestep_interval.c_str()) // Not parsed!
          ThrowErrorIllformed_("definitions->macros", "timestep_interval", "cycle_macro");

        // Shove this macro into our table.
        AddToTable(table, MakePPPrefix("amr", "cycle_macro", macro_name, "type"),
                                       MakePPEntry("period"));
        AddToTable(table, MakePPPrefix("amr", "cycle_macro", macro_name, "start"),
                                       MakePPEntry(start));
        AddToTable(table, MakePPPrefix("amr", "cycle_macro", macro_name, "stop"),
                                       MakePPEntry(stop));
        AddToTable(table, MakePPPrefix("amr", "cycle_macro", macro_name, "period"),
                                       MakePPEntry(timestep_interval));
      }
    }
    // List of all cycle macro names.
    AddToTable(table, MakePPPrefix("amr", "cycle_macros"), 
                                   MakePPEntry(cycle_macro_names));

    // FIXME: variable_macro not yet supported.
  }

  ParmParse::appendTable(table);
}

void InputConverterS::ParseExecutionControls_()
{
  list<ParmParse::PP_entry> table;
  bool found;
  DOMNode* echo_node = GetUniqueElementByTagsString_("echo_translated_input", found);
  if (found) {
    DOMElement* echo_elt = static_cast<DOMElement*>(echo_node);
    string echo_file_name = GetAttributeValueS_(echo_elt, "file_name", TYPE_NONE, true);
    string echo_file_format = GetAttributeValueS_(echo_elt, "format", TYPE_NONE, false);
    if (echo_file_format.empty()) {
      echo_file_format = "native";
    }
    if (echo_file_format != "native") {
      Errors::Message msg;
      msg << "An error occurred during parsing \"echo_translated_input\"\n";
      msg << "\"native\" (the default) is currently the only supported value for \"format\"\n";
      msg << "Please correct and try again.\n" ;
      Exceptions::amanzi_throw(msg);
    }
    else {
      echo_file_format = "structured_" + echo_file_format;
    }

    AddToTable(table, MakePPPrefix("translated_input", "file_name"), MakePPEntry(echo_file_name));
    AddToTable(table, MakePPPrefix("translated_input", "format"), MakePPEntry(echo_file_format));
  }

  DOMNode* controls_block = GetUniqueElementByTagsString_("execution_controls", found);
  if (found)
  {
    DOMNode* defaults = GetUniqueElementByTagsString_("execution_controls, execution_control_defaults", found);
    if (!found) {
      ThrowErrorMisschild_("execution_controls", "execution_control_defaults");
    }

    DOMElement* def = static_cast<DOMElement*>(defaults);

    map<string, string> default_vals;
    default_vals["max_cycles"] = "10000000"; // Default for the defaults

    // Reset defaults
    default_vals["init_dt"]          = GetAttributeValueS_(def, "init_dt",          TYPE_NUMERICAL, true);
    default_vals["max_dt"]           = GetAttributeValueS_(def, "max_dt",           TYPE_NUMERICAL, true);
    default_vals["reduction_factor"] = GetAttributeValueS_(def, "reduction_factor", TYPE_NUMERICAL, true);
    default_vals["increase_factor"]  = GetAttributeValueS_(def, "increase_factor",  TYPE_NUMERICAL, true);
    default_vals["mode"]             = GetAttributeValueS_(def, "mode",             TYPE_NUMERICAL, true);
    default_vals["method"]           = GetAttributeValueS_(def, "method",           TYPE_NUMERICAL, true);
    string this_max_cycles           = GetAttributeValueS_(def, "max_cycles",       TYPE_NUMERICAL, false);
    if (!this_max_cycles.empty()) {
      default_vals["max_cycles"] = this_max_cycles;
    }

    vector<DOMNode*> control_nodes = GetChildren_(controls_block, "execution_control", found);
    if (!found) {
      ThrowErrorMisschild_("execution_controls", "execution_control");
    }
    else {

      const size_t controls_size = 9;
      string control_names[controls_size] = {"start","end","init_dt","mode","max_dt","reduction_factor",
					     "increase_factor","method","max_cycles"};

      size_t ncn= control_nodes.size();
      vector<map<string, string> > exec_control(ncn);
      for (size_t i = 0; i < ncn; ++i)
      {
	exec_control[i] = default_vals;
        DOMElement* control_elt = static_cast<DOMElement*>(control_nodes[i]);
	for (size_t j=0; j<controls_size; ++j) {
	  string str = GetAttributeValueS_(control_elt, control_names[j].c_str(), TYPE_NONE, false);
	  if (!str.empty())  exec_control[i][control_names[j]] = str;
	}
      }

      for (int i=1; i<ncn; ++i) {
	for (map<string,string>::const_iterator it=exec_control[i].begin(); it!=exec_control[i].end(); ++it) {
	  if (it->first == "start") {
	    exec_control[i-1]["end"] = it->second;
	  }
	}
      }

      // If control is a time quantity, look up if labeled time, and convert to seconds
      for (int i=0; i<ncn; ++i) {
	for (map<string,string>::iterator it=exec_control[i].begin(); it!=exec_control[i].end(); ++it) {
	  if (it->first == "start"
	      || it->first == "end"
	      || it->first == "init_dt"
	      || it->first == "max_dt")
	  {
	    map<string, string>::const_iterator iter = labeled_times_.find(it->second);
	    if (iter != labeled_times_.end())
	      it->second = ConvertTimeToSeconds(iter->second);
	    else
	      it->second = ConvertTimeToSeconds(it->second);
	  }
	}
      }

      AddToTable(table, MakePPPrefix("max_step"), MakePPEntry(exec_control[ncn-1]["max_cycles"]));
      AddToTable(table, MakePPPrefix("stop_time"), MakePPEntry(exec_control[ncn-1]["end"]));
      AddToTable(table, MakePPPrefix("strt_time"), MakePPEntry(exec_control[0]["start"]));

      vector<string> ecnames(ncn);
      int ndigits = (int) (std::log10(std::max(size_t(1),ncn-1)) + .0001) + 1;
      for (int i=0; i<exec_control.size(); ++i) {
	ecnames[i] = BoxLib::Concatenate("exec_control_",i,ndigits);
	for (map<string,string>::const_iterator it=exec_control[i].begin(); it!=exec_control[i].end(); ++it) {
	  AddToTable(table, MakePPPrefix("exec_control",ecnames[i],it->first), MakePPEntry(it->second));
	}
      }
      AddToTable(table, MakePPPrefix("exec_controls"), MakePPEntry(ecnames));
    }
  }

  // Restart
  DOMNode* restart = GetUniqueElementByTagsString_("execution_controls, restart", found);
  if (found)
  {
    string restart_file = GetAttributeValueS_(static_cast<DOMElement*>(restart), "file");
    AddToTable(table, MakePPPrefix("amr", "restart"), MakePPEntry(restart_file));
  }

  // Verbosity level(s).
  DOMNode* verbosity = GetUniqueElementByTagsString_("execution_controls, verbosity", found);
  string level = "medium";
  if (found)
  {
    DOMElement* verb = static_cast<DOMElement*>(verbosity);
    level = GetAttributeValueS_(verb, "level");
    transform(level.begin(), level.end(), level.begin(), ::tolower);
  }

  // Each level of verbosity corresponds to several verbosity parameters 
  // for Amanzi-S.
  int prob_v(0), mg_v(0), cg_v(0), amr_v(0), diffuse_v(0), io_v(0), fab_v(0);
  if (level == "none") {
    prob_v = 0; mg_v = 0; cg_v = 0; amr_v = 0; diffuse_v = 0; io_v = 0; fab_v = 0;
  }
  else if (level == "low") {
    prob_v = 1; mg_v = 0; cg_v = 0; amr_v = 1;  diffuse_v = 0; io_v = 0; fab_v = 0;
  }
  else if (level == "medium") {
    prob_v = 1; mg_v = 0; cg_v = 0; amr_v = 2;  diffuse_v = 0; io_v = 0; fab_v = 0;
  }
  else if (level == "high") {
    prob_v = 2; mg_v = 0; cg_v = 0; amr_v = 3;  diffuse_v = 0; io_v = 0; fab_v = 0;
  }
  else if (level == "extreme") {
    prob_v = 3; mg_v = 2; cg_v = 2; amr_v = 3;  diffuse_v = 1; io_v = 1; fab_v = 1;
  }

  AddToTable(table, MakePPPrefix("prob",    "v"), MakePPEntry(prob_v));
  AddToTable(table, MakePPPrefix("mg",      "v"), MakePPEntry(mg_v));
  AddToTable(table, MakePPPrefix("cg",      "v"), MakePPEntry(cg_v));
  AddToTable(table, MakePPPrefix("amr",     "v"), MakePPEntry(amr_v));
  AddToTable(table, MakePPPrefix("diffuse", "v"), MakePPEntry(diffuse_v));
  AddToTable(table, MakePPPrefix("io",      "v"), MakePPEntry(io_v));
  AddToTable(table, MakePPPrefix("fab",     "v"), MakePPEntry(fab_v));

  ParmParse::appendTable(table);
}

void InputConverterS::ParseNumericalControls_(const string& flow_model)
{
  list<ParmParse::PP_entry> table;
  bool found;

  // AMR controls - These must be specified to backend code, declare/add whether or not specified in xml
  map<string,string> amr_controls;
  amr_controls["amr_levels"]                = "1";
  amr_controls["refinement_ratio"]          = "2";
  amr_controls["do_amr_subcycling"]         = "true";
  amr_controls["regrid_interval"]           = "1";
  amr_controls["blocking_factor"]           = "2";
  amr_controls["number_error_buffer_cells"] = "1";
  amr_controls["max_grid_size"]             = "64";

  // Common controls -- currently empty.
  GetUniqueElementByTagsString_("numerical_controls, common_controls", found);

#define CONTROLS_ARE_ATTRIBUTES // Otherwise, controls are children
#undef CONTROLS_ARE_ATTRIBUTES

  // Structured controls.
  DOMNode* structured_controls = GetUniqueElementByTagsString_("numerical_controls, structured_controls", found);
  if (found)
  {
    bool found;
    MemoryManager mm;

    // Time step controls.
    map<string,string> ts_controls;
    ts_controls["min_iterations"]                     = "10";
    ts_controls["max_iterations"]                     = "15";
    ts_controls["limit_iterations"]                   = "20";
    ts_controls["min_iterations_2"]                   = "2";
    ts_controls["time_step_increase_factor"]          = "1.6";
    ts_controls["time_step_increase_factor_2"]        = "10";
    ts_controls["max_consecutive_failures_1"]         = "3";
    ts_controls["time_step_retry_factor_1"]           = "0.2";
    ts_controls["max_consecutive_failures_2"]         = "4";
    ts_controls["time_step_retry_factor_2"]           = "0.01";
    ts_controls["time_step_retry_factor_f"]           = "0.001";
    ts_controls["max_consecutive_success"]            = "0";
    ts_controls["extra_time_step_increase"]           = "10";
    ts_controls["limit_function_evals"]               = "1000000";
    ts_controls["do_grid_sequence"]                   = "true";
    ts_controls["grid_sequence_new_level_dt_factor"]  = "1";
#ifdef CONTROLS_ARE_ATTRIBUTES
    vector<DOMNode*> tsc = GetChildren_(structured_controls, "str_time_step_controls", found);
    if (found) {
      for (size_t i = 0; i < tsc.size(); ++i)
      {
        DOMElement* control_elt = static_cast<DOMElement*>(tsc[i]);
	for (map<string,string>::iterator it=ts_controls.begin(); it!=ts_controls.end(); ++it) {
	  string str = GetAttributeValueS_(control_elt, it->first.c_str(), TYPE_NONE, false);
	  if (!str.empty())  it->second = str;
	}
      }
    }
#else
    DOMElement* ssc = GetChildByName_(structured_controls,"str_time_step_controls", found);
    if (found) {
      DOMNodeList* children = ssc->getChildNodes();
      for (int i=0; i<children->getLength(); ++i) {
	DOMNode* inode = children->item(i);
	if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
	char* cname = mm.transcode(inode->getNodeName());
	for (map<string,string>::iterator it=ts_controls.begin(); it!=ts_controls.end(); ++it) {
	  if (it->first == cname) {
	    it->second = mm.transcode(inode->getTextContent());
	  }
	}
      }
    }
#endif
    if (ts_controls.size() > 0) {
      for (map<string,string>::const_iterator it=ts_controls.begin(); it!=ts_controls.end(); ++it) {
	string s_parameter_name = "steady_" + it->first;
	AddToTable(table,MakePPPrefix("prob", (s_parameter_name).c_str() ), MakePPEntry(it->second));
      }
    }

    // Flow controls.
    map<string,string> flow_controls;
    flow_controls["petsc_options_file"]           = "";
    flow_controls["max_ls_iterations"]            = "10";
    flow_controls["ls_reduction_factor"]          = "0.1";
    flow_controls["min_ls_factor"]                = "1.e-8";
    flow_controls["ls_acceptance_factor"]         = "1.4";
    flow_controls["monitor_line_search"]          = "0";
    flow_controls["monitor_linear_solve"]         = "0";
    flow_controls["use_fd_jac"]                   = "true";
    flow_controls["perturbation_scale_for_J"]     = "1.e-8";
    flow_controls["use_dense_Jacobian"]           = "false";
    flow_controls["pressure_maxorder"]            = "3";
    flow_controls["scale_solution_before_solve"]  = "true";
    flow_controls["semi_analytic_J"]              = "false";
    flow_controls["atmospheric_pressure"]         = "101325";
    flow_controls["gravity"]                      = "9.807";
    flow_controls["gravity_dir"]                  = (dim_ == 2 ? "1" : "2");
    flow_controls["domain_thickness"]             = "1";

    if (flow_model == "saturated" || flow_model == "constant") {
      flow_controls["rel_perm_method"]            = "other-harmonic_average";
    }
    else {
      flow_controls["rel_perm_method"]            = "upwind-darcy_velocity";
    }

#ifdef CONTROLS_ARE_ATTRIBUTES
    vector<DOMNode*> flc = GetChildren_(structured_controls, "str_flow_controls", found);
    if (found) {
      for (size_t i = 0; i < flc.size(); ++i)
      {
        DOMElement* control_elt = static_cast<DOMElement*>(flc[i]);
	for (map<string,string>::iterator it=flow_controls.begin(); it!=flow_controls.end(); ++it) {
	  string str = GetAttributeValueS_(control_elt, it->first.c_str(), TYPE_NONE, false);
	  if (!str.empty())  it->second = str;
	}
      }
    }
#else
    DOMElement* stc = GetChildByName_(structured_controls,"str_flow_controls", found);
    if (found) {
      DOMNodeList* children = stc->getChildNodes();
      for (int i=0; i<children->getLength(); ++i) {
	DOMNode* inode = children->item(i);
	if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
	char* cname = mm.transcode(inode->getNodeName());
	for (map<string,string>::iterator it=flow_controls.begin(); it!=flow_controls.end(); ++it) {
	  if (it->first == cname) {
	    it->second = mm.transcode(inode->getTextContent());
	  }
	}
      }
    }
#endif

    if (flow_controls.size() > 0) {
      for (map<string,string>::const_iterator it=flow_controls.begin(); it!=flow_controls.end(); ++it) {
	string s_parameter_name;
	if ( (it->first=="petsc_options_file")
	     || (it->first=="gravity")
	     || (it->first=="gravity_dir")
	     || (it->first=="domain_thickness") )
	{
	  s_parameter_name = it->first;
	} else {
	  s_parameter_name = "richard_" + it->first;
	}
	AddToTable(table,MakePPPrefix("prob", (s_parameter_name).c_str() ), MakePPEntry(it->second));
      }
    }

    // Transport controls.
    map<string,string> transport_controls;
    transport_controls["max_n_subcycle_transport"] = "20";
    transport_controls["cfl"]                      = "1";
#ifdef CONTROLS_ARE_ATTRIBUTES
    vector<DOMNode*> trc = GetChildren_(structured_controls, "str_transport_controls", found);
    if (found) {
      for (size_t i = 0; i < trc.size(); ++i)
      {
        DOMElement* control_elt = static_cast<DOMElement*>(trc[i]);
	for (map<string,string>::iterator it=transport_controls.begin(); it!=transport_controls.end(); ++it) {
	  string str = GetAttributeValueS_(control_elt, it->first.c_str(), TYPE_NONE, false);
	  if (!str.empty())  it->second = str;
	}
      }
    }
#else
    DOMElement* trc = GetChildByName_(structured_controls,"str_transport_controls", found);
    if (found) {
      DOMNodeList* children = trc->getChildNodes();
      for (int i=0; i<children->getLength(); ++i) {
	DOMNode* inode = children->item(i);
	if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
	char* cname = mm.transcode(inode->getNodeName());
	for (map<string,string>::iterator it=transport_controls.begin(); it!=transport_controls.end(); ++it) {
	  if (it->first == cname) {
	    it->second = mm.transcode(inode->getTextContent());
	  }
	}
      }
    }
#endif
    if (transport_controls.size() > 0) {
      for (map<string,string>::const_iterator it=transport_controls.begin(); it!=transport_controls.end(); ++it) {
	string s_parameter_name = it->first;
	if (s_parameter_name == "cfl") {
	  AddToTable(table,MakePPPrefix("prob", (s_parameter_name).c_str() ), MakePPEntry(it->second));
	}
	else {
	  AddToTable(table,MakePPPrefix((s_parameter_name).c_str() ), MakePPEntry(it->second));
	}
      }
    }

    // AMR controls
    DOMElement* amrc = GetChildByName_(structured_controls,"str_amr_controls", found);
    if (found) {
      map<string,map<string,vector<string> > > rc_data; // rc_data[critName][parm] = val

      DOMNodeList* children = amrc->getChildNodes();
      for (int i=0; i<children->getLength(); ++i) {
	DOMNode* inode = children->item(i);
	if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
	char* cname = mm.transcode(inode->getNodeName());
	if (strcmp(cname, "refinement_indicator") == 0) {
	  DOMElement* ielt = static_cast<DOMElement*>(inode);
	  string rc_name = GetAttributeValueS_(ielt, "name");
	  DOMNodeList* rc_children = ielt->getChildNodes();

	  rc_data[rc_name]["max_refinement_level"].push_back("-1");
	  for (int j=0; j<rc_children->getLength(); ++j) {
	    DOMNode* jnode = rc_children->item(j);
	    if (jnode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
	    char* rc_ename = mm.transcode(jnode->getNodeName());

	    if (strcmp(rc_ename, "max_refinement_level") == 0)
	    {
	      rc_data[rc_name][rc_ename][0] = mm.transcode(jnode->getTextContent());
	    }
	    else if (strcmp(rc_ename, "field_name") == 0)
	    {
	      rc_data[rc_name][rc_ename].push_back(mm.transcode(jnode->getTextContent()));
	    }
	    else if (strcmp(rc_ename, "regions") == 0)
	    {
	      string region_list = mm.transcode(jnode->getTextContent());
	      vector<string> tokens = BoxLib::Tokenize(region_list,",");
	      if (tokens.size()<1) {
		Errors::Message msg;
		msg << "An error occurred during parsing\n";
		msg << "numerical_controls->structured_controls->str_amr_controls\n";
		msg << "Error in regions list \n";
		msg << "Please correct and try again.\n" ;
		Exceptions::amanzi_throw(msg);
	      }
	      for (int k=0; k<tokens.size(); ++k) {
		rc_data[rc_name][rc_ename].push_back(tokens[k]);
	      }
	    }
	    else if ( strcmp(rc_ename, "start_time") == 0
		 || strcmp(rc_ename, "end_time") == 0)
	    {
	      string thisTime = mm.transcode(jnode->getTextContent());
	      map<string, string>::const_iterator iter = labeled_times_.find(thisTime);
	      if (iter != labeled_times_.end())
		rc_data[rc_name][rc_ename].push_back(iter->second);
	      else
		rc_data[rc_name][rc_ename].push_back(thisTime);
	    }
	    else {
	      if (strcmp(rc_ename, "value_greater") == 0
		|| strcmp(rc_ename, "value_less") == 0
		|| strcmp(rc_ename, "adjacent_difference_greater") == 0
		  || strcmp(rc_ename, "inside_region") == 0)
	      {
		rc_data[rc_name][rc_ename].push_back(mm.transcode(jnode->getTextContent()));
	      }
	      else {
		Errors::Message msg;
		msg << "An error occurred during parsing\n";
		msg << "numerical_controls->structured_controls->str_amr_controls\n";
		msg << "Unknown option \"" << rc_ename << "\" \n";
		msg << "Please correct and try again.\n" ;
		Exceptions::amanzi_throw(msg);
	      }
	    }
	  }
	}
	else {
	  for (map<string,string>::iterator it=amr_controls.begin(); it!=amr_controls.end(); ++it) {
	    if (it->first == cname) {
	      it->second = mm.transcode(inode->getTextContent());
	    }
	  }
	}
      }

      if (rc_data.size() > 0) {
	vector<string> rc_names;;
	for (map<string,map<string,vector<string> > >::const_iterator it=rc_data.begin(); it!=rc_data.end(); ++it) {
	  rc_names.push_back(it->first);
	  const map<string,vector<string> >& rc_map = it->second;
	  for (map<string,vector<string> >::const_iterator it1=rc_map.begin(); it1!=rc_map.end(); ++it1) {
	      AddToTable(table,
			 MakePPPrefix("amr", (it->first).c_str(), (it1->first).c_str() ),
			 MakePPEntry(it1->second));
	    }
	}
	AddToTable(table, MakePPPrefix("amr", "refinement_indicators"), MakePPEntry(rc_names));
      }
    }
  }

  // Figure out max_level
  int max_level = atoi(amr_controls["amr_levels"].c_str()) - 1;
  std::stringstream is;
  is << max_level;
  AddToTable(table, MakePPPrefix("amr", "max_level"), MakePPEntry(is.str()));

  string strToAdd, valToAdd;
  for (map<string,string>::iterator it=amr_controls.begin(); it!=amr_controls.end(); ++it) {
    if (it->first != "amr_levels") {

      if (it->first == "refinement_ratio"
	  || it->first == "regrid_interval"
	  || it->first == "blocking_factor"
	  || it->first == "number_error_buffer_cells"
	  || it->first == "max_grid_size")
      {
	vector<string> tokens = BoxLib::Tokenize(it->second,", ");
	AddToTable(table, MakePPPrefix("amr", it->first), MakePPEntry(tokens));
      }
      else {
	AddToTable(table, MakePPPrefix("amr", it->first), MakePPEntry(it->second));
      }
    }
  }

  ParmParse::appendTable(table);
}

void InputConverterS::ParseMesh_()
{
  list<ParmParse::PP_entry> table;
  bool found;

  DOMElement* dimension = static_cast<DOMElement*>(GetUniqueElementByTagsString_("mesh, dimension", found));
  {
    string dim = XMLString::transcode(dimension->getTextContent());
    dim_ = atoi(dim.c_str());
    if (dim_ != BL_SPACEDIM) {
      Errors::Message msg;
      msg << "An error occurred during parsing mesh->dimension\n";
      msg << "MUST BE \"" << BL_SPACEDIM << ", consistent with compilation flags.\"\n";
      msg << "Please correct and try again.\n" ;
      Exceptions::amanzi_throw(msg);
    }
  }

  DOMElement* generate = static_cast<DOMElement*>(GetUniqueElementByTagsString_("mesh, generate", found));
  if (found)
  {
    bool found;
    DOMElement* number_of_cells = GetChildByName_(generate, "number_of_cells", found, true);
    nx_ = GetAttributeValueL_(number_of_cells, "nx");
    ny_ = GetAttributeValueL_(number_of_cells, "ny");
    vector<int> n(dim_);
    n[0] = nx_;
    n[1] = ny_;
    if (dim_ == 3)
    {
      nz_ = GetAttributeValueL_(number_of_cells, "nz");
      n[2] = nz_;
    }
    AddToTable(table, MakePPPrefix("amr", "n_cell"), MakePPEntry(n));

    // Stash min/max coordinates for our own porpoises.
    DOMElement* box = GetChildByName_(generate, "box", found, true);
    lo_coords_ = GetAttributeVectorD_(box, "low_coordinates", dim_, "m", true);
    if (lo_coords_.size() != dim_)
      ThrowErrorIllformed_("mesh->generate->box", "coordinate array", "low_coordinates");
    hi_coords_ = GetAttributeVectorD_(box, "high_coordinates", dim_, "m", true);
    if (hi_coords_.size() != dim_)
      ThrowErrorIllformed_("mesh->generate->box", "coordinate array", "high_coordinates");
    AddToTable(table, MakePPPrefix("geometry", "prob_lo"), MakePPEntry(lo_coords_));
    AddToTable(table, MakePPPrefix("geometry", "prob_hi"), MakePPEntry(hi_coords_));

    // Periodic boundaries are not supported.
    vector<int> is_periodic(dim_, 0);
    AddToTable(table, MakePPPrefix("geometry", "is_periodic"), MakePPEntry(is_periodic));

    // Coordinate system is 0, which is probably "Cartesian." -- MSD: Actually 0 == cartesian, so yes
    AddToTable(table, MakePPPrefix("geometry", "coord_sys"), MakePPEntry(0));
  }
  else
    ThrowErrorMisschild_("mesh", "generate", "mesh");

  ParmParse::appendTable(table);
}

static void MakeBox(list<ParmParse::PP_entry>& table,
		    const string&              name,
		    const string&              purpose,
		    const vector<double>&      lo,
		    const vector<double>&      hi)
{
  AddToTable(table, MakePPPrefix("geometry", name, "lo_coordinate"), MakePPEntry(lo));
  AddToTable(table, MakePPPrefix("geometry", name, "hi_coordinate"), MakePPEntry(hi));
  AddToTable(table, MakePPPrefix("geometry", name, "type"),          MakePPEntry("box"));
  AddToTable(table, MakePPPrefix("geometry", name, "purpose"),       MakePPEntry(purpose));
}
  
#include <PMAMR_Labels.H>
void InputConverterS::ParseRegions_()
{
  list<ParmParse::PP_entry> table;
  vector<string> region_names;
  
  // Determine the geometry tolerance geometry_eps. dim_, {x,y,z}{min/max}_
  // should all be available because they are computed in ParseMesh(), 
  // which should precede this method.
  double max_size = -FLT_MIN;
  for (int d = 0; d < dim_; ++d)
    max_size = max(max_size, hi_coords_[d] - lo_coords_[d]);
  double geometry_eps = 1e-6 * max_size; // FIXME: This factor is fixed.
  AddToTable(table, MakePPPrefix("geometry", "geometry_eps"), MakePPEntry(geometry_eps));

  // Create default regions
  string name;
  vector<double> lo(dim_), hi(dim_);
  name = "XLOBC"; lo=lo_coords_; hi=hi_coords_; hi[0]=lo[0]; MakeBox(table, name, "xlobc", lo, hi); region_names.push_back(name);
  name = "XHIBC"; lo=lo_coords_; hi=hi_coords_; lo[0]=hi[0]; MakeBox(table, name, "xhibc", lo, hi); region_names.push_back(name);
  name = "YLOBC"; lo=lo_coords_; hi=hi_coords_; hi[1]=lo[1]; MakeBox(table, name, "ylobc", lo, hi); region_names.push_back(name);
  name = "YHIBC"; lo=lo_coords_; hi=hi_coords_; lo[1]=hi[1]; MakeBox(table, name, "yhibc", lo, hi); region_names.push_back(name);    
  if (dim_ > 2) {
    name = "ZLOBC"; lo=lo_coords_; hi=hi_coords_; hi[2]=lo[2]; MakeBox(table, name, "zlobc", lo, hi); region_names.push_back(name);
    name = "ZHIBC"; lo=lo_coords_; hi=hi_coords_; lo[2]=hi[2]; MakeBox(table, name, "zhibc", lo, hi); region_names.push_back(name);    
  }
  // Leave lo,hi to define domain "All" below
  name = "All";   lo=lo_coords_; hi=hi_coords_;              MakeBox(table, name, "all", lo, hi); region_names.push_back(name);

  bool found;
  DOMNode* regions = GetUniqueElementByTagsString_("regions", found);
  if (found)
  {
    bool found;

    // FIXME: U supports an alternate form for region spec, for some mysterious reason.  Should we duplicate here?
    // box
    vector<DOMNode*> boxes = GetChildren_(regions, "box", found);
    for (size_t i = 0; i < boxes.size(); ++i)
    {
      DOMElement* box = static_cast<DOMElement*>(boxes[i]);
      string region_name = GetAttributeValueS_(box, "name");
      region_names.push_back(region_name);
      vector<double> lo_coords = GetAttributeVectorD_(box, "low_coordinates", dim_, "m", true);
      vector<double> hi_coords = GetAttributeVectorD_(box, "high_coordinates", dim_, "m", true);
      AddToTable(table, MakePPPrefix("geometry", region_name, "lo_coordinate"), MakePPEntry(lo_coords));
      AddToTable(table, MakePPPrefix("geometry", region_name, "hi_coordinate"), MakePPEntry(hi_coords));

      // Is this region a surface?
      string type = "box";
      string purpose = "all";
      for (int d = 0; d < dim_; ++d)
      {
	// Has no extent in direction d
	if (std::abs(hi_coords[d] - lo_coords[d]) < geometry_eps)
	{
	  bool is_plane = true;
	  for (int d1=0; d1<dim_; ++d1) {
	    if (d!=d1) {
	      // Has finite extent perpendicular to d
	      is_plane &= std::abs(hi_coords[d1] - lo_coords[d1]) > geometry_eps;
	    }
	  }

	  if (!is_plane) {
	    BoxLib::Abort("No support for box regions with zero extent in more than one dimension");
	  }
	  
	  type = "surface";

	  // Is this on the domain boundary?
	  if (lo_coords[d] == lo[d]) {
	    purpose = PMAMR::RpurposeDEF[d];
	  }
	  else if (hi_coords[d] == hi[d]) {
	    purpose = PMAMR::RpurposeDEF[d+3];
	  }
	}
      }
      AddToTable(table, MakePPPrefix("geometry", region_name, "type"), 
                                     MakePPEntry(type));

      AddToTable(table, MakePPPrefix("geometry", region_name, "purpose"), 
                                     MakePPEntry(purpose));
    }

    // FIXME: color functions (what files do we read from?)
    vector<DOMNode*> colors = GetChildren_(regions, "color", found);
    for (size_t i = 0; i < colors.size(); ++i)
    {
    }

    // point
    vector<DOMNode*> points = GetChildren_(regions, "point", found);
    for (size_t i = 0; i < points.size(); ++i)
    {
      DOMElement* point = static_cast<DOMElement*>(points[i]);
      string region_name = GetAttributeValueS_(point, "name");
      region_names.push_back(region_name);
      vector<double> coords = GetAttributeVectorD_(point, "coordinate", dim_, "m", true);
      AddToTable(table, MakePPPrefix("geometry", region_name, "coordinate"), MakePPEntry(coords));
      AddToTable(table, MakePPPrefix("geometry", region_name, "type"), MakePPEntry("point"));
      AddToTable(table, MakePPPrefix("geometry", region_name, "purpose"), MakePPEntry("all"));
    }

    // plane
    vector<DOMNode*> planes = GetChildren_(regions, "plane", found);
    for (size_t i = 0; i < planes.size(); ++i)
    {
      DOMElement* plane = static_cast<DOMElement*>(planes[i]);
      string region_name = GetAttributeValueS_(plane, "name");
      region_names.push_back(region_name);
      vector<double> location = GetAttributeVectorD_(plane, "location", dim_, "m", true);
      vector<double> normal = GetAttributeVectorD_(plane, "normal", dim_, "", true);

      // FIXME: We need to redo the orientation logic.
      vector<double> lo_coords;
      vector<double> hi_coords;
      
      AddToTable(table, MakePPPrefix("geometry", region_name, "lo_coordinate"), MakePPEntry(lo_coords));
      AddToTable(table, MakePPPrefix("geometry", region_name, "hi_coordinate"), MakePPEntry(hi_coords));
      AddToTable(table, MakePPPrefix("geometry", region_name, "type"), MakePPEntry("surface"));
      AddToTable(table, MakePPPrefix("geometry", region_name, "purpose"), MakePPEntry("all"));

      // Optional tolerance.
      string tolerance = GetAttributeValueS_(plane, "tolerance", TYPE_NUMERICAL, false);
      if (tolerance != "")
        AddToTable(table, MakePPPrefix("geometry", region_name, "tolerance"), MakePPEntry("all"));

    }

    // region - this appears to be valid only for AmanziU
    vector<DOMNode*> my_regions = GetChildren_(regions, "region", found);
    if (my_regions.size()>0) {
      Errors::Message msg;
      msg << "An error occurred during parsing regions\n";
      msg << "AmanziS does not support the region-type regions\n";
      msg << "Please correct and try again.\n" ;
      Exceptions::amanzi_throw(msg);
    }

    // Regions only available in 2D
    if (dim_ == 2)
    {
      // polygon 
      vector<DOMNode*> polygons = GetChildren_(regions, "polygon", found);
      for (size_t i = 0; i < polygons.size(); ++i)
      {
        bool found;
        DOMElement* polygon = static_cast<DOMElement*>(polygons[i]);
        string region_name = GetAttributeValueS_(polygon, "name");
        region_names.push_back(region_name);
        int num_points = GetAttributeValueL_(polygon, "num_points");
        vector<DOMNode*> points = GetChildren_(polygon, "point", found);

	if (!found || points.size() != num_points) {
	  std::cout << "points.size() != num_points " << points.size() << " " << num_points << std::endl;
	  ThrowErrorIllformed_("regions", "point", "polygon");
	}
        vector<double> v1, v2, v3;
        for (size_t j = 0; j < points.size(); ++j)
        {
          DOMElement* point = static_cast<DOMElement*>(points[j]);
          string coord_string = XMLString::transcode(point->getTextContent());
          vector<double> coords = MakeCoordinates_(coord_string);
          v1.push_back(coords[0]);
          v2.push_back(coords[1]);
          if (dim_ == 3)
            v3.push_back(coords[2]);
        }
        AddToTable(table, MakePPPrefix("geometry", region_name, "v1"), MakePPEntry(v1));
        AddToTable(table, MakePPPrefix("geometry", region_name, "v2"), MakePPEntry(v2));
        if (!v3.empty())
          AddToTable(table, MakePPPrefix("geometry", region_name, "v3"), MakePPEntry(v3));

        AddToTable(table, MakePPPrefix("geometry", region_name, "type"), MakePPEntry("polygon"));
        AddToTable(table, MakePPPrefix("geometry", region_name, "purpose"), MakePPEntry("all"));
      }


      // ellipse 
      vector<DOMNode*> ellipses = GetChildren_(regions, "ellipse", found);
      for (size_t i = 0; i < ellipses.size(); ++i)
      {
        DOMElement* ellipse = static_cast<DOMElement*>(ellipses[i]);
        string region_name = GetAttributeValueS_(ellipse, "name");
        region_names.push_back(region_name);

        string center_string = GetAttributeValueS_(ellipse, "center");
        vector<double> center = MakeCoordinates_(center_string);
        AddToTable(table, MakePPPrefix("geometry", region_name, "center"), MakePPEntry(center));

        string radius_string = GetAttributeValueS_(ellipse, "radius");
        vector<double> radius = MakeCoordinates_(radius_string);

        AddToTable(table, MakePPPrefix("geometry", region_name, "radius"), MakePPEntry(radius));
        AddToTable(table, MakePPPrefix("geometry", region_name, "type"), MakePPEntry("ellipse"));
        AddToTable(table, MakePPPrefix("geometry", region_name, "purpose"), MakePPEntry("all"));
      }
    }

    // 3D-only region types.
    else if (dim_ == 3)
    {
      // rotated_polygon 
      // FIXME
      vector<DOMNode*> rotated_polygons = GetChildren_(regions, "rotated_polygon", found);
      for (size_t i = 0; i < rotated_polygons.size(); ++i)
      {
      }

      // swept_polygon
      // FIXME
      vector<DOMNode*> swept_polygons = GetChildren_(regions, "swept_polygon", found);
      for (size_t i = 0; i < swept_polygons.size(); ++i)
      {
      }
    }

    // logical
    vector<DOMNode*> logicals = GetChildren_(regions, "logical", found);
    for (size_t i = 0; i < logicals.size(); ++i)
    {
      DOMElement* logical = static_cast<DOMElement*>(logicals[i]);
      string region_name = GetAttributeValueS_(logical, "name");
      string operation   = GetAttributeValueS_(logical, "operation");
      string region_list = GetAttributeValueS_(logical, "region_list");
      region_names.push_back(region_name);

      AddToTable(table,   MakePPPrefix("geometry", region_name, "type"),      MakePPEntry("logical"));
      AddToTable(table,   MakePPPrefix("geometry", region_name, "operation"), MakePPEntry(operation));
      AddToTable(table,   MakePPPrefix("geometry", region_name, "purpose"),   MakePPEntry("all"));
      if (operation == "complement") {
	AddToTable(table, MakePPPrefix("geometry", region_name, "region"),    MakePPEntry(region_list));
      }
      else {
	vector<string> tokens = BoxLib::Tokenize(region_list,",");
	AddToTable(table, MakePPPrefix("geometry", region_name, "regions"),   MakePPEntry(tokens));
      }
    }
  }

  // Record region names.
  if (!region_names.empty())
    AddToTable(table, MakePPPrefix("geometry", "regions"), MakePPEntry(region_names));

  ParmParse::appendTable(table);
}

void InputConverterS::ParseGeochemistry_()
{
  list<ParmParse::PP_entry> table;

  MemoryManager mm;
  DOMNode* node;
  DOMElement* element;

  // chemical engine
  bool flag;
  node = GetUniqueElementByTagsString_("process_kernels, chemistry", flag);
  string state = GetAttributeValueS_(static_cast<DOMElement*>(node), "state");
  if (!(state == "on" ^ state == "off")) {
    Errors::Message msg;
    msg << "process_kernal->chemistry->state must be \"on\" or \"off\".\n";
    msg << "Please correct and try again.\n";
    Exceptions::amanzi_throw(msg);
  }

  string model("Amanzi");

  if (state == "on") {
  
    string engine = GetAttributeValueS_(static_cast<DOMElement*>(node), "engine");

    // process engine
    if (engine ==  "amanzi") {

      model = "Amanzi";

      std::string bgdfilename, format("simple");
      node = GetUniqueElementByTagsString_("geochemistry, amanzi_chemistry, reaction_network", flag);
      if (flag) {
	element = static_cast<DOMElement*>(node);
	bgdfilename = GetAttributeValueS_(element, "file");
	format = GetAttributeValueS_(element, "format", TYPE_NONE, false, format);
      } else {
        Exceptions::amanzi_throw(Errors::Message("BGD file is no longer created. Aborting."));
	// bgdfilename = CreateBGDFile_(xmlfilename_, rank_, status);
      }

      AddToTable(table, MakePPPrefix("Chemistry", "Thermodynamic_Database_File"), MakePPEntry(bgdfilename));
      AddToTable(table, MakePPPrefix("Chemistry", "Thermodynamic_Database_Format"), MakePPEntry(format));

    } else {
      bool valid_engine(true);
      std::string file_location;

      model = "Alquimia";

      if (engine == "pflotran") {
	AddToTable(table, MakePPPrefix("Chemistry", "Engine"), MakePPEntry("PFloTran"));
	file_location = "process_kernels, chemistry";
      } else if (engine == "crunchflow") {
	AddToTable(table, MakePPPrefix("Chemistry", "Engine"), MakePPEntry("CrunchFlow"));
	file_location = "process_kernels, chemistry";
      } else if (engine == "pflotran+") {
	AddToTable(table, MakePPPrefix("Chemistry", "Engine"), MakePPEntry("PFloTran+"));
	file_location = "process_kernels, chemistry";
      } else if (engine == "crunchflow+") {
	AddToTable(table, MakePPPrefix("Chemistry", "Engine"), MakePPEntry("CrunchFlow+"));
	file_location = "process_kernels, chemistry";
 
      } else {
	valid_engine = false;
      }

      // Pass along chemistry engine info.
      if (valid_engine) {

	// Find the name of the engine-specific input file.
	node = GetUniqueElementByTagsString_(file_location, flag);

	if (flag) {
	  element = static_cast<DOMElement*>(node);
	  std::string inpfilename = 
	    (element->hasAttribute(mm.transcode("input_filename"))
	     ? GetAttributeValueS_(element, "input_filename")
	     : CreateINFile_(xmlfilename_, rank_));
	  AddToTable(table, MakePPPrefix("Chemistry", "Engine_Input_File"), MakePPEntry(inpfilename));
	} else {
	  Errors::Message msg;
	  msg << "Unique tag string \"" << file_location << "\" must exists.\n";
	  Exceptions::amanzi_throw(msg);
	}
      }

      bool confound;
      DOMNode* connode = GetUniqueElementByTagsString_("geochemistry, constraints", confound);
      if (confound)
      {
	bool gccsfound = false;
	vector<DOMNode*> gccs = GetChildren_(connode, "constraint", gccsfound);
	if (gccsfound)
	{
	  for (size_t g = 0; g < gccs.size(); ++g)
	  {
	    DOMElement* constraint = static_cast<DOMElement*>(gccs[g]);
	    string constraint_name = GetAttributeValueS_(constraint, "name");
	    if (!constraint_name.empty()) {
	      constraint_names_.push_back(constraint_name);
	    }
	  }
	}
      }
    }
  }
  else {
    model = "Off";
  }

  AddToTable(table, MakePPPrefix("prob", "chemistry_model"), MakePPEntry(model));

  ParmParse::appendTable(table);
}

void InputConverterS::ParseMaterials_(bool& do_tracer_diffusion)
{
  list<ParmParse::PP_entry> table;
  bool found;
  vector<string> material_names;

  DOMNode* materials = GetUniqueElementByTagsString_("materials", found);
  if (found)
  {
    bool found;
    vector<DOMNode*> mats = GetChildren_(materials, "material", found);
    for (size_t i = 0; i < mats.size(); ++i)
    {
      DOMElement* mat = static_cast<DOMElement*>(mats[i]);
      string mat_name = GetAttributeValueS_(mat, "name");
      material_names.push_back(mat_name);
      bool found;

      // Mechanical properties.
      DOMElement* mech_prop = GetChildByName_(mat, "mechanical_properties", found);
      if (found)
      {
        ParseMechProperty_(mech_prop, mat_name, "porosity", table, true);
        ParseMechProperty_(mech_prop, mat_name, "particle_density", table, false); // FIXME: Should be true for required!
        ParseMechProperty_(mech_prop, mat_name, "specific_storage", table, false); 
        ParseMechProperty_(mech_prop, mat_name, "specific_yield", table, false); 
        if (ParseMechProperty_(mech_prop, mat_name, "dispersion_tensor", table, false)) {
	  do_tracer_diffusion = true;
	}
        ParseMechProperty_(mech_prop, mat_name, "tortuosity", table, false); 
      }

      // Assigned regions.
      vector<string> assigned_regions = GetChildVectorS_(mat, "assigned_regions", found, true);
      AddToTable(table, MakePPPrefix("rock", mat_name, "regions"), 
                                     MakePPEntry(assigned_regions));

      // Permeability OR hydraulic conductivity.
      bool k_found, K_found;
      DOMElement* permeability = GetChildByName_(mat, "permeability", k_found, false);
      DOMElement* conductivity = GetChildByName_(mat, "hydraulic_conductivity", K_found, false);
      if (!k_found && !K_found)
      {
        Errors::Message msg;
        msg << "Neither permeability nor hydraulic_conductivity was found for material \"" << mat_name << "\".\n";
        msg << "Please correct and try again.\n";
        Exceptions::amanzi_throw(msg);
      }
      else if (k_found && K_found)
      {
        Errors::Message msg;
        msg << "Both permeability AND hydraulic_conductivity were found for material \"" << mat_name << "\".\n";
        msg << "Only one of these is allowed. Please correct and try again.\n";
        Exceptions::amanzi_throw(msg);
      }
      else if (k_found)
      {
        string x = GetAttributeValueS_(permeability, "x", TYPE_NUMERICAL, false);
        if (!x.empty())
        {
	  AddToTable(table, MakePPPrefix("rock", mat_name, "permeability", "horizontal", "vals"),
		     MakePPEntry(x));

          string y = GetAttributeValueS_(permeability, "y");
	  if (dim_ < 3) {
	    AddToTable(table, MakePPPrefix("rock", mat_name, "permeability", "vertical", "vals"),
		       MakePPEntry(y));
	  } else {
	    AddToTable(table, MakePPPrefix("rock", mat_name, "permeability", "horizontal", "vals"),
		       MakePPEntry(x));
	    AddToTable(table, MakePPPrefix("rock", mat_name, "permeability", "horizontal1", "vals"),
		       MakePPEntry(y));

	    string z = GetAttributeValueS_(permeability, "z");
	    AddToTable(table, MakePPPrefix("rock", mat_name, "permeability", "vertical", "vals"),
		       MakePPEntry(z));
	  }
        }
        else
        {
          string type = GetAttributeValueS_(permeability, "type");
          AddToTable(table, MakePPPrefix("rock", mat_name, "permeability", "type"),
                     MakePPEntry(type));
          if (type == "file")
          {
            string filename = GetAttributeValueS_(permeability, "filename");
            string attribute = GetAttributeValueS_(permeability, "attribute");
            AddToTable(table, MakePPPrefix("rock", mat_name, "permeability", "filename"),
                       MakePPEntry(filename));
            AddToTable(table, MakePPPrefix("rock", mat_name, "permeability", "attribute"),
                                           MakePPEntry(attribute));
          }
          else if (type == "gslib")
          {
            string parameter_file = GetAttributeValueS_(permeability, "parameter_file");
            string value = GetAttributeValueS_(permeability, "value");
            string data_file = GetAttributeValueS_(permeability, "data_file");
            AddToTable(table, MakePPPrefix("rock", mat_name, "permeability", "parameter_file"),
                       MakePPEntry(parameter_file));
            AddToTable(table, MakePPPrefix("rock", mat_name, "permeability", "value"),
                       MakePPEntry(value));
            AddToTable(table, MakePPPrefix("rock", mat_name, "permeability", "data_file"),
                       MakePPEntry(data_file));
          }
          else
            ThrowErrorIllformed_("materials", "type", "permeability");
        }
      }
      else
      {
        string x = GetAttributeValueS_(conductivity, "x", TYPE_NUMERICAL, false);
        if (!x.empty())
        {
          AddToTable(table, MakePPPrefix("rock", mat_name, "hydraulic_conductivity", "horizontal", "vals"),
                     MakePPEntry(x));
          string y = GetAttributeValueS_(conductivity, "y");
	  if (dim_ < 3) {
	    AddToTable(table, MakePPPrefix("rock", mat_name, "hydraulic_conductivity", "vertical", "vals"),
		       MakePPEntry(y));
	  } else {
	    AddToTable(table, MakePPPrefix("rock", mat_name, "hydraulic_conductivity", "horizontal1", "vals"),
		       MakePPEntry(y));
	    string z = GetAttributeValueS_(conductivity, "z");
	    AddToTable(table, MakePPPrefix("rock", mat_name, "hydraulic_conductivity", "vertical", "vals"),
		       MakePPEntry(z));
	  }
        }
        else
        {
          string type = GetAttributeValueS_(permeability, "type");
          AddToTable(table, MakePPPrefix("rock", mat_name, "hydraulic_conductivity", "type"),
                     MakePPEntry(type));
          if (type == "gslib")
          {
            string parameter_file = GetAttributeValueS_(conductivity, "parameter_file");
            string value = GetAttributeValueS_(conductivity, "value");
            string data_file = GetAttributeValueS_(conductivity, "data_file");
            AddToTable(table, MakePPPrefix("rock", mat_name, "hydraulic_conductivity", "parameter_file"),
                       MakePPEntry(parameter_file));
            AddToTable(table, MakePPPrefix("rock", mat_name, "hydraulic_conductivity", "value"),
                       MakePPEntry(value));
            AddToTable(table, MakePPPrefix("rock", mat_name, "hydraulic_conductivity", "data_file"),
                       MakePPEntry(data_file));
          }
          else
            ThrowErrorIllformed_("materials", "type", "hydraulic_conductivity");
        }
      }

      // Capillary pressure model.
      DOMElement* cap_pressure = GetChildByName_(mat, "cap_pressure", found, false);
      if (found)
      {
        bool found;
        string model = GetAttributeValueS_(cap_pressure, "model");
        if (model == "van_genuchten")
          AddToTable(table, MakePPPrefix("rock", mat_name, "cpl", "type"), MakePPEntry("VanGenuchten"));
        else if (model == "brooks_corey")
          AddToTable(table, MakePPPrefix("rock", mat_name, "cpl", "type"), MakePPEntry("BrooksCorey"));
        else if (model != "none")
          ThrowErrorIllformed_("materials", "type", "cap_pressure");

        if ((model == "van_genuchten") || (model == "brooks_corey"))
        {
          DOMElement* parameters = GetChildByName_(cap_pressure, "parameters", found, true);
          string alpha = GetAttributeValueS_(parameters, "alpha");
          string sr = GetAttributeValueS_(parameters, "sr");
          AddToTable(table, MakePPPrefix("rock", mat_name, "cpl", "alpha"), MakePPEntry(alpha));
          AddToTable(table, MakePPPrefix("rock", mat_name, "cpl", "Sr"), MakePPEntry(sr));
	  string m, lambda;
	  if (model == "van_genuchten") {
	    m = GetAttributeValueS_(parameters, "m");
	    AddToTable(table, MakePPPrefix("rock", mat_name, "cpl", "m"), MakePPEntry(m));
	  } else {
	    lambda = GetAttributeValueS_(parameters, "lambda");
	    AddToTable(table, MakePPPrefix("rock", mat_name, "cpl", "lambda"), MakePPEntry(lambda));
	  }
          string optional_krel_smoothing_interval = GetAttributeValueS_(parameters, "optional_krel_smoothing_interval", TYPE_NUMERICAL, false);
          if (!optional_krel_smoothing_interval.empty())
          {
            AddToTable(table, MakePPPrefix("rock", mat_name, "cpl", "Kr_smoothing_max_pcap"),
                       MakePPEntry(optional_krel_smoothing_interval));
          }
        }

        // FIXME: Something about a WRM plot file??
      }

      // Relative permeability.
      DOMElement* rel_perm = GetChildByName_(mat, "rel_perm", found, false);
      if (found)
      {
        bool found;
        string model = GetAttributeValueS_(rel_perm, "model");
        if (model == "mualem")
          AddToTable(table, MakePPPrefix("rock", mat_name, "Kr_model"), MakePPEntry("Mualem"));
        else if (model == "burdine")
          AddToTable(table, MakePPPrefix("rock", mat_name, "Kr_model"), MakePPEntry("Burdine"));
        else if (model != "none")
          ThrowErrorIllformed_("materials", "type", "rel_perm");

        if (model == "mualem")
        {
          // We stick in a default "ell" value, since ell doesn't appear 
          // in the v2.x input spec.
          double Kr_ell = 0.5;
          AddToTable(table, MakePPPrefix("rock", mat_name, "Kr_ell"), MakePPEntry(Kr_ell));
        }
        else if (model == "burdine")
        {
          // We stick in a default "ell" value, since ell doesn't appear 
          // in the v2.x input spec.
          double Kr_ell = 2.0;
          AddToTable(table, MakePPPrefix("rock", mat_name, "Kr_ell"), MakePPEntry(Kr_ell));

          // There's also an "exp" parameter.
          string Kr_exp = GetChildValueS_(cap_pressure, "exp", found, true);
          AddToTable(table, MakePPPrefix("rock", mat_name, "Kr_exp"), MakePPEntry(Kr_exp));
        }

      }

      // Sorption isotherms.
      DOMElement* sorption_isotherms = GetChildByName_(mat, "sorption_isotherms", found, false);
      if (found)
      {
        // Look for solutes.
        bool found;
        vector<DOMNode*> solutes = GetChildren_(sorption_isotherms, "primary", found, true);
        for (size_t i = 0; i < solutes.size(); ++i)
        {
          DOMElement* solute = static_cast<DOMElement*>(solutes[i]);
          string solute_name = GetAttributeValueS_(solute, "name");

          bool found;
          DOMElement* kd_model = GetChildByName_(solute, "kd_model", found, false);
          if (found)
          {
            // Search for kd, b, or n.
            string kd = GetAttributeValueS_(kd_model, "kd", TYPE_NUMERICAL, false);
            if (!kd.empty())
            {
              AddToTable(table, MakePPPrefix("rock", mat_name, "sorption_isotherms", solute_name, "Kd"), 
                         MakePPEntry(kd));
            }
            else
            {
              string b = GetAttributeValueS_(kd_model, "b", TYPE_NUMERICAL, false);
              if (!b.empty())
              {
                AddToTable(table, MakePPPrefix("rock", mat_name, "sorption_isotherms", solute_name, "Langmuir b"), 
                           MakePPEntry(b));
              }
              else
              {
                string n = GetAttributeValueS_(kd_model, "n", TYPE_NUMERICAL, false);
                if (!n.empty())
                {
                  AddToTable(table, MakePPPrefix("rock", mat_name, "sorption_isotherms", solute_name, "Freundlich n"), 
                             MakePPEntry(n));
                }
                else
                  ThrowErrorIllformed_("materials->sorption_isotherms", "kd_model", solute_name);
              }
            }
          }
        }
      }
    }
    AddToTable(table, MakePPPrefix("rock", "rocks"), MakePPEntry(material_names));
  }

  ParmParse::appendTable(table);
}

  void InputConverterS::ParseProcessKernels_(string& flow_model, bool& do_tracer_diffusion)
{
  list<ParmParse::PP_entry> table;
  bool found;

  // Flow, transport, and chemistry must all be present.
  DOMElement* flow = static_cast<DOMElement*>(GetUniqueElementByTagsString_("process_kernels, flow", found));
  if (!found)
    ThrowErrorMisschild_("process_kernels", "flow", "process_kernels");
  DOMElement* transport = static_cast<DOMElement*>(GetUniqueElementByTagsString_("process_kernels, transport", found));
  if (!found)
    ThrowErrorMisschild_("process_kernels", "transport", "process_kernels");
  GetUniqueElementByTagsString_("process_kernels, chemistry", found);
  if (!found)
    ThrowErrorMisschild_("process_kernels", "chemistry", "process_kernels");

  // Flow model.
  string flow_state = GetAttributeValueS_(flow, "state");
  AddToTable(table, MakePPPrefix("prob", "flow_state"), MakePPEntry(flow_state));
  flow_model = "unused";
  if (flow_state == "on")
  {
    flow_model = GetAttributeValueS_(flow, "model");
    AddToTable(table, MakePPPrefix("prob", "flow_model"), MakePPEntry(flow_model));
  }

  // Transport model.
  string transport_state = GetAttributeValueS_(transport, "state");
  AddToTable(table, MakePPPrefix("prob", "do_tracer_advection"), MakePPEntry((transport_state == "on")));

  AddToTable(table, MakePPPrefix("prob", "do_tracer_diffusion"), MakePPEntry(do_tracer_diffusion));
  // FIXME: What else here?

  // Chemistry model. -- processed by ParseGeochemistry

  ParmParse::appendTable(table);
}

void InputConverterS::ParseSoluteNames_()
{
  bool found;
  DOMElement* liquid_phase = static_cast<DOMElement*>(GetUniqueElementByTagsString_("phases, liquid_phase", found));
  if (found)
  {
    bool found;
    DOMElement* dissolved_comps = GetChildByName_(liquid_phase, "dissolved_components", found);
    if (found)
    {
      DOMElement* solutes = GetChildByName_(dissolved_comps, "primaries", found, true);
      vector<DOMNode*> sols = GetChildren_(solutes, "primary", found);
      for (size_t i = 0; i < sols.size(); ++i)
      {
        MemoryManager mm;
        string sol_name = TrimString_(mm.transcode(sols[i]->getTextContent()));

        // Record the name for later.
        solutes_.push_back(sol_name);
      }
    }
  }
  solute_names_parsed_ = true;
}

void InputConverterS::ParsePhases_(bool& do_tracer_diffusion)
{
  list<ParmParse::PP_entry> table;

  vector<string> phase_names;
  bool found;
  DOMElement* liquid_phase = static_cast<DOMElement*>(GetUniqueElementByTagsString_("phases, liquid_phase", found));
  liquid_density_ = -1;
  liquid_viscosity_ = -1;
  if (found)
  {
    string name = GetAttributeValueS_(liquid_phase, "name");
    phase_names.push_back(name);
    bool found;
    string viscosity = GetChildValueS_(liquid_phase, "viscosity", found, true);
    AddToTable(table, MakePPPrefix("phase", name, "viscosity"), MakePPEntry(viscosity));
    liquid_viscosity_ = atof(viscosity.c_str());
    string density = GetChildValueS_(liquid_phase, "density", found, true);
    liquid_density_ = atof(density.c_str());
    AddToTable(table, MakePPPrefix("phase", name, "density"), MakePPEntry(density));
    // FIXME: Not supported by structured
    //string eos = GetChildValueS_(liquid_phase, "eos", found, false);
    //if (found)
    //  AddToTable(table, MakePPPrefix("phase", name, "eos"), MakePPEntry(eos)); // FIXME
    DOMElement* dissolved_comps = GetChildByName_(liquid_phase, "dissolved_components", found);
    if (found)
    {
      DOMElement* solutes = GetChildByName_(dissolved_comps, "primaries", found, true);
      vector<DOMNode*> sols = GetChildren_(solutes, "primary", found);
      for (size_t i = 0; i < sols.size(); ++i)
      {
        MemoryManager mm;
        string sol_name = TrimString_(mm.transcode(sols[i]->getTextContent()));

        // Record the coefficient of diffusion.
        DOMElement* solute = static_cast<DOMElement*>(sols[i]);
        string diff_coeff = GetAttributeValueS_(solute, "coefficient_of_diffusion", TYPE_NONE, false);
	if (diff_coeff != "") {
	  AddToTable(table, MakePPPrefix("tracer", sol_name, "molecularDiffusivity"), MakePPEntry(diff_coeff));
	  if (atof(diff_coeff.c_str()) != 0) {
	    do_tracer_diffusion = true;
	  }
	}
      }
    }
    // Assume we have a single component with the same 
    // name as the liquid phase.
    vector<string> components(1, name);
    AddToTable(table, MakePPPrefix("phase", name, "comps"), MakePPEntry(components));
    
    // Zero diffusivity by default.
    AddToTable(table, MakePPPrefix("phase", name, "diffusivity"), MakePPEntry(0.0));
  }

  AddToTable(table, MakePPPrefix("tracer", "tracers"), MakePPEntry(solutes_));

  GetUniqueElementByTagsString_("phases, solid_phase", found);
  if (found) {
    // Currently, Amanzi-S doesn't think of things this way.
  }

  GetUniqueElementByTagsString_("phases, gas_phase", found);
  AddToTable(table, MakePPPrefix("phase", "phases"), MakePPEntry(phase_names));

  ParmParse::appendTable(table);
}

void InputConverterS::ParseInitialConditions_()
{
  list<ParmParse::PP_entry> table;
  bool found;
  DOMNode* initial_conditions = GetUniqueElementByTagsString_("initial_conditions", found);
  if (found)
  {
    bool ifound;
    vector<DOMNode*> ics = GetChildren_(initial_conditions, "initial_condition", ifound);
    if (ifound)
    {
      vector<string> ic_names;
      for (size_t i = 0; i < ics.size(); ++i)
      {
        DOMElement* ic = static_cast<DOMElement*>(ics[i]);
        string ic_name = GetAttributeValueS_(ic, "name");
        ic_names.push_back(ic_name);

        // phase/comp
        bool pfound = false;
        DOMElement* lp = GetChildByName_(ic, "liquid_phase", pfound);
        if (pfound) {
	  bool sfound = false;
          bool cfound = false;
          DOMElement* lc = GetChildByName_(lp, "liquid_component", cfound);
          if (cfound) {
            vector<DOMNode*> nodes;
            string ic_type_labels[5] = {
              "uniform_pressure",
              "linear_pressure",
              "uniform_saturation",
              "linear_saturation",
              "velocity"};
            bool tfound = false;
            for (int i=0; i<5 && !tfound; ++i) {
              nodes = GetChildren_(lc, ic_type_labels[i], tfound, false);
              if (tfound) {
                vector<string> values, vel;
                vector<double> gvalues, rvalues, rpoints;
                for (size_t j = 0; j < nodes.size(); ++j) {
                  DOMElement* elt = static_cast<DOMElement*>(nodes[j]);
                  if (ic_type_labels[i]=="velocity") {
                    vel.push_back(GetAttributeValueS_(elt, "x"));
                    vel.push_back(GetAttributeValueS_(elt, "y"));
                    if (dim_ > 2) {
                      vel.push_back(GetAttributeValueS_(elt, "z"));
                    }
                  } else {
                    values.push_back(GetAttributeValueS_(elt, "value"));
                  }

                  // Get extra info, if required
                  if (ic_type_labels[i]=="linear_pressure"
                      || ic_type_labels[i]=="linear_saturation") {

                    gvalues = GetAttributeVectorD_(elt, "gradient", dim_, "Pa/m", true);
                    if (gvalues.size() != dim_)
                      ThrowErrorIllformed_("initial_conditions->linear_pressure", "coordinate array", "gradient");

                    rpoints = GetAttributeVectorD_(elt, "reference_coord", dim_, "m", true);
                    if (rpoints.size() != dim_)
                      ThrowErrorIllformed_("initial_conditions->linear_pressure", "coordinate array", "reference_coord");
                  }
                }

                AddToTable(table, MakePPPrefix("comp", "ics", ic_name, "type"),MakePPEntry(ic_type_labels[i]));

                if (values.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "ics", ic_name, "val"),MakePPEntry(values));
                if (vel.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "ics", ic_name, "vel"),MakePPEntry(vel));
                if (gvalues.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "ics", ic_name, "grad"),MakePPEntry(gvalues));
                if (rvalues.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "ics", ic_name, "vals"),MakePPEntry(rvalues));
                if (rpoints.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "ics", ic_name, "loc"),MakePPEntry(rpoints));
              }
            }	    
          }
          else {
            Errors::Message msg;
            msg << "\"liquid_component\" not present in \"liquid_phase\" initial_condition \""
              << ic_name << "\".\n";
            msg << "Please correct and try again.\n";
            Exceptions::amanzi_throw(msg);
          }

	  // solute ICs
	  bool scfound = false;
	  DOMElement* sc = GetChildByName_(lp, "solute_component", scfound);
	  if (scfound)
	  {
	    sfound = true;

	    vector<string> conc(solutes_.size());
	    for (int s=0; s<solutes_.size(); ++s) {
	      conc[s] = "0.0";
	    }
	    vector<DOMNode*> sols = GetChildren_(sc, "uniform_conc", found);

	    if (sols.size() > 0) {
	      for (size_t si=0; si<sols.size(); ++si)
	      {
		DOMElement* ss = static_cast<DOMElement*>(sols[si]);
		string this_spec = GetAttributeValueS_(ss, "name");
		string this_conc = GetAttributeValueS_(ss, "value");
		bool found = false;
		for (int s=0; s<solutes_.size() && !found; ++s) {
		  if (this_spec == solutes_[s]) {
		    conc[s] = this_conc;
		    found = true;
		  }
		}
		if (!found) {
		  Errors::Message msg;
		  msg << "\"uniform_conc\" section names nonexisting solute in initial_condition \"" << ic_name << "\".\n";
		  msg << "Please correct and try again.\n";
		  Exceptions::amanzi_throw(msg);	  
		}
	      }
	      for (int s=0; s<solutes_.size(); ++s) {
		AddToTable(table, MakePPPrefix("tracer", solutes_[s], "ic", ic_name, "val"), MakePPEntry(conc[s]));
		AddToTable(table, MakePPPrefix("tracer", solutes_[s], "ic", ic_name, "type"), MakePPEntry("concentration"));
	      }
	    }
	  }

	  // Sniff out geochemical conditions, if any.
	  bool gcfound = false;
	  DOMElement* gc = GetChildByName_(lp, "geochemistry_component", gcfound);
	  if (gcfound)
	  {
	    sfound = true;

	    bool gfound;
	    DOMElement* geocon = GetChildByName_(gc, "constraint", gfound);
	    if (gfound)
	    {
	      string constraint_name = GetAttributeValueS_(geocon, "name");
	      if (!constraint_name.empty()) {
		bool known_cond = false;
		for (int c=0; c<constraint_names_.size() && !known_cond; ++c) {
		  if (constraint_names_[c] == constraint_name) known_cond = true;
		}
		if (known_cond) {
		  for (int s = 0; s<solutes_.size(); ++s) {
		    AddToTable(table, MakePPPrefix("tracer", solutes_[s], "ic", ic_name, "geochemical_condition"), MakePPEntry(constraint_name));
		    AddToTable(table, MakePPPrefix("tracer", solutes_[s], "ic", ic_name, "type"), MakePPEntry("concentration"));
		  }
		}
		else {
		  Errors::Message msg;
		  msg << "Undeclared constraint name \""  << constraint_name << "\" in initial_condition \"" << ic_name << "\".\n";
		  msg << "Please correct and try again.\n";
		  Exceptions::amanzi_throw(msg);	  
		}
	      }
	      else {
		Errors::Message msg;
		msg << "No \"name\" attibute specified for geochemical constraint in initial_condition \"" << ic_name << "\".\n";
		msg << "Please correct and try again.\n";
		Exceptions::amanzi_throw(msg);	  
	      }
	    }
	  }

	  // Assigned regions.
	  bool rfound = false;
	  vector<string> assigned_regions = GetChildVectorS_(ic, "assigned_regions", rfound, true);
	  if (rfound) {
	    AddToTable(table, MakePPPrefix("comp", "ics", ic_name, "regions"), MakePPEntry(assigned_regions));
	    if (sfound) {
	      for (size_t i = 0; i < solutes_.size(); ++i) {
		AddToTable(table, MakePPPrefix("tracer", solutes_[i], "ic", ic_name, "regions"), MakePPEntry(assigned_regions));
	      }
	    }
	  } else {
	    Errors::Message msg;
	    msg << "\"assigned_regions\" not present in initial_condition \"" << ic_name << "\".\n";
	    msg << "Please correct and try again.\n";
	    Exceptions::amanzi_throw(msg);
	  }

        }
        else {
          Errors::Message msg;
          msg << "\"liquid_phase\" not present in initial condition \"" << ic_name << "\".\n";
          msg << "Please correct and try again.\n";
          Exceptions::amanzi_throw(msg);
        }
      }

      // List initial conditions where needed.
      AddToTable(table, MakePPPrefix("comp", "ic_labels"), MakePPEntry(ic_names));
      for (size_t i = 0; i < solutes_.size(); ++i)
        AddToTable(table, MakePPPrefix("tracer", solutes_[i], "tinits"), MakePPEntry(ic_names));
    }
  }
  ParmParse::appendTable(table);
}

void InputConverterS::ParseBoundaryConditions_()
{
  list<ParmParse::PP_entry> table;
  bool found;
  DOMNode* boundary_conditions = GetUniqueElementByTagsString_("boundary_conditions", found);
  if (found)
  {
    bool bfound;
    vector<DOMNode*> bcs = GetChildren_(boundary_conditions, "boundary_condition", bfound);
    if (bfound)
    {
      vector<string> bc_names;
      vector<string> tbc_names;
      for (size_t i = 0; i < bcs.size(); ++i)
      {
        DOMElement* bc = static_cast<DOMElement*>(bcs[i]);
        string bc_name = GetAttributeValueS_(bc, "name");
        bc_names.push_back(bc_name);

        // phase/comp
        bool pfound = false;
        DOMElement* lp = GetChildByName_(bc, "liquid_phase", pfound);
        if (pfound) {
	  bool sfound = false;
          bool cfound = false;
          DOMElement* lc = GetChildByName_(lp, "liquid_component", cfound);
          if (cfound) {
            vector<DOMNode*> nodes;
            string bc_type_labels[10] = {"inward_mass_flux",
              "outward_mass_flux",
              "inward_volumetric_flux",
              "outward_volumetric_flux",
              "uniform_pressure",
              "linear_pressure",
              "seepage_face",
              "hydrostatic",
              "linear_hydrostatic",
              "no_flow"};
            bool tfound = false;
            for (int i=0; i<10 && !tfound; ++i) {
              nodes = GetChildren_(lc, bc_type_labels[i], tfound, false);
              if (tfound) {
                vector<string> functions, starts, values, imf;
                vector<string> gvalues, rvalues, rpoints, csys, rwths, sms;
                for (size_t j = 0; j < nodes.size(); ++j) {
                  DOMElement* elt = static_cast<DOMElement*>(nodes[j]);

                  functions.push_back(GetAttributeValueS_(elt, "function"));

                  if (bc_type_labels[i]!="linear_pressure"
                      && bc_type_labels[i]!="seepage_face"
                      && bc_type_labels[i]!="linear_hydrostatic")
                  {
		    string this_value = GetAttributeValueS_(elt, "value");
		    if (bc_type_labels[i]=="inward_mass_flux"
			|| bc_type_labels[i]=="outward_mass_flux") {
		      map<string, string>::const_iterator iter = labeled_area_mass_fluxes_.find(this_value);
		      if (iter != labeled_area_mass_fluxes_.end())
			values.push_back(iter->second);
		      else
			values.push_back(this_value);
		    }
		    else {
		      values.push_back(this_value);
		    }

		    string this_start = GetAttributeValueS_(elt, "start");
		    map<string, string>::const_iterator iter = labeled_times_.find(this_start);
		    if (iter != labeled_times_.end())
		      starts.push_back(ConvertTimeToSeconds(iter->second));
		    else
		      starts.push_back(ConvertTimeToSeconds(this_start));
                  }

                  // Get extra info, if required
                  if (bc_type_labels[i]=="linear_pressure"
                      || bc_type_labels[i]=="linear_hydrostatic") {
                    gvalues.push_back(GetAttributeValueS_(elt, "gradient_value"));
                    rpoints.push_back(GetAttributeValueS_(elt, "reference_point"));
                  }

                  if (bc_type_labels[i]=="linear_pressure") {
                    rvalues.push_back(GetAttributeValueS_(elt, "reference_value"));
                  }

                  if (bc_type_labels[i]=="linear_hydrostatic") {
                    rwths.push_back(GetAttributeValueS_(elt, "reference_water_table_height"));
                    sms.push_back(GetAttributeValueS_(elt, "submodel"));
                  }

                  if (bc_type_labels[i]=="hydrostatic") {
		    csys.push_back(GetAttributeValueS_(elt, "coordinate_system", TYPE_NONE, false, "absolute"));
		    string sms_ = GetAttributeValueS_(elt, "submodel", TYPE_NONE, false);
		    if (!sms_.empty()) {
		      sms.push_back(sms_);
		    }
                  }

                  if (bc_type_labels[i]=="seepage_face") {
		    Errors::Message msg;
		    msg << "\"seepage_face\" boundary_condition not supported by structured\".\n";
		    Exceptions::amanzi_throw(msg);
                    //imf.push_back(GetAttributeValueS_(elt, "inward_mass_flux"));
                  }
                }

                AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "type"),MakePPEntry(bc_type_labels[i]));

                if (functions.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "forms"),MakePPEntry(functions));
                if (starts.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "times"),MakePPEntry(starts));
                if (values.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "vals"),MakePPEntry(values));
                if (imf.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "inward_mass_flux"),MakePPEntry(imf));
                if (gvalues.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "grad"),MakePPEntry(gvalues));
                if (rvalues.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "vals"),MakePPEntry(rvalues));
                if (rpoints.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "loc"),MakePPEntry(rpoints));
                if (csys.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "coord_sys"),MakePPEntry(csys));
                if (rwths.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "wt"),MakePPEntry(functions));
                if (sms.size()>0)
                  AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "submodel"),MakePPEntry(functions));
              }
            }	    
          }
          else {
            Errors::Message msg;
            msg << "\"liquid_component\" not present in \"liquid_phase\" boundary_condition \""
              << bc_name << "\".\n";
            msg << "Please correct and try again.\n";
            Exceptions::amanzi_throw(msg);
          }

	  // solute BCs
	  bool scfound = false;
	  DOMElement* sc = GetChildByName_(lp, "solute_component", scfound);
	  if (scfound)
	  {
	    sfound = true;

	    vector<vector<string> > conc(solutes_.size());
	    vector<vector<string> > ctime(solutes_.size());
	    vector<vector<string> > cfunc(solutes_.size());
	    vector<DOMNode*> sols = GetChildren_(sc, "aqueous_conc", found);

	    if (sols.size() > 0)
	    {
	      for (size_t si=0; si<sols.size(); ++si)
	      {
		DOMElement* ss = static_cast<DOMElement*>(sols[si]);
		string this_spec = GetAttributeValueS_(ss, "name");
		string this_conc = GetAttributeValueS_(ss, "value");
		string this_func = GetAttributeValueS_(ss, "function");
		string this_time = GetAttributeValueS_(ss, "start");
		bool found = false;
		for (int s=0; s<solutes_.size() && !found; ++s) {
		  if (this_spec == solutes_[s]) {

		    map<string, string>::const_iterator iterV = labeled_numbers_.find(this_conc);
		    if (iterV != labeled_numbers_.end())
		      conc[s].push_back(iterV->second);
		    else
		      conc[s].push_back(this_conc);

		    map<string, string>::const_iterator iterT = labeled_times_.find(this_time);
		    if (iterT != labeled_times_.end())
		      ctime[s].push_back(ConvertTimeToSeconds(iterT->second));
		    else
		      ctime[s].push_back(ConvertTimeToSeconds(this_time));

		    cfunc[s].push_back(this_func);
		    found = true;
		  }
		}
		if (!found) {
		  Errors::Message msg;
		  msg << "\"aqueous_conc\" section names nonexisting solute in boundary_condition \"" << bc_name << "\".\n";
		  msg << "Please correct and try again.\n";
		  Exceptions::amanzi_throw(msg);	  
		}
	      }
	      for (int s=0; s<solutes_.size(); ++s) {
		if (conc[s].size() == 0) {
		  conc[s].push_back("0.0");
		}
	      }
	      for (int s=0; s<solutes_.size(); ++s) {
		AddToTable(table, MakePPPrefix("tracer", solutes_[s], "bc", bc_name, "vals"), MakePPEntry(conc[s]));
		AddToTable(table, MakePPPrefix("tracer", solutes_[s], "bc", bc_name, "type"), MakePPEntry("concentration"));
		if (conc[s].size() > 1) {
		  AddToTable(table, MakePPPrefix("tracer", solutes_[s], "bc", bc_name, "times"), MakePPEntry(ctime[s]));
		  AddToTable(table, MakePPPrefix("tracer", solutes_[s], "bc", bc_name, "forms"), MakePPEntry(cfunc[s]));
		}
	      }
	    }
	  }

	  // Sniff out geochemical conditions, if any.
	  bool gcfound = false;
	  DOMElement* gc = GetChildByName_(lp, "geochemistry_component", gcfound);
	  if (gcfound)
	  {
	    sfound = true;

	    bool gfound;

	    vector<DOMNode*> geocons = GetChildren_(gc, "constraint", gfound);
	    if (gfound)
	    {
	      vector<vector<string> > gicdata;
	      for (int g=0; g<geocons.size(); ++g) {
		DOMElement* geocon = static_cast<DOMElement*>(geocons[g]);
		string constraint_name = GetAttributeValueS_(geocon, "name");
		if (!constraint_name.empty()) {
		  bool known_cond = false;
		  for (int c=0; c<constraint_names_.size() && !known_cond; ++c) {
		    if (constraint_names_[c] == constraint_name) known_cond = true;
		  }
		  if (known_cond) {
		    vector<string> thisgic;
		    thisgic.push_back(constraint_name);
		    string condition_time = GetAttributeValueS_(geocon, "start");
		    if (!condition_time.empty()) {
		      thisgic.push_back(condition_time);
		    }
		    string condition_form = GetAttributeValueS_(geocon, "function");
		    if (!condition_form.empty()) {
		      if (condition_form != "constant") {
			Errors::Message msg;
			msg << "Unsupport non-constant constraint interval specified in boundary_condition \"" << bc_name << "\".\n";
			msg << "Please correct and try again.\n";
			Exceptions::amanzi_throw(msg);	  
		      }
		    }
		    gicdata.push_back(thisgic);
		  }
		  else {
		    Errors::Message msg;
		    msg << "Undeclared constraint name \""  << constraint_name << "\" in boundary_condition \"" << bc_name << "\".\n";
		    msg << "Please correct and try again.\n";
		    Exceptions::amanzi_throw(msg);	  
		  }
		}
		else {
		  Errors::Message msg;
		  msg << "No \"name\" attibute specified for geochemical constraint in boundary_condition \"" << bc_name << "\".\n";
		  msg << "Please correct and try again.\n";
		  Exceptions::amanzi_throw(msg);	  
		}
	      }

	      // Translate any labeled times to real values
	      for (int g=0; g<gicdata.size(); ++g) {
		if (gicdata[g].size() > 1) {
		  const string& this_time = gicdata[g][1];
		  map<string, string>::const_iterator iterT = labeled_times_.find(this_time);
		  if (iterT != labeled_times_.end())
		    gicdata[g][1] = ConvertTimeToSeconds(iterT->second);
		  else
		    gicdata[g][1] = ConvertTimeToSeconds(this_time);
		}
	      }

	      // Check that times provided if more than one constraint supplied
	      if (gicdata.size() > 1) {
		for (int g=0; g<gicdata.size(); ++g) {
		  if (gicdata[g].size() < 2) {
		    Errors::Message msg;
		    msg << "When using multiple geochemical constraints on a boundary, the \"start\" attribule is required for boundary_condition \"" << bc_name << "\".\n";
		    msg << "Please correct and try again.\n";
		    Exceptions::amanzi_throw(msg);	  
		  }
		}
	      }

	      if (gicdata.size() > 1) {
		  Errors::Message msg;
		  msg << "Multiple time intervals not yet supported in AmanziS for geochemical constraint boundary_condition \"" << bc_name << "\".\n";
		  msg << "Please correct and try again.\n";
		  Exceptions::amanzi_throw(msg);	  
	      }
	      else {
		for (int s = 0; s<solutes_.size(); ++s) {
		  AddToTable(table, MakePPPrefix("tracer", solutes_[s], "bc", bc_name, "type"), MakePPEntry("concentration"));
		}
		if (gicdata.size() == 1) {
		  for (int s = 0; s<solutes_.size(); ++s) {
		    AddToTable(table, MakePPPrefix("tracer", solutes_[s], "bc", bc_name, "geochemical_conditions"), MakePPEntry(gicdata[0][0]));
		  }
		}
		else {
		  vector<string> gicnames, gictimes, gicforms;
		  for (int g=0; g<gicdata.size(); ++g) {
		    gicnames.push_back(gicdata[g][0]);
		    gictimes.push_back(gicdata[g][1]);
		    gicnames.push_back("constant");
		  }
		  for (int s = 0; s<solutes_.size(); ++s) {
		    AddToTable(table, MakePPPrefix("tracer", solutes_[s], "bc", bc_name, "geochemical_conditions"), MakePPEntry(gicnames));
		    AddToTable(table, MakePPPrefix("tracer", solutes_[s], "bc", bc_name, "times"), MakePPEntry(gictimes));
		    AddToTable(table, MakePPPrefix("tracer", solutes_[s], "bc", bc_name, "forms"), MakePPEntry(gicforms));
		  }
		}
	      }
	    }
	  }

	  // Assigned regions.
	  bool rfound = false;
	  vector<string> assigned_regions = GetChildVectorS_(bc, "assigned_regions", rfound, true);
	  if (rfound) {
	    AddToTable(table, MakePPPrefix("comp", "bcs", bc_name, "regions"), MakePPEntry(assigned_regions));
	    if (sfound) {
	      tbc_names.push_back(bc_name);
	      for (size_t i = 0; i < solutes_.size(); ++i) {
		AddToTable(table, MakePPPrefix("tracer", solutes_[i], "bc", bc_name, "regions"), MakePPEntry(assigned_regions));
	      }
	    }
	  } else {
	    Errors::Message msg;
	    msg << "\"assigned_regions\" not present in boundary_condition \"" << bc_name << "\".\n";
	    msg << "Please correct and try again.\n";
	    Exceptions::amanzi_throw(msg);
	  }

        }
        else {
          Errors::Message msg;
          msg << "\"liquid_phase\" not present in boundary_condition \"" << bc_name << "\".\n";
          msg << "Please correct and try again.\n";
          Exceptions::amanzi_throw(msg);
        }
      }

      // List boundary conditions where needed.
      AddToTable(table, MakePPPrefix("comp", "bc_labels"), MakePPEntry(bc_names));
      for (size_t i = 0; i < solutes_.size(); ++i)
        AddToTable(table, MakePPPrefix("tracer", solutes_[i], "tbcs"), MakePPEntry(tbc_names));
    }
  }
  ParmParse::appendTable(table);
}

void InputConverterS::ParseOutput_()
{
  list<ParmParse::PP_entry> table;
  bool found;

  // Visualization files.
  DOMNode* vis = GetUniqueElementByTagsString_("output, vis", found);
  if (found)
  {
    bool found;
    string base_filename = GetChildValueS_(vis, "base_filename", found, true);
    string num_digits = GetChildValueS_(vis, "num_digits", found, true);
    vector<string> cycle_macros = GetChildVectorS_(vis, "cycle_macros", found, false);
    vector<string> time_macros = GetChildVectorS_(vis, "time_macros", found, false);

    AddToTable(table, MakePPPrefix("amr", "plot_file"), MakePPEntry(base_filename));
    AddToTable(table, MakePPPrefix("amr", "plot_file_digits"), MakePPEntry(num_digits));
    AddToTable(table, MakePPPrefix("amr", "viz_cycle_macros"), MakePPEntry(cycle_macros));
    AddToTable(table, MakePPPrefix("amr", "viz_time_macros"), MakePPEntry(time_macros));
  }

  // Checkpoint files.
  DOMNode* checkpoint = GetUniqueElementByTagsString_("output, checkpoint", found);
  if (found)
  {
    bool found;
    string base_filename = GetChildValueS_(checkpoint, "base_filename", found, true);
    string num_digits = GetChildValueS_(checkpoint, "num_digits", found, true);
    vector<string> cycle_macros = GetChildVectorS_(checkpoint, "cycle_macros", found, false);

    AddToTable(table, MakePPPrefix("amr", "check_file"), MakePPEntry(base_filename));
    AddToTable(table, MakePPPrefix("amr", "chk_file_digits"), MakePPEntry(num_digits));
    AddToTable(table, MakePPPrefix("amr", "chk_cycle_macros"), MakePPEntry(cycle_macros));
  }

  // Observations.
  DOMNode* observations = GetUniqueElementByTagsString_("output, observations", found);
  if (found)
  {
    bool found;
    string filename = GetChildValueS_(observations, "filename", found, true);
    // string num_digits = GetChildValueS_(checkpoint, "liquid_phase", found, true);
    // FIXME
  }

  // Walkabouts. 
  DOMNode* walkabout = GetUniqueElementByTagsString_("output, walkabout", found);
  if (found)
  {
    bool found;
    string base_filename = GetChildValueS_(walkabout, "base_filename", found, true);
    string num_digits = GetChildValueS_(walkabout, "num_digits", found, true);
    vector<string> cycle_macros = GetChildVectorS_(walkabout, "cycle_macros", found, true);

// FIXME
//    AddToTable(table, MakePPPrefix("amr", "check_file"), MakePPEntry(base_filename));
//    AddToTable(table, MakePPPrefix("amr", "chk_file_digits"), MakePPEntry(num_digits));
//    AddToTable(table, MakePPPrefix("amr", "chk_cycle_macros"), MakePPEntry(cycle_macros));
  }

  ParmParse::appendTable(table);
}

void InputConverterS::ParseMisc_()
{
  // FIXME: Not yet supported.
}

void InputConverterS::Translate(int rank) 
{
  rank_ = rank;
  solute_names_parsed_ = false;;
  ParseUnits_();
  ParseDefinitions_();
  ParseExecutionControls_();
  ParseMesh_();
  ParseRegions_();
  ParseSoluteNames_();
  ParseGeochemistry_();
  bool do_tracer_diffusion = false;
  ParseMaterials_(do_tracer_diffusion);
  ParsePhases_(do_tracer_diffusion);
  string flow_model;
  ParseProcessKernels_(flow_model,do_tracer_diffusion);
  ParseInitialConditions_();
  ParseBoundaryConditions_();
  ParseNumericalControls_(flow_model);
  ParseOutput_();
  ParseMisc_();
}

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
