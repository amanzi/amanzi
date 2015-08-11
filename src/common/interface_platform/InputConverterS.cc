/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson (jnjohnson@lbl.gov)
*/

#include "BoxLib.H"
#include "ParmParse.H"
#include "InputConverterS.hh"

namespace Amanzi {
namespace AmanziInput {

// For saving a bit of typing.
using namespace std;
using namespace xercesc;

// Internally-used stuff.
namespace {

typedef list<ParmParse::PP_entry>::iterator list_iterator;
typedef list<ParmParse::PP_entry>::const_iterator const_list_iterator;

// Construct a ParmParse prefix name from sets of strings.
string MakePPPrefix(const string& s1)
{
  return s1;
}

string MakePPPrefix(const string& s1, const string& s2)
{
  return s1 + string(".") + s2;
}

string MakePPPrefix(const string& s1, const string& s2, const string& s3)
{
  return s1 + string(".") + s2 + string(".") + s3;
}

string MakePPPrefix(const string& s1, const string& s2, const string& s3, const string& s4)
{
  return s1 + string(".") + s2 + string(".") + s3 + string(".") + s4;
}

// Construct a ParmParse entry from sets of strings/values.
list<string> MakePPEntry(const string& s1)
{
  list<string> pp;
  pp.push_back(s1);
  return pp;
}

list<string> MakePPEntry(double d)
{
  list<string> pp;
  stringstream s;
  s << d << '\0';
  pp.push_back(s.str());
  return pp;
}

list<string> MakePPEntry(int i)
{
  list<string> pp;
  stringstream s;
  s << i << '\0';
  pp.push_back(s.str());
  return pp;
}

list<string> MakePPEntry(long i)
{
  list<string> pp;
  stringstream s;
  s << i << '\0';
  pp.push_back(s.str());
  return pp;
}

list<string> MakePPEntry(const string& s1, const string& s2)
{
  list<string> pp;
  pp.push_back(s1);
  pp.push_back(s2);
  return pp;
}

list<string> MakePPEntry(const vector<string>& ss)
{
  list<string> pp;
  for (size_t i = 0; i < ss.size(); ++i)
    pp.push_back(ss[i]);
  return pp;
}

list<string> MakePPEntry(const vector<double>& ds)
{
  list<string> pp;
  for (size_t i = 0; i < ds.size(); ++i)
  {
    stringstream s;
    s << ds[i] << '\0';
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
    s << is[i] << '\0';
    pp.push_back(s.str());
  }
  return pp;
}

list<string> MakePPEntry(const vector<long>& is)
{
  list<string> pp;
  for (size_t i = 0; i < is.size(); ++i)
  {
    stringstream s;
    s << is[i] << '\0';
    pp.push_back(s.str());
  }
  return pp;
}

#if 0
static std::stack<std::string> prefix;

static std::string
buildPrefixedName(const std::string& name)
{
  std::string result = name;

  std::stack<std::string> pre_tmp = prefix;
  while ( ! pre_tmp.empty() )
    {
      result = pre_tmp.top() + "." + result;
      pre_tmp.pop();
    }
  return result;
}

void
buildTable(xercesc::DOMDocument* doc,
           std::list<ParmParse::PP_entry>& tab)
{    
  for (Teuchos::ParameterList::ConstIterator i=params.begin(); i!=params.end(); ++i)
    {
      const std::string& name = params.name(i);
      const Teuchos::ParameterEntry& val = params.getEntry(name);

      if (val.isList() )
        {
          prefix.push(name);
          bldTable(params.sublist(name), tab);          
        }
      else
        {
            // FIXME: It is unfortunate that Teuchos only stores the converted
            //  data, and this the Teuchos::toString methods apply an arbitrary 
            //  formatting rules buried in the bowels of trilinos...

            std::stringstream ppStr;
            std::ios::fmtflags oflags = ppStr.flags();
            ppStr.setf(std::ios::floatfield, std::ios::scientific);
            int old_prec = ppStr.precision(15);

            std::string prefixed_name = buildPrefixedName(name);
            std::list<std::string> ppStrList;

            Teuchos::ParameterEntry* entry = params.getEntryPtr(name);
            if (entry->isType<double>()) {
                double val = entry->getValue<double>(&val);
                ppStr << val; 
                ppStrList.push_back(ppStr.str());
            }
            else if (entry->isType<float>()) {
                float val = entry->getValue<float>(&val);
                ppStr << val;
                ppStrList.push_back(ppStr.str());
            }
            else if (entry->isType<short>()) {
                short val = entry->getValue<short>(&val);
                ppStr << val;
                ppStrList.push_back(ppStr.str());
            }
            else if (entry->isType<int>()) {
                int val = entry->getValue<int>(&val);
                ppStr << val;
                ppStrList.push_back(ppStr.str());
            }
            else if (entry->isType<bool>()) {
                bool val = entry->getValue<bool>(&val);
                ppStr << val;
                ppStrList.push_back(ppStr.str());
            }
            else if (entry->isType<std::string>()) {
                std::string val = entry->getValue<std::string>(&val);
                ppStr << val;
                ppStrList.push_back(ppStr.str());
            }
            else if (entry->isType<Teuchos::Array<int> >()) {
                Teuchos::Array<int> val = entry->getValue<Teuchos::Array<int> >(&val);
                for (int i=0; i<val.size(); ++i) {
                    ppStr.str(""); ppStr << val[i]; ppStrList.push_back(ppStr.str());
                }
            }
            else if (entry->isType<Teuchos::Array<short> >()) {
                Teuchos::Array<short> val = entry->getValue<Teuchos::Array<short> >(&val);
                for (int i=0; i<val.size(); ++i) {
                    ppStr.str(""); ppStr << val[i]; ppStrList.push_back(ppStr.str());
                }
            }
            else if (entry->isType<Teuchos::Array<float> >()) {
                Teuchos::Array<float> val = entry->getValue<Teuchos::Array<float> >(&val);
                for (int i=0; i<val.size(); ++i) {
                    ppStr.str(""); ppStr << val[i]; ppStrList.push_back(ppStr.str());
                }
            }
            else if (entry->isType<Teuchos::Array<double> >()) {
                Teuchos::Array<double> val = entry->getValue<Teuchos::Array<double> >(&val);
                for (int i=0; i<val.size(); ++i) {
                    ppStr.str(""); ppStr << val[i]; ppStrList.push_back(ppStr.str());
                }
            }
            else if (entry->isType<Teuchos::Array<std::string> >()) {
                Teuchos::Array<std::string> val = entry->getValue<Teuchos::Array<std::string> >(&val);
                for (int i=0; i<val.size(); ++i) {
                    ppStr.str(""); ppStr << val[i]; ppStrList.push_back(ppStr.str());
                }
            }
            else {
                BoxLib::Abort("Type is not supported");
            }

            tab.push_back(ParmParse::PP_entry(prefixed_name,ppStrList));

        }
    }
  if ( ! prefix.empty() )
    prefix.pop();
}
#endif

}

InputConverterS::InputConverterS():
  InputConverter()
{
}

InputConverterS::~InputConverterS()
{
}

void InputConverterS::ParseUnits()
{
  // Units are supposedly not supported yet. I guess we'll find out!
}

void InputConverterS::ParseDefinitions()
{
  std::list<ParmParse::PP_entry> table;
  bool found;
  DOMNode* macros = getUniqueElementByTagNames_("definitions", "macros", found);
  if (found)
  {
    bool found;
    vector<DOMNode*> time_macros = getChildren_(macros, "time_macro", found);
    if (found)
    {
      for (size_t i = 0; i < time_macros.size(); ++i)
      {
        DOMElement* time_macro = static_cast<DOMElement*>(time_macros[i]);
        string macro_name = GetAttributeValueS_(time_macro, "name");
        vector<string> times;

        // Before we look for specific times, check for other stuff.
        bool found;
        string start = GetChildValueS_(time_macro, "start", found);
        if (found)
        {
          // We've got one of the interval-based time macros.
          // Get the other required elements.
          string timestep_interval = GetChildValueS_(time_macro, "timestep_interval", found, true);
          string stop = GetChildValueS_(time_macro, "stop", found, true);
          
          char* endptr = const_cast<char*>(start.c_str());
          double t1 = strtod(start.c_str(), &endptr);
          if (endptr == start.c_str()) // Not parsed!
            ThrowErrorIllformed_("definitions->macros", "start", "time_macro");
          endptr = const_cast<char*>(stop.c_str());
          double t2 = strtod(stop.c_str(), &endptr);
          if (endptr == stop.c_str()) // Not parsed!
            ThrowErrorIllformed_("definitions->macros", "stop", "time_macro");
          endptr = const_cast<char*>(timestep_interval.c_str());
          double interval = strtod(timestep_interval.c_str(), &endptr);
          if (endptr == timestep_interval.c_str()) // Not parsed!
            ThrowErrorIllformed_("definitions->macros", "timestep_interval", "time_macro");

          // Shove this macro into our table.
          table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "time_macros", macro_name, "type"),
                                              MakePPEntry("period")));
          table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "time_macros", macro_name, "start"),
                                              MakePPEntry(start)));
          table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "time_macros", macro_name, "stop"),
                                              MakePPEntry(stop)));
          table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "time_macros", macro_name, "period"),
                                              MakePPEntry(interval)));
        }
        else
        {
          // We're just looking for times.
          bool found;
          vector<DOMNode*> time_nodes = getChildren_(time_macro, "time", found, true);
          vector<string> times;
          for (size_t j = 0; j < time_nodes.size(); ++j)
            times.push_back(XMLString::transcode(time_nodes[j]->getTextContent()));

          table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "time_macros", macro_name, "type"),
                                              MakePPEntry("times")));
          table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "time_macros", macro_name, "times"),
                                              MakePPEntry(times)));
        }
      }
    }

    vector<DOMNode*> cycle_macros = getChildren_(macros, "cycle_macro", found);
    if (found)
    {
      for (size_t i = 0; i < cycle_macros.size(); ++i)
      {
        DOMElement* cycle_macro = static_cast<DOMElement*>(cycle_macros[i]);
        string macro_name = GetAttributeValueS_(cycle_macro, "name");
        vector<string> cycles;

        bool found;
        string start = GetChildValueS_(cycle_macro, "start", found, true);
        string timestep_interval = GetChildValueS_(cycle_macro, "timestep_interval", found, true);
        string stop = GetChildValueS_(cycle_macro, "stop", found, true);
         
        // Verify the entries.
        char* endptr = const_cast<char*>(start.c_str());
        long c1 = strtol(start.c_str(), &endptr, 10);
        if (endptr == start.c_str()) // Not parsed!
          ThrowErrorIllformed_("definitions->macros", "start", "cycle_macro");
        endptr = const_cast<char*>(stop.c_str());
        long c2 = strtol(stop.c_str(), &endptr, 10);
        if (endptr == stop.c_str()) // Not parsed!
          ThrowErrorIllformed_("definitions->macros", "stop", "cycle_macro");
        endptr = const_cast<char*>(timestep_interval.c_str());
        long interval = strtol(timestep_interval.c_str(), &endptr, 10);
        if (endptr == timestep_interval.c_str()) // Not parsed!
          ThrowErrorIllformed_("definitions->macros", "timestep_interval", "cycle_macro");

        // Shove this macro into our table.
        table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "cycle_macro", macro_name, "type"),
                                            MakePPEntry("period")));
        table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "cycle_macro", macro_name, "start"),
                                            MakePPEntry(start)));
        table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "cycle_macro", macro_name, "stop"),
                                            MakePPEntry(stop)));
        table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "cycle_macro", macro_name, "period"),
                                            MakePPEntry(interval)));
      }
    }

    // FIXME: variable_macro not yet supported.
  }

  if (!table.empty())
    ParmParse::appendTable(table);
}

void InputConverterS::ParseExecutionControls()
{
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
}

void InputConverterS::ParseNumericalControls()
{
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
}

void InputConverterS::ParseMesh()
{
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
}

void InputConverterS::ParseRegions()
{
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
}

void InputConverterS::ParseGeochemistry()
{
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
}

void InputConverterS::ParseMaterials()
{
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
}

void InputConverterS::ParseProcessKernels()
{
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
}

void InputConverterS::ParsePhases()
{
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
}

void InputConverterS::ParseInitialConditions()
{
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
}

void InputConverterS::ParseBoundaryConditions()
{
  std::list<ParmParse::PP_entry> table;
  bool found;
  DOMNode* boundary_conditions = getUniqueElementByTagNames_(doc_, "boundary_conditions", found);
  if (found)
  {
    bool found;
    vector<DOMNode*> bcs = getChildren_(boundary_conditions, "boundary_condition", found);
    if (found)
    {
      for (size_t i = 0; i < bcs.size(); ++i)
      {
        DOMElement* bc = static_cast<DOMElement*>(bcs[i]);
        string bc_name = GetAttributeValueS_(bc, "name");
        vector<string> times;
      }
    }
  }
  if (!table.empty())
    ParmParse::appendTable(table);
}

void InputConverterS::ParseOutput()
{
  std::list<ParmParse::PP_entry> table;
  bool found;

  // Visualization files.
  DOMNode* vis = getUniqueElementByTagNames_("output", "vis", found);
  if (found)
  {
    bool found;
    string base_filename = GetChildValueS_(vis, "base_filename", found, true);
    string num_digits = GetChildValueS_(vis, "num_digits", found, true);
    vector<string> cycle_macros = GetChildVectorS_(vis, "cycle_macros", found, false);
    vector<string> time_macros = GetChildVectorS_(vis, "time_macros", found, false);

    table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "plot_file"),
                                        MakePPEntry(base_filename)));
    table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "plot_file_digits"),
                                        MakePPEntry(num_digits)));
    table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "viz_cycle_macros"),
                                        MakePPEntry(cycle_macros)));
    table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "viz_time_macros"),
                                        MakePPEntry(time_macros)));
  }

  // Checkpoint files.
  DOMNode* checkpoint = getUniqueElementByTagNames_("output", "checkpoint", found);
  if (found)
  {
    bool found;
    string base_filename = GetChildValueS_(checkpoint, "base_filename", found, true);
    string num_digits = GetChildValueS_(checkpoint, "num_digits", found, true);
    vector<string> cycle_macros = GetChildVectorS_(checkpoint, "cycle_macros", found, true);

    table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "check_file"),
                                        MakePPEntry(base_filename)));
    table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "chk_file_digits"),
                                        MakePPEntry(num_digits)));
    table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "chk_cycle_macros"),
                                        MakePPEntry(cycle_macros)));
  }

  // Observations.
  DOMNode* observations = getUniqueElementByTagNames_("output", "observations", found);
  if (found)
  {
    bool found;
    string filename = GetChildValueS_(checkpoint, "filename", found, true);
//    string num_digits = GetChildValueS_(checkpoint, "liquid_phase", found, true);

    // FIXME
  }

  // Walkabouts. 
  DOMNode* walkabout = getUniqueElementByTagNames_("output", "walkabout", found);
  if (found)
  {
    bool found;
    string base_filename = GetChildValueS_(walkabout, "base_filename", found, true);
    string num_digits = GetChildValueS_(walkabout, "num_digits", found, true);
    vector<string> cycle_macros = GetChildVectorS_(walkabout, "cycle_macros", found, true);

// FIXME
//    table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "check_file"),
//                                        MakePPEntry(base_filename)));
//    table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "chk_file_digits"),
//                                        MakePPEntry(num_digits)));
//    table.push_back(ParmParse::PP_entry(MakePPPrefix("amr", "chk_cycle_macros"),
//                                        MakePPEntry(cycle_macros)));
  }

  if (!table.empty())
    ParmParse::appendTable(table);
}

void InputConverterS::ParseMisc()
{
  // FIXME: Not yet supported.
}

void InputConverterS::Translate() 
{
  ParseUnits();
  ParseDefinitions();
  ParseExecutionControls();
  ParseNumericalControls();
  ParseMesh();
  ParseRegions();
  ParseGeochemistry();
  ParseMaterials();
  ParseProcessKernels();
  ParsePhases();
  ParseInitialConditions();
  ParseBoundaryConditions();
  ParseOutput();
  ParseMisc();
}

}  // namespace AmanziInput
}  // namespace Amanzi

