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

// Internally-used stuff.
namespace {

typedef std::list<ParmParse::PP_entry>::iterator list_iterator;
typedef std::list<ParmParse::PP_entry>::const_iterator const_list_iterator;

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
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
}

void InputConverterS::ParseDefinitions()
{
  std::list<ParmParse::PP_entry> table;
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
  ParmParse::appendTable(table);
}

void InputConverterS::ParseOutput()
{
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
}

void InputConverterS::ParseMisc()
{
  std::list<ParmParse::PP_entry> table;
  ParmParse::appendTable(table);
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

