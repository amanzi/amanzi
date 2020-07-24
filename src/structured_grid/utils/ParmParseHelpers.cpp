#include <BoxLib.H>
#include <ParmParseHelpers.H>
#include <ParallelDescriptor.H>
#include <sstream>
#include <strstream>

#include "Teuchos_StrUtils.hpp"
#include "Utility.H"

typedef std::list<ParmParse::PP_entry>::iterator list_iterator;
typedef std::list<ParmParse::PP_entry>::const_iterator const_list_iterator;

static std::list<std::string>
split_string_to_list(const std::string& str)
{
  std::vector<std::string> tokens = BoxLib::Tokenize(str," ,");
  std::list<std::string> token_list;
  for (int i=0; i<tokens.size(); ++i)
    {
      token_list.push_back(tokens[i]);
    }
  return token_list;
}

static std::stack<std::string> prefix;

static std::string
build_prefixed_name(const std::string& name)
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
bldTable (Teuchos::ParameterList& params,
	  std::list<ParmParse::PP_entry>& tab)
{    
  for (Teuchos::ParameterList::ConstIterator i=params.begin(); i!=params.end(); ++i)
    {
      const std::string& name = params.name(i);
      const Teuchos::ParameterEntry& entry = params.getEntry(name);

      if (entry.isList() )
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

            std::string prefixed_name = build_prefixed_name(name);
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

static void
print_table (const std::string& pfx, const ParmParse::Table& table)
{

  for ( const_list_iterator li = table.begin(); li != table.end(); ++li )
    {
      if ( li->m_table )
	{
          std::cout << "Record " << li->m_name << std::endl;
          std::string prefix;
          if (pfx!=std::string())
            prefix = pfx + "::";
          prefix += li->m_name;
          print_table(prefix, *li->m_table);

	}
      else
	{
          std::string prefix;
          if (pfx!=std::string())
            prefix = pfx + "::";
          std::cout << prefix << li->m_name << "(nvals = " << li->m_vals.size() << ") " << " :: [";
          int n = li->m_vals.size();
          for ( int i = 0; i < n; i++ )
            {
              std::cout << li->m_vals[i];
              if ( i < n-1 ) std::cout << ", ";
            }
          std::cout << "]" << '\n';
          
	}
    }
}

void
BoxLib::Initialize_ParmParse(Teuchos::ParameterList& params)
{
  std::list<ParmParse::PP_entry> table;
  bldTable(params,table);
  ParmParse::appendTable(table);
}


