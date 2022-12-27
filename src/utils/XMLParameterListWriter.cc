/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Brendt
*/

/*
  Utils

*/

#include <string>
#include "XMLParameterListWriter.hh"

namespace Teuchos {

Amanzi_XMLParameterListWriter::Amanzi_XMLParameterListWriter()
{
  precision_ = 17; // round(52 * log10(2)) + 1 = 17 digits
}


XMLObject
Amanzi_XMLParameterListWriter::toXML(const ParameterList& p) const
{
  XMLObject rtn("ParameterList");

  for (auto i = p.begin(); i != p.end(); ++i) {
    const ParameterEntry& val = p.entry(i);
    const std::string& name = p.name(i);
    XMLObject child = toXML(val);
    child.addAttribute("name", name);
    rtn.addChild(child);
  }

  return rtn;
}


XMLObject
Amanzi_XMLParameterListWriter::toXML(const ParameterEntry& entry) const
{
  if (entry.isList()) { return toXML(getValue<ParameterList>(entry)); }

  XMLObject rtn("Parameter");
  std::string type;
  std::string value;

  if (entry.isType<int>()) {
    type = "int";
    value = toString(any_cast<int>(entry.getAny(false)));
  } else if (entry.isType<short>()) {
    type = "short";
    value = toString(any_cast<short>(entry.getAny(false)));
  } else if (entry.isType<double>()) {
    type = "double";
    value = toString(any_cast<double>(entry.getAny(false)));
  } else if (entry.isType<float>()) {
    type = "float";
    value = toString(any_cast<float>(entry.getAny(false)));
  } else if (entry.isType<std::string>()) {
    type = "string";
    value = toString(any_cast<std::string>(entry.getAny(false)));
  } else if (entry.isType<char>()) {
    type = "char";
    value = toString(any_cast<char>(entry.getAny(false)));
  } else if (entry.isType<bool>()) {
    type = "bool";
    value = toString(any_cast<bool>(entry.getAny(false)));
  } else if (entry.isType<Array<int>>()) {
    const Array<int>& a = any_cast<Array<int>>(entry.getAny(false));
    type = "Array(int)";
    value = a.toString();
  } else if (entry.isType<Array<short>>()) {
    const Array<short>& a = any_cast<Array<short>>(entry.getAny(false));
    type = "Array(short)";
    value = a.toString();
  } else if (entry.isType<Array<float>>()) {
    const Array<float>& a = any_cast<Array<float>>(entry.getAny(false));
    type = "Array(float)";
    value = a.toString();
  } else if (entry.isType<Array<double>>()) {
    const Array<double>& a = any_cast<Array<double>>(entry.getAny(false));
    type = "Array(double)";
    value = Amanzi_toString(a);
  } else if (entry.isType<Array<std::string>>()) {
    const Array<std::string>& a = any_cast<Array<std::string>>(entry.getAny(false));
    type = "Array(string)";
    value = a.toString();
  } else {
    type = "any";
    std::ostringstream ss;
    ss << entry;
    value = ss.str().c_str();
  }

  rtn.addAttribute("type", type);
  rtn.addAttribute("value", value);

  if (entry.isDefault()) { rtn.addAttribute("isDefault", "true"); }

  if (entry.isUsed()) { rtn.addAttribute("isUsed", "true"); }

  return rtn;
}

} // namespace Teuchos
