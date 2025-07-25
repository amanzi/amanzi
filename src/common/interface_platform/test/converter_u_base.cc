/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  This is the input component of the Amanzi code.

*/

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"
#include "XMLParameterListWriter.hh"
#include "Epetra_MpiComm.h"

// Amanzi
#include "InputConverterU.hh"


bool
ComparePLists(const Teuchos::ParameterList& plist1,
              const Teuchos::ParameterList& plist2,
              std::string& name)
{
  for (auto it = plist1.begin(); it != plist1.end(); ++it) {
    name = plist1.name(it);
    if (plist1.isSublist(name)) {
      if (!plist2.isSublist(name) ) return false;

      std::string name_tmp;
      bool flag = ComparePLists(plist1.sublist(name), plist2.sublist(name), name_tmp);
      if (!flag) {
        name.append("->");
        name.append(name_tmp);
        return false;
      }

    } else if (plist1.isParameter(name)) {
      if (!plist2.isParameter(name) ) return false;
      const Teuchos::ParameterEntry& e1 = plist1.getEntry(name);
      const Teuchos::ParameterEntry& e2 = plist2.getEntry(name);

      if (e1.isType<double>()) {
        if (!e2.isType<double>() ) return false;
        double v1 = plist1.get<double>(name);
        double v2 = plist2.get<double>(name);
        if (fabs(v1 - v2) > fabs(v1) * 1e-12) return false;
      }

      if (e1.isType<std::string>()) {
        if (!e2.isType<std::string>() ) return false;
        std::string v1 = plist1.get<std::string>(name);
        std::string v2 = plist2.get<std::string>(name);
        if (v1 != v2) return false;
      }

      if (e1.isType<int>()) {
        if (!e2.isType<int>() ) return false;
        int v1 = plist1.get<int>(name);
        int v2 = plist2.get<int>(name);
        if (v1 != v2) return false;
      }

      if (e1.isType<Teuchos::Array<std::string>>()) {
        if (!e2.isType<Teuchos::Array<std::string>>() ) return false;
        std::vector<std::string> v1 = plist1.get<Teuchos::Array<std::string>>(name).toVector();
        std::vector<std::string> v2 = plist2.get<Teuchos::Array<std::string>>(name).toVector();
        if (v1.size() != v2.size()) return false;
        for (int i = 0; i < v1.size(); ++i)
          if (v1[i] != v2[i]) return false;
      }

      if (e1.isType<Teuchos::Array<int>>()) {
        if (!e2.isType<Teuchos::Array<int>>() ) return false;
        std::vector<int> v1 = plist1.get<Teuchos::Array<int>>(name).toVector();
        std::vector<int> v2 = plist2.get<Teuchos::Array<int>>(name).toVector();
        if (v1.size() != v2.size()) return false;
        for (int i = 0; i < v1.size(); ++i)
          if (v1[i] != v2[i]) return false;
      }

      if (e1.isType<Teuchos::Array<double>>()) {
        if (!e2.isType<Teuchos::Array<double>>() ) return false;
        std::vector<double> v1 = plist1.get<Teuchos::Array<double>>(name).toVector();
        std::vector<double> v2 = plist2.get<Teuchos::Array<double>>(name).toVector();
        if (v1.size() != v2.size()) return false;
        for (int i = 0; i < v1.size(); ++i)
          if (fabs(v1[i] - v2[i]) > fabs(v1[i]) * 1e-12) return false;
      }
    }
  }
  return true;
}


/* **************************************************************** */
TEST(CONVERTER_BASE)
{
  using namespace Teuchos;
  using namespace Amanzi;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  // read parameter list
  for (int i = 1; i < 11; i++) {
    std::stringstream xmlFileName, id;
    id << std::setw(2) << std::setfill('0') << i;
    xmlFileName << "test/converter_u_test" << id.str() << ".xml";

    if (MyPID == 0)
      std::cout << std::endl
                << std::endl
                << "Test " << xmlFileName.str() << ": convert xml" << std::endl
                << "---------------------------------------------" << std::endl;

    Amanzi::AmanziInput::InputConverterU converter(xmlFileName.str());
    Teuchos::ParameterList new_xml;
    try {
      new_xml = converter.Translate(0, 1);

      Teuchos::Amanzi_XMLParameterListWriter XMLWriter;
      Teuchos::XMLObject XMLobj = XMLWriter.toXML(new_xml);

      std::stringstream ss;
      ss << "test" << id.str() << "_native.xml";
      std::ofstream xmlfile;
      xmlfile.open(ss.str().c_str());
      xmlfile << XMLobj;

      // development
      Teuchos::RCP<Teuchos::ParameterList> old_xml;
      xmlFileName.str("");
      xmlFileName << "test/converter_u_validate" << id.str() << ".xml";

      std::cout << std::endl
                << "Successful translation. Validating the result..." << std::endl
                << "    gold file: " << xmlFileName.str() << std::endl
                << "    native file: " << ss.str() << std::endl
                << std::endl;

      old_xml = Teuchos::getParametersFromXmlFile(xmlFileName.str());
      new_xml.validateParameters(*old_xml);
      old_xml->validateParameters(new_xml);

      // data validation
      std::string name;
      bool flag = ComparePLists(new_xml, *old_xml, name);
      if (!flag) {
        std::cout << "Test:" << i << ", error at \"" << name << "\"." << std::endl;
        CHECK(false);
        break;
      } else {
        std::cout << "Success" << std::endl;
      }
    } catch (std::exception& e) {
      std::cout << "Test:" << i << " threw exception:" << std::endl;
      std::cout << "  \"" << e.what() << "\"" << std::endl;
      CHECK(false);
      break;
    }
  }
}
