/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Erin Barker (original version)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>

// TPLs
#include "boost/bind.hpp"
#include "boost/algorithm/string.hpp"

#define  BOOST_FILESYTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"

#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/DOMLSParserImpl.hpp>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>

#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Amanzi's
#include "ErrorHandler.hpp"
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputConverterU.hh"


namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE


/* ******************************************************************
* Translate unstructured mesh. Introduces global parameter dim_.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateMesh_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating unstructured mesh" << std::endl;
  }

  MemoryManager mm;
  DOMNodeList *node_list, *children;
  DOMNode* node;
  DOMElement* element;

  bool flag, read(false), generate(false);
  std::string framework, verify;
  Errors::Message msg;
  Teuchos::ParameterList mesh_list;
    
  // read in new stuff
  node_list = doc_->getElementsByTagName(mm.transcode("mesh"));
  if (node_list->getLength() == 0)
      ThrowErrorIllformed_("mesh", "element", "framework");

  DOMNode* inode = node_list->item(0);
  element = static_cast<DOMElement*>(inode);
  framework = GetAttributeValueS_(element, "framework");

  // Define global parameter dim_ = the space dimension.
  node = GetUniqueElementByTagsString_(inode, "dimension", flag);
  if (flag) {
    char* tmp = mm.transcode(node->getTextContent());
    dim_ = std::strtol(tmp, NULL, 10);
  } 

  if (!flag || dim_ <= 0)
      ThrowErrorIllformed_("mesh", "dimension", "dimension");

  // Now we can properly parse the generate/read list.
  children = inode->getChildNodes();

  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      char* tagname = mm.transcode(inode->getNodeName());   

      // A structured mesh is generated.
      if (strcmp(tagname, "generate") == 0) {
        generate = true;
        mesh_rectangular_ = "true";
        node = GetUniqueElementByTagsString_(inode, "number_of_cells", flag);
        if (!flag) 
            ThrowErrorIllformed_("mesh", "number_of_cells", "generate");
        element = static_cast<DOMElement*>(node);

        std::vector<int> ncells; 
        int nx = GetAttributeValueL_(element, "nx");
        if (nx > 0) ncells.push_back(nx);
        int ny = GetAttributeValueL_(element, "ny", TYPE_NUMERICAL, false, 0);
        if (ny > 0) ncells.push_back(ny); 
        int nz = GetAttributeValueL_(element, "nz", TYPE_NUMERICAL, false, 0);
        if (nz > 0) ncells.push_back(nz); 

        if (ncells.size() != dim_) 
            ThrowErrorIllformed_("mesh", "number_of_cells", "generate");

        // get Box
        node = GetUniqueElementByTagsString_(inode, "box", flag);
        if (!flag) 
            ThrowErrorIllformed_("mesh", "box", "generate");
        element = static_cast<DOMElement*>(node);

        std::string tmp = GetAttributeValueS_(element, "low_coordinates");
        std::vector<double> low = MakeCoordinates_(tmp);
        if (low.size() != dim_)
            ThrowErrorIllformed_("mesh", "low_coordinates", "generate");

        tmp = GetAttributeValueS_(element, "high_coordinates");
        std::vector<double> high = MakeCoordinates_(tmp);
        if (high.size() != dim_)
            ThrowErrorIllformed_("mesh", "high_coordinates", "generate");

        mesh_list.set<Teuchos::Array<int> >("Number of Cells", ncells);
        mesh_list.set<Teuchos::Array<double> >("Domain Low Coordinate", low);
        mesh_list.set<Teuchos::Array<double> >("Domain High Coordinate", high);
      }

      // Un unstructured mesh will be read from a file.
      else if (strcmp(tagname, "read") == 0) {
        bool flag1, flag2, flag3;

        node = GetUniqueElementByTagsString_(inode, "format", flag1);
        if (flag1) {
          std::string format = GetTextContentS_(node, "exodus ii, exodus II, Exodus II, Exodus ii, H5M, h5m");

          if (boost::iequals(format, "exodus ii")) {
            mesh_list.set<std::string>("Format", "Exodus II");
          } else if (boost::iequals(format, "h5m")) {
            mesh_list.set<std::string>("Format", "H5M");
          }
        }

        node = GetUniqueElementByTagsString_(inode, "file", flag2);
        if (flag2) {
          std::string filename = TrimString_(mm.transcode(node->getTextContent()));
          flag2 = (filename.size() > 0);
          if (flag2) {
            if (num_proc_ > 1) {
              std::string par_filename(filename);
              par_filename.replace(par_filename.size() - 4, 4, ".par");

              // attach the right extensions as required by Nemesis file naming conventions
              // in which files are named as mymesh.par.N.r where N = numproc and r is rank
              int ndigits = (int)floor(log10(num_proc_)) + 1;
              std::string fmt = boost::str(boost::format("%%s.%%d.%%0%dd") % ndigits);
              std::string tmp = boost::str(boost::format(fmt) % par_filename % num_proc_ % rank_);
              boost::filesystem::path p(tmp);

              if (boost::filesystem::exists(p)) filename = par_filename;
            }
            mesh_list.set<std::string>("File", filename);
          } 
        }
        node = GetUniqueElementByTagsString_(inode, "verify", flag3);
        if (flag3) {
          verify = mm.transcode(node->getTextContent());
        }
        read = flag1 && flag2;
      }
    }
  }

  if (generate || read) {
    Teuchos::ParameterList& tmp_list = out_list.sublist("Unstructured").sublist("Expert");
    if (strcmp(framework.c_str(), "mstk") == 0) {
      tmp_list.set<std::string>("Framework", "MSTK");
    } else if (strcmp(framework.c_str() ,"moab") == 0) {
      tmp_list.set<std::string>("Framework", "MOAB");
    } else if (strcmp(framework.c_str(), "simple") == 0) {
      tmp_list.set<std::string>("Framework", "Simple");
    } else if (strcmp(framework.c_str(), "stk::mesh") == 0) {
      tmp_list.set<std::string>("Framework", "stk::mesh");
    } else {
      msg << "Amanzi::InputConverter: an error occurred during parsing mesh.\n"
          << "  Unknown framework \"" << framework << "\".\n";
      Exceptions::amanzi_throw(msg); 
    }
    if (strcmp(verify.c_str(), "true") == 0) {
      tmp_list.set<bool>("Verify Mesh", (strcmp(verify.c_str(), "true") == 0));
    }
  }

  if (generate) {
    out_list.sublist("Unstructured").sublist("Generate Mesh") = mesh_list;
  } else if (read) {
    out_list.sublist("Unstructured").sublist("Read Mesh File") = mesh_list;
  } else {
    msg << "Amanzi::InputConverter: an error occurred during parsing mesh.\n";
    Exceptions::amanzi_throw(msg);
  }

  return out_list;
}


/* ******************************************************************
* Convert regions.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateRegions_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Translating regions" << std::endl;

  MemoryManager mm;
  DOMNodeList* node_list;
  DOMNode *node, *node_attr;
  DOMElement* element;
  DOMNamedNodeMap* attr_map;

  char *tagname, *node_name;
  char *text_content, *text_content2;
  std::string reg_name, text;

  // get regions node
  node_list = doc_->getElementsByTagName(mm.transcode("regions"));
  node = node_list->item(0);

  // new options: comment, region, box, point
  DOMNodeList* childern = node->getChildNodes();
  for (int i = 0; i < childern->getLength(); i++) {
    DOMNode* inode = childern->item(i);

    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      tagname = mm.transcode(inode->getNodeName());
      bool have_name(false);
      
      DOMElement* reg_elem;
      // If region is under a region tag, get get the name
      // and set child element as region element.
      if (strcmp(tagname, "region") == 0) {
        reg_name = GetAttributeValueS_(static_cast<DOMElement*>(inode), "name");
        have_name = true;

        // loop over children to get region element
        DOMNodeList* kids = inode->getChildNodes();
        for (int j = 0; j < kids->getLength(); j++) {
          DOMNode* jnode = kids->item(j);

          if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
            node_name = mm.transcode(jnode->getNodeName());
            if (strcmp(node_name, "comments") != 0) 
              reg_elem = static_cast<DOMElement*>(jnode);
          }
        }
      } else {
        // else set the current element as region element
        reg_elem = static_cast<DOMElement*>(inode);
      }
      
      // get reg_elem type
      node_name = mm.transcode(reg_elem->getNodeName());
      
      // get name if needed
      if (!have_name) {
        if (strcmp(node_name, "comments") != 0) {
          reg_name = GetAttributeValueS_(reg_elem, "name");
        }
      }
      
      // loop over attributes of region element
      if (strcmp(node_name, "box") == 0) {
        tree_["regions"].push_back(reg_name);
        
        std::vector<double> low = GetAttributeVector_(reg_elem, "low_coordinates");
        out_list.sublist(reg_name).sublist("Region: Box").set<Teuchos::Array<double> >("Low Coordinate", low);
        
        std::vector<double> high = GetAttributeVector_(reg_elem, "high_coordinates");
        out_list.sublist(reg_name).sublist("Region: Box").set<Teuchos::Array<double> >("High Coordinate", high);
      }

      else if (strcmp(node_name, "plane") == 0) {
        tree_["regions"].push_back(reg_name);
        
        text = GetAttributeValueS_(reg_elem, "location");
        std::vector<double> loc = MakeCoordinates_(text);
        out_list.sublist(reg_name).sublist("Region: Plane").set<Teuchos::Array<double> >("Location", loc);
        
        text = GetAttributeValueS_(reg_elem, "normal");
        std::vector<double> dir = MakeCoordinates_(text);
        out_list.sublist(reg_name).sublist("Region: Plane").set<Teuchos::Array<double> >("Direction", dir);
      }

      else if (strcmp(node_name,"region_file") == 0) {
        tree_["regions"].push_back(reg_name);

        Teuchos::ParameterList rfPL;
        text = GetAttributeValueS_(reg_elem, "name");
        rfPL.set<std::string>("File", text);
        
        text = GetAttributeValueS_(reg_elem, "type");
        if (strcmp(text.c_str(), "color") == 0) {
          int value = GetAttributeValueD_(reg_elem, "label");
          rfPL.set<int>("Value", value);
          out_list.sublist(reg_name).sublist("Region: Color Function") = rfPL;
        }
        else if (strcmp(text.c_str(), "labeled set") == 0) {
          std::string value = GetAttributeValueS_(reg_elem, "label");
          rfPL.set<std::string>("Label", value);
          
          value = GetAttributeValueS_(reg_elem, "format");
          if (strcmp(value.c_str(), "exodus ii") == 0) {
            rfPL.set<std::string>("Format", "Exodus II");
          }
          
          value = GetAttributeValueS_(reg_elem, "entity");
          rfPL.set<std::string>("Entity", value);
          
          out_list.sublist(reg_name).sublist("Region: Labeled Set") = rfPL;
        }
      }

      else if (strcmp(node_name, "point") == 0) {
        tree_["regions"].push_back(reg_name);
        std::vector<double> coord = GetAttributeVector_(reg_elem, "coordinate");
        out_list.sublist(reg_name).sublist("Region: Point").set<Teuchos::Array<double> >("Coordinate", coord);
      }

      else if (strcmp(node_name,"polygonal_surface") == 0) {
        tree_["regions"].push_back(reg_name);
        
        std::vector<double> point, points;
        DOMNodeList* point_list = reg_elem->getElementsByTagName(mm.transcode("point"));
        int num_points = point_list->getLength();

        for (int j = 0; j < num_points; j++) {
          DOMNode* jnode = point_list->item(j);
          if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
            point = MakeCoordinates_(mm.transcode(jnode->getTextContent()));
            points.insert(points.end(), point.begin(), point.end());
          }
        }
        // get expert parameters
        if (reg_elem->hasAttribute(XMLString::transcode("tolerance"))) {
          text_content2 = mm.transcode(reg_elem->getAttribute(mm.transcode("tolerance")));
          out_list.sublist(reg_name).sublist("Region: Polygon").sublist("Expert Parameters")
              .set<double>("Tolerance", std::strtod(text_content2, NULL));
        }
        out_list.sublist(reg_name).sublist("Region: Polygon")
            .set<Teuchos::Array<double> >("Points", points)
            .set<int>("Number of points", num_points);
      }

      else if (strcmp(node_name, "logical") == 0) {
        tree_["regions"].push_back(reg_name);
       
        bool haveOp(false), haveRL(false);
        Teuchos::Array<Teuchos::Array<double> > points;
        DOMNodeList* gkids = reg_elem->getChildNodes();

        for (int j = 0; j < gkids->getLength(); j++) {
          DOMNode* jnode = gkids->item(j);
          if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
            node_name = mm.transcode(jnode->getNodeName());
            // deal with operation
            if (strcmp(node_name, "operation") == 0) {
              text_content2 = mm.transcode(jnode->getTextContent());
              if (strcmp(text_content2,"union") == 0) {
                out_list.sublist(reg_name).sublist("Region: Logical").set<std::string>("Operation", "Union");
              }
              else if (strcmp(text_content2,"intersection") == 0) {
                out_list.sublist(reg_name).sublist("Region: Logical").set<std::string>("Operation", "Intersection");
              }
              else if (strcmp(text_content2,"subtraction") == 0) {
                out_list.sublist(reg_name).sublist("Region: Logical").set<std::string>("Operation","Subtraction");
              }
              else if (strcmp(text_content2,"complement") == 0) {
                out_list.sublist(reg_name).sublist("Region: Logical").set<std::string>("Operation", "Complement");
              }
              else {
                ThrowErrorIllformed_("regions", "element", "operation", "union, intersection, subtraction, or complement");
              }
              haveOp = true;
            }
            // deal with region list
            else if (strcmp(node_name, "region_list") == 0) {
              text_content2 = mm.transcode(jnode->getTextContent());
              Teuchos::Array<std::string> regs = CharToStrings_(text_content2);
              out_list.sublist(reg_name).sublist("Region: Logical").set<Teuchos::Array<std::string> >("Regions", regs);
              haveRL = true;
            }
          }
        }
        if (!haveOp) {
          ThrowErrorMissattr_("Regions", "element", "operation", "logical");
        }
        if (!haveRL) {
          ThrowErrorMissattr_("Regions", "element", "region_list", "logical");
        }
      }

      else if (strcmp(node_name, "boundary") == 0) {
        tree_["regions"].push_back(reg_name);
        std::string type = GetAttributeValueS_(reg_elem, "entity");
        out_list.sublist(reg_name).sublist("Region: Boundary").set<std::string>("entity", type);
      }
    }
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "found " << tree_["regions"].size() << " regions plus region \"All\"." << std::endl;
  }

  out_list.sublist("All") = CreateRegionAll_();
  return out_list;
}


/* ******************************************************************
* Create global region.
****************************************************************** */
Teuchos::ParameterList InputConverterU::CreateRegionAll_()
{
  Teuchos::ParameterList out_list;
  Teuchos::ParameterList& all = out_list.sublist("Region: Box");

  std::vector<double> low(2, -1e99), high(2, 1e99);

  if (dim_ == 3) {
    low.push_back(-1e99);
    high.push_back(1e99);
  }

  all.set<Teuchos::Array<double> >("Low Coordinate", low);
  all.set<Teuchos::Array<double> >("High Coordinate", high);

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
