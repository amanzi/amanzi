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

  bool generate(true), read(false), all_good(false);
  std::string framework;
  std::stringstream helper;
  Errors::Message msg;
  Teuchos::ParameterList mesh_list;
    
  // read in new stuff
  node_list = doc_->getElementsByTagName(mm.transcode("mesh"));

  // read the attribute to set the framework sublist
  if (node_list->getLength() > 0) {
    node = node_list->item(0);
    element = static_cast<DOMElement*>(node);
    framework = GetAttributeValueS_(element, "framework");

    // Define global parameter dim_ = the space dimension.
    children = node->getChildNodes();
    all_good = false;

    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
        char* tagname = mm.transcode(inode->getNodeName());
        if (strcmp(tagname, "dimension") == 0) {
          char* tmp = mm.transcode(inode->getTextContent());
          if (strlen(tmp) > 0) {
            dim_ = std::strtol(tmp, NULL, 10);
            all_good = true;
          }
        }
      }
    }

    if (!all_good) {
      ThrowErrorIllformed_("mesh", "element", "dimension");
    }

    // Now we can properly parse the generate/read list.
    all_good = false;
    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
        char* tagname = mm.transcode(inode->getNodeName());   

        // A structured mesh is generated.
        if (strcmp(tagname,"generate") == 0) {
          all_good = true;
          generate = true;
          read = false;
          DOMElement* element_gen = static_cast<DOMElement*>(inode);

          node_list = element_gen->getElementsByTagName(mm.transcode("number_of_cells"));
          node = node_list->item(0);
          DOMElement* element_node = static_cast<DOMElement*>(node);
          DOMNamedNodeMap *attr_map = node->getAttributes();

          Teuchos::Array<int> ncells; 
          DOMNode* node_attr;
          char* attr_name;
          char* tmp;

          // make sure number of attributes equals dimension
          if (attr_map->getLength() == dim_) {
            // loop over attributes to get nx, ny, nz as needed
            for (int j = 0; j < attr_map->getLength(); j++) {
              node_attr = attr_map->item(j);
              attr_name = mm.transcode(node_attr->getNodeName());

              if (attr_name) {
                tmp = mm.transcode(node_attr->getNodeValue());
                if (strlen(tmp) > 0) {
                  ncells.append(std::strtol(tmp, NULL, 10));
                } else {
                  all_good = false;
                  helper << "number_of_cells " << attr_name;
                }
              } else {
               all_good = false;
               helper << "number_of_cells " << attr_name;
              }
            }
            mesh_list.set<Teuchos::Array<int> >("Number of Cells", ncells);
          } else {
            helper << "number_of_cells";
            all_good = false;
          }

          // get Box - generalize
          node_list = element_gen->getElementsByTagName(mm.transcode("box"));
          node = node_list->item(0);
          element_node = static_cast<DOMElement*>(node);

          tmp = mm.transcode(element_node->getAttribute(mm.transcode("low_coordinates")));
          if (strlen(tmp) > 0) {
            // translate to array
            std::vector<double> low = MakeCoordinates_(tmp);
            mesh_list.set<Teuchos::Array<double> >("Domain Low Coordinate", low);
            if (low.size() != dim_) {
              helper << "low_coordinates";
              all_good = false;
            }
          } else {
            helper << "low_coordinates";
            all_good = false;
          }

          tmp = mm.transcode(element_node->getAttribute(mm.transcode("high_coordinates")));
          if (strlen(tmp) > 0) {
            // translate to array
            std::vector<double> high = MakeCoordinates_(tmp);
            mesh_list.set<Teuchos::Array<double> >("Domain High Coordinate", high);
            if (high.size() != dim_) {
              helper << "high_coordinates";
              all_good = false;
            }
          } else {
            helper << "high_coordinates";
            all_good = false;
          }
        }

        // Un unstructured mesh will be read from a file.
        else if (strcmp(tagname, "read") == 0) {
          read = true;
          generate = false;
          bool flag1, flag2;

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
          all_good = flag1 && flag2;
        }
      }
    }

    if (!all_good) {
      ThrowErrorIllformed_("mesh", helper.str(), "generate/read");
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
    }

    if (generate) {
      out_list.sublist("Unstructured").sublist("Generate Mesh") = mesh_list;
    } else if (read) {
      out_list.sublist("Unstructured").sublist("Read Mesh File") = mesh_list;
    } else {
      msg << "Amanzi::InputConverter: an error occurred during parsing mesh.\n";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    ThrowErrorIllformed_("mesh", "element", "framework");
  }

  return out_list;
}


/* ******************************************************************
* Convert regions.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateRegions_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating regions" << std::endl;
  }

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
        
        // if attribute 'num_points' exists, get it
        int num_points(-1);
        int pt_cnt(0);
        if (reg_elem->hasAttribute(mm.transcode("num_points"))) {
          text_content2 = mm.transcode(reg_elem->getAttribute(mm.transcode("num_points")));
          std::string str(text_content2);
          boost::algorithm::trim(str);
          num_points = std::strtol(text_content2, NULL, 10);
          out_list.sublist(reg_name).sublist("Region: Polygon").set<int>("Number of points", num_points);
        }
        // get verticies (add count them)
        std::vector<double> points;
        DOMNodeList* gkids = reg_elem->getChildNodes();
        for (int j = 0; j < gkids->getLength(); j++) {
          DOMNode* jnode = gkids->item(j);
          if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
            node_name = mm.transcode(jnode->getNodeName());
            if (strcmp(node_name, "point") == 0) {
              text_content2 = mm.transcode(jnode->getTextContent());
              std::vector<double> point = MakeCoordinates_(text_content2);
              for (std::vector<double>::iterator pt = point.begin(); pt != point.end(); ++pt) {
                points.push_back(*pt);
              }
              pt_cnt++;
            }
          }
        }
        out_list.sublist(reg_name).sublist("Region: Polygon").set<Teuchos::Array<double> >("Points", points);
        if (!out_list.sublist(reg_name).sublist("Region: Polygon").isParameter("Number of points")) {
          out_list.sublist(reg_name).sublist("Region: Polygon").set<int>("Number of points", pt_cnt);
        }
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
    }
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
