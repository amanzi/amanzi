/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#define  BOOST_FILESYTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"

// TPLs
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/DOMLSParserImpl.hpp>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>

#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Amanzi's
#include "ErrorHandler.hpp"
#include "InputConverterU.hh"


namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

const Teuchos::Array<std::string> meshfile_strings = 
   Teuchos::tuple<std::string>("exodus ii", "exodus II", "Exodus II", "Exodus ii", "H5M", "h5m");

/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputConverterU::Translate()
{
  Teuchos::ParameterList out_list;
  
  out_list.sublist("Mesh") = TranslateMesh_();
  out_list.sublist("Domain").set<int>("Spatial Dimension", dim_);
  out_list.sublist("Regions") = TranslateRegions_();
  
  return out_list;
}


/* ******************************************************************
* Translate unstructured mesh. Introduces global parameter dim_.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateMesh_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating unstructured mesh" << std::endl;
  }

  bool generate(true), read(false), all_good(false);
  char *framework;
  std::stringstream helper;
  Errors::Message msg;
  Teuchos::ParameterList mesh_list;
    
  Teuchos::RCP<Teuchos::StringValidator> meshfile_validator = Teuchos::rcp(new Teuchos::StringValidator(meshfile_strings));
  
  // read in new stuff
  XMLCh* tag = XMLString::transcode("mesh");
  DOMNodeList* node_list = doc_->getElementsByTagName(tag);
  XMLString::release(&tag);

  // read the attribute to set the framework sublist
  if (node_list->getLength() > 0) {
    DOMNode* node_mesh = node_list->item(0);
    DOMElement* element_mesh = static_cast<DOMElement*>(node_mesh);

    if (element_mesh->hasAttribute(XMLString::transcode("framework"))) {
      framework = XMLString::transcode(element_mesh->getAttribute(XMLString::transcode("framework")));
    } else { 
      msg << "Amanzi::InputConverter: An error occurred during parsing mesh.\n"
          << "Framework was missing or ill-formed.\n" 
          << "Use default framework='mstk' if unsure.\n";
      Exceptions::amanzi_throw(msg);
    }

    // Define global parameter dim_ = the space dimension.
    DOMNodeList* children = node_mesh->getChildNodes();
    all_good = false;

    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
	char* tagname = XMLString::transcode(inode->getNodeName());
	if (strcmp(tagname, "dimension") == 0) {
	  char* tmp = XMLString::transcode(inode->getTextContent());
	  if (strlen(tmp) > 0) {
	      dim_ = atoi(tmp);
	      all_good = true;
	  }
          XMLString::release(&tmp);
	}
        XMLString::release(&tagname);
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
	char* tagname = XMLString::transcode(inode->getNodeName());   

        // A structured mesh is generated.
	if (strcmp(tagname,"generate") == 0) {
	  all_good = true;
	  generate = true;
	  read = false;
	  DOMElement* element_gen = static_cast<DOMElement*>(inode);

	  // get Number of Cells
	  DOMNodeList* node_list = element_gen->getElementsByTagName( XMLString::transcode("number_of_cells"));
	  DOMNode* node = node_list->item(0);
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
	      attr_name =XMLString::transcode(node_attr->getNodeName());

	      if (attr_name) {
	        tmp = XMLString::transcode(node_attr->getNodeValue());
	        if (strlen(tmp) > 0) {
	          ncells.append(atoi(tmp));
		} else {
		  all_good = false;
	          helper << "number_of_cells " << attr_name;
	        }
                XMLString::release(&tmp);
	      } else {
		all_good = false;
	        helper << "number_of_cells " << attr_name;
	      }
	      XMLString::release(&attr_name);
	    }
            mesh_list.set<Teuchos::Array<int> >("Number of Cells", ncells);
	  } else {
	    helper << "number_of_cells";
	    all_good = false;
	  }

	  // get Box - generalize
	  node_list = element_gen->getElementsByTagName(XMLString::transcode("box"));
	  node = node_list->item(0);
	  element_node = static_cast<DOMElement*>(node);
	  tmp = XMLString::transcode(element_node->getAttribute(XMLString::transcode("low_coordinates")));

	  if (strlen(tmp) > 0) {
	    // translate to array
	    Teuchos::Array<double> low = MakeCoordinates_(tmp);
            mesh_list.set<Teuchos::Array<double> >("Domain Low Coordinate", low);
	    if (low.length() != dim_) {
	      helper << "low_coordinates";
	      all_good = false;
	    }
	  } else {
	    helper << "low_coordinates";
	    all_good = false;
	  }
	  XMLString::release(&tmp);

	  tmp = XMLString::transcode(element_node->getAttribute( XMLString::transcode("high_coordinates")));
	  if (strlen(tmp) > 0) {
	    // translate to array
	    Teuchos::Array<double> high = MakeCoordinates_(tmp);
            mesh_list.set<Teuchos::Array<double> >("Domain High Coordinate", high);
	    if (high.length() != dim_) {
	      helper << "high_coordinates";
	      all_good = false;
	    }
	  }
          else {
	    helper << "high_coordinates";
	    all_good = false;
	  }
	  XMLString::release(&tmp);
	}

        // Un unstructured mesh will be read from a file.
	else if (strcmp(tagname, "read") == 0) {
	  read = true;
	  generate = false;
	  bool goodtype = false;
	  bool goodname = false;
	  DOMElement* elementRead = static_cast<DOMElement*>(inode);

	  char* value = XMLString::transcode(elementRead->getElementsByTagName(
              XMLString::transcode("format"))->item(0)->getTextContent());
          std::string format(TrimString_(value));

          if (boost::iequals(format, "exodus ii")) {
            mesh_list.set<std::string>("Format", "Exodus II");
            goodtype = true;
	  } else if (boost::iequals(format, "h5m")) {
            mesh_list.set<std::string>("Format", "H5M");
            goodtype = true;
          } else {
            mesh_list.set<std::string>("Format", format, "Format of meshfile", meshfile_validator);
          }

	  char* filename = XMLString::transcode(elementRead->getElementsByTagName(
              XMLString::transcode("file"))->item(0)->getTextContent());
	  if (strlen(filename) > 0) {
            mesh_list.set<std::string>("File", TrimString_(filename));
	    goodname = true;
	  }
	  XMLString::release(&value);
	  XMLString::release(&filename);
	  if (goodtype && goodname) all_good = true;
	}
      }
    }

    if (!all_good) {
      ThrowErrorIllformed_("mesh", helper.str(), "generate/read");
    }
    
    if (generate || read) {
      if (strcmp(framework,"mstk") == 0) {
        out_list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","MSTK");
      } else if (strcmp(framework,"moab") == 0) {
        out_list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","MOAB");
      } else if (strcmp(framework,"simple") == 0) {
        out_list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","Simple");
      } else if (strcmp(framework,"stk::mesh") == 0) {
        out_list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","stk::mesh");
      } else {
        msg << "Amanzi::InputConverter: An error occurred during parsing mesh - "
            << "unknown framework=" << framework << ". See the schema for acceptable types.\n";
        Exceptions::amanzi_throw(msg); 
      }
    }

    if (generate) {
      out_list.sublist("Unstructured").sublist("Generate Mesh").sublist("Uniform Structured") = mesh_list;
    } else if (read) {
      out_list.sublist("Unstructured").sublist("Read Mesh File") = mesh_list;
    } else {
      // bad mesh, again if validated shouldn't need this
      msg << "Amanzi::InputConverter: An error occurred during parsing mesh.\n";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    // no mesh sub-elements
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

  DOMNodeList* node_list;
  DOMNode* nodeTmp;
  DOMNode* node_attr;
  DOMNamedNodeMap* attr_map;
  char* tag_name;
  char* node_name;
  char* reg_name;
  char* text_content;
  char* text_content2;
  char* char_array;

  // get regions node
  node_list = doc_->getElementsByTagName(XMLString::transcode("regions"));
  DOMNode* node_rgn = node_list->item(0);
  DOMElement* elementRgn = static_cast<DOMElement*>(node_rgn);

  // new options: comment, region, box, point
  DOMNodeList* childern = node_rgn->getChildNodes();
  for (int i = 0; i < childern->getLength(); i++) {
    DOMNode* inode = childern->item(i);

    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      tag_name = XMLString::transcode(inode->getNodeName());
      bool have_name = false;
      
      DOMElement* reg_elem;
      // If region is under a region tag, get get the name
      // and set child element as region element.
      if (strcmp(tag_name,"region") == 0) {
        attr_map = inode->getAttributes();
        node_attr = attr_map->getNamedItem(XMLString::transcode("name"));
        if (node_attr) {
          reg_name = XMLString::transcode(node_attr->getNodeValue());
          have_name = true;
        } else {
          ThrowErrorMissattr_("Regions", "attribute", "name", "region");
        }

        // loop over children to get region element
        DOMNodeList* kids = inode->getChildNodes();
        for (int j = 0; j < kids->getLength(); j++) {
          DOMNode* jnode = kids->item(j);

          if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
            node_name = XMLString::transcode(jnode->getNodeName());
            if (strcmp(node_name,"comments") != 0) 
                reg_elem =static_cast<DOMElement*>(jnode);
            XMLString::release(&node_name);
          }
        }
      } else {
        // else set the current element as region element
        reg_elem =static_cast<DOMElement*>(inode);
      }
      
      // get reg_elem type
      node_name = XMLString::transcode(reg_elem->getNodeName());
      
      // get name if needed
      if (!have_name) {
        if (reg_elem->hasAttribute(XMLString::transcode("name"))) {
          reg_name = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("name")));
        } else {
          if (strcmp(node_name, "comments") != 0)
              ThrowErrorMissattr_("Regions", "attribute", "name", node_name);
        }
      }
      
      // loop over attributes of region element
      if (strcmp(node_name, "box") == 0) {
        tree_["regions"].push_back(reg_name);
        
        if (reg_elem->hasAttribute(XMLString::transcode("low_coordinates"))) {
          text_content2 = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("low_coordinates")));
        } else {
          ThrowErrorMissattr_("Regions", "attribute", "low_coordinates", "box");
        }
        Teuchos::Array<double> low = MakeCoordinates_(text_content2);
        out_list.sublist(reg_name).sublist("Region: Box").set<Teuchos::Array<double> >("Low Coordinate", low);
        XMLString::release(&text_content2);
        
        if (reg_elem->hasAttribute(XMLString::transcode("high_coordinates"))) {
          text_content2 = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("high_coordinates")));
        } else {
          ThrowErrorMissattr_("Regions","attribute","high_coordinates","box");
        }
        Teuchos::Array<double> high = MakeCoordinates_(text_content2);
        out_list.sublist(reg_name).sublist("Region: Box").set<Teuchos::Array<double> >("High Coordinate",high);
        XMLString::release(&text_content2);
      }

      else if (strcmp(node_name, "plane") == 0) {
        tree_["regions"].push_back(reg_name);
        
        if (reg_elem->hasAttribute(XMLString::transcode("location"))) {
          text_content2 = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("location")));
        } else {
          ThrowErrorMissattr_("Regions","attribute","location","plane");
        }
        Teuchos::Array<double> loc = MakeCoordinates_(text_content2);
        out_list.sublist(reg_name).sublist("Region: Plane").set<Teuchos::Array<double> >("Location", loc);
        
        if (reg_elem->hasAttribute(XMLString::transcode("normal"))) {
          text_content2 = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("normal")));
        } else {
          ThrowErrorMissattr_("Regions","attribute","normal","plane");
        }
        Teuchos::Array<double> dir = MakeCoordinates_(text_content2);
        out_list.sublist(reg_name).sublist("Region: Plane").set<Teuchos::Array<double> >("Direction", dir);
        XMLString::release(&text_content2);
      }

      else if (strcmp(node_name,"region_file") == 0) {
        tree_["regions"].push_back(reg_name);

        Teuchos::ParameterList rfPL;
        
        if (reg_elem->hasAttribute(XMLString::transcode("name"))) {
          text_content2 = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("name")));
        } else {
          ThrowErrorMissattr_("Regions", "attribute", "name", "region_file");
        }
        rfPL.set<std::string>("File", TrimString_(text_content2));
        XMLString::release(&text_content2);
        
        if (reg_elem->hasAttribute(XMLString::transcode("type"))) {
          text_content2 = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("type")));
        } else {
          ThrowErrorMissattr_("Regions","attribute","type","region_file");
        }
        if  (strcmp(text_content2,"color") == 0){
          char* value;
          if (reg_elem->hasAttribute(XMLString::transcode("label"))) {
            value = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("label")));
          } else {
            ThrowErrorMissattr_("Regions","attribute","label","color");
          }
          rfPL.set<int>("Value",atoi(value));
          XMLString::release(&value);
          out_list.sublist(reg_name).sublist("Region: Color Function") = rfPL;
        }

        else if (strcmp(text_content2, "labeled set") == 0) {
          char* value;
          if (reg_elem->hasAttribute(XMLString::transcode("label"))) {
            value = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("label")));
          } else {
            ThrowErrorMissattr_("Regions","attribute","label","labeled set");
          }
          rfPL.set<std::string>("Label", TrimString_(value));
          XMLString::release(&value);
          
          if (reg_elem->hasAttribute(XMLString::transcode("format"))) {
            value = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("format")));
          } else {
            ThrowErrorMissattr_("Regions","attribute","format","labeled set");
          }
          if  (strcmp(value,"exodus ii") == 0){
            rfPL.set<std::string>("Format","Exodus II");
          }
          XMLString::release(&value);
          
          if (reg_elem->hasAttribute(XMLString::transcode("entity"))) {
            value = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("entity")));
          } else {
            ThrowErrorMissattr_("Regions", "attribute", "entity", "labeled set");
          }
          rfPL.set<std::string>("Entity", TrimString_(value));
          XMLString::release(&value);
          
          out_list.sublist(reg_name).sublist("Region: Labeled Set") = rfPL;
        }
        XMLString::release(&text_content2);
      }

      else if (strcmp(node_name,"point") == 0) {
        tree_["regions"].push_back(reg_name);
        
        if (reg_elem->hasAttribute(XMLString::transcode("coordinate"))) {
          text_content2 = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("coordinate")));
        } else {
          ThrowErrorMissattr_("Regions", "attribute", "coordinate", "point");
        }
        
        Teuchos::Array<double> coord = MakeCoordinates_(text_content2);
        out_list.sublist(reg_name).sublist("Region: Point").set<Teuchos::Array<double> >("Coordinate", coord);
        XMLString::release(&text_content2);
      }

      else if (strcmp(node_name,"polygonal_surface") == 0) {
        tree_["regions"].push_back(reg_name);
        
        // if attribute 'num_points' exists, get it
        int num_points(-1);
        int pt_cnt(0);
        if (reg_elem->hasAttribute(XMLString::transcode("num_points"))) {
          text_content2 = XMLString::transcode(reg_elem->getAttribute(XMLString::transcode("num_points")));
          std::string str(text_content2);
          boost::algorithm::trim(str);
          num_points = atoi(text_content2);
          //list.sublist(reg_name).sublist("Region: Polygonal Surface").set<int>("Number of points",num_points);
          out_list.sublist(reg_name).sublist("Region: Polygon").set<int>("Number of points",num_points);
          XMLString::release(&text_content2);
        }
        // get verticies (add count them)
        Teuchos::Array<double> points;
        DOMNodeList* gkids = reg_elem->getChildNodes();
        for (int j=0; j<gkids->getLength(); j++) {
          DOMNode* curGKid = gkids->item(j);
          if (DOMNode::ELEMENT_NODE == curGKid->getNodeType()) {
            node_name = XMLString::transcode(curGKid->getNodeName());
            if (strcmp(node_name,"point") == 0){
              text_content2 = XMLString::transcode(curGKid->getTextContent());
              Teuchos::Array<double> point = MakeCoordinates_(text_content2);
              for (Teuchos::Array<double>::iterator pt = point.begin(); pt != point.end(); ++pt) {
                points.append(*pt);
              }
              pt_cnt++;
              XMLString::release(&text_content2);
            }
          }
        }
        //list.sublist(reg_name).sublist("Region: Polygonal Surface").set<Teuchos::Array<double> >("Points",points);
        out_list.sublist(reg_name).sublist("Region: Polygon").set<Teuchos::Array<double> >("Points",points);
        // check that 'num_points' with current count
        //if (!list.sublist(reg_name).sublist("Region: Polygonal Surface").isParameter("Number of points")) {
        //  list.sublist(reg_name).sublist("Region: Polygonal Surface").set<int>("Number of points",pt_cnt);
        //}
        if (!out_list.sublist(reg_name).sublist("Region: Polygon").isParameter("Number of points")) {
          out_list.sublist(reg_name).sublist("Region: Polygon").set<int>("Number of points",pt_cnt);
        }
      }

      else if (strcmp(node_name, "logical") == 0) {
        tree_["regions"].push_back(reg_name);
       
        // Loop over childern
	bool haveOp = false;
	bool haveRL = false;
	Teuchos::Array<Teuchos::Array<double> > points;
        DOMNodeList* gkids = reg_elem->getChildNodes();

        for (int j = 0; j < gkids->getLength(); j++) {
          DOMNode* curGKid = gkids->item(j);
          if (DOMNode::ELEMENT_NODE == curGKid->getNodeType()) {
            node_name = XMLString::transcode(curGKid->getNodeName());
	    // deal with operation
            if (strcmp(node_name,"operation") == 0){
              text_content2 = XMLString::transcode(curGKid->getTextContent());
              if ( strcmp(text_content2,"union") == 0) {
                out_list.sublist(reg_name).sublist("Region: Logical").set<std::string>("Operation","Union");
              }
              else if (strcmp(text_content2,"intersection") == 0) {
                out_list.sublist(reg_name).sublist("Region: Logical").set<std::string>("Operation","Intersection");
              }
              else if (strcmp(text_content2,"subtraction") == 0) {
                out_list.sublist(reg_name).sublist("Region: Logical").set<std::string>("Operation","Subtraction");
              }
              else if (strcmp(text_content2,"complement") == 0) {
                out_list.sublist(reg_name).sublist("Region: Logical").set<std::string>("Operation","Complement");
              }
	      else {
		// error about missing or unsupported operation
		ThrowErrorIllformed_("Regions","element","operation","union, intersection, subtraction, or complement");
	      }
	      haveOp = true;
              XMLString::release(&text_content2);
            }
	    // deal with region list
	    else if (strcmp(node_name,"region_list") == 0){
              text_content2 = XMLString::transcode(curGKid->getTextContent());
	      Teuchos::Array<std::string> regs; // KLN = make_regions_list(text_content2);
              out_list.sublist(reg_name).sublist("Region: Logical").set<Teuchos::Array<std::string> >("Regions", regs);
	      haveRL = true;
              XMLString::release(&text_content2);
	    }
          }
        }
	if (!haveOp) {
	  ThrowErrorMissattr_("Regions","element","operation","logical");
	}
	if (!haveRL) {
	  ThrowErrorMissattr_("Regions","element","region_list","logical");
	}
      }
      
      XMLString::release(&node_name);
      
      XMLString::release(&tag_name);
      if (have_name) XMLString::release(&reg_name);
    }
  }

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
