/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Erin Barker (original version)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Input Converter

*/

#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <climits>

// TPLs
#define BOOST_FILESYTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/format.hpp"

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
Teuchos::ParameterList
InputConverterU::TranslateMesh_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating unstructured mesh" << std::endl;
  }

  MemoryManager mm;
  DOMNodeList *node_list, *children;
  DOMNode* node;

  bool flag, read(false), generate(false);
  std::string framework, verify;
  std::string partitioner = "NOT_SPECIFIED";
  Errors::Message msg;
  Teuchos::ParameterList mesh_list;

  // read in new stuff
  node_list = doc_->getElementsByTagName(mm.transcode("mesh"));
  if (node_list->getLength() == 0) ThrowErrorIllformed_("mesh", "element", "framework");

  DOMNode* inode = node_list->item(0);
  framework = GetAttributeValueS_(inode, "framework");

  // Define global parameter dim_ = the space dimension.
  node = GetUniqueElementByTagsString_(inode, "dimension", flag);
  if (flag) {
    char* tmp = mm.transcode(node->getTextContent());
    dim_ = std::strtol(tmp, NULL, 10);
  }

  if (!flag || dim_ <= 0) ThrowErrorIllformed_("mesh", "dimension", "dimension");

  // Define the parameter partitioner_
  node = GetUniqueElementByTagsString_(inode, "partitioner", flag);
  if (flag) { partitioner = mm.transcode(node->getTextContent()); }

  // Now we can properly parse the generate/read list.
  children = inode->getChildNodes();

  for (int i = 0; i < children->getLength(); i++) {
    inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      char* tagname = mm.transcode(inode->getNodeName());

      // A structured mesh is generated.
      if (strcmp(tagname, "generate") == 0) {
        generate = true;
        mesh_rectangular_ = "true";
        node = GetUniqueElementByTagsString_(inode, "number_of_cells", flag);
        if (!flag) ThrowErrorIllformed_("mesh", "number_of_cells", "generate");

        std::vector<int> ncells;
        int nx = GetAttributeValueL_(node, "nx");
        if (nx > 0) ncells.push_back(nx);
        int ny = GetAttributeValueL_(node, "ny", TYPE_NUMERICAL, 0, INT_MAX, false, 0);
        if (ny > 0) ncells.push_back(ny);
        int nz = GetAttributeValueL_(node, "nz", TYPE_NUMERICAL, 0, INT_MAX, false, 0);
        if (nz > 0) ncells.push_back(nz);

        if (ncells.size() != dim_) ThrowErrorIllformed_("mesh", "number_of_cells", "generate");

        // get Box
        node = GetUniqueElementByTagsString_(inode, "box", flag);
        if (!flag) ThrowErrorIllformed_("mesh", "box", "generate");

        std::vector<double> low = GetAttributeVectorD_(node, "low_coordinates", dim_, "m");
        if (low.size() != dim_) ThrowErrorIllformed_("mesh", "low_coordinates", "generate");

        std::vector<double> high = GetAttributeVectorD_(node, "high_coordinates", dim_, "m");
        if (high.size() != dim_) ThrowErrorIllformed_("mesh", "high_coordinates", "generate");

        mesh_list.set<Teuchos::Array<int>>("number of cells", ncells);
        mesh_list.set<Teuchos::Array<double>>("domain low coordinate", low);
        mesh_list.set<Teuchos::Array<double>>("domain high coordinate", high);
      }

      // Un unstructured mesh will be read from a file.
      else if (strcmp(tagname, "read") == 0) {
        bool flag1, flag2, flag3;

        node = GetUniqueElementByTagsString_(inode, "format", flag1);
        if (flag1) {
          std::string format =
            GetTextContentS_(node, "exodus ii, exodus II, Exodus II, Exodus ii, H5M, h5m");

          if (format == "exodus ii") {
            mesh_list.set<std::string>("format", "Exodus II");
          } else if (format == "h5m") {
            mesh_list.set<std::string>("format", "H5M");
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
            mesh_list.set<std::string>("file", filename);
          }
        }
        node = GetUniqueElementByTagsString_(inode, "verify", flag3);
        if (flag3) { verify = mm.transcode(node->getTextContent()); }
        read = flag1 && flag2;
      }
    }
  }

  if (generate || read) {
    Teuchos::ParameterList& tmp_list = out_list.sublist("unstructured").sublist("expert");
    if (strcmp(framework.c_str(), "mstk") == 0) {
      tmp_list.set<std::string>("framework", "MSTK");
    } else if (strcmp(framework.c_str(), "moab") == 0) {
      tmp_list.set<std::string>("framework", "MOAB");
    } else if (strcmp(framework.c_str(), "simple") == 0) {
      tmp_list.set<std::string>("framework", "Simple");
    } else if (strcmp(framework.c_str(), "stk::mesh") == 0) {
      tmp_list.set<std::string>("framework", "stk::mesh");
    } else {
      msg << "Amanzi::InputConverter: an error occurred during parsing mesh.\n"
          << "  Unknown framework \"" << framework << "\".\n";
      Exceptions::amanzi_throw(msg);
    }
    if (strcmp(verify.c_str(), "true") == 0) {
      tmp_list.set<bool>("verify mesh", (strcmp(verify.c_str(), "true") == 0));
    }
    if (partitioner != "") tmp_list.set<std::string>("partitioner", partitioner);
  }

  if (generate) {
    out_list.sublist("unstructured").sublist("generate mesh") = mesh_list;
  } else if (read) {
    out_list.sublist("unstructured").sublist("read mesh file") = mesh_list;
  } else {
    msg << "Amanzi::InputConverter: an error occurred during parsing mesh.\n";
    Exceptions::amanzi_throw(msg);
  }

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");
  return out_list;
}


/* ******************************************************************
* Convert regions.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateRegions_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) *vo_->os() << "Translating regions" << std::endl;

  MemoryManager mm;

  DOMNodeList* node_list;
  DOMNode* node;

  char *tagname, *node_name;
  char* text_content2;
  std::string reg_name, text;

  // get regions node
  node_list = doc_->getElementsByTagName(mm.transcode("regions"));
  node = node_list->item(0);

  // new options: comment, region, box, point
  DOMNode* inode = node->getFirstChild();
  while (inode != NULL) {
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      tagname = mm.transcode(inode->getNodeName());
      bool have_name(false);

      DOMElement* reg_elem;
      // If region is under a region tag, get the name and set
      // child element as region element.
      if (strcmp(tagname, "region") == 0) {
        reg_name = GetAttributeValueS_(static_cast<DOMElement*>(inode), "name");
        have_name = true;

        // loop over children to get region element
        DOMNodeList* kids = inode->getChildNodes();
        for (int j = 0; j < kids->getLength(); j++) {
          DOMNode* jnode = kids->item(j);

          if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
            node_name = mm.transcode(jnode->getNodeName());
            if (strcmp(node_name, "comments") != 0) reg_elem = static_cast<DOMElement*>(jnode);
          }
        }
      } else {
        // else set the current element as region element
        reg_elem = static_cast<DOMElement*>(inode);
      }

      // get reg_elem type
      node_name = mm.transcode(reg_elem->getNodeName());
      region_type_[reg_name] = 0;

      // get name if needed
      if (!have_name) {
        if (strcmp(node_name, "comments") != 0) {
          reg_name = GetAttributeValueS_(reg_elem, "name");
        }
      }

      // loop over attributes of region element
      if (strcmp(node_name, "box") == 0) {
        tree_["regions"].push_back(reg_name);

        std::vector<double> low = GetAttributeVectorD_(reg_elem, "low_coordinates", dim_, "m");
        std::vector<double> high = GetAttributeVectorD_(reg_elem, "high_coordinates", dim_, "m");
        out_list.sublist(reg_name)
          .sublist("region: box")
          .set<Teuchos::Array<double>>("low coordinate", low)
          .set<Teuchos::Array<double>>("high coordinate", high);
      }

      else if (strcmp(node_name, "plane") == 0) {
        tree_["regions"].push_back(reg_name);

        std::vector<double> loc = GetAttributeVectorD_(reg_elem, "location", dim_, "m");
        std::vector<double> dir = GetAttributeVectorD_(reg_elem, "normal", dim_, "m");

        out_list.sublist(reg_name)
          .sublist("region: plane")
          .set<Teuchos::Array<double>>("point", loc)
          .set<Teuchos::Array<double>>("normal", dir);
      }

      else if (strcmp(node_name, "cylinder") == 0) {
        tree_["regions"].push_back(reg_name);

        std::vector<double> loc = GetAttributeVectorD_(reg_elem, "location", dim_, "m");
        std::vector<double> dir = GetAttributeVectorD_(reg_elem, "axis", dim_, "m");
        double rad = GetAttributeValueD_(reg_elem, "radius", TYPE_NUMERICAL, 0.0, DVAL_MAX, "-");

        out_list.sublist(reg_name)
          .sublist("region: cylinder")
          .set<Teuchos::Array<double>>("point", loc)
          .set<Teuchos::Array<double>>("axis", dir)
          .set<double>("radius", rad);
      }

      else if (strcmp(node_name, "halfspace") == 0) {
        tree_["regions"].push_back(reg_name);

        std::vector<double> loc = GetAttributeVectorD_(reg_elem, "location", dim_, "m");
        std::vector<double> dir = GetAttributeVectorD_(reg_elem, "normal", dim_, "m");

        out_list.sublist(reg_name)
          .sublist("region: halfspace")
          .set<Teuchos::Array<double>>("point", loc)
          .set<Teuchos::Array<double>>("normal", dir);
      }

      else if (strcmp(node_name, "region_file") == 0) {
        tree_["regions"].push_back(reg_name);

        Teuchos::ParameterList rfPL;
        text = GetAttributeValueS_(reg_elem, "name");
        rfPL.set<std::string>("file", text);

        text = GetAttributeValueS_(reg_elem, "type");
        if (strcmp(text.c_str(), "color") == 0) {
          int value = GetAttributeValueL_(reg_elem, "label");
          rfPL.set<int>("value", value);
          out_list.sublist(reg_name).sublist("region: color function") = rfPL;
        } else if (strcmp(text.c_str(), "labeled set") == 0) {
          std::string value = GetAttributeValueS_(reg_elem, "label");
          rfPL.set<std::string>("label", value);

          value = GetAttributeValueS_(reg_elem, "format");
          if (strcmp(value.c_str(), "exodus ii") == 0) {
            rfPL.set<std::string>("format", "Exodus II");
          }

          value = GetAttributeValueS_(reg_elem, "entity");
          rfPL.set<std::string>("entity", value);

          out_list.sublist(reg_name).sublist("region: labeled set") = rfPL;
        }
      }

      else if (strcmp(node_name, "point") == 0) {
        tree_["regions"].push_back(reg_name);
        std::vector<double> coord = GetAttributeVectorD_(reg_elem, "coordinate", dim_, "m");
        out_list.sublist(reg_name)
          .sublist("region: point")
          .set<Teuchos::Array<double>>("coordinate", coord);
      }

      else if (strcmp(node_name, "polygonal_surface") == 0) {
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
        if (reg_elem->hasAttribute(mm.transcode("tolerance"))) {
          text_content2 = mm.transcode(reg_elem->getAttribute(mm.transcode("tolerance")));
          out_list.sublist(reg_name)
            .sublist("region: polygon")
            .sublist("expert parameters")
            .set<double>("tolerance", std::strtod(text_content2, NULL));
        }
        out_list.sublist(reg_name)
          .sublist("region: polygon")
          .set<Teuchos::Array<double>>("points", points)
          .set<int>("number of points", num_points);
      }

      else if (strcmp(node_name, "logical") == 0) {
        tree_["regions"].push_back(reg_name);

        bool haveOp(false), haveRL(false);
        Teuchos::Array<Teuchos::Array<double>> points;

        DOMNodeList* kids = inode->getChildNodes();
        for (int j = 0; j < kids->getLength(); j++) {
          DOMNode* jnode = kids->item(j);

          if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
            char* elem_name = mm.transcode(jnode->getNodeName());
            // deal with operation
            if (strcmp(elem_name, "operation") == 0) {
              text_content2 = mm.transcode(jnode->getTextContent());
              if (strcmp(text_content2, "union") == 0) {
                out_list.sublist(reg_name)
                  .sublist("region: logical")
                  .set<std::string>("operation", "union");
              } else if (strcmp(text_content2, "intersection") == 0) {
                out_list.sublist(reg_name)
                  .sublist("region: logical")
                  .set<std::string>("operation", "intersect");
              } else if (strcmp(text_content2, "subtraction") == 0) {
                out_list.sublist(reg_name)
                  .sublist("region: logical")
                  .set<std::string>("operation", "subtract");
              } else if (strcmp(text_content2, "complement") == 0) {
                out_list.sublist(reg_name)
                  .sublist("region: logical")
                  .set<std::string>("operation", "complement");
              } else {
                ThrowErrorIllformed_(
                  "regions", "element", "operation", "union, intersect, subtract, or complement");
              }
              haveOp = true;
            }
            // deal with region list
            else if (strcmp(elem_name, "region_list") == 0) {
              text_content2 = mm.transcode(jnode->getTextContent());
              Teuchos::Array<std::string> regs = CharToStrings_(text_content2);
              out_list.sublist(reg_name)
                .sublist("region: logical")
                .set<Teuchos::Array<std::string>>("regions", regs);
              haveRL = true;
            }
          }
        }

        if (!haveOp) { ThrowErrorMissing_("regions", "element", "operation", "logical"); }
        if (!haveRL) { ThrowErrorMissing_("regions", "element", "region_list", "logical"); }
      }

      else if (strcmp(node_name, "boundary") == 0) {
        tree_["regions"].push_back(reg_name);
        std::string type = GetAttributeValueS_(reg_elem, "entity");
        out_list.sublist(reg_name).sublist("region: boundary").set<std::string>("entity", type);
      }

      else if (strcmp(node_name, "box_volume_fractions") == 0) {
        tree_["regions"].push_back(reg_name);
        region_type_[reg_name] = 1;

        std::vector<double> low = GetAttributeVectorD_(reg_elem, "corner_coordinates", dim_, "m");
        std::vector<double> high =
          GetAttributeVectorD_(reg_elem, "opposite_corner_coordinates", dim_, "m");
        std::vector<double> normals = GetAttributeVectorD_(reg_elem, "normals", dim_, "", false);

        out_list.sublist(reg_name)
          .sublist("region: box volume fractions")
          .set<Teuchos::Array<double>>("corner coordinate", low)
          .set<Teuchos::Array<double>>("opposite corner coordinate", high);

        if (normals.size() > 0)
          out_list.sublist(reg_name)
            .sublist("region: box volume fractions")
            .set<Teuchos::Array<double>>("normals", normals);
      } else if (strcmp(node_name, "line_segment") == 0) {
        tree_["regions"].push_back(reg_name);
        std::vector<double> p1 = GetAttributeVectorD_(reg_elem, "end_coordinates", dim_, "m");
        std::vector<double> p2 =
          GetAttributeVectorD_(reg_elem, "opposite_end_coordinates", dim_, "m");
        out_list.sublist(reg_name)
          .sublist("region: line segment")
          .set<Teuchos::Array<double>>("end coordinate", p1)
          .set<Teuchos::Array<double>>("opposite end coordinate", p2);
      }
    }

    inode = inode->getNextSibling();
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "found " << tree_["regions"].size() << " regions plus region \"All\"."
               << std::endl;
  }

  out_list.sublist("All") = CreateRegionAll_();
  return out_list;
}


/* ******************************************************************
* Create global region.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::CreateRegionAll_()
{
  Teuchos::ParameterList out_list;
  Teuchos::ParameterList& all = out_list.sublist("region: box");

  std::vector<double> low(2, -1e99), high(2, 1e99);

  if (dim_ == 3) {
    low.push_back(-1e99);
    high.push_back(1e99);
  }

  all.set<Teuchos::Array<double>>("low coordinate", low);
  all.set<Teuchos::Array<double>>("high coordinate", high);

  return out_list;
}

} // namespace AmanziInput
} // namespace Amanzi
