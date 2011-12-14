/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   RegionFactory.cc
 * @author Rao Garimella
 * @date   
 * 
 * @brief  
 * 
 * 
 */

#include "Region.hh"
#include "BoxRegion.hh"
#include "PlaneRegion.hh"
#include "LabeledSetRegion.hh"
#include "PointRegion.hh"
#include "RegionFactory.hh"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


// Create region from XML specification
// NOTE: Can we get the region name also from reg_spec

Amanzi::AmanziGeometry::RegionPtr 
Amanzi::AmanziGeometry::RegionFactory(const std::string reg_name,
                                      const unsigned int reg_id, 
                                      const Teuchos::ParameterList& reg_params)
{
  
  // There should be only one item below the region name
  // which indicates the shape of the
  // region. Unfortunately, there is nothing to prevent
  // there from being multiple sublists specifying the shape
  
  // The right way may be to have a keyword for SHAPE in the XML spec
  
  Teuchos::ParameterList::ConstIterator k = reg_params.begin();

  std::string shape = reg_params.name(k);

  if (shape == "Region: Box") 
    {

      Teuchos::ParameterList box_params = reg_params.sublist(shape);

      Teuchos::Array<double> p0_vec = box_params.get< Teuchos::Array<double> >("Low Coordinate");
        
      Teuchos::Array<double> p1_vec = box_params.get< Teuchos::Array<double> >("High Coordinate");

      Point p0, p1;

      if (p0_vec.size() == 3) 
        {
          p0.init(3);
          p0.set(p0_vec[0], p0_vec[1], p0_vec[2]);
          p1.init(3);
          p1.set(p1_vec[0], p1_vec[1], p1_vec[2]);
        }
      else if (p0_vec.size() == 2)
        {
          p0.init(2);
          p0.set(p0_vec[0], p0_vec[1]);
          p1.init(2);
          p1.set(p1_vec[0], p1_vec[1]);          
        }

      RegionPtr regptr = new BoxRegion(reg_name, reg_id, p0, p1);
      return regptr;
    }
  else if (shape == "Region: Plane")
    {
      Teuchos::ParameterList plane_params = reg_params.sublist(shape);

      Teuchos::Array<double> p_vec = plane_params.get< Teuchos::Array<double> >("Location");
        
      Teuchos::Array<double> n_vec = plane_params.get< Teuchos::Array<double> >("Direction");

      Point p, n;

      if (p_vec.size() == 3) 
        {
          p.init(3);
          p.set(p_vec[0], p_vec[1], p_vec[2]);
          n.init(3);
          n.set(n_vec[0], n_vec[1], n_vec[2]);
        }
      else if (p_vec.size() == 2)
        {
          p.init(2);
          p.set(p_vec[0], p_vec[1]);
          n.init(2);
          n.set(n_vec[0], n_vec[1]);          
        }

      RegionPtr regptr = new PlaneRegion(reg_name, reg_id, p, n);
      return regptr;
    }
  else if (shape == "Region: Labeled Set")
    {
      Teuchos::ParameterList labeledset_params = reg_params.sublist(shape);

      std::string file = labeledset_params.get<std::string>("File");
      std::string format = labeledset_params.get<std::string>("Format");
      std::string name = labeledset_params.get<std::string>("Label");
      std::string entity_str = labeledset_params.get<std::string>("Entity");

      if (entity_str == "Cell" || entity_str == "cell" || entity_str == "CELL")
        entity_str = "CELL";
      else if (entity_str == "Face" || entity_str == "face" || entity_str == "FACE")
        entity_str = "FACE";
      else if (entity_str == "Node" || entity_str == "node" || entity_str == "NODE")
        entity_str = "NODE";

      RegionPtr regptr = new LabeledSetRegion(reg_name, reg_id, entity_str, file, format, name);
      return regptr;
    }
  else if (shape == "Region: Point")
    {
      Teuchos::ParameterList point_params = reg_params.sublist(shape);

      Teuchos::Array<double> p_vec = point_params.get< Teuchos::Array<double> >("Coordinate");
        
      Point pnt, n;

      if (p_vec.size() == 3) 
        {
          pnt.init(3);
          pnt.set(p_vec[0], p_vec[1], p_vec[2]);
        }
      else if (p_vec.size() == 2)
        {
          pnt.init(2);
          pnt.set(p_vec[0], p_vec[1]);
        }

      RegionPtr regptr = new PointRegion(reg_name, reg_id, pnt);
      return regptr;
    }
  else 
    {
      std::cerr << "ERROR: Cannot process region with shape " << shape << std::endl;
      throw std::exception();
    }

}


