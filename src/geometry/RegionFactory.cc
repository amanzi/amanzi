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

#include <iostream>
#include <sstream>

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Region.hh"
#include "BoxRegion.hh"
#include "PlaneRegion.hh"
#include "LabeledSetRegion.hh"
#include "ColorFunctionRegion.hh"
#include "PointRegion.hh"
#include "LogicalRegion.hh"
#include "PolygonRegion.hh"

#include "dbc.hh"
#include "errors.hh"

#include "RegionFactory.hh"


// Create region from XML specification
// NOTE: Can we get the region name also from reg_spec

Amanzi::AmanziGeometry::RegionPtr 
Amanzi::AmanziGeometry::RegionFactory(const std::string reg_name,
                                      const unsigned int reg_id, 
                                      const Teuchos::ParameterList& reg_params,
                                      const int space_dim,
                                      const Epetra_MpiComm *comm)
{
  std::stringstream sstream;

  // There should be only one item below the region name
  // which indicates the shape of the
  // region. Unfortunately, there is nothing to prevent
  // there from being multiple sublists specifying the shape
  
  // The right way may be to have a keyword for SHAPE in the XML spec
  
  Teuchos::ParameterList::ConstIterator k = reg_params.begin();

  std::string shape = reg_params.name(k);

  LifeCycleType lifecycle=PERMANENT;
  if (reg_params.isParameter("Lifecycle")) {

    std::string lifecycle_str = reg_params.get<std::string>("Lifecycle");

    if (lifecycle_str != "Permanent" || lifecycle_str != "Temporary") {
      std::cerr << "Lifecycle can only be Temporary or Permanent. Reset to Permanent" << std::endl;
      lifecycle = PERMANENT;
    }

  }

  if (shape == "Region: Box") 
    {

      Teuchos::ParameterList box_params = reg_params.sublist(shape);

      Teuchos::Array<double> p0_vec = box_params.get< Teuchos::Array<double> >("Low Coordinate");
        
      Teuchos::Array<double> p1_vec = box_params.get< Teuchos::Array<double> >("High Coordinate");

      int dim = p0_vec.size();
      Point p0(dim), p1(dim);

      if (dim == 3) 
        {
          // if (space_dim != 3) {
          //   sstream << "Box " << reg_name << 
          //     " specified using 3D coordinates but problem is " << 
          //     space_dim << " dimensional. Check input!" << std::endl;
          //   Errors::Message mesg(sstream.str());
          //   amanzi_throw(mesg);
          // }

          p0.set(p0_vec[0], p0_vec[1], p0_vec[2]);
          p1.set(p1_vec[0], p1_vec[1], p1_vec[2]);
        }
      else if (dim == 2)
        {
          // if (space_dim != 2) {
          //   sstream << "Box " << reg_name << 
          //     " specified using 2D coordinates but problem is " << 
          //     space_dim << " dimensional. Check input!" << std::endl;
          //   Errors::Message mesg(sstream.str());
          //   amanzi_throw(mesg);
          // }

          p0.set(p0_vec[0], p0_vec[1]);
          p1.set(p1_vec[0], p1_vec[1]);          
        }

      try {
        RegionPtr regptr = new BoxRegion(reg_name, reg_id, p0, p1, lifecycle);        

        // Verify that we have a usable box

        if (comm->MyPID() == 0) {
          int ndeg=0;
          if (((BoxRegionPtr) regptr)->is_degenerate(&ndeg) && ndeg > 1) {
            std::cerr << "\n" << "Box region \"" << reg_name << 
              "\" is degenerate in 2 or more directions" << std::endl;
            std::cerr << 
            "This means it is a line or point in 3D, or it is a point in 2D" <<
              std::endl;
            std::cerr << 
            "Can only ask for nodes (not cells or faces) on this region" << 
              std::endl << std::endl;
          }
        }

        return regptr;
      }
      catch (Errors::Message mesg) {
        mesg << "\n" << "Cannot create region of type Box";
        Exceptions::amanzi_throw(mesg);
      }
    }
  else if (shape == "Region: Plane")
    {
      Teuchos::ParameterList plane_params = reg_params.sublist(shape);

      Teuchos::Array<double> p_vec = plane_params.get< Teuchos::Array<double> >("Location");
        
      Teuchos::Array<double> n_vec = plane_params.get< Teuchos::Array<double> >("Direction");

      int dim = p_vec.size();
      Point p(dim), n(dim);

      if (dim == 3) 
        {
          if (space_dim != 3) {
            sstream << "Plane " << reg_name << 
              " specified using 3D coordinates but problem is " << 
              space_dim << " dimensional. Check input!" << std::endl;
            Errors::Message mesg(sstream.str());
            amanzi_throw(mesg);
          }

          p.set(p_vec[0], p_vec[1], p_vec[2]);
          n.set(n_vec[0], n_vec[1], n_vec[2]);
        }
      else if (dim == 2)
        {
          if (space_dim != 2) {
            sstream << "Plane " << reg_name << 
              " specified using 2D coordinates but problem is " << 
              space_dim << " dimensional. Check input!" << std::endl;
            Errors::Message mesg(sstream.str());
            amanzi_throw(mesg);
          }

          p.set(p_vec[0], p_vec[1]);
          n.set(n_vec[0], n_vec[1]);          
        }

      try {
        RegionPtr regptr = new PlaneRegion(reg_name, reg_id, p, n, lifecycle);
        return regptr;
      }
      catch (Errors::Message mesg) {
        mesg << "\n" << "Cannot create region of type Plane";
        Exceptions::amanzi_throw(mesg);
      }
    }
  else if (shape == "Region: Polygon")
    {
      Teuchos::ParameterList poly_params = reg_params.sublist(shape);

      int num_points = poly_params.get<int>("Number of points");
        
      Teuchos::Array<double> pvec = poly_params.get< Teuchos::Array<double> >("Points");

      if (pvec.size()%num_points != 0) {
        Errors::Message mesg("Incorrect number of values specified for polygon point specification");
        amanzi_throw(mesg);
      }
      int dim = pvec.size()/num_points;

      if (dim == 3 && space_dim != 3) {
        sstream << "Polygon " << reg_name << 
          " specified using 3D coordinates but problem is " << 
          space_dim << " dimensional. Check input!" << std::endl;
        Errors::Message mesg(sstream.str());
        amanzi_throw(mesg);
      }
      else if (dim == 2 && space_dim != 2) {
        sstream << "Polygon " << reg_name << 
          " specified using 3D coordinates but problem is " << 
          space_dim << " dimensional. Check input!" << std::endl;
        Errors::Message mesg(sstream.str());
        amanzi_throw(mesg);
      }



      std::vector<Point> points;
      Point pnt(dim);
      for (int i = 0; i < num_points; i++) {
        if (dim == 2)
          pnt.set(pvec[i*dim],pvec[i*dim+1]);
        else if (dim == 3)
          pnt.set(pvec[i*dim],pvec[i*dim+1],pvec[i*dim+2]);
        points.push_back(pnt);
      }

      try {
        RegionPtr regptr = new PolygonRegion(reg_name, reg_id, num_points, points, lifecycle);
        return regptr;
      }
      catch (Errors::Message mesg) {
        mesg << "\n" << "Cannot create region of type Polygon";
        Exceptions::amanzi_throw(mesg);
      }
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

      try {
        RegionPtr regptr = new LabeledSetRegion(reg_name, reg_id, entity_str, file, format, name, lifecycle);
        return regptr;
      }
      catch (Errors::Message mesg) {
        mesg << "\n" << "Cannot create region of type LabeledSet";
        Exceptions::amanzi_throw(mesg);
      }

    }
  else if (shape == "Region: Color Function")
    {
      Teuchos::ParameterList colorfunc_params = reg_params.sublist(shape);

      std::string file = colorfunc_params.get<std::string>("File");
      int value = colorfunc_params.get<int>("Value");

      try {
        RegionPtr regptr = new ColorFunctionRegion(reg_name, reg_id, file, value, comm, lifecycle);
        return regptr;
      }
      catch (Errors::Message mesg) {
        mesg << "\n" << "Cannot create region of type Color Function";
        Exceptions::amanzi_throw(mesg);
      }

    }
  else if (shape == "Region: Point")
    {
      Teuchos::ParameterList point_params = reg_params.sublist(shape);

      Teuchos::Array<double> p_vec = point_params.get< Teuchos::Array<double> >("Coordinate");
        
      int dim = p_vec.size();
      Point pnt(dim), n(dim);

      if (dim == 3) {
          if (space_dim != 3) {
            sstream << "Point " << reg_name << 
              " specified using 3D coordinates but problem is " << 
              space_dim << " dimensional. Check input!" << std::endl;
            Errors::Message mesg(sstream.str());
            amanzi_throw(mesg);
          }

        pnt.set(p_vec[0], p_vec[1], p_vec[2]);
      }
      else if (dim == 2) {
          if (space_dim != 2) {
            sstream << "Plane " << reg_name << 
              " specified using 3D coordinates but problem is " << 
              space_dim << " dimensional. Check input!" << std::endl;
            Errors::Message mesg(sstream.str());
            amanzi_throw(mesg);
          }

        pnt.set(p_vec[0], p_vec[1]);
      }

      try {
        RegionPtr regptr = new PointRegion(reg_name, reg_id, pnt, lifecycle);
        return regptr;
      }
      catch (Errors::Message mesg) {
        mesg << "\n" << "Cannot create region of type Point";
        Exceptions::amanzi_throw(mesg);
      }
    }
  else if (shape == "Region: Logical")
    {
      Teuchos::ParameterList logical_params = reg_params.sublist(shape);
      std::string opstr = logical_params.get<std::string>("Operation");
      Teuchos::Array<std::string> region_names = 
        logical_params.get< Teuchos::Array<std::string> >("Regions");

      // No idea how to directly convert Teuchos::Array to std::vector
      // do it element by element

      std::vector<std::string> region_names1;
      int nelem = region_names.size();
      for (int i = 0; i < nelem; i++)
        region_names1.push_back(region_names[i]);      

      try {
        RegionPtr regptr = new LogicalRegion(reg_name, reg_id, opstr, region_names1, lifecycle);
        return regptr;
      }
      catch (Errors::Message mesg) {
        mesg << "\n" << "Cannot create region of type Logical";
        Exceptions::amanzi_throw(mesg);
      }
    }
  else 
    {
      Errors::Message mesg("ERROR: Cannot process region with given shape ");
      Exceptions::amanzi_throw(mesg);
    }

}


