/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/**
 * @file   GenerationSpec.cc
 * @author William A. Perkins
 * @date Wed Sep 28 10:20:33 2011
 *
 * @brief
 *
 *
 */


#include <Teuchos_ParameterList.hpp>

#include "GenerationSpec.hh"

namespace Amanzi {
namespace AmanziMesh {

// -------------------------------------------------------------
//  class GenerationSpec
// -------------------------------------------------------------

// -------------------------------------------------------------
// GenerationSpec:: constructors / destructor
// -------------------------------------------------------------
GenerationSpec::GenerationSpec(const Teuchos::ParameterList& parameter_list)
    : domain_(Teuchos::null), nx_(0), ny_(0), nz_(0), blocks_()
{
  parse_(parameter_list);
}

GenerationSpec::~GenerationSpec(void)
{
  // empty
}

// -------------------------------------------------------------
// GenerationSpec::parse_
// -------------------------------------------------------------
void
GenerationSpec::parse_(const Teuchos::ParameterList& parameter_list)
{
  // read the parameters from the parameter list
  Teuchos::Array<int> ncells = parameter_list.get< Teuchos::Array<int> >("number of cells");
  Teuchos::Array<double> low_corner = parameter_list.get< Teuchos::Array<double> >("domain low coordinate");
  Teuchos::Array<double> high_corner = parameter_list.get< Teuchos::Array<double> >("domain high coordinate");

  unsigned int dimension = ncells.size();

  nx_ = ncells[0];
  ny_ = ncells[1];
  if (dimension == 3) nz_ = ncells[2];

  AmanziGeometry::Point p0(dimension), p1(dimension);

  p0.set(&(low_corner[0]));
  p1.set(&(high_corner[0]));

  if (parameter_list.isParameter("partitioner")) {
    std::string partitioner_str = parameter_list.get<std::string>("partitioner");
    if (partitioner_str == "METIS" || partitioner_str == "metis")
      partitioner_ = Partitioner_kind::METIS;
    else if (partitioner_str == "ZOLTAN_GRAPH" || partitioner_str == "zoltan_graph")
      partitioner_ = Partitioner_kind::ZOLTAN_GRAPH;
    else if (partitioner_str == "ZOLTAN_RCB" || partitioner_str == "zoltan_rcb")
      partitioner_ = Partitioner_kind::ZOLTAN_RCB;
  }
  
  domain_ = Teuchos::rcp(new AmanziGeometry::RegionBox("GenDomain", 0, p0, p1));
}

} // end namespace AmanziMesh
} // end namespace Amanzi
