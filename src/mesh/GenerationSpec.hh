/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

#pragma once

#include "Region.hh"
#include "RegionBox.hh"
#include "MeshDefs.hh"

namespace Teuchos {
  class ParameterList;
}

namespace Amanzi {
namespace AmanziMesh {

// -------------------------------------------------------------
//  class GenerationSpec
// -------------------------------------------------------------
/// Encapsulates parsing of mesh generation specifications
/**
 * The sole purpose of this class is to encapsulate the input
 * specifications for a generated unstructured mesh.  This should be
 * the only place in Amanzi where the parameter list for the
 * generation of an unstructured mesh is parsed.
 *
 */
class GenerationSpec {
 public:

  /// Constructor that uses a parameter list
  GenerationSpec(const Teuchos::ParameterList& parameter_list);

  /// Destructor
  ~GenerationSpec(void);

  /// Get the overall spatial domain
  const AmanziGeometry::RegionBox& domain(void) const
  { return *domain_; }

  /// Get the mesh dimensions in the x-direction
  const unsigned int& xcells(void) const { return nx_; }

  /// Get the mesh dimensions in the y-direction
  const unsigned int& ycells(void) const { return ny_; }

  /// Get the mesh dimensions in the z-direction
  const unsigned int& zcells(void) const { return nz_; }

  /// Get access to the list of blocks or zones
  AmanziGeometry::RegionVector::const_iterator block_begin(void) const
  { return blocks_.begin(); }

  /// Get access to the list of blocks or zones
  AmanziGeometry::RegionVector::const_iterator block_end(void) const
  { return blocks_.end(); }

  Partitioner_type partitioner() const { return partitioner_; }

 protected:

  /// overall mesh domain  FIXME: We already have a domain

  Teuchos::RCP<AmanziGeometry::RegionBox> domain_;

  unsigned int nx_;                     /**< number of cells in the x-direction */
  unsigned int ny_;                     /**< number of cells in the y-direction */
  unsigned int nz_;                     /**< number of cells in the y-direction */

  AmanziGeometry::RegionVector blocks_; /**< list of mesh subdomains */

  Partitioner_type partitioner_;

  /// fill attributes from specified list
  void parse_(const Teuchos::ParameterList &parameter_list);

 private:

  /// private, undefined copy constructor to avoid unwanted copies

  GenerationSpec(const GenerationSpec& old);

};


} // end namespace AmanziMesh
} // end namespace Amanzi


