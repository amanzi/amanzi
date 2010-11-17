/**
 * @file   HexMeshGenerator.hh
 * @author William A. Perkins
 * @date Sun Nov 14 18:58:57 2010
 * 
 * @brief  Declaration of the HexMeshGenerator class
 * 
 * 
 */
#ifndef _HexMeshGenerator_hh_
#define _HexMeshGenerator_hh_

#include <vector>
#include <map>

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

#include "Data.hh"


namespace Mesh_data
{
// -------------------------------------------------------------
//  class HexMeshGenerator
// -------------------------------------------------------------

/// A class to make a hexahedral mesh and put in a Mesh_data::Data instance
/**
 * This class is used to generate a hex mesh (in the form of a
 * Mesh_data::Data instance) in parallel.  The mesh is divided evenly,
 * by elements, amongst the involved processors.
 * 
 */
class HexMeshGenerator
{
protected:

  const Epetra_Comm *comm_;      /**< The parallel environment */

  const unsigned int ni_;       /**< Number of elements in the i-direction */
  const unsigned int nj_;       /**< Number of elements in the j-direction */
  const unsigned int nk_;       /**< Number of elements in the k-direction */

  const unsigned int ncell_;    /**< Number of global elemets */
  const unsigned int nvert_;    /**< Number of global vertexes */

  const double xorig_;          /**< The mesh origin x-coordinate */
  const double yorig_;          /**< The mesh origin y-coordinate */
  const double zorig_;          /**< The mesh origin z-coordinate */

  const double dx_;              /**< The element size in the x-direction */
  const double dy_;              /**< The element size in the y-direction */
  const double dz_;              /**< The element size in the z-direction */

  unsigned int cell0_;          /**< The (global, 0-based) index of first cell to generate */
  unsigned int cell1_;          /**< The (global, 0-based) index of last cell to generate */

  /// The global vertex ids used in this instance
  std::vector<unsigned int> vertex_gidx_;

  /// A map from global vertex id to local vertex id
  std::map<unsigned int, unsigned int> vertex_idxmap_;

  /// Get a global vertex index from direction indexes
  unsigned int global_vertex(const unsigned int& i, const unsigned int& j, const unsigned int& k) const;
  
  /// Get vertex direction indexes from global index
  void global_rvertex(const unsigned int& index, 
                      unsigned int& i, unsigned int& j, unsigned int& k) const;
  
  /// Get a global cell index from direction indexes
  unsigned int global_cell(const unsigned int& i, const unsigned int& j, const unsigned int& k) const;

  /// Get cell direction indexes from a global index
  void global_rcell(const unsigned int& index,
                    unsigned int& i, unsigned int& j, unsigned int& k) const;

  /// create the elements (in one block)
  void generate_the_elements(std::vector<int>& connectivity_map);

  /// create side sets for all sides of the domain
  std::vector<Side_set *> generate_the_sidesets(void);

  /// compute the vertex coordinates
  Coordinates<double>* generate_the_coordinates(void);

private:

  /// Private, unimplemented default constructor
  HexMeshGenerator();

  /// Private, unimplemented copy constructor to avoid copies
  HexMeshGenerator(const HexMeshGenerator& old);

public:

  /// (Collective) Default constructor.
  HexMeshGenerator(const Epetra_Comm *comm, 
                   const unsigned int& ni, const unsigned int& nj, const unsigned int& nk,
                   const double& xorigin = 0.0, 
                   const double& yorigin = 0.0, 
                   const double& zorigin = 0.0, 
                   const double& xdelta = 1.0, 
                   const double& ydelta = 1.0, 
                   const double& zdelta = 1.0);

  /// (Local) Destructor
  ~HexMeshGenerator(void);

  /// Get the number of cells
  unsigned int cells(void) const { return ncell_; }

  /// Get the number of local cells 
  unsigned int mycells(void) const { return cell1_ - cell0_ + 1; }

  /// Get the number of vertexes 
  unsigned int vertexes(void) const { return nvert_; }

  /// Get the number of local vertexes 
  unsigned int myvertexes(void) const { return vertex_gidx_.size(); }

  /// (Collective) Generate the (local local part of the) mesh and fill a Data instance
  Data *generate(void);

  /// (Collective) Generate a Epetra_Map for the cells
  Epetra_Map *cellmap(void);

  /// (Collective) Generate a Epetra_Map for the vertexes
  Epetra_Map *vertexmap(void);

}; // close class HexMeshGenerator

} // close namespace Mesh_data

#endif
