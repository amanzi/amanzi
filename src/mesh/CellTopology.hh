/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//
// Static variables containing face and node numbering conventions for cell
// types listed in MeshDefs.hh
//
// Knowing the above cell type, one can look into the following arrays
// to determine how many faces are expected in a cell, how many nodes
// are expected in a particular face in the cell and which nodes of
// the cell_get_nodes result define that face in the cell
//

#ifndef AMANZI_CELL_TOPOLOGY_H_
#define AMANZI_CELL_TOPOLOGY_H_


// Node and face numbering conventions for standard 3D cell types
namespace Amanzi {
namespace AmanziMesh {
namespace Topology {

// expected number of faces for standard cells
// 0 indicates that it is not possible to specify the number for this cell type
//
static const int nface_std[9] = {0,3,4,0,4,5,5,6,0};

// Expected number of nodes for each face of standard cells
// 0 indicates that it is not possible to specify the number for this cell type
//
static const int nfnodes_std[9][6] = {
  {0,0,0,0,0,0},         // UNKNOWN
  {2,2,2,0,0,0},         // TRI
  {2,2,2,2,0,0},         // QUAD
  {0,0,0,0,0,0},         // POLYGON
  {3,3,3,3,0,0},         // TET
  {4,4,4,3,3,0},         // PRISM
  {3,3,3,3,4,0},         // PYRAMID
  {4,4,4,4,4,4}          // HEX
};

// Expected cell-face-node pattern for standard cells
// -1 indicates that it is not possible to specify the number for this cell type
// or that it is a filler for extra nodes in the array
//
static const int fnodes_std[9][6][4] = {
  { // UNKNOWN
    {-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}
  },
  {
    // TRI
    {0,1,-1,-1},{1,2,-1,-1},{2,0,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}
  },
  {
    // QUAD
    {0,1,-1,-1},{1,2,-1,-1},{2,3,-1,-1},{3,0,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}
  },
  { // POLYGON
    {-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}
  },
  {
    // TET
    {0,1,3,-1},{1,2,3,-1},{2,0,3,-1},{0,2,1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}
  },
  {
    // PRISM
    {0,1,4,3},{1,2,5,4},{2,0,3,5},{0,2,1,-1},{3,4,5,-1},{-1,-1,-1,-1}
  },
  {
    // PYRAMID
    {0,1,4,-1},{1,2,4,-1},{2,3,4,-1},{3,0,4,-1},{0, 3, 2, 1},{-1,-1,-1,-1}
  },
  {
    // HEX
    {0,1,5,4},{1,2,6,5},{2,3,7,6},{0,4,7,3},{0,3,2,1},{4,5,6,7}
  }
};



// LEFT THESE IN HERE BUT WE SHOULD LOOK CAREFULLY IF WE NEED IT

// Outward oriented faces of the reference hexahedron.
static const int HexFaceVert[6][4] = {{0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {3,0,4,7}, {0,3,2,1}, {4,5,6,7}};

// Indices of the three faces adjacent to each corner of the reference hexahedron.
static const int HexCornerFace[8][3] = {{0,3,4}, {0,1,4}, {1,2,4}, {2,3,4}, {0,3,5}, {0,1,5}, {1,2,5}, {2,3,5}};

// Corner tetrahedra (first 8) followed by the internal tetrahedra (last 2).
static const int HexTetVert[10][4] = {{0,1,3,4}, {1,2,0,5}, {2,3,1,6},
                                      {3,0,2,7}, {4,7,5,0}, {5,4,6,1},
                                      {6,5,7,2}, {7,6,4,3}, {0,2,7,5},
                                      {1,3,4,6}
};

} // namespace Topology
} // namespace AmanziMesh
} // namespace Amanzi

#endif
