#ifndef __CELL_TOPOLOGY_H__
#define __CELL_TOPOLOGY_H__

//      7---------6
//     /|        /|
//    / |       / |            5
//   /  |      /  |            | 2
//  4---------5   |            |/
//  |   |     |   |       3----+----1
//  |   3-----|---2           /|
//  |  /      |  /           0 |
//  | /       | /              4
//  |/        |/
//  0---------1

namespace cell_topology {

// Outward oriented faces of the reference hexahedron.
static int HexFaceVert[6][4] = {{0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {3,0,4,7}, {0,3,2,1}, {4,5,6,7}};

// Indices of the three faces adjacent to each corner of the reference hexahedron.
static int HexCornerFace[8][3] = {{0,3,4}, {0,1,4}, {1,2,4}, {2,3,4}, {0,3,5}, {0,1,5}, {1,2,5}, {2,3,5}};  

}

#endif
