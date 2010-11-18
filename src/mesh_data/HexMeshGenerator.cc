/**
 * @file   HexMeshGenerator.cc
 * @author William A. Perkins
 * @date Thu Nov 18 10:05:06 2010
 * 
 * @brief  Implementation of the HexMeshGenerator class
 * 
 * 
 */

#include <algorithm>
#include <boost/lambda/lambda.hpp>
namespace bl = boost::lambda;

#include "HexMeshGenerator.hh"

namespace Mesh_data
{

// -------------------------------------------------------------
//  class HexMeshGenerator
// -------------------------------------------------------------

// -------------------------------------------------------------
// HexMeshGenerator:: constructors / destructor
// -------------------------------------------------------------
/** 
 * All the constructor does is decide how many cells each process owns
 * 
 * @param comm parallel environment
 * @param ni number of cells in the i-direction
 * @param nj number of cells in the j-direction
 * @param nk number of cells in the k-direction
 * @param xorigin x-coordinate of mesh origin
 * @param yorigin y-coordinate of mesh origin
 * @param zorigin z-coordinate of mesh origin
 * @param xdelta cell size in the x-direction
 * @param ydelta cell size in the y-direction
 * @param zdelta cell size in the z-direction
 */
HexMeshGenerator::HexMeshGenerator(const Epetra_Comm *comm, 
                                   const unsigned int& ni, 
                                   const unsigned int& nj, 
                                   const unsigned int& nk,
                                   const double& xorigin, 
                                   const double& yorigin, 
                                   const double& zorigin, 
                                   const double& xdelta, 
                                   const double& ydelta, 
                                   const double& zdelta)
  : comm_(comm), 
    ni_(ni), nj_(nj), nk_(nk),
    ncell_(ni_*nj_*nk_), nvert_((ni_+1)*(nj_+1)*(nk_+1)),
    xorig_(xorigin), yorig_(yorigin), zorig_(zorigin),
    dx_(xdelta), dy_(ydelta), dz_(zdelta),
    cell0_(0), cell1_(0)
{
  ASSERT(ni_ > 0);
  ASSERT(nj_ > 0);
  ASSERT(nk_ > 0);
  const int p_size(comm_->NumProc());
  ASSERT(ncell_ > p_size);      // require at least 1 cell per process
  const int p_rank(comm_->MyPID());
  unsigned int proccell(ncell_/p_size);
  cell0_ = p_rank * proccell;
  cell1_ = (p_rank + 1) * proccell - 1;

  // make sure all of the cells are owned by someone
  if ( (p_rank == p_size - 1)) {
    cell1_ = ncell_ - 1;
  }

}

HexMeshGenerator::~HexMeshGenerator(void)
{
  // empty
}

// -------------------------------------------------------------
// HexMeshGenerator::global_vertex
// -------------------------------------------------------------
/** 
 * Determines the global vertex index from the directional
 * indexes
 * 
 * @param i vertex index in the i-direction (0-based)
 * @param j vertex index in the j-direction (0-based)
 * @param k vertex index in the k-direction (0-based)
 * 
 * @return global vertex index (0-based)
 */
unsigned int
HexMeshGenerator::global_vertex(const unsigned int& i, 
                                const unsigned int& j, 
                                const unsigned int& k) const
{
  ASSERT(i >= 0 && i < ni_ + 1);
  ASSERT(j >= 0 && j < nj_ + 1);
  ASSERT(k >= 0 && k < nk_ + 1);
  unsigned int result;
  result = (i + j*(ni_+1) + k*(ni_+1)*(nj_+1)); // 0-based
  return result;
}

// -------------------------------------------------------------
// HexMeshGenerator::global_rvertex
// -------------------------------------------------------------
/** 
 * Does the opposite of ::global_vertex.
 * 
 * @param index global vertex index (0-based)
 * @param i vertex index in the i-direction (0-based)
 * @param j vertex index in the j-direction (0-based)
 * @param k vertex index in the k-direction (0-based)
 */
void 
HexMeshGenerator::global_rvertex(const unsigned int& index, 
                                 unsigned int& i, unsigned int& j, unsigned int& k) const
{
  ASSERT (index >= 0 && index <= nvert_ - 1); // index is 0-based
  unsigned int tmp(index);
  k = tmp / ((ni_+1)*(nj_+1));
  tmp -= k*(nj_+1)*(nj_+1);
  j = tmp / (ni_+1);
  tmp -= j*(ni_+1);
  i = tmp;
  return;
}

// -------------------------------------------------------------
// HexMeshGenerator::global_cell
// -------------------------------------------------------------
/** 
 * Compute a global cell index from the directional indexes
 * 
 * @param i cell index in the i-direction (0-based)
 * @param j cell index in the j-direction (0-based)
 * @param k cell index in the k-direction (0-based)
 * 
 * @return global cell index (0-based)
 */
unsigned int
HexMeshGenerator::global_cell(const unsigned int& i, 
                              const unsigned int& j, 
                              const unsigned int& k) const
{
  ASSERT(i >= 0 && i < ni_);
  ASSERT(j >= 0 && j < nj_);
  ASSERT(k >= 0 && k < nk_);
  unsigned int result;
  result = (i + j*(ni_) + k*(ni_)*(nj_)); // 0-based
  return result;
}

// -------------------------------------------------------------
// HexMeshGenerator::global_rcell
// -------------------------------------------------------------
/** 
 * Does the opposite of ::global_cell.
 * 
 * @param index global cell index (0-based)
 * @param i cell index in the i-direction (0-based)
 * @param j cell index in the j-direction (0-based)
 * @param k cell index in the k-direction (0-based)
 */
void
HexMeshGenerator::global_rcell(const unsigned int& index,
                    unsigned int& i, unsigned int& j, unsigned int& k) const
{
  ASSERT (index >= 0 && index < ncell_); // index is 1-based
  unsigned int tmp(index);
  k = tmp / ((ni_)*(nj_));
  tmp -= k*(nj_)*(nj_);
  j = tmp / (ni_);
  tmp -= j*(ni_);
  i = tmp;
  return;
}  

// -------------------------------------------------------------
// HexMeshGenerator::generate_the_elements
// -------------------------------------------------------------
/** 
 * this also generates this->vertex_gidx_ and this->vertex_gidx_

 * 
 * @param connectivity_map cell connectivity to vertexes (local 0-based vertex indexes)
 * @param involved_global_vertexes an ordered list of the global indexes (0-based) of all the vertexes involved in @connectivity_map
 * @param vertex_idxmap a map for looking up local vertex indexes with global vertex indexes
 */


void
HexMeshGenerator::generate_the_elements(std::vector<int>& connectivity_map)
{
  static const unsigned int nvcell(8);
  int num_elements(cell1_-cell0_+1);

  ASSERT(num_elements > 0);

  connectivity_map.clear();
  connectivity_map.resize(num_elements*nvcell);

  std::vector<int>::iterator c(connectivity_map.begin());
  for (unsigned int g = cell0_; g <= cell1_; g++) {
    unsigned int i, j, k;
    global_rcell(g, i, j, k);   // g is 1-based; i, j, & k are 0-based

    *(c+0) = i + j*(ni_+1) + k*(ni_+1)*(nj_+1);
    *(c+1) = (i+1) + j*(ni_+1) + k*(ni_+1)*(nj_+1);
    *(c+2) = (i+1) + (j+1)*(ni_+1) + k*(ni_+1)*(nj_+1);
    *(c+3) = i + (j+1)*(ni_+1) + k*(ni_+1)*(nj_+1);
    *(c+4) = i + j*(ni_+1) + (k+1)*(ni_+1)*(nj_+1);
    *(c+5) = (i+1) + j*(ni_+1) + (k+1)*(ni_+1)*(nj_+1);
    *(c+6) = (i+1) + (j+1)*(ni_+1) + (k+1)*(ni_+1)*(nj_+1);
    *(c+7) = i + (j+1)*(ni_+1) + (k+1)*(ni_+1)*(nj_+1);
    c += nvcell;
  }

  // collect up the global indexes of the vertexes used

  { 
    std::vector<int> tmp(connectivity_map);
    std::sort(tmp.begin(), tmp.end());
    vertex_gidx_.clear();
    vertex_gidx_.reserve(tmp.size());
    std::unique_copy(tmp.begin(), tmp.end(), std::back_inserter(vertex_gidx_));
  }

  // relate the global index back to the local index (a std::map may
  // not be a good thing here, a large enough, if it can be made,
  // std::vector would work the same)

  vertex_idxmap_.clear();
  unsigned int lidx(0);
  for (std::vector<unsigned int>::iterator i = vertex_gidx_.begin();
       i != vertex_gidx_.end(); i++, lidx++) {
    vertex_idxmap_[*i] = lidx;
  }

  // replace the global indexes in the connectivity vector with local
  // indexes 

  for (std::vector<int>::iterator i = connectivity_map.begin();
       i != connectivity_map.end(); i++) {
    unsigned int gidx(*i);
    unsigned int lidx(vertex_idxmap_[gidx]);
    *i = lidx;
  }

  return;
}

// -------------------------------------------------------------
// HexMeshGenerator::generate_the_sidesets
// -------------------------------------------------------------

/** 
 * Cell indexes are 0-based.  Side indexes are 1-based:
 * 
 *  - West = side 5 (i = 0)
 *  - East = side 3
 *  - South = side 2 (j = 0)
 *  - North = side 4
 *  - Bottom = side 1 (k = 0)
 *  - Top = side 6
 * 
 * @return 
 */
std::vector<Side_set *>
HexMeshGenerator::generate_the_sidesets(void)
{
  std::vector<Side_set *> result;

  static const int SIX(6);
  for (int side = 1; side <= SIX; side++) {
    int imin(0), imax(ni_);
    int jmin(0), jmax(nj_);
    int kmin(0), kmax(nk_);

    std::string ssname("Bogus");

    switch (side) {
    case (1):
      ssname = "Bottom";
      kmax = 1;
      break;
    case (2):
      ssname = "South";
      jmax = 1;
      break;
    case (3):
      ssname = "East";
      imin = imax - 1;
      break;
    case (4):
      ssname = "North";
      jmin = jmax - 1;
      break;
    case (5):
      ssname = "West";
      imax = 1;
      break;
    case (6):
      ssname = "Top";
      kmin = kmax - 1;
      break;
    }

    std::vector<int> clist;
    std::vector<int> slist;

    for (unsigned int i = imin; i < imax; i++) {
      for (unsigned int j = jmin; j < jmax; j++) {
        for (unsigned int k = kmin; k < kmax; k++) {
          unsigned cellidx(global_cell(i, j, k));
          if (cellidx >= cell0_ && cellidx <= cell1_) {
            clist.push_back(cellidx);
            slist.push_back(side);
          }
        }
      }
    }

    Side_set* ss = Side_set::build_from(side, clist, slist, ssname);
    result.push_back(ss);
  }
  return result;
}

// -------------------------------------------------------------
// HexMeshGenerator::generate_the_coordinates
// -------------------------------------------------------------
/** 
 * Computes a set of coordinates for the vertexes specified by global
 * index
 * 
 * @return the coordinates
 */
Coordinates<double>*
HexMeshGenerator::generate_the_coordinates()
{
  ASSERT(vertex_gidx_.size() > 0);
  Coordinates<double>* result(new Coordinates<double>(vertex_gidx_.size(), 3));
  
  for (std::vector<unsigned int>::const_iterator v = vertex_gidx_.begin();
       v != vertex_gidx_.end(); v++) {
    unsigned int gidx(*v);
    unsigned int lidx(vertex_idxmap_[gidx]);
    unsigned int i, j, k;
    global_rvertex(gidx, i, j, k);

    (*result)(lidx, 0) = xorig_ + i*dx_;
    (*result)(lidx, 1) = yorig_ + j*dy_;
    (*result)(lidx, 2) = zorig_ + k*dz_;
  }

  return result;
}


// -------------------------------------------------------------
// HexMeshGenerator::generate
// -------------------------------------------------------------
Data *
HexMeshGenerator::generate(void)
{
  const std::vector<int> one(1, 1);

  Element_block *blk;
  std::vector<double> attribute;
  {
    // Just to be safe, put this in a block, so connectivity is not
    // used after it's destroyed by Element_block::build_from()

    std::vector<int> connectivity;
    attribute.clear();

    generate_the_elements(connectivity);


    blk = Element_block::build_from(0, "Generated Elements", 
                                    cell1_ - cell0_ + 1, HEX,
                                    connectivity,
                                    attribute);
    attribute.clear();
  }

  // Build a list of nodes for the node set.  All nodes involved in
  // this process's cells are included.

  Node_set *nodes;
  { 
    // Just to be safe, put this in a block, so nodelist is not
    // used after it's destroyed by Node_set::build_from()

    std::vector<int> nodelist(vertex_gidx_.size());
    {
      int i(0);
      for (std::vector<int>::iterator n = nodelist.begin();
           n != nodelist.end(); n++, i++) *n = i;
    }

    // the temporary string seems necessary for some reason

    attribute.clear();
    std::string n("Generated Nodes");
    nodes = Node_set::build_from(0, nodelist, attribute, n);
    attribute.clear();
  }

  // generate side sets

  std::vector<Side_set *> side_sets(generate_the_sidesets());
  std::vector<int> ssids(side_sets.size());
  for (int i = 0; i < ssids.size(); i++) ssids[i] = side_sets[i]->id();

  // compute coordinates

  Coordinates<double> *coords(generate_the_coordinates());

  // assemble mesh parameters

  Parameters* params(new Parameters("Generated", 3, 
                                    vertex_gidx_.size(), 
                                    cell1_ - cell0_ + 1,
                                    1, 1, ssids.size(),
                                    one, one, ssids));


  // finally assemble the mesh data

  std::vector<Element_block *> tmpe; 
  tmpe.push_back(blk);

  std::vector<Node_set *> tmpn;
  tmpn.push_back(nodes);

  Data *result = Data::build_from(params, coords, tmpe, side_sets, tmpn);

  // FIMME: Why can't I delete these things?  
  // Because the Data instance takes ownership of the pointers

  // delete blk;
  // delete nodes;
  // delete coords;
  // delete params;

  return result;
}

// -------------------------------------------------------------
// HexMeshGenerator::cellmap
// -------------------------------------------------------------
/** 
 * Creates a parallel map relating local cell index to global index.
 * The map's global indexes are 0-based unless @c onebased is true.
 * The created map should not overlap with other processors, so it
 * does represent cell ownership.
 * 
 * @param onebased if true, generate a map with 1-based global indexes
 * 
 * @return map relating local to global cell indexes
 */
Epetra_Map*
HexMeshGenerator::cellmap(bool onebased)
{
  std::vector<int> myidx(cell1_-cell0_ + 1);
  int cidx(cell0_);
  if (onebased) cidx += 1;

  std::for_each(myidx.begin(), myidx.end(), bl::_1 = bl::var(cidx)++);

  Epetra_Map *result(new Epetra_Map(ncell_, myidx.size(), &myidx[0], 
                                    (onebased ? 1 : 0), *comm_));

  return result;
}

// -------------------------------------------------------------
// HexMeshGenerator::vertexmap
// -------------------------------------------------------------
/** 
 * Creates a parallel map relating local vertex index to global index.
 * The map's global indexes are 0-based unless @c onebased is true.
 * The created map will have to overlap with other processors, so it
 * should not be used for vertex ownership.
 * 
 * @param onebased if true, generate an Epetra_Map with 1-based global indexes
 * 
 * @return map relating local to global vertex indexes
 */
Epetra_Map *
HexMeshGenerator::vertexmap(bool onebased)
{
  std::vector<int> myidx(vertex_gidx_.size());
  std::copy(vertex_gidx_.begin(), vertex_gidx_.end(), myidx.begin());

  if (onebased)
    std::for_each(myidx.begin(), myidx.end(), bl::_1 += 1);
  
  Epetra_Map *result(new Epetra_Map(-1, vertex_gidx_.size(), &myidx[0], 
                                    (onebased ? 1 : 0), *comm_));

  return result;
}

} // close namespace Mesh_data



