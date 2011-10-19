/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   HexMeshGenerator.cc
 * @author William A. Perkins
 * @date Mon Aug  8 13:14:41 2011
 * 
 * @brief  Implementation of the HexMeshGenerator class
 * 
 * 
 */

#include <algorithm>
#include <boost/format.hpp>
#include <boost/lambda/lambda.hpp>
namespace bl = boost::lambda;

#include "HexMeshGenerator.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace Data {

// -------------------------------------------------------------
//  class HexMeshGenerator
// -------------------------------------------------------------

const unsigned int HexMeshGenerator::nvcell(8);

// -------------------------------------------------------------
// HexMeshGenerator:: constructors / destructor
// -------------------------------------------------------------
/** 
 * All the constructor does is decide how many cells each process owns
 *
 * Also, initializes the blocks from the geometric model regions
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
    cell_gidx_(), blocks_(),
    cell_idxmap_(), vertex_gidx_(), vertex_idxmap_()
{
  ASSERT(ni_ > 0);
  ASSERT(nj_ > 0);
  ASSERT(nk_ > 0);
  const int p_size(comm_->NumProc());
  ASSERT(ncell_ > p_size);      // require at least 1 cell per process
  const int p_rank(comm_->MyPID());
  unsigned int proccell(ncell_/p_size);
  int cell0_ = p_rank * proccell;
  int cell1_ = (p_rank + 1) * proccell - 1;

  // make sure all of the cells are owned by someone
  if ( (p_rank == p_size - 1)) {
    cell1_ = ncell_ - 1;
  }

  cell_gidx_.resize(cell1_ - cell0_ + 1);
  std::for_each(cell_gidx_.begin(), cell_gidx_.end(), bl::_1 = bl::var(cell0_)++);
  ASSERT(cell_gidx_.back() == cell1_);
  
  int lid(0);
  for (std::vector<unsigned int>::iterator c = cell_gidx_.begin();
       c != cell_gidx_.end(); c++, lid++) {
    cell_idxmap_[*c] = lid;
  }

  // Put the default block in the block list

  Block b;
  b.region = NULL;
  blocks_.push_back(b);
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
  ASSERT (index >= 0 && index < nvert_); // index is 0-based
  
  unsigned int tmp(index);
  k = tmp / ((ni_+1)*(nj_+1));
  tmp -= k*(ni_+1)*(nj_+1);
  j = tmp / (ni_+1);
  tmp -= j*(ni_+1);
  i = tmp;

  // std::cerr << "vertex " << index << ": " << i << ", " << j << ", " << k << std::endl;
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
  tmp -= k*(ni_)*(nj_);
  j = tmp / (ni_);
  tmp -= j*(ni_);
  i = tmp;
  return;
}  

// -------------------------------------------------------------
// HexMeshGenerator::generate_the_elements
// -------------------------------------------------------------
/** 
 * Generates the connectivity for @e all (local) elements.  During the
 * process, each cell index and connectivity is assigned to a block.
 * This also generates ::vertex_gidx_ and ::vertex_idxmap_
 * 
 */
void
HexMeshGenerator::generate_the_elements_(void)
{
  static const unsigned int nvcell(8);
  int num_elements(cell_gidx_.size());

  ASSERT(num_elements > 0);

  std::vector<int> connectivity_map(num_elements*nvcell);

  std::vector<int>::iterator c(connectivity_map.begin());
  for (std::vector<unsigned int>::iterator g = cell_gidx_.begin();
       g != cell_gidx_.end(); g++) {
    unsigned int i, j, k;
    unsigned int gidx(*g);
    global_rcell(gidx, i, j, k);   // g is 1-based; i, j, & k are 0-based

    // std::cerr << "cell " << gidx << ": " << i << ", " << j << ", " << k << std::endl;

    // This numbering does not match ExodusII
    // *(c+0) = global_vertex(i  , j  , k  );
    // *(c+1) = global_vertex(i+1, j  , k  );
    // *(c+2) = global_vertex(i+1, j+1, k  );
    // *(c+3) = global_vertex(i  , j+1, k  );
    // *(c+4) = global_vertex(i  , j  , k+1);
    // *(c+5) = global_vertex(i+1, j  , k+1);
    // *(c+6) = global_vertex(i+1, j+1, k+1);
    // *(c+7) = global_vertex(i  , j+1, k+1);

    // This numbering matches ExodusII
    *(c+0) = global_vertex(i  , j  , k  );
    *(c+1) = global_vertex(i  , j  , k+1);
    *(c+2) = global_vertex(i+1, j  , k+1);
    *(c+3) = global_vertex(i+1, j  , k  );
    *(c+4) = global_vertex(i  , j+1, k  );
    *(c+5) = global_vertex(i  , j+1, k+1);
    *(c+6) = global_vertex(i+1, j+1, k+1);
    *(c+7) = global_vertex(i+1, j+1, k  );

    // assign the global index to a block

    AmanziGeometry::Point p(xorig_ + (static_cast<double>(i) + 0.5)*dx_,
                            yorig_ + (static_cast<double>(j) + 0.5)*dy_,
                            zorig_ + (static_cast<double>(k) + 0.5)*dz_);
    
    std::vector<Block>::iterator r;
    for (r = blocks_.begin(); r != blocks_.end(); ++r) {
      if (r->region) {
        int id = r->region->id();
        if (r->region->inside(p)) {
          break;
        }
      }
    }
    if (r == blocks_.end()) {
      r = blocks_.begin();
    }
    r->gidx.push_back(gidx);
    std::copy(c, c+nvcell, std::back_inserter(r->connectivity));
    
    // std::cerr << comm_->MyPID() << ": "
    //           << "cell " << gidx << ": " << i << ", " << j << ", " << k 
    //           << ": " << p << ": block " << r->id << std::endl;

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

  for (std::vector<Block>::iterator r = blocks_.begin(); 
       r != blocks_.end(); ++r) {
    for (std::vector<int>::iterator i = r->connectivity.begin();
         i != r->connectivity.end(); i++) {
      unsigned int gidx(*i);
      unsigned int lidx(vertex_idxmap_[gidx]);
      *i = lidx;
    }
  }

  return;
}


// -------------------------------------------------------------
// HexMeshGenerator::generate_the_sidesets
// -------------------------------------------------------------

/** 
    Now these will be generated on demand based on the input specification
 * @return 
 */
std::vector<Side_set *>
HexMeshGenerator::generate_the_sidesets(void)
{
  std::map<unsigned int, bool> mine;
  std::vector<Side_set *> result;

  // Note: result contains pointers to allocate memory that must be
  // deleted elsewhere

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

  std::vector<double> attribute;
  attribute.clear();

  // Build the blocks

  generate_the_elements_();


  std::vector<Element_block *> tmpe; 
  std::vector<int> blkids;
  
  // Cell sets will be generated on demand at the Mesh_STK level based
  // on the input spec. Nothing to do here

  // Node sets will be generated on demand at the Mesh_STK level based
  // on the input spec. Nothing to do here

  Node_set *nodes;

  // generate side sets. This is also a dummy call since the side sets
  // will be generated on demand at the Mesh_STK level based on the
  // input spec

  std::vector<Side_set *> side_sets(generate_the_sidesets());
  std::vector<int> ssids(side_sets.size());
  for (int i = 0; i < ssids.size(); i++) ssids[i] = side_sets[i]->id();

  // compute coordinates

  Coordinates<double> *coords(generate_the_coordinates());

  // assemble mesh parameters

  Parameters* params(new Parameters("Generated", 3, 
                                    vertex_gidx_.size(), 
                                    cell_gidx_.size(),
                                    tmpe.size(), 1, ssids.size(),
                                    blkids, one, ssids));


  // finally assemble the mesh data

  std::vector<Node_set *> tmpn;
  // FIXME: delete  tmpn.push_back(nodes);

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
  std::vector<int> myidx;
  myidx.reserve(cell_gidx_.size());
  for (std::vector<Block>::const_iterator b = blocks_.begin();
       b != blocks_.end(); ++b) {
    std::copy(b->gidx.begin(), b->gidx.end(), 
              std::back_inserter(myidx));
  }
  if (onebased)
    std::for_each(myidx.begin(), myidx.end(), bl::_1 += 1);

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

} // namespace Data
} // namespace Mesh
} // namespace Amanzi



