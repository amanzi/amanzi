/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _Layout_H_
#define _Layout_H_


#include <winstd.H>

#include <Amr.H>
#include <BCRec.H>

#include <MFTower.H>

#ifdef BL_USE_PETSC
#include <petsc.h>

#define CHKPETSC(n) CHKERRABORT(ParallelDescriptor::Communicator(),n);

#endif


struct IVSMap
  : std::map<IntVect,Stencil,IntVect::Compare>
{
  // Really a typedef, but compiler was unhappy in some cases
};

typedef IVSMap::const_iterator IVScit;

struct Node
{
  enum typeEnum{INIT=0, COVERED=1, VALID=2};
  Node()
    : iv(IntVect(D_DECL(-1,-1,-1))), level(-1), type(Node::INIT) {}
  Node(const IntVect& idx, int lev, typeEnum typ = INIT)
    : iv(idx), level(lev), type(typ) {}
  inline bool operator< (const Node& rhs) const
    {
      if (level < rhs.level) return true;
      if ((level == rhs.level) && IntVect::Compare()(iv,rhs.iv)) return true;
      return false;
    }
  inline bool operator!= (const Node& rhs) const
    {
      return ((*this) < rhs || rhs < (*this));
    }
  IntVect iv;
  int level;
  typeEnum type;
};


class Stencil
  : public std::map<Node,Real>
{
public:
  Stencil() {}
  typedef std::map<Node,Real>::const_iterator const_iterator;
  typedef std::map<Node,Real>::iterator iterator;
  Stencil& operator*=(Real val);
  Stencil& operator+=(const Stencil& rhs);
  Stencil& operator-=(const Stencil& rhs);
};

const Stencil operator*(Real val, const Stencil& lhs);
const Stencil operator*(const Stencil& lhs, Real val);
const Stencil operator+(const Stencil& lhs, const Stencil& rhs);
const Stencil operator-(const Stencil& lhs, const Stencil& rhs);

struct Layout
{
  typedef BaseFab<Node> NodeFab;
  typedef FabArray<NodeFab> MultiNodeFab;
  typedef BaseFab<int> IntFab;
  typedef FabArray<IntFab> MultiIntFab;

  Layout(Amr* parent = 0,
         int  nLevs = -1);

  Layout(const Array<IntVect>&  refRatio_array,
         const Array<BoxArray>& grid_array,
         const Array<Geometry>& geom_array,
         int                    n_levs);

  ~Layout();
  void Clear();
  void SetParent(Amr* new_parent);
  void Build();
  void SetNodeIds(BaseFab<int>& idFab, int lev, int grid) const;
  int NumberOfGlobalNodeIds() const {return nNodes_global;}
  int NumberOfLocalNodeIds() const {return nNodes_local;}
  int NumLevels() const {return nLevs;}
  const Array<BoxArray>& GridArray() const {return gridArray;}
  const Array<Geometry>& GeomArray() const {return geomArray;}
  const DistributionMapping& DistributionMap(int lev) const {return nodeIds[lev].DistributionMap();}
  const PArray<Layout::MultiNodeFab>& Nodes() const {return nodes;}
  const PArray<Layout::MultiNodeFab>& CrseNodes() const {return crseNodes;}
  const PArray<Layout::MultiIntFab>& NodeIds() const {return nodeIds;}
  const Array<BoxArray>& BndryCells() const {return bndryCells;}


  static BoxArray GetBndryCells (const BoxArray& ba,
                                 const IntVect&  ngrow,
                                 const Geometry& geom);


#ifdef BL_USE_PETSC
  Mat& Jacobian();
  Vec& JRowScale();
  PetscErrorCode MFTowerToVec(Vec& V,
                              const MFTower& mft,
                              int comp); // non-const because we track all created Vecs and destroy them in ~Layout

  PetscErrorCode VecToMFTower(MFTower& mft,
                              const Vec& V,
                              int        comp) const;
#endif

  const Array<IntVect>& RefRatio() const {return refRatio;}

  bool IsCompatible(const MFTower& mft) const;
  bool IsCompatible(Amr* parent, int  nLevs = -1);

  const MultiFab& Volume(int lev) const;
  const MultiFab& Area(int lev, int dir) const;

protected:
  void DestroyPetscStructures();
  void BuildMetrics();
  void SetGridsFromParent();

#ifdef BL_USE_PETSC
    Mat J_mat;
  Vec JRowScale_vec;
  Array<Vec*> vecs_I_created;
  Array<Mat*> mats_I_created;
#endif

  PArray<MultiNodeFab> nodes;
  PArray<MultiNodeFab> crseNodes;
  PArray<MultiIntFab> nodeIds;
  PArray<MultiIntFab> crseIds;
  Array<BoxArray> bndryCells;

  Array<Geometry> geomArray;
  Array<BoxArray> gridArray;
  Array<IntVect> refRatio;
  int nGrow, nLevs;
  bool initialized;
  Amr* parent;
  int nNodes_global, nNodes_local;

  PArray<MultiFab> volume;
  Array<PArray<MultiFab> > area;
};


std::ostream& operator<< (std::ostream&  os, const Stencil& a);
std::ostream& operator<< (std::ostream&  os, const Node& node);
std::ostream& operator<< (std::ostream&  os, const Layout::IntFab& ifab);
std::ostream& operator<< (std::ostream&  os, const Layout::NodeFab& nfab);
std::ostream& operator<< (std::ostream&  os, const Layout::MultiNodeFab& mnf);

#endif
