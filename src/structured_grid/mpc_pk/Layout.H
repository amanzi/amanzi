#ifndef _MLConn_
#define _MLConn_


#include <winstd.H>

#ifndef WIN32
#include <unistd.h>
#endif

#include <Amr.H>

#ifdef BL_USE_PETSC
#include <petscmat.h>
#define CHKPETSC(n) CHKERRABORT(ParallelDescriptor::Communicator(),n);
#endif

struct Node
{
    enum typeEnum{INIT=0, COVERED=1, VALID=2};
    Node()
        : iv(IntVect(D_DECL(-1,-1,-1))), level(-1), grid(-1), type(Node::INIT) {}
    Node(const IntVect& idx, int lev, int grd, typeEnum typ = INIT)
        : iv(idx), level(lev), grid(grd), type(typ) {}
    inline bool operator< (const Node& rhs) const
        {
            if (level < rhs.level) return true;
            if ((level == rhs.level) && iv.lexLT(rhs.iv)) return true;
            return false;
        }
    inline bool operator!= (const Node& rhs) const
        {
            return ((*this) < rhs || rhs < (*this));
        }
    IntVect iv;
    int level;
    int grid;
    typeEnum type;
};

std::ostream& operator<< (std::ostream&  os, const Node& node);

struct MFTower
{
    MFTower();
    MFTower(const PArray<MultiFab>& pamf,
            const Array<IntVect>&   ratio);

    MFTower(const Array<BoxArray>& bat,
            const Array<IntVect>&  ratio,
            int                    nComp=1,
            int                    nGrow=0);

    void define(const Array<BoxArray>& bat,
                const Array<IntVect>&  ratio,
                int                    nComp=1,
                int                    nGrow=0);

    int size() const {return mft.size();}
    IntVect refRatio(int lev) const {return ref_ratio[lev];} 
    MultiFab& operator[](int i) {return mft[i];}
    const MultiFab& operator[](int i) const {return mft[i];}

    bool IsCompatible(const MFTower& rhs) const;
    void AXPY(const MFTower& rhs,
              Real           p = 1);
    Real norm() const;

protected:
    Array<IntVect> ref_ratio;
    PArray<MultiFab> mft;
};

#ifdef BL_USE_PETSC
struct Layout
{
    Layout(Amr* parent = 0);
    void Clear();
    void SetParent(Amr* new_parent);
    void Rebuild();
    void SetNodeIds(BaseFab<int>& idFab, int lev, int grid) const;
    int NumberOfGlobalNodeIds() const {return nNodes_global;}
    int NumberOfLocalNodeIds() const {return nNodes_local;}

#ifdef BL_USE_PETSC
    Mat& Jacobian();
    PetscErrorCode MFTowerToVec(Vec& V,
                                const MFTower& mft,
                                int comp) const;

    PetscErrorCode VecToMFTower(MFTower&   mft,
                                const Vec& V,
                                int        comp) const;
#endif
    const Array<IntVect>& RefRatio() const {return refRatio;}

    bool IsCompatible(const MFTower& mft) const;
    void BuildMFTower(MFTower& mft,
                      int      nCompMF,
                      int      nGrowMF) const;


    typedef BaseFab<Node> NodeFab;
    typedef FabArray<NodeFab> MultiNodeFab; 
    typedef BaseFab<int> IntFab;
    typedef FabArray<IntFab> MultiIntFab; 

protected:
#ifdef BL_USE_PETSC
    void AllocateJ();
    Mat* J;
#endif

    PArray<MultiNodeFab> nodes;
    PArray<MultiIntFab> nodeIds;
    Array<Geometry> geomArray;
    Array<BoxArray> gridArray;
    Array<IntVect> refRatio;
    int nGrow, nLevs;
    bool initialized;
    Amr* parent;
    int nNodes_global, nNodes_local;

};


std::ostream& operator<< (std::ostream&  os, const Layout::IntFab& ifab);
#endif

#endif
