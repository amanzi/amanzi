#ifndef _MLConn_
#define _MLConn_


#include <winstd.H>

#ifndef WIN32
#include <unistd.h>
#endif

#include <Amr.H>
#include <BCRec.H>

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
    MFTower() {}
    virtual ~MFTower() {}
    const IntVect& refRatio(int lev) const {return ref_ratio[lev];} 
    int NumLevels() const {return ref_ratio.size()+1;}

    virtual bool IsCompatible(const MFTower& rhs) const = 0;
    virtual Real norm() const = 0;

protected:
    Array<IntVect> ref_ratio;
};

struct CCMFTower
    : public MFTower
{
    CCMFTower();
    ~CCMFTower();
    CCMFTower(const PArray<MultiFab>& pamf,
              const Array<IntVect>&   ratio);

    CCMFTower(const Array<BoxArray>& bat,
              const Array<IntVect>&  ratio,
              int                    nComp=1,
              int                    nGrow=0);

    void define(const Array<BoxArray>& bat,
                const Array<IntVect>&  ratio,
                int                    nComp=1,
                int                    nGrow=0);

    MultiFab& operator[](int i) {return mft[i];}
    const MultiFab& operator[](int i) const {return mft[i];}

    virtual bool IsCompatible(const MFTower& rhs) const;
    void AXPY(const CCMFTower& rhs,
              Real             p = 1);
    virtual Real norm() const;

protected:
    PArray<MultiFab> mft;
};

struct ECMFTower
    : public MFTower
{
    ECMFTower();
    virtual ~ECMFTower();
    ECMFTower(const Array<PArray<MultiFab> >& pamf,
              const Array<IntVect>&           ratio);

    ECMFTower(const Array<BoxArray>& bat,
              const Array<IntVect>&  ratio,
              int                    nComp=1,
              int                    nGrow=0);

    void define(const Array<BoxArray>& bat,
                const Array<IntVect>&  ratio,
                int                    nComp=1,
                int                    nGrow=0);

    PArray<MultiFab>& operator[](int i) {return mft[i];}
    const PArray<MultiFab>& operator[](int i) const {return mft[i];}

    virtual bool IsCompatible(const MFTower& rhs) const;
    void AXPY(const ECMFTower& rhs,
              Real             p = 1);
    virtual Real norm() const;
    
protected:
    Array<PArray<MultiFab> > mft;
};

#ifdef BL_USE_PETSC
class Stencil
    : public std::map<IntVect,Real,IntVect::Compare>
{
public:
    Stencil() {}
    typedef std::map<IntVect,Real>::const_iterator const_iterator;
    typedef std::map<IntVect,Real>::iterator iterator;
    Stencil& operator*=(Real val);
    Stencil& operator+=(const Stencil& rhs);
    Stencil& operator-=(const Stencil& rhs);
};

std::ostream& operator<< (std::ostream&  os, const Stencil& a);

struct Layout
{
    Layout(Amr* parent = 0);
    ~Layout();
    void Clear();
    void SetParent(Amr* new_parent);
    void Rebuild();
    void SetNodeIds(BaseFab<int>& idFab, int lev, int grid) const;
    int NumberOfGlobalNodeIds() const {return nNodes_global;}
    int NumberOfLocalNodeIds() const {return nNodes_local;}

    Mat& Jacobian();
    Vec& JRowScale();
    PetscErrorCode CCMFTowerToVec(Vec& V,
                                  const CCMFTower& mft,
                                  int comp); // non-const because we track all created Vecs and destroy them in ~Layout

    PetscErrorCode VecToCCMFTower(CCMFTower& mft,
                                  const Vec& V,
                                  int        comp) const;

    const Array<IntVect>& RefRatio() const {return refRatio;}

    bool IsCompatible(const CCMFTower& mft) const;
    void BuildCCMFTower(CCMFTower& mft,
                        int        nCompMF,
                        int        nGrowMF) const;

    void DoCoarseFineParallelInterp(CCMFTower& mft,
                                    int        sComp,
                                    int        nComp) const;

    void DoOp(const CCMFTower&    mft,
              int                 sComp,
              int                 nComp,
              const Array<BCRec>& bc,
              const ECMFTower&    beta);

    typedef BaseFab<Node> NodeFab;
    typedef FabArray<NodeFab> MultiNodeFab; 
    typedef BaseFab<int> IntFab;
    typedef FabArray<IntFab> MultiIntFab; 

protected:
    void DestroyPetscStructures();
    Mat J_mat;
    Vec JRowScale_vec;
    Array<Vec*> vecs_I_created;
    Array<Mat*> mats_I_created;

    PArray<MultiNodeFab> nodes;
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

    Array<Array<std::map<IntVect,Stencil,IntVect::Compare> > > bndryStencil; //[lev][dir][iv] = stencil to perform tangential interp
};


std::ostream& operator<< (std::ostream&  os, const Layout::IntFab& ifab);
std::ostream& operator<< (std::ostream&  os, const Layout::NodeFab& nfab);
std::ostream& operator<< (std::ostream&  os, const Layout::MultiNodeFab& mnf);
#endif

#endif
