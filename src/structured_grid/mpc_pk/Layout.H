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

std::ostream& operator<< (std::ostream&  os, const Node& node);

class Layout; // Forward declaration


struct MFTower
{
    MFTower(const Layout&    layout,
            const IndexType& t = IndexType(IntVect::TheZeroVector()),
            int              nComp = 1,
            int              nGrow = 0);

    MFTower(Layout&           layout,
            PArray<MultiFab>& pamf);

    int NumLevels() const {return nLevs;}
    int NComp() const {return nComp;}
    int NGrow() const {return nGrow;}
    MultiFab& operator[](int i) {return mft[i];}
    const MultiFab& operator[](int i) const {return mft[i];}
    const Array<BoxArray>& GridArray() const {return gridArray;}
    const Array<Geometry>& GeomArray() const {return geomArray;}
    const Array<IntVect>& RefRatio() const {return refRatio;}
    bool IsCompatible(const MFTower& rhs) const;

    Real norm() const;
    void AXPY(const MFTower& rhs,
              Real           p = 1);

protected:
    void define_alloc();
    void define_noalloc(PArray<MultiFab>& pamf);
    const Layout& layout;
    const Array<BoxArray>& gridArray;
    const Array<Geometry>& geomArray;
    const Array<IntVect>& refRatio;
    int nLevs, nComp, nGrow;
    PArray<MultiFab> mft;
    IndexType iType;
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

std::ostream& operator<< (std::ostream&  os, const Stencil& a);

struct Layout
{
    typedef std::map<IntVect,Stencil,IntVect::Compare> IVSMap;
    typedef IVSMap::const_iterator IVScit;
    typedef BaseFab<Node> NodeFab;
    typedef FabArray<NodeFab> MultiNodeFab; 
    typedef BaseFab<int> IntFab;
    typedef FabArray<IntFab> MultiIntFab; 

    Layout(Amr* parent = 0);
    ~Layout();
    void Clear();
    void SetParent(Amr* new_parent);
    void Rebuild();
    void SetNodeIds(BaseFab<int>& idFab, int lev, int grid) const;
    int NumberOfGlobalNodeIds() const {return nNodes_global;}
    int NumberOfLocalNodeIds() const {return nNodes_local;}
    int NumLevels() const {return nLevs;}
    const Array<BoxArray>& GridArray() const {return gridArray;}
    const Array<Geometry>& GeomArray() const {return geomArray;}
    const DistributionMapping& DistributionMap(int lev) const {return nodeIds[lev].DistributionMap();}
    const PArray<Layout::MultiNodeFab>& Nodes() const {return nodes;}
    const PArray<Layout::MultiNodeFab>& CrseNodes() const {return crseNodes;}
    const Array<BoxArray>& BndryCells() const {return bndryCells;}

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

protected:
    void DestroyPetscStructures();
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
};


std::ostream& operator<< (std::ostream&  os, const Layout::IntFab& ifab);
std::ostream& operator<< (std::ostream&  os, const Layout::NodeFab& nfab);
std::ostream& operator<< (std::ostream&  os, const Layout::MultiNodeFab& mnf);

struct MFTGrowFill
{
    MFTGrowFill(const Layout& layout);
    int NumLevels() const {return layout.NumLevels();}
    const Layout& GetLayout() const {return layout;}
    const Array<BoxArray>& GridArray() const {return layout.GridArray();}
    const Array<Geometry>& GeomArray() const {return layout.GeomArray();}
    void BuildStencil(const BCRec& bc,
                      int maxorder);

    void DoCoarseFineParallelInterp(MFTower& mft,
                                    int      sComp,
                                    int      nComp) const;
    void FillGrowCells(MFTower& mft,
                       int      sComp,
                       int      nComp) const;
    
    void CCtoECgrad(PArray<MFTower>& mft_ec,
                    const MFTower&   mft_cc,
                    Real             mult,
                    int              sComp,
                    int              dComp,
                    int              nComp) const;

    void ECtoCCdiv(MFTower&               mft_cc,
                   const PArray<MFTower>& mft_ec,
                   Real                   mult,
                   int                    sComp,
                   int                    dComp,
                   int                    nComp) const;
    
protected:
    void BuildCFParallelInterpStencil();
    
    int nLevs;
    const Layout& layout;
    const Array<BoxArray>& gridArray;
    const Array<Geometry>& geomArray;
    const Array<IntVect>& refRatio;

    Array<Array<Layout::IVSMap> > parallelInterpStencil; //[lev][dir][iv] = stencil to perform parallel bc Interp
    Array<Array<Layout::IVSMap> > perpInterpStencil; //[lev][dir][iv] = stencil to perform perpendicular bc Interp
};

class RichardOp
{
public:
    RichardOp(const MFTGrowFill& mftgrow);

    const MFTower& Coefficients(int dir) const {return coefs[dir];}
    MFTower& Coefficients(int dir) {return coefs[dir];}

    void Residual(MFTower&       Rmft,
                  const MFTower& Pnew_mft,
                  const MFTower& Pold_mft,
                  const Real     dt);

protected:
    int nLevs;
    const MFTGrowFill& mftgrow;
    const Layout& layout;
    const Array<BoxArray>& gridArray;
    const Array<Geometry>& geomArray;
    const Array<IntVect>& refRatio;

    PArray<MFTower> coefs;
    MFTower* lambda;     // Managed in/externally.  lambda_allocated must be true prior to many of the function calls
    MFTower* saturation; // Managed in/externally.  saturation_allocated must be true prior to many of the function calls
    bool lambda_allocated, saturation_allocated;
};

#endif
