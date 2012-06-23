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
              Real           p = 1,
              int            sComp = 0,
              int            dComp = 0,
              int            nComp = 1,
              int            nGrow = 0);

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
                      int          maxorder);

    void DoCoarseFineParallelInterp(MFTower& mft,
                                    int      sComp=1,
                                    int      nComp=1) const;
    void FillGrowCells(MFTower& mft,
                       int      sComp=0,
                       int      nComp=1) const;
    
    void FillGrowCellsSimple(MFTower& mft,
                             int      sComp=0,
                             int      nComp=1) const;
    
    void CCtoECgrad(PArray<MFTower>& mft_ec,
                    const MFTower&   mft_cc,
                    Real             mult=1,
                    int              sComp=0,
                    int              dComp=0,
                    int              nComp=1) const;

    void ECtoCCdiv(MFTower&               mft_cc,
                   const PArray<MFTower>& mft_ec,
                   Real                   mult=1,
                   int                    sComp=0,
                   int                    dComp=0,
                   int                    nComp=1) const;
    
    struct MyHash // Unique characteristics
    {
        MyHash()
            : maxorder(-1) {}
        MyHash(const BCRec& _bc, int maxorder);
        BCRec bc;
        int maxorder;
    };

protected:
    void BuildCFParallelInterpStencil();


    MyHash myhash;
    int nLevs;
    const Layout& layout;
    const Array<BoxArray>& gridArray;
    const Array<Geometry>& geomArray;
    const Array<IntVect>& refRatio;

    Array<Array<Layout::IVSMap> > parallelInterpStencil; //[lev][dir][iv] = stencil to perform parallel bc Interp
    Array<Array<Layout::IVSMap> > perpInterpStencil; //[lev][dir][iv] = stencil to perform perpendicular bc Interp
};

inline bool operator==(const MFTGrowFill::MyHash& lhs, const MFTGrowFill::MyHash& rhs);
inline bool operator!=(const MFTGrowFill::MyHash& lhs, const MFTGrowFill::MyHash& rhs);

class PMAmr;
class MFTGrowFill;
#include <PMAmr.H>
#include <PorousMedia.H>

class CalcCoefficients
{
public:
    CalcCoefficients(MFTGrowFill& mftgrow,
                     PMAmr&       pmamr);

    void operator()(PArray<MFTower>&       coefficients,
                    MFTower&               pressure,
                    MFTower&               saturation,
                    MFTower&               lambda,
                    const PArray<MFTower>& DarcyVelocity,
                    Real                   mult = 1,
                    int                    presComp = 0,
                    int                    satComp = 0,
                    int                    lamComp = 0,
                    int                    darcComp = 0,
                    int                    coefComp = 0,
                    int                    nComp = 1);
protected:
    int nLevs;
    MFTGrowFill& mftgrow;
    const Layout& layout;
    PMAmr& pmamr;
};

struct RichardContext
{
    RichardContext(PMAmr&           pmamr,
                   MFTGrowFill&     mftgfill,
                   MFTower&         pressure_old,
                   MFTower&         saturation,
                   MFTower&         lambda,
                   PArray<MFTower>& darcyVelocity,
                   const BCRec&     pressure_bc,
                   int              pressure_maxorder);
    PMAmr& PMAMR() {return pmamr;}
    CalcCoefficients& CalcCoefs() {return *calcCoefs;}
    MFTower& PressureOld() {return pressure_old;}
    MFTower& Saturation() {return saturation;}
    MFTower& Lambda() {return lambda;}
    PArray<MFTower>& DarcyVelocity() {return darcyVelocity;}

protected:
    MFTower& pressure_old;
    MFTower& saturation;
    MFTower& lambda;
    PArray<MFTower>& darcyVelocity;

    PMAmr& pmamr;
    Layout& layout;
    MFTGrowFill& mftgfill;
    CalcCoefficients* calcCoefs;
};

class RichardOp
{
public:
    RichardOp(RichardContext& richardContext);
    ~RichardOp();

    const Array<BoxArray>& gridArray() const {return layout.GridArray();}
    const Array<Geometry>& geomArray() const {return layout.GeomArray();}
    const Array<IntVect>& refRatio() const {return layout.RefRatio();}

    void SetBC(const BCRec& bc, 
               int          maxorder); // Note: triggers build of mftgrow stencil

    void Residual(MFTower&               residual,
                  MFTower&               pressure,
                  const Real             dt);

    void DarcyVelocity(PArray<MFTower>& velocity_ec,
                       MFTower&         pressure_cc,
                       int              sComp = 0,
                       int              dComp = 0,
                       int              nComp=1) const;


    const MFTower& Coefficients(int dir) const {return coefs[dir];}
    MFTower& Coefficients(int dir) {return coefs[dir];}

protected:
    void AllocLambda();
    void AllocSaturation();

    int nLevs;
    MFTGrowFill* mftgrow;
    PMAmr& pmamr;
    const Layout& layout;

    PArray<MFTower> coefs;
    RichardContext& richardContext;
    bool bc_initialized;
};

#endif
