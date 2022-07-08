#ifndef _MFTower_
#define _MFTower_


#include <winstd.H>
#include <MultiFab.H>
#include <PArray.H>
#include <Array.H>
#include <Geometry.H>
#include <BCRec.H>

struct Layout;
class Stencil;
struct IVSMap;

struct MFTower
{
    MFTower(const Layout&    layout,
            const IndexType& t = IndexType(IntVect::TheZeroVector()),
            int              nComp = 1,
            int              nGrow = 0,
	    int              numLevs = -1);

    MFTower(Layout&           layout,
            PArray<MultiFab>& pamf,
	    int               numLevs = -1);

    int NComp() const;
    int NGrow() const;
    int NumLevels() const;
    MultiFab& operator[](int i);
    const MultiFab& operator[](int i) const;
    bool IsCompatible(const MFTower& rhs) const;
    Real norm(int numLevs) const;
    const Layout& GetLayout() const;
    const IndexType& ixType() const;

    void SetValCovered(Real value);


    void Write(const std::string& fileName,
	       const std::string& names,
	       Real               time) const;

    static void WriteSet(const std::string&          fileName,
			 const Array<MFTower*> mfta,
			 const Array<std::string>&   names,
			 Real                        time);
 
    static void CCtoECgrad(PArray<MFTower>& mft_ec,
                           const MFTower&   mft_cc,
                           Real             mult,
                           int              sComp,
                           int              dComp,
                           int              nComp,
			   int              numLevs=-1);
    
    static void CCtoECavg(PArray<MFTower>& mft_ec,
                          const MFTower&   mft_cc,
                          Real             mult,
                          int              sComp,
                          int              dComp,
                          int              nComp,
                          int              do_harmonic=1,
                          int              numLevs=-1);
    
    static void ECtoCCdiv(MFTower&               mft_cc,
                          const PArray<MFTower>& mft_ec,
                          const Array<Real>&     mult,
                          int                    sComp,
                          int                    dComp,
                          int                    nComp,
			  int                    numLevs=-1);
    
    static void AverageDown(MFTower& mft,
                            int      sComp,
                            int      nComp,
			    int      numLevs=-1);
    
    void SetVal(Real     val,
		int      sComp=0,
		int      nComp=1,
		int      numLevs=-1);
  
    static IndexType EC[BL_SPACEDIM], CC, NC;

protected:
    void define_alloc();
    void define_noalloc(PArray<MultiFab>& pamf);
    const Layout& layout;
    int nLevs, nComp, nGrow;
    PArray<MultiFab> mft;
    IndexType iType;
};


struct MFTFillPatch
{
    MFTFillPatch(Layout& layout);
    int NumLevels() const;
    Layout& GetLayout();

    int BuildStencil(const BCRec& bc,
                     int          maxorder); // Returns 0 if successful, returns 1 if bc not properly defined

    void DoCoarseFineParallelInterp(MFTower& mft,
                                    int      sComp=1,
                                    int      nComp=1,
				    int      numLevs=-1) const;
    void FillGrowCells(MFTower& mft,
                       int      sComp=0,
                       int      nComp=1,
                       bool     do_piecewise_constant = false,
		       int      numLevs=-1) const;

    void FillGrowCellsSimple(MFTower& mft,
                             int      sComp=0,
                             int      nComp=1,
			     int      numLevs=-1) const;

    const Array<Array<IVSMap> >& GrowCellStencil() const {return perpInterpStencil;} //[lev][dir][iv] = stencil to fill grow cells

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
    Layout& layout;
    Array<Array<IVSMap> > parallelInterpStencil; //[lev][dir][iv] = stencil to perform parallel bc Interp
    Array<Array<IVSMap> > perpInterpStencil; //[lev][dir][iv] = stencil to perform perpendicular bc Interp
};

inline bool operator==(const MFTFillPatch::MyHash& lhs, const MFTFillPatch::MyHash& rhs);
inline bool operator!=(const MFTFillPatch::MyHash& lhs, const MFTFillPatch::MyHash& rhs);

#endif
