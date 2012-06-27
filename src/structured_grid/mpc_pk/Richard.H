#ifndef _Richard_H_
#define _Richard_H_


#include <winstd.H>

#include <Layout.H>

class MFTower;
class MFTFillPatch;

#ifdef BL_USE_PETSC
#include <petsc.h>
#endif

class PMAmr;
#include <PMAmr.H>
#include <PorousMedia.H>

class CalcCoefficients
{
public:
    CalcCoefficients(MFTFillPatch& mftfp,
                     PMAmr&        pmamr);

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
    MFTFillPatch& mftfp;
    const Layout& layout;
    PMAmr& pmamr;
};

struct RichardContext
{
    RichardContext(PMAmr&           pmamr,
                   MFTFillPatch&    mftfp,
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
    const Array<Real>& Density() {return density;}
    const Array<Real>& Gravity() {return gravity;}
    MFTFillPatch& FillPatch() {return mftfp;}
    const Layout& GetLayout() {return mftfp.GetLayout();}

protected:
    MFTower& pressure_old;
    MFTower& saturation;
    MFTower& lambda;
    PArray<MFTower>& darcyVelocity;
    Array<Real> gravity;
    const Array<Real>& density;

    PMAmr& pmamr;
    MFTFillPatch& mftfp;
    CalcCoefficients* calcCoefs;
};

class RichardOp
{
public:
    RichardOp(RichardContext& richardContext);
    ~RichardOp();

    const Layout& GetLayout() const {return richardContext.GetLayout();}
    const Array<BoxArray>& GridArray() const {return GetLayout().GridArray();}
    const Array<Geometry>& GeomArray() const {return GetLayout().GeomArray();}
    const Array<IntVect>& RefRatio() const {return GetLayout().RefRatio();}

    void SetBC(const BCRec& bc, 
               int          maxorder); // Note: triggers build of mftfp stencil

    void Residual(MFTower&               residual,
                  MFTower&               pressure,
                  const Real             dt);

    void DarcyVelocity(PArray<MFTower>& velocity,
                       MFTower&         pressure,
                       int              sComp = 0,
                       int              dComp = 0,
                       int              nComp=1) const;


    const MFTower& Coefficients(int dir) const {return coefs[dir];}
    MFTower& Coefficients(int dir) {return coefs[dir];}

    void Write(const std::string& fileName) const;

protected:

    // Compute Y = (Y + a)*bX, where a,b are per component, and X,Y are vectors and the * is element-wisevoid 
    void YpambX(MFTower&           Y,
                const MFTower&     X,
                const Array<Real>& a,
                const Array<Real>& b,
                int                sComp,
                int                dComp,
                int                nComp,
                int                nGrow) const;
    void AllocLambda();
    void AllocSaturation();

    PMAmr& pmamr;
    PArray<MFTower> coefs;
    RichardContext& richardContext;
    bool bc_initialized;
};

#endif
