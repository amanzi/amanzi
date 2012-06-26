#include <Richard.H>
#include <RICHARD_F.H>

#include <Utility.H>

#ifdef BL_USE_PETSC
#include <petscmat.h>
#define CHKPETSC(n) CHKERRABORT(ParallelDescriptor::Communicator(),n);
#endif


CalcCoefficients::CalcCoefficients(MFTFillPatch& _mftfp,
                                   PMAmr&        _pmamr)
    : mftfp(_mftfp),
      pmamr(_pmamr),
      layout(_mftfp.GetLayout()),
      nLevs(_mftfp.GetLayout().NumLevels())
{
}


void
CalcCoefficients::operator()(PArray<MFTower>&       coefficients,
                             MFTower&               pressure,
                             MFTower&               saturation,
                             MFTower&               lambda,
                             const PArray<MFTower>& DarcyVelocity,
                             Real                   mult,
                             int                    presComp,
                             int                    satComp,
                             int                    lamComp,
                             int                    darcComp,
                             int                    coefComp,
                             int                    nComp)
{
    BL_ASSERT(layout.IsCompatible(pressure));

    // This will set physbc and c-f grow cells without interpolation (ie, "piecewise-constant interp")
    // For c-f: copy the underlaying coarse cell value
    // For phys: the Dirichlet value (should already be there on entry to this routine)
    mftfp.FillGrowCellsSimple(pressure,presComp,nComp);

    // Assumes lev=0 here corresponds to Amr.level=0
    for (int lev=0; lev<nLevs; ++lev) {
        PorousMedia* pmp = dynamic_cast<PorousMedia*>(&(pmamr.getLevel(lev)));
        if (!pmp) {
            BoxLib::Abort("Bad cast in CalcCoefficients::operator()");
        }
        MultiFab& pressureLev = pressure[lev];
        MultiFab& saturationLev = saturation[lev];
        MultiFab& lambdaLev = saturation[lev];
        pmp->calcInvPressure(saturationLev,pressureLev); // FIXME: Writes/reads only to comp=0
        pmp->calcLambda(&lambdaLev,saturationLev); // FIXME: Writes/reads only to comp=0
    }

    for (int lev=0; lev<nLevs; ++lev) {
        MultiFab& lamlev = lambda[lev];
        for (MFIter mfi(lamlev); mfi.isValid(); ++mfi) {
            FArrayBox& lamfab = lamlev[mfi];
            const Box& vccbox = mfi.validbox();
            for (int d=0; d<BL_SPACEDIM; ++d) {            
                FArrayBox& coefab = coefficients[d][lev][mfi];
                const FArrayBox& darfab = DarcyVelocity[d][lev][mfi];
                BL_ASSERT(coefComp+nComp<=coefab.nComp());
                BL_ASSERT(darcComp+nComp<=darfab.nComp());
                BL_ASSERT(Box(vccbox).surroundingNodes(d).contains(coefab.box()));
                BL_ASSERT(Box(vccbox).surroundingNodes(d).contains(darfab.box()));

                FORT_CC_TO_EC_UPW(coefab.dataPtr(coefComp),ARLIM(coefab.loVect()), ARLIM(coefab.hiVect()),
                                  lamfab.dataPtr(lamComp), ARLIM(lamfab.loVect()), ARLIM(lamfab.hiVect()),
                                  darfab.dataPtr(darcComp),ARLIM(darfab.loVect()), ARLIM(darfab.hiVect()),
                                  vccbox.loVect(), vccbox.hiVect(), &mult, &d, &nComp);
            }
        }
    }
}

RichardContext::RichardContext(PMAmr&           _pmamr,
                               MFTFillPatch&     _mftfp,
                               MFTower&         _pressure_old,
                               MFTower&         _saturation,
                               MFTower&         _lambda,
                               PArray<MFTower>& _darcyVelocity,
                               const BCRec&     pressure_bc,
                               int              pressure_maxorder)
    : pmamr(_pmamr), 
      mftfp(_mftfp),
      pressure_old(_pressure_old),
      saturation(_saturation),
      lambda(_lambda),
      darcyVelocity(_darcyVelocity),
      density(PorousMedia::Density()),
      gravity(BL_SPACEDIM,0)
{
    gravity[BL_SPACEDIM-1] = PorousMedia::getGravity();
    mftfp.BuildStencil(pressure_bc, pressure_maxorder);
    calcCoefs = new CalcCoefficients(mftfp,pmamr);
}

RichardOp::RichardOp(RichardContext& _richardContext)
    : richardContext(_richardContext),
      pmamr(_richardContext.PMAMR()),
      coefs(BL_SPACEDIM,PArrayManage),
      bc_initialized(false)
{
    Array<IndexType> itype(BL_SPACEDIM);
    const Layout& layout = GetLayout(); 
    for (int d=0; d<BL_SPACEDIM; ++d) {
        coefs.set(d, new MFTower(layout,MFTower::EC[d]));
    }
}

RichardOp::~RichardOp()
{
}

// Compute Y = (Y + a)*bX, where a,b are per component, and X,Y are vectors and the * is element-wise
void 
RichardOp::YpambX(MFTower&           Y,
                  const MFTower&     X,
                  const Array<Real>& a,
                  const Array<Real>& b,
                  int                sComp,
                  int                dComp,
                  int                nComp,
                  int                nGrow) const
{
    const Layout& layout = GetLayout(); 
    BL_ASSERT(layout.IsCompatible(X));
    BL_ASSERT(layout.IsCompatible(Y));
    BL_ASSERT(X.NComp()>=sComp+nComp);
    BL_ASSERT(Y.NComp()>=dComp+nComp);
    BL_ASSERT(X.NGrow()>=nGrow);
    BL_ASSERT(Y.NGrow()>=nGrow);
    BL_ASSERT(a.size()==0 || a.size()>=nComp);
    BL_ASSERT(b.size()==0 || b.size()>=nComp);

    const Array<BoxArray>& gridArray = GridArray();
    const Array<Geometry>& geomArray = GeomArray();
    const Array<IntVect>& refRatio = RefRatio();

    FArrayBox fab;
    int nLevs = layout.NumLevels();
    for (int lev=0; lev<nLevs; ++lev)
    {
        MultiFab& Ylev = Y[lev];
        const MultiFab& Xlev = X[lev];
        const BoxArray& ba = Xlev.boxArray();
        BoxArray fba;
        if (lev<nLevs-1) {
            fba = BoxArray(ba).coarsen(refRatio[lev]);
        }

        for (MFIter mfi(Y[lev]); mfi.isValid(); ++mfi)
        {
            const Box& vbox = mfi.validbox();
            Box gbox = Box(vbox).grow(nGrow);
            FArrayBox& Yfab = Ylev[mfi];
            const FArrayBox& Xfab = Xlev[mfi];

            if (a.size()) {
                for (int n=0; n<nComp; ++n) {
                    Yfab.plus(a[n],dComp+n,1);
                }
            }

            fab.resize(gbox,nComp);
            fab.copy(Xfab,sComp,0,nComp);

            // Zero out parts of X covered by fine grid (to ensure valid data)
            if (lev<nLevs-1) {
                std::vector< std::pair<int,Box> > isects = fba.intersections(gbox);
                for (int i = 0; i < isects.size(); i++)
                {
                    fab.setVal(0,isects[i].second,0,nComp);
                }
            }

            if (b.size()) {
                for (int n=0; n<nComp; ++n) {
                    fab.mult(b[n],n,1);
                }
            }
            Yfab.plus(fab,0,dComp,nComp);
        }
    }
}



//
// Flux = H(Grad(P) - rho.g.zhat)
//
void
RichardOp::DarcyVelocity(PArray<MFTower>& velocity,
                         MFTower&         pressure,
                         int              sComp,
                         int              dComp,
                         int              nComp) const
{
    // Note: Assumes pressure[sComp] corresponds to the density[0]
    BL_ASSERT(pressure.NComp()>=sComp+nComp);
    BL_ASSERT(nComp<=richardContext.Density().size());
    Array<Array<Real> > rhog(BL_SPACEDIM,Array<Real>(nComp));
    for (int d=0; d<BL_SPACEDIM; ++d) {
        for (int n=0; n<nComp; ++n) {
            rhog[d][n] = richardContext.Density()[n] * richardContext.Gravity()[d];
        }
    }
    Array<Real> ones(nComp,1);
    Real mult = 1;

    MFTFillPatch& mftfp = richardContext.FillPatch();
    mftfp.FillGrowCells(pressure,sComp,nComp);
    mftfp.CCtoECgrad(velocity,pressure,mult,sComp,dComp,nComp);

    // Compute vel = H(Grad(p) - rho.g)
    for (int d=0; d<BL_SPACEDIM; ++d) {
        BL_ASSERT(velocity[d].NComp()>dComp);
        YpambX(velocity[d],coefs[d],rhog[d],ones,0,0,1,0);
    }
}

void
RichardOp::Residual(MFTower&               residual,
                    MFTower&               pressure,
                    const Real             dt)
{
    // Grab refs to temporary data
    MFTower& pressure_old = richardContext.PressureOld();
    MFTower& saturation = richardContext.Saturation();
    MFTower& lambda = richardContext.Lambda();
    PArray<MFTower>& darcyVelocity = richardContext.DarcyVelocity();
    CalcCoefficients& calcCoefs = richardContext.CalcCoefs();
    
    // Update coefficients
    calcCoefs(coefs,pressure,saturation,lambda,darcyVelocity);

    // Compute Darcy velocity
    DarcyVelocity(darcyVelocity,pressure);

    // Compute the divergence
    MFTFillPatch& mftfp = richardContext.FillPatch();
    mftfp.ECtoCCdiv(residual,darcyVelocity);
}

void
RichardOp::Write(const std::string& fileName) const
{
    std::string FullPath = fileName;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
        FullPath += '/';
    
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
            BoxLib::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
    
    const Layout& layout = richardContext.GetLayout();
    int nLevs = layout.NumLevels();
    for (int lev=0; lev<nLevs; ++lev) {
        std::string FullPath = fileName;
        if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
        {
            FullPath += '/';
        }
        std::string LevelPath = FullPath + "RichardOp_Level_";
        LevelPath = BoxLib::Concatenate(LevelPath,lev,1);
        
        if (ParallelDescriptor::IOProcessor()) {
            if (!BoxLib::UtilCreateDirectory(LevelPath, 0755)) {
                BoxLib::CreateDirectoryFailed(LevelPath);
            }
        }
        ParallelDescriptor::Barrier();

        LevelPath += '/';
        for (int d=0; d<BL_SPACEDIM; ++d) {
            std::string DirLevelPathCoefs = LevelPath + "coefs_";
            DirLevelPathCoefs = BoxLib::Concatenate(DirLevelPathCoefs,d,1);
            VisMF::Write(coefs[d][lev],DirLevelPathCoefs);
        }
        std::string LevelPathPressureOld = LevelPath + "pressure_old";
        VisMF::Write(richardContext.PressureOld()[lev],LevelPathPressureOld);
        
        std::string LevelPathLambda = LevelPath + "lambda";
        VisMF::Write(richardContext.Lambda()[lev],LevelPathLambda);
        
        std::string LevelPathSaturation = LevelPath + "saturation";
        VisMF::Write(richardContext.Saturation()[lev],LevelPathSaturation);
        
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::string MiscPath = FullPath + "miscData.dat";
        std::ofstream os; os.open(MiscPath.c_str());
        os << "Number of Levels: " << nLevs << '\n';
        os << "Refinement Ratios: ";
        for (int lev=1; lev<nLevs; ++lev) {
            os << layout.RefRatio()[lev-1] << " ";
        }
        os << '\n';
        os << "Grid Arrays: " << '\n';
        for (int lev=0; lev<nLevs; ++lev) {
            os << layout.GridArray()[lev] << '\n';
        }
        os.close();
    }
}
