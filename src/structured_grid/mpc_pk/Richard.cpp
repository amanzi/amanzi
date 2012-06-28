#include <Richard.H>
#include <RICHARD_F.H>

#include <Utility.H>

#ifdef BL_USE_PETSC
#include <petscmat.h>
#define CHKPETSC(n) CHKERRABORT(ParallelDescriptor::Communicator(),n);
#endif

#include <ostream>
std::ostream& operator<< (std::ostream&  os, const MultiFab& mf)
{
  int myproc = ParallelDescriptor::MyProc();
    os << "MultiFab" << '\n';
    for (int i=0; i<mf.size(); ++i) {
      const FArrayBox& fab = mf[i];
      if (mf.DistributionMap()[i]==myproc) {
	const Box& box = fab.box();
	os << "Box: " << box << '\n';
	for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
	  {
	    os << iv << " ";
	    for (int n=0; n<fab.nComp(); ++n) {
	      os << fab(iv,n) << " ";
	    }
	    os << '\n';
	  } 
      }
      ParallelDescriptor::Barrier();
    }
    return os;
}


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
                             MFTower&               lambda)
{
    BL_ASSERT(layout.IsCompatible(pressure));
    BL_ASSERT(pressure.NGrow()>=1);
    BL_ASSERT(lambda.NGrow()>=pressure.NGrow());
    BL_ASSERT(saturation.NGrow()>=pressure.NGrow());

    mftfp.FillGrowCells(pressure);

    // Assumes lev=0 here corresponds to Amr.level=0
    for (int lev=0; lev<nLevs; ++lev) {
        PorousMedia* pmp = dynamic_cast<PorousMedia*>(&(pmamr.getLevel(lev)));
        if (!pmp) {
            BoxLib::Abort("Bad cast in CalcCoefficients::operator()");
        }
        MultiFab& pressureLev = pressure[lev];
        MultiFab& saturationLev = saturation[lev];
        MultiFab& lambdaLev = lambda[lev];
        pmp->calcInvPressure(saturationLev,pressureLev); // FIXME: Writes/reads only to comp=0
        pmp->calcLambda(&lambdaLev,saturationLev); // FIXME: Writes/reads only to comp=0
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

void
RichardOp::CenterToEdgeUpwind(PArray<MFTower>&       mfte,
                              MFTower&               mftc,
                              const PArray<MFTower>& sgn,
                              int                    nComp) const
{
    int nLevs = richardContext.GetLayout().NumLevels();
    for (int lev=0; lev<nLevs; ++lev) {
        MultiFab& clev = mftc[lev];
        BL_ASSERT(nComp<=clev.nComp());
        for (MFIter mfi(clev); mfi.isValid(); ++mfi) {
            FArrayBox& cfab = clev[mfi];
            const Box& vccbox = mfi.validbox();
            for (int d=0; d<BL_SPACEDIM; ++d) {            
                FArrayBox& efab = mfte[d][lev][mfi];
                const FArrayBox& sgnfab = sgn[d][lev][mfi];
                BL_ASSERT(nComp<=efab.nComp());
                BL_ASSERT(nComp<=sgnfab.nComp());
                BL_ASSERT(Box(vccbox).surroundingNodes(d).contains(efab.box()));
                BL_ASSERT(Box(vccbox).surroundingNodes(d).contains(sgnfab.box()));
                
                FORT_RICHARD_CTE_UPW(efab.dataPtr(), ARLIM(efab.loVect()), ARLIM(efab.hiVect()),
                                     cfab.dataPtr(), ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                                     sgnfab.dataPtr(),ARLIM(sgnfab.loVect()), ARLIM(sgnfab.hiVect()),
                                     vccbox.loVect(), vccbox.hiVect(), &d, &nComp);
            }
        }
    }
}

void 
RichardOp::YmultX(MFTower&           Y,
                  const MFTower&     X,
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

            FORT_RICHARD_YMULTX(Yfab.dataPtr(sComp),ARLIM(Yfab.loVect()), ARLIM(Yfab.hiVect()),
                                Xfab.dataPtr(dComp),ARLIM(Xfab.loVect()), ARLIM(Xfab.hiVect()),
                                gbox.loVect(), gbox.hiVect(), &nComp);
        }
    }
}



//
// Compute mfte[dir][comp] = Grad(mftc[comp]) + a[dir][comp]
//
void
RichardOp::CCtoECgradAdd(PArray<MFTower>& mfte,
                         const MFTower&   mftc,
                         const FArrayBox& a,
                         int              sComp,
                         int              dComp,
                         int              nComp) const
{
    const Layout& layout = richardContext.GetLayout();
    for (int d=0; d<BL_SPACEDIM; ++d) {            
        BL_ASSERT(layout.IsCompatible(mfte[d]));
    }
    BL_ASSERT(layout.IsCompatible(mftc));

    const Array<Geometry>& geomArray = layout.GeomArray();
    int nLevs = layout.NumLevels();
    for (int lev=0; lev<nLevs; ++lev) {
        const MultiFab& mfc = mftc[lev];
        BL_ASSERT(mfc.nGrow()>=1);
        BL_ASSERT(sComp+nComp<=mfc.nComp());
        const Real* dx = geomArray[lev].CellSize();

        for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
            const FArrayBox& cfab = mfc[mfi];
            const Box& vcbox = mfi.validbox();
            
            for (int d=0; d<BL_SPACEDIM; ++d) {            
                FArrayBox& efab = mfte[d][lev][mfi];
                BL_ASSERT(dComp+nComp<=efab.nComp());
                BL_ASSERT(Box(vcbox).surroundingNodes(d).contains(efab.box()));
                efab.setVal(0);
                FORT_RICHARD_GXPA(efab.dataPtr(dComp),ARLIM(efab.loVect()), ARLIM(efab.hiVect()),
                                  cfab.dataPtr(sComp),ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                                  vcbox.loVect(),vcbox.hiVect(),dx,a.dataPtr(d),&d,&nComp);
            }
        }
    }
}

void
RichardContext::FillPatch(MFTower& mft,
                          bool     do_piecewise_constant)
{
    mftfp(mft,do_piecewise_constant);
}

void
RichardOp::FillPatch(MFTower& mft,
                     bool do_piecewise_constant)
{
    richardContext.FillPatch(mft,do_piecewise_constant);
}

void
RichardOp::Residual(MFTower&   residual,
                    MFTower&   pressure,
                    const Real dt)
{
    // Grab refs to temporary data
    MFTower& pressure_old = richardContext.PressureOld();
    MFTower& saturation = richardContext.Saturation();
    MFTower& lambda = richardContext.Lambda();
    PArray<MFTower>& velocity = richardContext.DarcyVelocity();
    CalcCoefficients& calcCoefs = richardContext.CalcCoefs();
    int nComp = 1;
    int sComp = 0; 
    int dComp = 0;
    Box abox(IntVect::TheZeroVector(),(nComp-1)*BoxLib::BASISV(0));
    FArrayBox a(abox,BL_SPACEDIM); // Make a funny box for a to simplify passing to Fortran
    for (int d=0; d<BL_SPACEDIM; ++d) {
        Real* ap = a.dataPtr(d);
        for (int n=0; n<nComp; ++n) {
            ap[n] = richardContext.Density()[n] * richardContext.Gravity()[d];
        }
    }
    
    // Update cell-centered coefficients
    calcCoefs(coefs,pressure,saturation,lambda);

    // Get  Grad(p) + rho.g
    FillPatch(pressure);
    CCtoECgradAdd(velocity,pressure,a,sComp,dComp,nComp);

    // Get edge-centered coefficients based on the sign of (Grad(p) - rho.g)
    CenterToEdgeUpwind(coefs,lambda,velocity,nComp);

    // Get Darcy flux = H * (Grad(p) - rho.g)
    for (int d=0; d<BL_SPACEDIM; ++d) {
        BL_ASSERT(velocity[d].NComp()>dComp);
        YmultX(velocity[d],coefs[d]);
    }

    // Average down Darcy flux to coarse faces

    // Get the divergence of the Darcy Flux
    MFTower::ECtoCCdiv(residual,velocity);
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
