#include <winstd.H>

#include <PressBndry.H>
#include <MacOperator.H>
#include <MacOpMacDrivers.H>
#include <MACOPERATOR_F.H>
#include <CGSolver.H>
#include <MultiGrid.H>
#include <DataServices.H>

#ifdef MG_USE_HYPRE
#include <HypreABec.H>
#endif

#ifdef MG_USE_FBOXLIB
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

MacOperator::MacOperator (const BndryData& mgb,
                          const Real*      h)
    :
    ABecLaplacian(mgb,h)
{}

MacOperator::~MacOperator () {}

//
// Define the meaning of gradient for the multigrid object.
//

void
MacOperator::setCoefficients (MultiFab*    area,
                              MultiFab*    mac_coef,
                              const Real*  dx,
		              const BCRec& phys_bc,
		              const Box&   domain)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    const BoxArray& ba = gbox[0];

    //
    // First set scalar coeficients.
    //
    setScalars(0.0,1.0);
    //
    // Don't need to set a because alpha is set to zero.
    //
    const int n_grow = 0;

    D_TERM(MultiFab bxcoef(area[0].boxArray(),area[0].nComp(),n_grow);,
           MultiFab bycoef(area[1].boxArray(),area[1].nComp(),n_grow);,
           MultiFab bzcoef(area[2].boxArray(),area[2].nComp(),n_grow););
    D_TERM(bxcoef.setVal(0);,
           bycoef.setVal(0);,
           bzcoef.setVal(0););

    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();

    for (MFIter mfi(mac_coef[0]); mfi.isValid(); ++mfi)
    {
        const Box& grd       = ba[mfi.index()];
        const int* lo        = grd.loVect();
        const int* hi        = grd.hiVect();
        FArrayBox& bx        = bxcoef[mfi];
        FArrayBox& by        = bycoef[mfi];
        const FArrayBox& lx  = mac_coef[0][mfi];
        const FArrayBox& ly  = mac_coef[1][mfi];
        const FArrayBox& ax  = area[0][mfi];
        const FArrayBox& ay  = area[1][mfi];

        DEF_LIMITS(bx,bx_dat,bxlo,bxhi);
        DEF_LIMITS(by,by_dat,bylo,byhi);
        DEF_CLIMITS(lx,lx_dat,lxlo,lxhi);
        DEF_CLIMITS(ly,ly_dat,lylo,lyhi);
        DEF_CLIMITS(ax,ax_dat,axlo,axhi);
        DEF_CLIMITS(ay,ay_dat,aylo,ayhi);

        const int* domlo    = domain.loVect();
        const int* domhi    = domain.hiVect();

#if (BL_SPACEDIM == 2)
        FORT_MACCOEF(bx_dat,ARLIM(bxlo),ARLIM(bxhi),
                     by_dat,ARLIM(bylo),ARLIM(byhi),
                     lx_dat,ARLIM(lxlo),ARLIM(lxhi),
                     ly_dat,ARLIM(lylo),ARLIM(lyhi),
                     ax_dat,ARLIM(axlo),ARLIM(axhi),
                     ay_dat,ARLIM(aylo),ARLIM(ayhi),
                     lo,hi,domlo,domhi,dx,lo_bc,hi_bc);
#endif
#if (BL_SPACEDIM == 3)
        FArrayBox& bz       = bzcoef[mfi];
	const FArrayBox& lz = mac_coef[2][mfi];
        const FArrayBox& az = area[2][mfi];

	DEF_CLIMITS(lz,lz_dat,lzlo,lzhi);
        DEF_CLIMITS(az,az_dat,azlo,azhi);
        DEF_LIMITS(bz,bz_dat,bzlo,bzhi);

        FORT_MACCOEF(bx_dat,ARLIM(bxlo),ARLIM(bxhi),
                     by_dat,ARLIM(bylo),ARLIM(byhi),
                     bz_dat,ARLIM(bzlo),ARLIM(bzhi),
                     lx_dat,ARLIM(lxlo),ARLIM(lxhi),
                     ly_dat,ARLIM(lylo),ARLIM(lyhi),
                     lz_dat,ARLIM(lzlo),ARLIM(lzhi),
                     ax_dat,ARLIM(axlo),ARLIM(axhi),
                     ay_dat,ARLIM(aylo),ARLIM(ayhi),
                     az_dat,ARLIM(azlo),ARLIM(azhi),
                     lo,hi,domlo,domhi,dx,lo_bc,hi_bc);
#endif
    }
    D_TERM(bCoefficients(bxcoef,0);,
           bCoefficients(bycoef,1);,
           bCoefficients(bzcoef,2););
}

//
// This function creates the initial rhs for use in the mac multgrid solve.
//

void
MacOperator::defRHS (MultiFab*    area,
                     MultiFab&    volume,
                     MultiFab&    Rhs,
                     MultiFab*    vel,
		     MultiFab*    RhoD,
                     const BCRec& phys_bc,
                     const Box&   domain,
                     Real         scale)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    const BoxArray& ba = gbox[0];
    BL_ASSERT(Rhs.boxArray() == ba);

    const int* domlo    = domain.loVect();
    const int* domhi    = domain.hiVect();

    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();

    for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
        BL_ASSERT(ba[Rhsmfi.index()] == Rhsmfi.validbox());

        const Box& grd       = Rhsmfi.validbox();
        const int* lo        = grd.loVect();
        const int* hi        = grd.hiVect();
        const FArrayBox& ax  = area[0][Rhsmfi];
        const FArrayBox& ay  = area[1][Rhsmfi];
        const FArrayBox& vol = volume[Rhsmfi];
        const FArrayBox& ux  = vel[0][Rhsmfi];
        const FArrayBox& uy  = vel[1][Rhsmfi];
        const FArrayBox& rx  = RhoD[0][Rhsmfi];
        const FArrayBox& ry  = RhoD[1][Rhsmfi];
        FArrayBox& rhs       = Rhs[Rhsmfi];

        DEF_CLIMITS(ux,ux_dat,uxlo,uxhi);
        DEF_CLIMITS(uy,uy_dat,uylo,uyhi);
        DEF_CLIMITS(rx,rx_dat,rxlo,rxhi);
        DEF_CLIMITS(ry,ry_dat,rylo,ryhi);
        DEF_CLIMITS(ax,ax_dat,axlo,axhi);
        DEF_CLIMITS(ay,ay_dat,aylo,ayhi);
        DEF_CLIMITS(vol,vol_dat,vlo,vhi);
        DEF_LIMITS(rhs,rhs_dat,rlo,rhi);
      
#if (BL_SPACEDIM == 2)
        FORT_MACRHS(ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
		    rx_dat,ARLIM(rxlo),ARLIM(rxhi),
                    ry_dat,ARLIM(rylo),ARLIM(ryhi),
                    ax_dat,ARLIM(axlo),ARLIM(axhi),
                    ay_dat,ARLIM(aylo),ARLIM(ayhi),
                    vol_dat,ARLIM(vlo),ARLIM(vhi), 
                    rhs_dat,ARLIM(rlo),ARLIM(rhi),
                    lo,hi,domlo,domhi,lo_bc,hi_bc,&scale);
#endif
#if (BL_SPACEDIM == 3)
        const FArrayBox& az = area[2][Rhsmfi];
        DEF_CLIMITS(az,az_dat,azlo,azhi);

        const FArrayBox& uz = vel[2][Rhsmfi];
        DEF_CLIMITS(uz,uz_dat,uzlo,uzhi);

        const FArrayBox& rz = RhoD[2][Rhsmfi];
        DEF_CLIMITS(rz,rz_dat,rzlo,rzhi);

        FORT_MACRHS(ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
                    uz_dat,ARLIM(uzlo),ARLIM(uzhi),
		    rx_dat,ARLIM(rxlo),ARLIM(rxhi),
                    ry_dat,ARLIM(rylo),ARLIM(ryhi),
                    rz_dat,ARLIM(rzlo),ARLIM(rzhi),
                    ax_dat,ARLIM(axlo),ARLIM(axhi),
                    ay_dat,ARLIM(aylo),ARLIM(ayhi),
                    az_dat,ARLIM(azlo),ARLIM(azhi),
                    vol_dat,ARLIM(vlo),ARLIM(vhi),
                    rhs_dat,ARLIM(rlo),ARLIM(rhi),
                    lo,hi,domlo,domhi,lo_bc,hi_bc,&scale);

#endif
    }
    Rhs.mult(-1.0,Rhs.nGrow());
}

//
// Apply the mac pressure gradient to a velocity field.
// init, means that velocities are initialized here.
//

void
mac_vel_update (int              init,
                D_DECL(FArrayBox& ux,
                       FArrayBox& uy,
                       FArrayBox& uz),
                D_DECL(const FArrayBox& lx,
                       const FArrayBox& ly,
                       const FArrayBox& lz),
		const FArrayBox& phi,
                const BCRec&     phys_bc,
                const Box&       grd,
                const Real*      dx,
                const Box&       domain,
                Real             scale)
{
    const int* lo    = grd.loVect();
    const int* hi    = grd.hiVect();

    const int* domlo = domain.loVect();
    const int* domhi = domain.hiVect();

    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();

    DEF_LIMITS(ux,ux_dat,uxlo,uxhi);
    DEF_LIMITS(uy,uy_dat,uylo,uyhi);
    DEF_CLIMITS(lx,lx_dat,lxlo,lxhi);
    DEF_CLIMITS(ly,ly_dat,lylo,lyhi);
    DEF_CLIMITS(phi,phi_dat,p_lo,p_hi);

    
#if (BL_SPACEDIM == 2)
    FORT_MACUPDATE(&init,
                   ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
		   lx_dat,ARLIM(lxlo),ARLIM(lxhi),
                   ly_dat,ARLIM(lylo),ARLIM(lyhi),
                   phi_dat,ARLIM(p_lo),ARLIM(p_hi),
                   lo,hi,domlo,domhi,dx,lo_bc,hi_bc,&scale);
#endif
#if (BL_SPACEDIM == 3)
    DEF_LIMITS(uz,uz_dat,uzlo,uzhi);
    DEF_CLIMITS(lz,lz_dat,lzlo,lzhi);

    FORT_MACUPDATE(&init,
                   ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
                   uz_dat,ARLIM(uzlo),ARLIM(uzhi),
		   lx_dat,ARLIM(lxlo),ARLIM(lxhi),
                   ly_dat,ARLIM(lylo),ARLIM(lyhi), 
		   lz_dat,ARLIM(lzlo),ARLIM(lzhi),
                   phi_dat,ARLIM(p_lo),ARLIM(p_hi),
                   lo,hi,domlo,domhi,dx,lo_bc,hi_bc,&scale);
#endif
}

void
mac_vel_sync (D_DECL(FArrayBox& ux,
		     FArrayBox& uy,
		     FArrayBox& uz),
	      D_DECL(const FArrayBox& rhodx,
		     const FArrayBox& rhody,
		     const FArrayBox& rhodz),
	      D_DECL(const FArrayBox& areax,
		     const FArrayBox& areay,
		     const FArrayBox& areaz),
	      D_DECL(const FArrayBox& lx,
		     const FArrayBox& ly,
		     const FArrayBox& lz),
	      const FArrayBox& phi,
	      const Box&       grd,
	      const Real*      dx,
	      Real             scale)
{
    const int* lo    = grd.loVect();
    const int* hi    = grd.hiVect();
    
    DEF_LIMITS(ux,ux_dat,uxlo,uxhi);
    DEF_LIMITS(uy,uy_dat,uylo,uyhi);
    DEF_CLIMITS(rhodx,rx_dat,rxlo,rxhi);
    DEF_CLIMITS(rhody,ry_dat,rylo,ryhi);
    DEF_CLIMITS(areax,ax_dat,axlo,axhi);
    DEF_CLIMITS(areay,ay_dat,aylo,ayhi);
    DEF_CLIMITS(lx,lx_dat,lxlo,lxhi);
    DEF_CLIMITS(ly,ly_dat,lylo,lyhi);
    DEF_CLIMITS(phi,phi_dat,p_lo,p_hi);
    
#if (BL_SPACEDIM == 2)
    FORT_MAC_SYNC( ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
		   rx_dat,ARLIM(rxlo),ARLIM(rxhi),
                   ry_dat,ARLIM(rylo),ARLIM(ryhi),
		   ax_dat,ARLIM(axlo),ARLIM(axhi),
                   ay_dat,ARLIM(aylo),ARLIM(ayhi),
		   lx_dat,ARLIM(lxlo),ARLIM(lxhi),
                   ly_dat,ARLIM(lylo),ARLIM(lyhi),
                   phi_dat,ARLIM(p_lo),ARLIM(p_hi),
                   lo,hi,dx,&scale);
#else

    DEF_LIMITS(uz,uz_dat,uzlo,uzhi);
    DEF_CLIMITS(rhodz,rz_dat,rzlo,rzhi);
    DEF_CLIMITS(areaz,az_dat,azlo,azhi);
    DEF_CLIMITS(lz,lz_dat,lzlo,lzhi);

    FORT_MAC_SYNC( ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
                   uz_dat,ARLIM(uzlo),ARLIM(uzhi),
		   rx_dat,ARLIM(rxlo),ARLIM(rxhi),
                   ry_dat,ARLIM(rylo),ARLIM(ryhi),
		   rz_dat,ARLIM(rzlo),ARLIM(rzhi),
		   ax_dat,ARLIM(axlo),ARLIM(axhi),
                   ay_dat,ARLIM(aylo),ARLIM(ayhi),
		   az_dat,ARLIM(azlo),ARLIM(azhi),
		   lx_dat,ARLIM(lxlo),ARLIM(lxhi),
                   ly_dat,ARLIM(lylo),ARLIM(lyhi),
		   lz_dat,ARLIM(lzlo),ARLIM(lzhi),
                   phi_dat,ARLIM(p_lo),ARLIM(p_hi),
                   lo,hi,dx,&scale);
#endif
}

void
mac_vel_norml1 (Real* normval,
                D_DECL(FArrayBox& ux,
                       FArrayBox& uy,
                       FArrayBox& uz),
                const Box&       grd,
                const Real*      dx)
{
    const int* lo    = grd.loVect();
    const int* hi    = grd.hiVect();

 
    DEF_LIMITS(ux,ux_dat,uxlo,uxhi);
    DEF_LIMITS(uy,uy_dat,uylo,uyhi);
    
#if (BL_SPACEDIM == 2)
    FORT_MACNORML1(normval,ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
                   lo,hi,dx);
#endif
#if (BL_SPACEDIM == 3)
    DEF_LIMITS(uz,uz_dat,uzlo,uzhi);
    
    FORT_MACNORML1(normval,
	           ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
                   uz_dat,ARLIM(uzlo),ARLIM(uzhi),
                   lo,hi,dx);
#endif
}

//
// Apply the mac pressure gradient to the divergent mac velocities.
// The resultant velocity field is nondivergent.
//

void
MacOperator::velUpdate (MultiFab*       Vel,
                        MultiFab&       Phi,
                        const MultiFab* mac_coef,
                        const BCRec&    phys_bc,
                        const Real*     dx,
                        const Box&      domain,
			int             level,
                        Real            scale)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    const BoxArray& ba = gbox[0];

    //
    // Set bndry data in ghost zones.
    //
    int apply_lev = 0;//((level == 0))? 0 : level;
    applyBC(Phi,0,1,apply_lev);
    for (MFIter Phimfi(Phi); Phimfi.isValid(); ++Phimfi)
    {
        BL_ASSERT(ba[Phimfi.index()] == Phimfi.validbox());

        const Box& grd = Phimfi.validbox();

        mac_vel_update(0, 
                       D_DECL(Vel[0][Phimfi],Vel[1][Phimfi],Vel[2][Phimfi]),
		       D_DECL(mac_coef[0][Phimfi],mac_coef[1][Phimfi],mac_coef[2][Phimfi]),
                       Phi[Phimfi],
                       phys_bc, grd, dx, domain, scale );
    }
}

//
// Multiply by volume*rhs_scale since reflux step (which computed rhs)
// divided by volume.
//

void
MacOperator::syncRhs (const MultiFab& Volume,
                      MultiFab&       Rhs,
		      MultiFab*       u_mac_sync,
		      MultiFab*       rhoD_sync,
		      MultiFab*       area,
                      Real            rhs_scale,
		      const BCRec&    phys_bc,
		      const Box&      domain)
{
    const BoxArray& ba = gbox[0];

    const int* domlo    = domain.loVect();
    const int* domhi    = domain.hiVect();

    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();

    for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
        BL_ASSERT(ba[Rhsmfi.index()] == Rhsmfi.validbox());

        const Box& grd       = Rhsmfi.validbox();
        const int* lo        = grd.loVect();
        const int* hi        = grd.hiVect();
        FArrayBox& rhs       = Rhs[Rhsmfi];
        const FArrayBox& vol = Volume[Rhsmfi];
        const FArrayBox& ux  = u_mac_sync[0][Rhsmfi];
        const FArrayBox& uy  = u_mac_sync[1][Rhsmfi];
        const FArrayBox& rx  = rhoD_sync[0][Rhsmfi];
        const FArrayBox& ry  = rhoD_sync[1][Rhsmfi];
        const FArrayBox& ax  = area[0][Rhsmfi];
        const FArrayBox& ay  = area[1][Rhsmfi];

        DEF_CLIMITS(vol,vol_dat,vlo,vhi);
	DEF_CLIMITS(ux,ux_dat,uxlo,uxhi);
        DEF_CLIMITS(uy,uy_dat,uylo,uyhi);
	DEF_CLIMITS(rx,rx_dat,rxlo,rxhi);
        DEF_CLIMITS(ry,ry_dat,rylo,ryhi);
	DEF_CLIMITS(ax,ax_dat,axlo,axhi);
        DEF_CLIMITS(ay,ay_dat,aylo,ayhi);
	DEF_LIMITS(rhs,rhs_dat,rlo,rhi);

#if (BL_SPACEDIM == 2)
        FORT_MACSYNCRHS(rhs_dat,ARLIM(rlo),ARLIM(rhi),
			ux_dat,ARLIM(uxlo),ARLIM(uxhi),
			uy_dat,ARLIM(uylo),ARLIM(uyhi),
			rx_dat,ARLIM(rxlo),ARLIM(rxhi),
			ry_dat,ARLIM(rylo),ARLIM(ryhi),
			ax_dat,ARLIM(axlo),ARLIM(axhi),
			ay_dat,ARLIM(aylo),ARLIM(ayhi),
			vol_dat,ARLIM(vlo),ARLIM(vhi),
			lo,hi,
			domlo,domhi,lo_bc,hi_bc,&rhs_scale);
	
#else 
	const FArrayBox& uz  = u_mac_sync[2][Rhsmfi];
        const FArrayBox& rz  = rhoD_sync[2][Rhsmfi];
        const FArrayBox& az  = area[2][Rhsmfi];
        DEF_CLIMITS(uz,uz_dat,uzlo,uzhi);
        DEF_CLIMITS(rz,rz_dat,rzlo,rzhi);
        DEF_CLIMITS(az,az_dat,azlo,azhi);

        FORT_MACSYNCRHS(rhs_dat,ARLIM(rlo),ARLIM(rhi),
			ux_dat,ARLIM(uxlo),ARLIM(uxhi),
			uy_dat,ARLIM(uylo),ARLIM(uyhi),
			uz_dat,ARLIM(uzlo),ARLIM(uzhi),
			rx_dat,ARLIM(rxlo),ARLIM(rxhi),
			ry_dat,ARLIM(rylo),ARLIM(ryhi),
			rz_dat,ARLIM(rzlo),ARLIM(rzhi),
			ax_dat,ARLIM(axlo),ARLIM(axhi),
			ay_dat,ARLIM(aylo),ARLIM(ayhi),
			az_dat,ARLIM(azlo),ARLIM(azhi),
			vol_dat,ARLIM(vlo),ARLIM(vhi),
			lo,hi,
			domlo,domhi,lo_bc,hi_bc,&rhs_scale);
#endif
    }
    
    //Rhs.mult(-1.0,Rhs.nGrow());    
}

//
// Driver functions follow.
//

//
// A driver function for computing a level MAC solve.
//

void
mac_level_driver (const PressBndry& mac_bndry,
		  const BCRec&    phys_bc,
                  const BoxArray& grids,
                  int             the_solver,
                  int             level,
                  const Real*     dx,
                  Real            mac_tol,
                  Real            mac_abs_tol,
                  Real            rhs_scale,
                  MultiFab*       area,
                  MultiFab&       volume,
                  MultiFab*       mac_coef,
                  MultiFab&       Rhs,
                  MultiFab*       u_mac,
		  MultiFab*       RhoD,
                  MultiFab*       mac_phi,
                  const Box&      domain,
		  int             verbose)
{

  BL_PROFILE("mac_level_driver");
  MacOperator mac_op(mac_bndry,dx);
  mac_op.setCoefficients(area,mac_coef,dx,phys_bc,domain);
  mac_op.defRHS(area,volume,Rhs,u_mac,RhoD,phys_bc,domain,rhs_scale);   
  mac_op.maxOrder(4);
  
  if (verbose > 1)
  {
    Real sum = 0.0;
    for (MFIter mfi(Rhs); mfi.isValid(); ++mfi)
      sum += Rhs[mfi].sum(mfi.validbox(),0);
  
    ParallelDescriptor::ReduceRealSum(sum);
    if (ParallelDescriptor::IOProcessor())
      std::cout << "Sum of RHS at level " << level 
  		<< " = " << sum << std::endl;
  }

  if (the_solver == 1 && mac_op.maxOrder() != 2)
    {
      BoxLib::Error("Can't use CGSolver with maxorder > 2");
    }
  //
  // Construct MultiGrid or CGSolver object and solve system.
  //
  if (the_solver == 1)
    {
      bool use_mg_precond = true;
      CGSolver mac_cg(mac_op,use_mg_precond);
      mac_cg.solve(*mac_phi,Rhs,mac_tol,mac_abs_tol);
    }
  else if (the_solver == 3 ) 
    {
#ifdef MG_USE_FBOXLIB
      std::vector<BoxArray> bav(1);
      bav[0] = mac_phi->boxArray();
      std::vector<DistributionMapping> dmv(1);
      dmv[0] = Rhs.DistributionMap();
      bool nodal = false;
      std::vector<Geometry> geom(1);
      geom[0] = mac_bndry.getGeom();

      int mg_bc[2*BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i )
      { 
        if ( geom[0].isPeriodic(i) )
        {
          mg_bc[i*2 + 0] = 0;
          mg_bc[i*2 + 1] = 0;
        }
      else
        {
	  mg_bc[i*2 + 0] =((phys_bc.lo(i)== EXT_DIR))? MGT_BC_DIR : MGT_BC_NEU;
          mg_bc[i*2 + 1] =((phys_bc.hi(i)== EXT_DIR))? MGT_BC_DIR : MGT_BC_NEU;
        }
      }

      Array< Array<Real> > xa(1);
      Array< Array<Real> > xb(1);
      
      xa[0].resize(BL_SPACEDIM);
      xb[0].resize(BL_SPACEDIM);
      
      if (level == 0) {
	for ( int i = 0; i < BL_SPACEDIM; ++i ) {
	  xa[0][i] = 0.;
	  xb[0][i] = 0.;
	}
      } else {
	// a bit of a problem here.  parent not defined.
	const Real* dx_crse = geom[0].CellSize();
	for ( int i = 0; i < BL_SPACEDIM; ++i ) {
	  xa[0][i] = 0.5*dx_crse[i];
	  xb[0][i] = 0.5*dx_crse[i];
	}
      }
      
      MGT_Solver mgt_solver(geom, mg_bc, bav, dmv, nodal);
      mgt_solver.set_maxorder(4);

      const MultiFab* aa_p[1]; 
      aa_p[0] = &(mac_op.aCoefficients());
      const MultiFab* bb_p[1][BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i )
      {
          bb_p[0][i] = &(mac_op.bCoefficients(i));
      }
      mgt_solver.set_mac_coefficients(aa_p, bb_p, xa, xb);

      MultiFab* mac_phi_p[1];
      MultiFab* Rhs_p[1];
      mac_phi_p[0] = mac_phi;
      Rhs_p[0] = &Rhs;

      Real final_resnorm;
      mgt_solver.solve(mac_phi_p, Rhs_p, mac_tol, mac_abs_tol, mac_bndry, final_resnorm);
#else
      BoxLib::Error("mac_level_driver::mg_cpp not in this build");
#endif
    }
  else
    {
      MultiGrid mac_mg(mac_op);
      mac_mg.solve(*mac_phi,Rhs,mac_tol,mac_abs_tol);
    }

  Real vel_update_scale = -1.0;
  mac_op.velUpdate(u_mac,*mac_phi,mac_coef,phys_bc,dx,domain,level,vel_update_scale);
}

//
// A driver function for computing a sync MAC solve.
//

void
mac_sync_driver (const PressBndry& mac_bndry,
	         const BCRec&    phys_bc,
                 const BoxArray& grids,
                 int             the_solver,
                 int             level, 
                 const Real*     dx,
                 Real            mac_sync_tol,
                 Real            mac_abs_tol,
                 Real            rhs_scale,
                 MultiFab*       area,
                 MultiFab&       volume,
                 MultiFab&       Rhs,
                 MultiFab*       mac_coef,
                 MultiFab*       mac_sync_phi,
		 MultiFab*       u_mac_sync,
		 MultiFab*       rhoD_sync,
                 const Box&      domain,
		 int             verbose)
{
  BL_PROFILE("MacOperator::mac_sync_driver()");
  MacOperator mac_op(mac_bndry,dx);
  mac_op.maxOrder(4);
  mac_op.setCoefficients(area,mac_coef, dx, phys_bc,domain);
  mac_op.syncRhs(volume,Rhs,u_mac_sync,rhoD_sync,area,rhs_scale,
		 phys_bc,domain);
 
  bool allNeumann = true;
  std::vector<Geometry> fgeom(1);
  fgeom[0] = mac_bndry.getGeom();
  for ( int i = 0; i < BL_SPACEDIM; ++i )
  {
    if ((~fgeom[0].isPeriodic(i)) &&
	(phys_bc.lo(i) == EXT_DIR || phys_bc.hi(i) == EXT_DIR))
        allNeumann = false;
  }
  
  if (level == 0 && allNeumann) 
  {
    Real sum = 0.0;
    Real npts = 0.0;	
    for (MFIter mfi(Rhs); mfi.isValid(); ++mfi)
    {
      sum += Rhs[mfi].sum(mfi.validbox(),0);
      npts += mfi.validbox().d_numPts();
    }

    ParallelDescriptor::ReduceRealSum(sum);	
    ParallelDescriptor::ReduceRealSum(npts);	

    if (verbose > 1)
    {
      if (ParallelDescriptor::IOProcessor())
	std::cout << "Before adjustment: sum of RHS at level " << level 
		  << " = " << sum << std::endl;
    }

    if (std::fabs(sum) > 1.e-16) Rhs.plus(-sum/npts,0);

    if (verbose > 1)
    {
      sum = 0.0;
      for (MFIter mfi(Rhs); mfi.isValid(); ++mfi)
	sum += Rhs[mfi].sum(mfi.validbox(),0);
      ParallelDescriptor::ReduceRealSum(sum);
      
      if (ParallelDescriptor::IOProcessor())
	std::cout << "After adjustment: sum of RHS at level " << level 
		  << " = " << sum << std::endl;
    }
  }
  
  if (the_solver == 1 && mac_op.maxOrder() != 2)
    {
      BoxLib::Error("Can't use CGSolver with maxorder > 2");
    }
  //
  // Now construct MultiGrid or CGSolver object to solve system.
  //
  if (the_solver == 1)
    {
      bool use_mg_precond = true;
      CGSolver mac_cg(mac_op,use_mg_precond);
      mac_cg.solve(*mac_sync_phi,Rhs,mac_sync_tol,mac_abs_tol);
    }
  else if ( the_solver == 2 )
    {
#ifdef MG_USE_HYPRE
      HypreABec hp(mac_sync_phi->boxArray(), mac_bndry, dx, 0, false);
      hp.setScalars(mac_op.get_alpha(), mac_op.get_beta());
      hp.aCoefficients(mac_op.aCoefficients());
      for ( int i = 0; i < BL_SPACEDIM; ++i )
        {
	  hp.bCoefficients(mac_op.bCoefficients(i), i);
        }
      hp.setup_solver(mac_sync_tol, mac_abs_tol, 50);
      hp.solve(*mac_sync_phi, Rhs, true);
      hp.clear_solver();
#else
      BoxLib::Error("mac_sync_driver: HypreABec not in this build");
#endif
    }
  else if (the_solver == 3 )
    {
#ifdef MG_USE_FBOXLIB
      std::vector<BoxArray> bav(1);
      bav[0] = mac_sync_phi->boxArray();
      std::vector<DistributionMapping> dmv(1);
      dmv[0] = Rhs.DistributionMap();
      bool nodal = false;
      std::vector<Geometry> geom(1);
      geom[0] = mac_bndry.getGeom();

      int mg_bc[2*BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i )
      { 
        if ( geom[0].isPeriodic(i) )
        {
          mg_bc[i*2 + 0] = 0;
          mg_bc[i*2 + 1] = 0;
        }
	else
        {
	  mg_bc[i*2 + 0] =((phys_bc.lo(i)== EXT_DIR))? MGT_BC_DIR : MGT_BC_NEU;
          mg_bc[i*2 + 1] =((phys_bc.hi(i)== EXT_DIR))? MGT_BC_DIR : MGT_BC_NEU;
        }
      }

      Array< Array<Real> > xa(1);
      Array< Array<Real> > xb(1);
      
      xa[0].resize(BL_SPACEDIM);
      xb[0].resize(BL_SPACEDIM);
      
      if (level == 0) {
	for ( int i = 0; i < BL_SPACEDIM; ++i ) {
	  xa[0][i] = 0.;
	  xb[0][i] = 0.;
	}
      } else {
	const Real* dx_crse   = geom[0].CellSize();
	for ( int i = 0; i < BL_SPACEDIM; ++i ) {
	  xa[0][i] = 0.5*dx_crse[i];
	  xb[0][i] = 0.5*dx_crse[i];
	}
      }
      MGT_Solver mgt_solver(geom, mg_bc, bav, dmv, nodal);

      const MultiFab* aa_p[1];
      aa_p[0] = &(mac_op.aCoefficients());
      const MultiFab* bb_p[1][BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i )
      {
          bb_p[0][i] = &(mac_op.bCoefficients(i));
      }
      mgt_solver.set_mac_coefficients(aa_p, bb_p, xa, xb);

      MultiFab* mac_phi_p[1];
      MultiFab* Rhs_p[1];
      mac_phi_p[0] = mac_sync_phi;
      Rhs_p[0] = &Rhs;

      Real final_resnorm;
      mgt_solver.solve(mac_phi_p, Rhs_p, mac_sync_tol, mac_abs_tol, mac_bndry, final_resnorm);
#else
      BoxLib::Error("mac_sync_driver::mg_cpp not in this build");
#endif
    }
  else
    {
      MultiGrid mac_mg(mac_op);
      mac_mg.solve(*mac_sync_phi,Rhs,mac_sync_tol,mac_abs_tol);
    }
  
  int mac_op_lev = 0;
  mac_op.applyBC(*mac_sync_phi,0,1,mac_op_lev);
}
