
#ifndef _MACOPMACDRIVERS_H_

//====================================================================
// This header file contains declarations for the MacOperator mac projection
// driver routines
//====================================================================


// the driver functions for the MAC solves
void mac_level_driver( const MacBndry &mac_bndry, 
                       const BCRec& phys_bc,
                       const BoxArray &grids,
                       int use_cg_solve, int level, 
                       const Real *dx, 
                       Real mac_tol, Real mac_abs_tol, Real rhs_scale,
                       MultiFab *area,  MultiFab &volume,
                       MultiFab *mac_coef,     MultiFab &Rhs,
                       MultiFab *u_mac, MultiFab *RhoD, 
		       MultiFab *mac_phi, const Box& domain,
		       int verbose);

// a driver function for computing a sync MAC solve
void mac_sync_driver( const MacBndry &mac_bndry, const BCRec& phys_bc,
                      const BoxArray &grids,
                      int use_cg_solve, int level, 
                      const Real *dx, 
                      Real mac_sync_tol, Real mac_abs_tol, Real rhs_scale,
                      MultiFab *area,  MultiFab &volume,
                      MultiFab &Rhs,   MultiFab *mac_coef,
                      MultiFab *mac_sync_phi, MultiFab* u_mac_sync,
		      MultiFab* rhoD_sync, const Box& domain,
		      int verbose);

// apply the mac pressure gradient to a velocity field
// init, means that velocities are initialized here
void mac_vel_update( int init,
                     FArrayBox &uvel,
                     FArrayBox &vvel,
#if (BL_SPACEDIM == 3 )
                     FArrayBox &wvel,
#endif
		     const FArrayBox &lx,
                     const FArrayBox &ly,
#if (BL_SPACEDIM == 3 )
                     const FArrayBox &lz,
#endif
                     const FArrayBox &phi,
                     const BCRec& phys_bc,
                     const Box &grd, int level, int n,
                     const Real *dx, const Box& domain, Real scale );

void mac_vel_norml1( Real* normval,
                     FArrayBox &uvel,
                     FArrayBox &vvel,
#if (BL_SPACEDIM == 3 )
                     FArrayBox &wvel,
#endif
                     const Box &grd, 
                     const Real *dx);
#endif

void mac_vel_sync   (FArrayBox &uvel,
                     FArrayBox &vvel,
#if (BL_SPACEDIM == 3 )
                     FArrayBox &wvel,
#endif
		     const FArrayBox &rhodx,
                     const FArrayBox &rhody,
#if (BL_SPACEDIM == 3 )
                     const FArrayBox &rhodz,
#endif
		     const FArrayBox &areax,
                     const FArrayBox &areay,
#if (BL_SPACEDIM == 3 )
                     const FArrayBox &areaz,
#endif 
		     const FArrayBox &lx,
                     const FArrayBox &ly,
#if (BL_SPACEDIM == 3 )
                     const FArrayBox &lz,
#endif
		     const FArrayBox &phi,
                     const Box &grd,
		     const Real *dx, 
		     Real scale );
