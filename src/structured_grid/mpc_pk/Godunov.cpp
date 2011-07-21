//
// $Id: Godunov.cpp,v 1.41 2011-07-07 17:30:49 gpau Exp $
//

//
// Godunov is the object which calculates advective terms.
//
#include <winstd.H>

#include <LO_BCTYPES.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <FArrayBox.H>
#include <Godunov.H>
#include <GODUNOV_F.H>

#include <algorithm>

#define GEOM_GROW 1
#define HYP_GROW  3
#define XVEL 0
#define YVEL 1
#define ZVEL 2

int Godunov::verbose               = 0;
int Godunov::slope_order           = 4;
int Godunov::use_forces_in_trans   = 0;
const int use_unlimited_slopes_DEF = 0;

//
// Construct the Godunov Object.
//

Godunov::Godunov ()
  :
  max_1d(0)
{
  read_params();
    
  ZeroScratch();
  SetScratch(128);
}

//
// Size the 1D workspace explicitly.
//

Godunov::Godunov (int max_size)
  :
  max_1d(max_size)
{
  read_params();
  ZeroScratch();
  SetScratch(max_size);
}

//
// Read parameters from input file and command line.
//

void
Godunov::read_params ()
{
  //
  // Read parameters from input file and command line.
  //
  ParmParse pp("godunov");

  pp.query("v",verbose);

  pp.query("slope_order",slope_order);
#if (BL_SPACEDIM==2)
  BL_ASSERT(slope_order==1 || slope_order==2 || slope_order==4);
#else
  BL_ASSERT(slope_order==1 || slope_order==4);
#endif
  pp.query("use_forces_in_trans",use_forces_in_trans);
  int use_unlimited_slopes=use_unlimited_slopes_DEF;
  pp.query("use_unlimited_slopes",use_unlimited_slopes);

  FORT_SET_PARAMS(slope_order,use_unlimited_slopes);
}

//
// Set 1d scratch space as empty.
//

void
Godunov::ZeroScratch ()
{
  D_TERM(stxlo=0;,  stylo=0;,  stzlo=0;);
  D_TERM(stxhi=0;,  styhi=0;,  stzhi=0;);
  D_TERM(slxscr=0;, slyscr=0;, slzscr=0;);
}

//
// Initialize 1d scratch space to a bogus value.
//

void
Godunov::SetBogusScratch ()
{
#ifndef NDEBUG
  const Real bogus_value = 1.e200;

  for (int i = 0 ; i < scr_size ; i++)
    {
      D_TERM(stxlo[i]=bogus_value;,
	     stylo[i]=bogus_value;,
	     stzlo[i]=bogus_value;);

      D_TERM(stxhi[i]=bogus_value;,
	     styhi[i]=bogus_value;,
	     stzhi[i]=bogus_value;);

      D_TERM(slxscr[i]=bogus_value;,
	     slyscr[i]=bogus_value;,
	     slzscr[i]=bogus_value;);
    }
#endif /*NDEBUG*/
}

//
// Set 1d scratch space.
//

void
Godunov::SetScratch (int max_size)
{
  //
  // Set sizing parameters.
  //
  if (max_size <= max_1d)
    return;
  else
    max_1d = std::max(max_1d,max_size);
  scr_size = (max_size+2*HYP_GROW)*6;
  if (verbose > 1 && ParallelDescriptor::IOProcessor())
    {
      std::cout << "Set scratch size: "
		<< "max_size = " << max_size 
		<< "scr_size = " << scr_size << std::endl;
    }
  //
  // Get rid of the old scratch space.
  //
  RemScratch();
  //
  // Construct arrays.
  //
  D_TERM(stxlo  = new Real[scr_size];,
	 stylo  = new Real[scr_size];,
	 stzlo  = new Real[scr_size];);

  D_TERM(stxhi  = new Real[scr_size];,
	 styhi  = new Real[scr_size];,
	 stzhi  = new Real[scr_size];);

  D_TERM(slxscr = new Real[scr_size];,
	 slyscr = new Real[scr_size];,
	 slzscr = new Real[scr_size];);

}

//
// Remove 1D scratch space.
//

void
Godunov::RemScratch ()
{

  D_TERM(delete [] stxlo;,  delete [] stylo;,  delete [] stzlo;);
  D_TERM(delete [] stxhi;,  delete [] styhi;,  delete [] stzhi;);
  D_TERM(delete [] slxscr;, delete [] slyscr;, delete [] slzscr;);

}

//
// Destructor destroys work arrays.
//

Godunov::~Godunov ()
{
  RemScratch();
}

//
// Set up the farrayboxes for computing edge states, This function
// returns, the Farrayboxes where the fluxes are stored as these
// are used in PorousMedia.
//
// It also computes the transverse advective velocities.
//
// The amount of workspace needed in FArrayBox work is currently 2*SDIM+1.
//

void
Godunov::Setup (const Box& grd,
                FArrayBox& xflux,
                FArrayBox& yflux,
#if (BL_SPACEDIM == 3)
                FArrayBox& zflux,
#endif
                const int  nscal,
		const int  model) 
{
  //
  // Compute the edge boxes.
  //
  D_TERM(xflux_bx = grd; xflux_bx.surroundingNodes(0);,
	 yflux_bx = grd; yflux_bx.surroundingNodes(1);,
	 zflux_bx = grd; zflux_bx.surroundingNodes(2););
  //
  // Create storage for fluxes.
  //
  D_TERM(xflux.resize(xflux_bx,nscal);,
	 yflux.resize(yflux_bx,nscal);,
	 zflux.resize(zflux_bx,nscal););
    
  //
  // Ensure 1D scratch space is large enough.
  //
  SetScratch(BoxLib::grow(grd,HYP_GROW).longside());
  //
  // Create the advective velocities and FAB workspace for GODUNOV Box.
  int work_sz;
  work_bx = BoxLib::grow(grd,1);
  if (model == -1)
    work.resize(work_bx,2*BL_SPACEDIM+1);
  else if (model == 0 || model == 1)
    {
      work_sz = BL_SPACEDIM*nscal;
      work_rmn.resize(work_bx,work_sz);
    }
  else if (model == 2)
    {
      work_sz = 3*BL_SPACEDIM*nscal+3*BL_SPACEDIM*nscal*nscal;
      work_rmn.resize(work_bx,work_sz);
    }
  else
    { 
      work_sz = 3*BL_SPACEDIM*nscal+3*BL_SPACEDIM*nscal*nscal;
      work_rmn.resize(work_bx,work_sz);
      work_bx = BoxLib::grow(grd,3);
      utmp.resize(work_bx,nscal-1);  
      utmpn.resize(work_bx,nscal-1);
    }

  SetBogusScratch();
}


void
Godunov::Setup_tracer (const Box&       grd,
		       FArrayBox&       xflux,
		       FArrayBox&       yflux,
#if (BL_SPACEDIM == 3)
		       FArrayBox&       zflux,
#endif
		       const int        nscal) 
{
  //
  // Compute the edge boxes.
  //
  D_TERM(xflux_bx = grd; xflux_bx.surroundingNodes(0);,
	 yflux_bx = grd; yflux_bx.surroundingNodes(1);,
	 zflux_bx = grd; zflux_bx.surroundingNodes(2););
  //
  // Create storage for fluxes.
  //
  D_TERM(xflux.resize(xflux_bx,nscal);,
	 yflux.resize(yflux_bx,nscal);,
	 zflux.resize(zflux_bx,nscal););

  D_TERM(xflux_bx.grow(1);,
	 yflux_bx.grow(1);,
	 zflux_bx.grow(1););

  D_TERM(xlo.resize(xflux_bx,nscal); xlo.setVal(0.);,
	 ylo.resize(yflux_bx,nscal); ylo.setVal(0.);,
	 zlo.resize(zflux_bx,nscal); zlo.setVal(0.););

  D_TERM(xhi.resize(xflux_bx,nscal); xhi.setVal(0.);,
	 yhi.resize(yflux_bx,nscal); yhi.setVal(0.);,
	 zhi.resize(zflux_bx,nscal); zhi.setVal(0.););

  //
  // Ensure 1D scratch space is large enough.
  //
  SetScratch(BoxLib::grow(grd,HYP_GROW).longside());
  //
  // Create the advective velocities and FAB workspace for GODUNOV Box.
  work_bx = BoxLib::grow(grd,1);
  work.resize(work_bx,BL_SPACEDIM*nscal);

  utmp.clear();
  utmpn.clear();
  work_rmn.clear();

  SetBogusScratch();
}

//
// Advection functions follow.
//

//
// Compute the edge states using the advective transverse velocities
// The amount of workspace needed in FArrayBox work is currently 2*SDIM+1.
//

void
Godunov::edge_states (const Box&  grd,
                      const Real* dx,
                      Real        dt,
                      int         velpred,
                      FArrayBox&  uedge,
                      FArrayBox&  stx,
                      FArrayBox&  vedge,
                      FArrayBox&  sty,
#if (BL_SPACEDIM == 3)               
                      FArrayBox&  wedge,
                      FArrayBox&  stz,
#endif
                      FArrayBox&  S,
                      FArrayBox&  tforces,
                      FArrayBox&  divu,
                      int         fab_ind,
                      int         state_ind,
                      const int*  bc,
                      int         iconserv)
{
  //
  // Error block.
  //
  BL_ASSERT(S.box().contains(work_bx));

  BL_ASSERT(S.nComp()       >= fab_ind    );

  BL_ASSERT(uedge.nComp()   >= 1          );
  BL_ASSERT(stx.nComp()     >= 1          );

  BL_ASSERT(vedge.nComp()   >= 1          );
  BL_ASSERT(sty.nComp()     >= 1          );
#if (BL_SPACEDIM == 3)
  BL_ASSERT(wedge.nComp()   >= 1          );
  BL_ASSERT(stz.nComp()     >= 1          );
#endif    
  //
  // Create the bounds and pointers.
  //
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();
  const int *s_lo       = S.loVect();
  const int *s_hi       = S.hiVect();
  const int *ww_lo      = work.loVect();
  const int *ww_hi      = work.hiVect();
  const Real *s_dat     = S.dataPtr(fab_ind);
  const Real *tfr_dat   = tforces.dataPtr(fab_ind);
  const Real *divu_dat  = divu.dataPtr();
  //
  // Set work space to bogus values.
  //
  SetBogusScratch();

#if (BL_SPACEDIM == 3)
  const Real *xhi_dat   = work.dataPtr(0);
  const Real *yhi_dat   = work.dataPtr(0);
  const Real *zhi_dat   = work.dataPtr(0);
  const Real *xlo_dat   = work.dataPtr(1);
  const Real *ylo_dat   = work.dataPtr(2);
  const Real *zlo_dat   = work.dataPtr(3);
  const Real *slx_dat   = work.dataPtr(4);
  const Real *sly_dat   = work.dataPtr(5);
  const Real *slz_dat   = work.dataPtr(6);
#else
  const Real *xhi_dat   = work.dataPtr(0);
  const Real *yhi_dat   = work.dataPtr(0);
  const Real *xlo_dat   = work.dataPtr(1);
  const Real *ylo_dat   = work.dataPtr(2);
  const Real *slx_dat   = work.dataPtr(3);
  const Real *sly_dat   = work.dataPtr(4);
#endif
  //
  // C component indices starts from 0, Fortran from 1
  //
  int fort_ind = state_ind+1;  
      
  FORT_ESTATE_FPU(s_dat, tfr_dat, divu_dat, ARLIM(s_lo), ARLIM(s_hi),
                    
		  xlo_dat, xhi_dat, slx_dat,
		  slxscr, stxlo, stxhi,
		  uedge.dataPtr(), ARLIM(uedge.loVect()), ARLIM(uedge.hiVect()),
		  stx.dataPtr(),   ARLIM(stx.loVect()),   ARLIM(stx.hiVect()),
                    
		  ylo_dat, yhi_dat, sly_dat,
		  slyscr, stylo, styhi,
		  vedge.dataPtr(), ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
		  sty.dataPtr(),   ARLIM(sty.loVect()),   ARLIM(sty.hiVect()),
#if (BL_SPACEDIM == 3)
		  zlo_dat, zhi_dat, slz_dat,
		  slzscr, stzlo, stzhi,
		  wedge.dataPtr(), ARLIM(wedge.loVect()), ARLIM(wedge.hiVect()),
		  stz.dataPtr(),   ARLIM(stz.loVect()),   ARLIM(stz.hiVect()),
#endif
		  ARLIM(ww_lo), ARLIM(ww_hi),
		  bc, lo, hi, &dt, dx, &fort_ind,
		  &use_forces_in_trans, &iconserv);
}

void
Godunov::edge_states_lin (const Box&  grd,
			  const Real* dx,
			  Real        dt,
			  FArrayBox&  uedge,
			  FArrayBox&  stx,
			  FArrayBox&  vedge,
			  FArrayBox&  sty,
#if (BL_SPACEDIM == 3)               
			  FArrayBox&  wedge,
			  FArrayBox&  stz,
#endif
			  FArrayBox&  S,
			  FArrayBox&  S_new,
			  FArrayBox&  tforces,
			  FArrayBox&  rock_phi,
			  int         fab_ind,
			  int         state_ind,
			  const int*  bc,
			  int         nscal)
{

  // 
  // Compute the fluxes at edges.  
  // NB: The dimension of work is box + 1 ghost cell.
  //     This is thus valid for all non-state FArrayBox.
  // 
    
  BL_ASSERT(S.box().contains(work_bx));

  BL_ASSERT(S.nComp()       >= nscal      );
  BL_ASSERT(S_new.nComp()   >= nscal      );

  BL_ASSERT(uedge.nComp()   >= 1          );
  BL_ASSERT(stx.nComp()     >= nscal      );

  BL_ASSERT(vedge.nComp()   >= 1          );
  BL_ASSERT(sty.nComp()     >= nscal      );

#if (BL_SPACEDIM == 3)
  BL_ASSERT(wedge.nComp()   >= 1          );
  BL_ASSERT(stz.nComp()     >= nscal     );
#endif    

  //
  // Create the bounds and pointers.
  //
  const int *lo           = grd.loVect();
  const int *hi           = grd.hiVect();

  const Real *s_dat       = S.dataPtr(fab_ind);
  const int *s_lo         = S.loVect();
  const int *s_hi         = S.hiVect();

  const Real *sn_dat      = S_new.dataPtr(fab_ind);
  const int *sn_lo        = S_new.loVect();
  const int *sn_hi        = S_new.hiVect();

  const Real *tfr_dat     = tforces.dataPtr(fab_ind);
  const int *ww_lo        = work_rmn.loVect();
  const int *ww_hi        = work_rmn.hiVect();
  
  //
  // Set work space to bogus values.
  //
  SetBogusScratch();

#if (BL_SPACEDIM == 3)
  const Real *slx_dat   = work_rmn.dataPtr(0*nscal);
  const Real *sly_dat   = work_rmn.dataPtr(1*nscal);
  const Real *slz_dat   = work_rmn.dataPtr(2*nscal);
#else
  const Real *slx_dat   = work_rmn.dataPtr(0*nscal);
  const Real *sly_dat   = work_rmn.dataPtr(1*nscal);
#endif
  //
  // C component indices starts from 0, Fortran from 1
  //
  FORT_ESTATE_LIN(s_dat,  ARLIM(s_lo), ARLIM(s_hi),
		  sn_dat, ARLIM(sn_lo), ARLIM(sn_hi), 
		  tfr_dat, 
		  rock_phi.dataPtr(),ARLIM(rock_phi.loVect()),ARLIM(rock_phi.hiVect()),
		  slx_dat, slxscr,
		  uedge.dataPtr(),ARLIM(uedge.loVect()),ARLIM(uedge.hiVect()),
		  stx.dataPtr(),ARLIM(stx.loVect()),ARLIM(stx.hiVect()),
		  sly_dat, slyscr,
		  vedge.dataPtr(),ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
		  sty.dataPtr(),ARLIM(sty.loVect()),ARLIM(sty.hiVect()),
#if (BL_SPACEDIM == 3)
		  slz_dat, slzscr,
		  wedge.dataPtr(),ARLIM(wedge.loVect()),ARLIM(wedge.hiVect()),
		  stz.dataPtr(),ARLIM(stz.loVect()),ARLIM(stz.hiVect()),
#endif
		  ARLIM(ww_lo), ARLIM(ww_hi),
		  bc, lo, hi, &dt, dx, &nscal);
}

void
Godunov::edge_states_rmn (const Box&  grd,
			  const Real* dx,
			  Real        dt,
			  FArrayBox&  uedge,
			  FArrayBox&  stx,
			  FArrayBox&  kappax,
			  FArrayBox&  vedge,
			  FArrayBox&  sty,
			  FArrayBox&  kappay,
#if (BL_SPACEDIM == 3)               
			  FArrayBox&  wedge,
			  FArrayBox&  stz,
			  FArrayBox&  kappaz,
#endif
			  FArrayBox&  S,
			  FArrayBox&  S_new,
			  FArrayBox&  tforces,
			  FArrayBox&  divu,
			  FArrayBox&  rock_phi,
			  FArrayBox&  kappa,
			  FArrayBox&  lambda_cc,
			  FArrayBox&  dlambda_cc,
			  FArrayBox&  kr_coef,
			  const int   n_kr_coef,
			  int         fab_ind,
			  int         state_ind,
			  const int*  bc,
			  int         iconserv,
			  int         nscal)
{

  // 
  // Compute the fluxes at edges.  
  // NB: The dimension of work is box + 1 ghost cell.
  //     This is thus valid for all non-state FArrayBox.
  // 
    
  BL_ASSERT(S.box().contains(work_bx));

  BL_ASSERT(S.nComp()       >= nscal      );
  BL_ASSERT(S_new.nComp()   >= nscal      );

  BL_ASSERT(uedge.nComp()   >= 1          );
  BL_ASSERT(stx.nComp()     >= nscal      );

  BL_ASSERT(vedge.nComp()   >= 1          );
  BL_ASSERT(sty.nComp()     >= nscal      );

#if (BL_SPACEDIM == 3)
  BL_ASSERT(wedge.nComp()   >= 1          );
  BL_ASSERT(stz.nComp()     >= nscal     );
#endif    

  //
  // Create the bounds and pointers.
  //
  const int *lo           = grd.loVect();
  const int *hi           = grd.hiVect();

  const Real *s_dat       = S.dataPtr(fab_ind);
  const int *s_lo         = S.loVect();
  const int *s_hi         = S.hiVect();

  const Real *sn_dat      = S_new.dataPtr(fab_ind);
  const int *sn_lo        = S_new.loVect();
  const int *sn_hi        = S_new.hiVect();

  const Real *tfr_dat     = tforces.dataPtr(fab_ind);
  const Real *divu_dat    = divu.dataPtr();
  const Real *lbd_dat     = lambda_cc.dataPtr();
  const Real *dlbd_dat    = dlambda_cc.dataPtr();
  const int *ww_lo        = work_rmn.loVect();
  const int *ww_hi        = work_rmn.hiVect();
  
  const Real *kr_dat      = kr_coef.dataPtr();
  const int *kr_lo        = kr_coef.loVect();
  const int *kr_hi        = kr_coef.hiVect();
  //
  // Set work space to bogus values.
  //
  SetBogusScratch();

#if (BL_SPACEDIM == 3)
  const Real *slx_dat   = work_rmn.dataPtr(0*nscal);
  const Real *sly_dat   = work_rmn.dataPtr(1*nscal);
  const Real *slz_dat   = work_rmn.dataPtr(2*nscal);
  const Real *eigvx_dat = work_rmn.dataPtr(3*nscal);
  const Real *eigvy_dat = work_rmn.dataPtr(5*nscal);
  const Real *eigvz_dat = work_rmn.dataPtr(7*nscal);
  const Real *eiglx_dat = work_rmn.dataPtr(9*nscal);
  const Real *eigly_dat = work_rmn.dataPtr(9*nscal+nscal*nscal);
  const Real *eiglz_dat = work_rmn.dataPtr(9*nscal+2*nscal*nscal);
  const Real *eigrx_dat = work_rmn.dataPtr(9*nscal+3*nscal*nscal);
  const Real *eigry_dat = work_rmn.dataPtr(9*nscal+4*nscal*nscal);
  const Real *eigrz_dat = work_rmn.dataPtr(9*nscal+5*nscal*nscal);
#else
  const Real *slx_dat   = work_rmn.dataPtr(0*nscal);
  const Real *sly_dat   = work_rmn.dataPtr(1*nscal);
  const Real *eigvx_dat = work_rmn.dataPtr(2*nscal);
  const Real *eigvy_dat = work_rmn.dataPtr(4*nscal);
  const Real *eiglx_dat = work_rmn.dataPtr(6*nscal);
  const Real *eigly_dat = work_rmn.dataPtr(6*nscal+nscal*nscal);
  const Real *eigrx_dat = work_rmn.dataPtr(6*nscal+2*nscal*nscal);
  const Real *eigry_dat = work_rmn.dataPtr(6*nscal+3*nscal*nscal);
#endif
  //
  // C component indices starts from 0, Fortran from 1
  //
  FORT_ESTATE_RMN(s_dat,  ARLIM(s_lo), ARLIM(s_hi),
		  sn_dat, ARLIM(sn_lo), ARLIM(sn_hi), 
		  tfr_dat, divu_dat, lbd_dat, dlbd_dat,
		  rock_phi.dataPtr(),ARLIM(rock_phi.loVect()),ARLIM(rock_phi.hiVect()),
		  kappa.dataPtr(),ARLIM(kappa.loVect()),ARLIM(kappa.hiVect()),
		  slx_dat, slxscr,
		  uedge.dataPtr(),ARLIM(uedge.loVect()),ARLIM(uedge.hiVect()),
		  stx.dataPtr(),ARLIM(stx.loVect()),ARLIM(stx.hiVect()),
		  kappax.dataPtr(),ARLIM(kappax.loVect()),ARLIM(kappax.hiVect()),
		  eigvx_dat, eiglx_dat, eigrx_dat,
		  sly_dat, slyscr,
		  vedge.dataPtr(),ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
		  sty.dataPtr(),ARLIM(sty.loVect()),ARLIM(sty.hiVect()),
		  kappay.dataPtr(),ARLIM(kappay.loVect()),ARLIM(kappay.hiVect()),
		  eigvy_dat, eigly_dat, eigry_dat, 
#if (BL_SPACEDIM == 3)
		  slz_dat, slzscr,
		  wedge.dataPtr(),ARLIM(wedge.loVect()),ARLIM(wedge.hiVect()),
		  stz.dataPtr(),ARLIM(stz.loVect()),ARLIM(stz.hiVect()),
		  kappaz.dataPtr(),ARLIM(kappaz.loVect()),ARLIM(kappaz.hiVect()),
		  eigvz_dat, eiglz_dat, eigrz_dat,
#endif
		  ARLIM(ww_lo), ARLIM(ww_hi),
		  kr_dat, ARLIM(kr_lo), ARLIM(kr_hi), &n_kr_coef,
		  bc, lo, hi, &dt, dx, 
		  &use_forces_in_trans, &iconserv, &nscal);

}


void
Godunov::edge_states_cpl (const Box&  grd,
			  const Real* dx,
			  Real        dt,
			  FArrayBox&  uedge,
			  FArrayBox&  stx,
			  FArrayBox&  kappax,
			  FArrayBox&  lambdax,
			  FArrayBox&  vedge,
			  FArrayBox&  sty,
			  FArrayBox&  kappay,
			  FArrayBox&  lambday,
#if (BL_SPACEDIM == 3)               
			  FArrayBox&  wedge,
			  FArrayBox&  stz,
			  FArrayBox&  kappaz,
			  FArrayBox&  lambdaz,
#endif
			  FArrayBox&  S,
			  FArrayBox&  S_new,
			  FArrayBox&  tforces,
			  FArrayBox&  divu,
			  FArrayBox&  rock_phi,
			  FArrayBox&  kappa,
			  FArrayBox&  pc,
			  FArrayBox&  lambda_cc,
			  FArrayBox&  dlambda_cc,
			  FArrayBox&  kr_coef,
			  const int   n_kr_coef,
			  int         fab_ind,
			  int         state_ind,
			  const int*  bc,
			  int         iconserv,
			  int         nscal)
{

  // 
  // Compute the fluxes at edges.  
  // NB: The dimension of work is box + 1 ghost cell.
  //     This is thus valid for all non-state FArrayBox.
  // 
    
  BL_ASSERT(S.box().contains(work_bx));

  BL_ASSERT(S.nComp()       >= nscal      );
  BL_ASSERT(S_new.nComp()   >= nscal      );

  BL_ASSERT(uedge.nComp()   >= 1          );
  BL_ASSERT(stx.nComp()     >= nscal      );

  BL_ASSERT(vedge.nComp()   >= 1          );
  BL_ASSERT(sty.nComp()     >= nscal      );

#if (BL_SPACEDIM == 3)
  BL_ASSERT(wedge.nComp()   >= 1          );
  BL_ASSERT(stz.nComp()     >= nscal     );
#endif    

  //
  // Create the bounds and pointers.
  //
  const int *lo           = grd.loVect();
  const int *hi           = grd.hiVect();

  const Real *s_dat       = S.dataPtr(fab_ind);
  const int *s_lo         = S.loVect();
  const int *s_hi         = S.hiVect();

  const Real *sn_dat      = S_new.dataPtr(fab_ind);
  const int *sn_lo        = S_new.loVect();
  const int *sn_hi        = S_new.hiVect();

  const Real *pc_dat      = pc.dataPtr();
  const int *pc_lo        = pc.loVect();
  const int *pc_hi        = pc.hiVect();

  const Real *tfr_dat     = tforces.dataPtr(fab_ind);
  const Real *divu_dat    = divu.dataPtr();
  const Real *lbd_dat     = lambda_cc.dataPtr();
  const Real *dlbd_dat    = dlambda_cc.dataPtr();
  const int *ww_lo        = work_rmn.loVect();
  const int *ww_hi        = work_rmn.hiVect();
  
  const Real *kr_dat      = kr_coef.dataPtr();
  const int *kr_lo        = kr_coef.loVect();
  const int *kr_hi        = kr_coef.hiVect();
  //
  // Set work space to bogus values.
  //
  SetBogusScratch();

#if (BL_SPACEDIM == 3)
  const Real *slx_dat   = work_rmn.dataPtr(0*nscal);
  const Real *sly_dat   = work_rmn.dataPtr(1*nscal);
  const Real *slz_dat   = work_rmn.dataPtr(2*nscal);
  const Real *eigvx_dat = work_rmn.dataPtr(3*nscal);
  const Real *eigvy_dat = work_rmn.dataPtr(5*nscal);
  const Real *eigvz_dat = work_rmn.dataPtr(7*nscal);
  const Real *eiglx_dat = work_rmn.dataPtr(9*nscal);
  const Real *eigly_dat = work_rmn.dataPtr(9*nscal+nscal*nscal);
  const Real *eiglz_dat = work_rmn.dataPtr(9*nscal+2*nscal*nscal);
  const Real *eigrx_dat = work_rmn.dataPtr(9*nscal+3*nscal*nscal);
  const Real *eigry_dat = work_rmn.dataPtr(9*nscal+4*nscal*nscal);
  const Real *eigrz_dat = work_rmn.dataPtr(9*nscal+5*nscal*nscal);
#else
  const Real *slx_dat   = work_rmn.dataPtr(0*nscal);
  const Real *sly_dat   = work_rmn.dataPtr(1*nscal);
  const Real *eigvx_dat = work_rmn.dataPtr(2*nscal);
  const Real *eigvy_dat = work_rmn.dataPtr(4*nscal);
  const Real *eiglx_dat = work_rmn.dataPtr(6*nscal);
  const Real *eigly_dat = work_rmn.dataPtr(6*nscal+nscal*nscal);
  const Real *eigrx_dat = work_rmn.dataPtr(6*nscal+2*nscal*nscal);
  const Real *eigry_dat = work_rmn.dataPtr(6*nscal+3*nscal*nscal);
#endif
  //
  // C component indices starts from 0, Fortran from 1
  //
  FORT_ESTATE_CPL(s_dat,  ARLIM(s_lo),  ARLIM(s_hi),
		  sn_dat, ARLIM(sn_lo), ARLIM(sn_hi), 
		  pc_dat, ARLIM(pc_lo), ARLIM(pc_hi),
		  tfr_dat, divu_dat, lbd_dat, dlbd_dat,
		  rock_phi.dataPtr(),ARLIM(rock_phi.loVect()),ARLIM(rock_phi.hiVect()),
		  kappa.dataPtr(),ARLIM(kappa.loVect()),ARLIM(kappa.hiVect()),
		  slx_dat, slxscr,
		  uedge.dataPtr(),ARLIM(uedge.loVect()),ARLIM(uedge.hiVect()),
		  stx.dataPtr(),ARLIM(stx.loVect()),ARLIM(stx.hiVect()),
		  kappax.dataPtr(),ARLIM(kappax.loVect()),ARLIM(kappax.hiVect()),
		  lambdax.dataPtr(),ARLIM(lambdax.loVect()),ARLIM(lambdax.hiVect()),
		  eigvx_dat, eiglx_dat, eigrx_dat,
		  sly_dat, slyscr,
		  vedge.dataPtr(),ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
		  sty.dataPtr(),ARLIM(sty.loVect()),ARLIM(sty.hiVect()),
		  kappay.dataPtr(),ARLIM(kappay.loVect()),ARLIM(kappay.hiVect()),
		  lambday.dataPtr(),ARLIM(lambday.loVect()),ARLIM(lambday.hiVect()),
		  eigvy_dat, eigly_dat, eigry_dat, 
#if (BL_SPACEDIM == 3)
		  slz_dat, slzscr,
		  wedge.dataPtr(),ARLIM(wedge.loVect()),ARLIM(wedge.hiVect()),
		  stz.dataPtr(),ARLIM(stz.loVect()),ARLIM(stz.hiVect()),
		  kappaz.dataPtr(),ARLIM(kappaz.loVect()),ARLIM(kappaz.hiVect()),
		  lambdaz.dataPtr(),ARLIM(lambdaz.loVect()),ARLIM(lambdaz.hiVect()),
		  eigvz_dat, eiglz_dat, eigrz_dat,
#endif
		  ARLIM(ww_lo), ARLIM(ww_hi),
		  kr_dat, ARLIM(kr_lo), ARLIM(kr_hi), &n_kr_coef,
		  bc, lo, hi, &dt, dx, 
		  &use_forces_in_trans, &iconserv, &nscal);
}

void
Godunov::edge_states_pmr (const Box&  grd,
			  const Real* dx,
			  Real        dt,
			  FArrayBox&  uedge,
			  FArrayBox&  stx,
			  FArrayBox&  kappax,
			  FArrayBox&  vedge,
			  FArrayBox&  sty,
			  FArrayBox&  kappay,
#if (BL_SPACEDIM == 3)               
			  FArrayBox&  wedge,
			  FArrayBox&  stz,
			  FArrayBox&  kappaz,
#endif
			  FArrayBox&  S,
			  FArrayBox&  S_new,
			  FArrayBox&  tforces,
			  FArrayBox&  divu,
			  FArrayBox&  rock_phi,
			  FArrayBox&  kappa,
			  int         fab_ind,
			  int         state_ind,
			  const int*  bc,
			  int         iconserv,
			  int         nscal,
			  Real        gravity,
			  Real*       eigmax)
{
  //
  // Error block.
  //
  BL_ASSERT(S.box().contains(work_bx));
  BL_ASSERT(S_new.box().contains(work_bx));

  BL_ASSERT(S.nComp()       >= nscal      );
  BL_ASSERT(S_new.nComp()   >= nscal      );
  BL_ASSERT(tforces.nComp() >= fab_ind    );

  BL_ASSERT(uedge.nComp()   >= 1          );
  BL_ASSERT(stx.nComp()     >= nscal      );

  BL_ASSERT(vedge.nComp()   >= 1          );
  BL_ASSERT(sty.nComp()     >= nscal      );
#if (BL_SPACEDIM == 3)
  BL_ASSERT(wedge.nComp()   >= 1          );
  BL_ASSERT(stz.nComp()     >= nscal     );
#endif    
  //
  // Create the bounds and pointers.
  //
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();
  const int *s_lo       = S.loVect();
  const int *s_hi       = S.hiVect();
  const int *sn_lo      = S_new.loVect();
  const int *sn_hi      = S_new.hiVect();
  const int *ww_lo      = work_rmn.loVect();
  const int *ww_hi      = work_rmn.hiVect();
  const Real *s_dat     = S.dataPtr(fab_ind);
  const Real *sn_dat    = S_new.dataPtr(fab_ind);
  const Real *u_dat     = utmp.dataPtr(fab_ind);
  const Real *un_dat    = utmpn.dataPtr(fab_ind);
  const Real *tfr_dat   = tforces.dataPtr(fab_ind);
  const Real *divu_dat  = divu.dataPtr();
  //
  // Set work space to bogus values.
  //
  SetBogusScratch();

#if (BL_SPACEDIM == 3)
  const Real *slx_dat   = work_rmn.dataPtr(0*nscal);
  const Real *sly_dat   = work_rmn.dataPtr(1*nscal);
  const Real *slz_dat   = work_rmn.dataPtr(2*nscal);
  const Real *eigvx_dat = work_rmn.dataPtr(3*nscal);
  const Real *eigvy_dat = work_rmn.dataPtr(4*nscal);
  const Real *eigvz_dat = work_rmn.dataPtr(5*nscal);
  const Real *eiglx_dat = work_rmn.dataPtr(6*nscal);
  const Real *eigly_dat = work_rmn.dataPtr(6*nscal+nscal*nscal);
  const Real *eiglz_dat = work_rmn.dataPtr(6*nscal+2*nscal*nscal);
  const Real *eigrx_dat = work_rmn.dataPtr(6*nscal+3*nscal*nscal);
  const Real *eigry_dat = work_rmn.dataPtr(6*nscal+4*nscal*nscal);
  const Real *eigrz_dat = work_rmn.dataPtr(6*nscal+5*nscal*nscal);
  const Real *strcx_dat = work_rmn.dataPtr(6*nscal+6*nscal*nscal);
  const Real *strcy_dat = work_rmn.dataPtr(6*nscal+7*nscal*nscal);
  const Real *strcz_dat = work_rmn.dataPtr(6*nscal+8*nscal*nscal);
#else
  const Real *slx_dat   = work_rmn.dataPtr(0*nscal);
  const Real *sly_dat   = work_rmn.dataPtr(1*nscal);
  const Real *eigvx_dat = work_rmn.dataPtr(2*nscal);
  const Real *eigvy_dat = work_rmn.dataPtr(3*nscal);
  const Real *eiglx_dat = work_rmn.dataPtr(4*nscal);
  const Real *eigly_dat = work_rmn.dataPtr(4*nscal+nscal*nscal);
  const Real *eigrx_dat = work_rmn.dataPtr(4*nscal+2*nscal*nscal);
  const Real *eigry_dat = work_rmn.dataPtr(4*nscal+3*nscal*nscal);
  const Real *strcx_dat = work_rmn.dataPtr(4*nscal+4*nscal*nscal);
  const Real *strcy_dat = work_rmn.dataPtr(4*nscal+5*nscal*nscal);
#endif
  //
  // C component indices starts from 0, Fortran from 1
  //
  FORT_ESTATE_PMR(s_dat, u_dat, ARLIM(s_lo), ARLIM(s_hi),
		  sn_dat, un_dat, ARLIM(sn_lo), ARLIM(sn_hi), tfr_dat, divu_dat, 
		  rock_phi.dataPtr(),ARLIM(rock_phi.loVect()),ARLIM(rock_phi.hiVect()),
		  kappa.dataPtr(),ARLIM(kappa.loVect()),ARLIM(kappa.hiVect()),
		  slx_dat, slxscr,
		  uedge.dataPtr(),ARLIM(uedge.loVect()),ARLIM(uedge.hiVect()),
		  stx.dataPtr(),ARLIM(stx.loVect()),ARLIM(stx.hiVect()),
		  kappax.dataPtr(),ARLIM(kappax.loVect()),ARLIM(kappax.hiVect()),
		  strcx_dat,eigvx_dat, eiglx_dat, eigrx_dat, 
		  sly_dat, slyscr,
		  vedge.dataPtr(),ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
		  sty.dataPtr(),ARLIM(sty.loVect()),ARLIM(sty.hiVect()),
		  kappay.dataPtr(),ARLIM(kappay.loVect()),ARLIM(kappay.hiVect()),
		  strcy_dat,eigvy_dat, eigly_dat, eigry_dat, 
#if (BL_SPACEDIM == 3)
		  slz_dat, slzscr,
		  wedge.dataPtr(),ARLIM(wedge.loVect()),ARLIM(wedge.hiVect()),
		  stz.dataPtr(),ARLIM(stz.loVect()),ARLIM(stz.hiVect()),
		  kappaz.dataPtr(),ARLIM(kappaz.loVect()),ARLIM(kappaz.hiVect()),
		  strcz_dat,eigvz_dat, eiglz_dat, eigrz_dat, 
#endif
		  ARLIM(ww_lo), ARLIM(ww_hi),
		  bc, lo, hi, &dt, dx, 
		  &use_forces_in_trans, &iconserv, &gravity, eigmax);

}

void
Godunov::edge_states_tracer (const Box&  grd,
			     const Real* dx,
			     Real        dt,
			     FArrayBox&  uedge,
			     FArrayBox&  stx,
			     FArrayBox&  xlo,
			     FArrayBox&  xhi,
			     FArrayBox&  vedge,
			     FArrayBox&  sty,
			     FArrayBox&  ylo,
			     FArrayBox&  yhi,
#if (BL_SPACEDIM == 3)               
			     FArrayBox&  wedge,
			     FArrayBox&  stz,
			     FArrayBox&  zlo,
			     FArrayBox&  zhi,
#endif
			     FArrayBox&  S,
			     FArrayBox&  S_new,
			     FArrayBox&  St,
			     FArrayBox&  St_new,
			     FArrayBox&  rock_phi,
			     const int*  bc,
			     int         nscal)
{
    
  BL_ASSERT(S.box().contains(work_bx));

  BL_ASSERT(S.nComp()       >= nscal      );
  BL_ASSERT(S_new.nComp()   >= nscal      );

  BL_ASSERT(uedge.nComp()   >= 1          );
  BL_ASSERT(stx.nComp()     >= nscal      );

  BL_ASSERT(vedge.nComp()   >= 1          );
  BL_ASSERT(sty.nComp()     >= nscal      );

#if (BL_SPACEDIM == 3)
  BL_ASSERT(wedge.nComp()   >= 1          );
  BL_ASSERT(stz.nComp()     >= nscal     );
#endif    

  //
  // Create the bounds and pointers.
  //
  const int *lo           = grd.loVect();
  const int *hi           = grd.hiVect();
  const int *s_lo         = S.loVect();
  const int *s_hi         = S.hiVect();
  const int *sn_lo        = S_new.loVect();
  const int *sn_hi        = S_new.hiVect();
  const int *st_lo        = St.loVect();
  const int *st_hi        = St.hiVect();
  const int *stn_lo       = St_new.loVect();
  const int *stn_hi       = St_new.hiVect();
  const int *ww_lo        = work.loVect();
  const int *ww_hi        = work.hiVect();
  const Real *s_dat       = S.dataPtr();
  const Real *sn_dat      = S_new.dataPtr();
  const Real *st_dat      = St.dataPtr();
  const Real *stn_dat     = St_new.dataPtr();

#if (BL_SPACEDIM == 3)
  const Real *slx_dat   = work.dataPtr(0*nscal);
  const Real *sly_dat   = work.dataPtr(1*nscal);
  const Real *slz_dat   = work.dataPtr(2*nscal);
#else
  const Real *slx_dat   = work.dataPtr(0*nscal);
  const Real *sly_dat   = work.dataPtr(1*nscal);
#endif
  //
  // C component indices starts from 0, Fortran from 1
  //
  FORT_ESTATE_TRACER(s_dat, ARLIM(s_lo), ARLIM(s_hi),
		     sn_dat, ARLIM(sn_lo), ARLIM(sn_hi), 
		     st_dat, ARLIM(st_lo), ARLIM(st_hi),
		     stn_dat, ARLIM(stn_lo), ARLIM(stn_hi), 
		     slx_dat, slxscr,xlo.dataPtr(),xhi.dataPtr(),
		     uedge.dataPtr(),ARLIM(uedge.loVect()),ARLIM(uedge.hiVect()),
		     stx.dataPtr(),ARLIM(stx.loVect()),ARLIM(stx.hiVect()),
		     sly_dat, slyscr,ylo.dataPtr(),yhi.dataPtr(),
		     vedge.dataPtr(),ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
		     sty.dataPtr(),ARLIM(sty.loVect()),ARLIM(sty.hiVect()),
#if (BL_SPACEDIM == 3)
		     slz_dat, slzscr,zlo.dataPtr(),zhi.dataPtr(),
		     wedge.dataPtr(),ARLIM(wedge.loVect()),ARLIM(wedge.hiVect()),
		     stz.dataPtr(),ARLIM(stz.loVect()),ARLIM(stz.hiVect()),
#endif
		     rock_phi.dataPtr(),ARLIM(rock_phi.loVect()),ARLIM(rock_phi.hiVect()),
		     ARLIM(ww_lo), ARLIM(ww_hi),
		     bc, lo, hi, &dt, dx, &nscal);

}


void
Godunov::edge_sync_rmn (const Box&  grd,
			const Real* dx,
			Real        dt,
			FArrayBox&  uedge,
			FArrayBox&  stx,
			FArrayBox&  kappax,
			FArrayBox&  vedge,
			FArrayBox&  sty,
			FArrayBox&  kappay,
#if (BL_SPACEDIM == 3)               
			FArrayBox&  wedge,
			FArrayBox&  stz,
			FArrayBox&  kappaz,
#endif
			FArrayBox&  S,
			FArrayBox&  S_new,
			FArrayBox&  tforces,
			FArrayBox&  divu,
			FArrayBox&  rock_phi,
			FArrayBox&  kappa,
			FArrayBox&  lbd_cc,
			FArrayBox&  dlbd_cc,
			FArrayBox&  kr_coef,
			const int   n_kr_coef,
			int         fab_ind,
			int         state_ind,
			const int*  bc,
			int         iconserv,
			int         nscal)
{
    
  BL_ASSERT(S.box().contains(work_bx));

  BL_ASSERT(S.nComp()       >= nscal      );
  BL_ASSERT(S_new.nComp()   >= nscal      );

  BL_ASSERT(uedge.nComp()   >= 1          );
  BL_ASSERT(stx.nComp()     >= nscal      );

  BL_ASSERT(vedge.nComp()   >= 1          );
  BL_ASSERT(sty.nComp()     >= nscal      );

#if (BL_SPACEDIM == 3)
  BL_ASSERT(wedge.nComp()   >= 1          );
  BL_ASSERT(stz.nComp()     >= nscal     );
#endif    

  //
  // Create the bounds and pointers.
  //
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();
  const int *s_lo       = S.loVect();
  const int *s_hi       = S.hiVect();
  const int *sn_lo      = S_new.loVect();
  const int *sn_hi      = S_new.hiVect();
  const int *ww_lo      = work_rmn.loVect();
  const int *ww_hi      = work_rmn.hiVect();
  const Real *s_dat     = S.dataPtr(fab_ind);
  const Real *sn_dat    = S_new.dataPtr(fab_ind);
  const Real *tfr_dat   = tforces.dataPtr(fab_ind);
  const Real *divu_dat  = divu.dataPtr();
  const Real *lbd_dat   = lbd_cc.dataPtr();
  const Real *dlbd_dat  = dlbd_cc.dataPtr();

  const Real *kr_dat      = kr_coef.dataPtr();
  const int *kr_lo        = kr_coef.loVect();
  const int *kr_hi        = kr_coef.hiVect();

  //
  // Set work space to bogus values.
  //
  SetBogusScratch();

#if (BL_SPACEDIM == 3)
  const Real *slx_dat   = work_rmn.dataPtr(0*nscal);
  const Real *sly_dat   = work_rmn.dataPtr(1*nscal);
  const Real *slz_dat   = work_rmn.dataPtr(2*nscal);
  const Real *eigvx_dat = work_rmn.dataPtr(3*nscal);
  const Real *eigvy_dat = work_rmn.dataPtr(5*nscal);
  const Real *eigvz_dat = work_rmn.dataPtr(7*nscal);
  const Real *eiglx_dat = work_rmn.dataPtr(9*nscal);
  const Real *eigly_dat = work_rmn.dataPtr(9*nscal+nscal*nscal);
  const Real *eiglz_dat = work_rmn.dataPtr(9*nscal+2*nscal*nscal);
  const Real *eigrx_dat = work_rmn.dataPtr(9*nscal+3*nscal*nscal);
  const Real *eigry_dat = work_rmn.dataPtr(9*nscal+4*nscal*nscal);
  const Real *eigrz_dat = work_rmn.dataPtr(9*nscal+5*nscal*nscal);
#else
  const Real *slx_dat   = work_rmn.dataPtr(0*nscal);
  const Real *sly_dat   = work_rmn.dataPtr(1*nscal);
  const Real *eigvx_dat = work_rmn.dataPtr(2*nscal);
  const Real *eigvy_dat = work_rmn.dataPtr(4*nscal);
  const Real *eiglx_dat = work_rmn.dataPtr(6*nscal);
  const Real *eigly_dat = work_rmn.dataPtr(6*nscal+nscal*nscal);
  const Real *eigrx_dat = work_rmn.dataPtr(6*nscal+2*nscal*nscal);
  const Real *eigry_dat = work_rmn.dataPtr(6*nscal+3*nscal*nscal);
#endif
  //
  // C component indices starts from 0, Fortran from 1
  //
  FORT_SYNC_RMN(s_dat, ARLIM(s_lo), ARLIM(s_hi),
		sn_dat, ARLIM(sn_lo), ARLIM(sn_hi), 
		tfr_dat, divu_dat, lbd_dat, dlbd_dat,
		rock_phi.dataPtr(),ARLIM(rock_phi.loVect()),ARLIM(rock_phi.hiVect()),
		kappa.dataPtr(),ARLIM(kappa.loVect()),ARLIM(kappa.hiVect()),
		slx_dat, slxscr,
		uedge.dataPtr(),ARLIM(uedge.loVect()),ARLIM(uedge.hiVect()),
		stx.dataPtr(),ARLIM(stx.loVect()),ARLIM(stx.hiVect()),
		kappax.dataPtr(),ARLIM(kappax.loVect()),ARLIM(kappax.hiVect()),
		eigvx_dat, eiglx_dat, eigrx_dat,
		sly_dat, slyscr,
		vedge.dataPtr(),ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
		sty.dataPtr(),ARLIM(sty.loVect()),ARLIM(sty.hiVect()),
		kappay.dataPtr(),ARLIM(kappay.loVect()),ARLIM(kappay.hiVect()),
		eigvy_dat, eigly_dat, eigry_dat, 
#if (BL_SPACEDIM == 3)
		slz_dat, slzscr,
		wedge.dataPtr(),ARLIM(wedge.loVect()),ARLIM(wedge.hiVect()),
		stz.dataPtr(),ARLIM(stz.loVect()),ARLIM(stz.hiVect()),
		kappaz.dataPtr(),ARLIM(kappaz.loVect()),ARLIM(kappaz.hiVect()),
		eigvz_dat, eiglz_dat, eigrz_dat,
#endif
		ARLIM(ww_lo), ARLIM(ww_hi),
		kr_dat, ARLIM(kr_lo), ARLIM(kr_hi), &n_kr_coef,
		bc, lo, hi, &dt, dx, 
		&use_forces_in_trans, &iconserv, &nscal);

}

//
// Advect a state component.
// This routine assumes uad,vad,wad have been precomputed.
// FArrayBox work sized as in edge_states.
//

void
Godunov::AdvectState (const Box&  grd,
                      const Real* dx,
                      Real        dt, 
                      FArrayBox&  areax,
                      FArrayBox&  uedge,
                      FArrayBox&  xflux,  
                      FArrayBox&  areay,
                      FArrayBox&  vedge,
                      FArrayBox&  yflux,  
#if (BL_SPACEDIM == 3)                               
                      FArrayBox&  areaz,
                      FArrayBox&  wedge,
                      FArrayBox&  zflux,
#endif
                      FArrayBox&  S,
                      FArrayBox&  tforces,
                      FArrayBox&  divu,
                      int         fab_ind,
                      FArrayBox&  aofs,
                      int         aofs_ind,
                      int         iconserv,
                      int         state_ind,
                      const int*  bc,
                      FArrayBox&  vol)
{
  int velpred = 0;
  //
  // Compute edge states for an advected quantity.
  //
  edge_states(grd, dx, dt, velpred,
	      uedge, xflux,
	      vedge, yflux,
#if (BL_SPACEDIM == 3)             
	      wedge, zflux,
#endif
	      S, tforces, divu, fab_ind, state_ind, bc, 
	      iconserv);
  //
  // Compute the advective tendency.
  //
  ComputeAofs( grd,
	       areax, uedge, xflux,  
	       areay, vedge, yflux,  
#if (BL_SPACEDIM == 3)                             
	       areaz, wedge, zflux,
#endif                     
	       vol, aofs, aofs_ind, iconserv);
}

void
Godunov::AdvectStateLin (const Box&  grd,
			 const Real* dx,
			 Real        dt, 
			 FArrayBox&  areax,
			 FArrayBox&  uedge,
			 FArrayBox&  xflux,
			 FArrayBox&  areay,
			 FArrayBox&  vedge,
			 FArrayBox&  yflux, 
#if (BL_SPACEDIM == 3)                               
			 FArrayBox&  areaz,
			 FArrayBox&  wedge,
			 FArrayBox&  zflux,
#endif
			 FArrayBox&  S,
			 FArrayBox&  S_new,
			 FArrayBox&  tforces,
			 int         fab_ind,
			 FArrayBox&  aofs,
			 int         aofs_ind,
			 FArrayBox&  rock_phi,
			 int         state_ind,
			 const int*  bc,
			 FArrayBox&  vol,
			 int         nscal)
{
  //
  // Compute edge states for an advected quantity.
  //
  edge_states_lin(grd, dx, dt,
		  uedge, xflux, 
		  vedge, yflux, 
#if (BL_SPACEDIM == 3)             
		  wedge, zflux,
#endif
		  S, S_new, tforces,  rock_phi, 
		  fab_ind, state_ind, bc, nscal);
  //
  // Compute the advective tendency.
  //
  ComputeAofsRmn( grd,
		  areax, uedge, xflux,  
		  areay, vedge, yflux,  
#if (BL_SPACEDIM == 3)                             
		  areaz, wedge, zflux,
#endif                     
		  vol, aofs, aofs_ind, nscal);
    
}

void
Godunov::AdvectStateRmn (const Box&  grd,
			 const Real* dx,
			 Real        dt, 
			 FArrayBox&  areax,
			 FArrayBox&  uedge,
			 FArrayBox&  xflux,
			 FArrayBox&  kappax,
			 FArrayBox&  areay,
			 FArrayBox&  vedge,
			 FArrayBox&  yflux, 
			 FArrayBox&  kappay,
#if (BL_SPACEDIM == 3)                               
			 FArrayBox&  areaz,
			 FArrayBox&  wedge,
			 FArrayBox&  zflux,
			 FArrayBox&  kappaz,
#endif
			 FArrayBox&  S,
			 FArrayBox&  S_new,
			 FArrayBox&  tforces,
			 FArrayBox&  divu,
			 int         fab_ind,
			 FArrayBox&  aofs,
			 int         aofs_ind,
			 FArrayBox&  rock_phi,
			 FArrayBox&  kappa,
			 FArrayBox&  lambda_cc,
			 FArrayBox&  dlambda_cc,
			 FArrayBox&  kr_coef,
			 const int   n_kr_coef,
			 int         iconserv,
			 int         state_ind,
			 const int*  bc,
			 FArrayBox&  vol,
			 int         nscal)
{
  //
  // Compute edge states for an advected quantity.
  //
  edge_states_rmn(grd, dx, dt,
		  uedge, xflux, kappax,
		  vedge, yflux, kappay,
#if (BL_SPACEDIM == 3)             
		  wedge, zflux, kappaz,
#endif
		  S, S_new, tforces, divu, rock_phi, kappa,
		  lambda_cc, dlambda_cc, kr_coef, n_kr_coef,
		  fab_ind, state_ind, bc, 
		  iconserv, nscal);
  //
  // Compute the advective tendency.
  //
  ComputeAofsRmn( grd,
		  areax, uedge, xflux,  
		  areay, vedge, yflux,  
#if (BL_SPACEDIM == 3)                             
		  areaz, wedge, zflux,
#endif                     
		  vol, aofs, aofs_ind, nscal);

    
}

void
Godunov::AdvectStateCpl (const Box&  grd,
			 const Real* dx,
			 Real        dt, 
			 FArrayBox&  areax,
			 FArrayBox&  uedge,
			 FArrayBox&  xflux,
			 FArrayBox&  kappax,
			 FArrayBox&  lambdax,
			 FArrayBox&  areay,
			 FArrayBox&  vedge,
			 FArrayBox&  yflux, 
			 FArrayBox&  kappay,
			 FArrayBox&  lambday,
#if (BL_SPACEDIM == 3)                               
			 FArrayBox&  areaz,
			 FArrayBox&  wedge,
			 FArrayBox&  zflux,
			 FArrayBox&  kappaz,
			 FArrayBox&  lambdaz,
#endif
			 FArrayBox&  S,
			 FArrayBox&  S_new,
			 FArrayBox&  tforces,
			 FArrayBox&  divu,
			 int         fab_ind,
			 FArrayBox&  aofs,
			 int         aofs_ind,
			 FArrayBox&  rock_phi,
			 FArrayBox&  kappa,
			 FArrayBox&  pc,
			 FArrayBox&  lambda_cc,
			 FArrayBox&  dlambda_cc,
			 FArrayBox&  kr_coef,
			 const int   n_kr_coef,
			 int         iconserv,
			 int         state_ind,
			 const int*  bc,
			 FArrayBox&  vol,
			 int         nscal)
{
  //
  // Compute edge states for an advected quantity.
  //
  edge_states_cpl(grd, dx, dt,
		  uedge, xflux, kappax, lambdax,
		  vedge, yflux, kappay, lambday,
#if (BL_SPACEDIM == 3)             
		  wedge, zflux, kappaz, lambdaz,
#endif
		  S, S_new, tforces, divu, rock_phi, kappa, pc,
		  lambda_cc, dlambda_cc, kr_coef, n_kr_coef,
		  fab_ind, state_ind, bc, 
		  iconserv, nscal);
  //
  // Compute the advective tendency.
  //
  ComputeAofsRmn( grd,
		  areax, uedge, xflux,  
		  areay, vedge, yflux,  
#if (BL_SPACEDIM == 3)                             
		  areaz, wedge, zflux,
#endif                     
		  vol, aofs, aofs_ind, nscal);

    
}

void
Godunov::AdvectStatePmr (const Box&  grd,
			 const Real* dx,
			 Real        dt, 
			 FArrayBox&  areax,
			 FArrayBox&  uedge,
			 FArrayBox&  xflux,
			 FArrayBox&  kappax,
			 FArrayBox&  areay,
			 FArrayBox&  vedge,
			 FArrayBox&  yflux, 
			 FArrayBox&  kappay,
#if (BL_SPACEDIM == 3)                               
			 FArrayBox&  areaz,
			 FArrayBox&  wedge,
			 FArrayBox&  zflux,
			 FArrayBox&  kappaz,
#endif
			 FArrayBox&  S,
			 FArrayBox&  S_new,
			 FArrayBox&  tforces,
			 FArrayBox&  divu,
			 int         fab_ind,
			 FArrayBox&  aofs,
			 int         aofs_ind,
			 FArrayBox&  rock_phi,
			 FArrayBox&  kappa,
			 int         iconserv,
			 int         state_ind,
			 const int*  bc,
			 FArrayBox&  vol,
			 int         nscal,
			 Real        gravity,
                         Real*       eigmax)
//
// Compute the advective derivative from fluxes.
//
{
  //
  // Compute edge states for an advected quantity.
  //
  edge_states_pmr(grd, dx, dt,
		  uedge, xflux, kappax,
		  vedge, yflux, kappay,
#if (BL_SPACEDIM == 3)             
		  wedge, zflux, kappaz,
#endif
		  S, S_new, tforces, divu, rock_phi, kappa,
		  fab_ind, state_ind, bc, 
		  iconserv, nscal,gravity,eigmax);
  //
  // Compute the advective tendency.
  //
  ComputeAofsRmn( grd,
		  areax, uedge, xflux,  
		  areay, vedge, yflux,  
#if (BL_SPACEDIM == 3)                             
		  areaz, wedge, zflux,
#endif                     
		  vol, aofs, state_ind,nscal);

    
}

void
Godunov::AdvectTracer (const Box&  grd,
		       const Real* dx,
		       Real        dt, 
		       FArrayBox&  areax,
		       FArrayBox&  uedge,
		       FArrayBox&  xflux,
		       FArrayBox&  areay,
		       FArrayBox&  vedge,
		       FArrayBox&  yflux, 
#if (BL_SPACEDIM == 3)                               
		       FArrayBox&  areaz,
		       FArrayBox&  wedge,
		       FArrayBox&  zflux,
#endif
		       FArrayBox&  S,
		       FArrayBox&  S_new,
		       FArrayBox&  St,
		       FArrayBox&  St_new,
		       FArrayBox&  tforces,
		       FArrayBox&  divu,
		       int         fab_ind,
		       FArrayBox&  aofs,
		       int         aofs_ind,
		       FArrayBox&  rock_phi,
		       int         iconserv,
		       int         state_ind,
		       const int*  bc,
		       FArrayBox&  vol,
		       int         nscal)
{
  //
  // Compute edge states for an advected quantity.
  //
  edge_states_tracer(grd, dx, dt,
		     uedge, xflux, xlo, xhi,
		     vedge, yflux, ylo, yhi,
#if (BL_SPACEDIM == 3)             
		     wedge, zflux, zlo, zhi,
#endif
		     S, S_new, St, St_new, rock_phi, 
		     bc, nscal);
  //
  // Compute the advective tendency.
  //
  ComputeAofsRmn(grd,
		 areax, uedge, xflux,  
		 areay, vedge, yflux,  
#if (BL_SPACEDIM == 3)                             
		 areaz, wedge, zflux,
#endif                     
		 vol, aofs, aofs_ind,  nscal);
    
}


void
Godunov::Getdfdn(FArrayBox& S,
		 int state_ind,
		 int ncomps,
		 int dir,
		 int idx)
{
  const int *s_lo       = S.loVect();
  const int *s_hi       = S.hiVect();
  const int *ww_lo      = work_rmn.loVect();
  const int *ww_hi      = work_rmn.hiVect();
  const Real *s_dat     = S.dataPtr(state_ind);

#if (BL_SPACEDIM == 3)
  const Real *eigvx_dat = work_rmn.dataPtr(3*ncomps);
  const Real *eigvy_dat = work_rmn.dataPtr(5*ncomps);
  const Real *eigvz_dat = work_rmn.dataPtr(7*ncomps);
#else
  const Real *eigvx_dat = work_rmn.dataPtr(2*ncomps);
  const Real *eigvy_dat = work_rmn.dataPtr(4*ncomps);
#endif
    
  if (dir == 0) 
    FORT_GET_EIG(s_dat,ARLIM(s_lo),ARLIM(s_hi),eigvx_dat,ARLIM(ww_lo),ARLIM(ww_hi),&ncomps, &idx);
  if (dir == 1) 
    FORT_GET_EIG(s_dat,ARLIM(s_lo),ARLIM(s_hi),eigvy_dat,ARLIM(ww_lo),ARLIM(ww_hi),&ncomps, &idx);
#if (BL_SPACEDIM == 3)
  if (dir == 2) 
    FORT_GET_EIG(s_dat,ARLIM(s_lo),ARLIM(s_hi),eigvz_dat,ARLIM(ww_lo),ARLIM(ww_hi),&ncomps, &idx);
#endif
}

void
Godunov::Getsandc(FArrayBox& S,
		  int state_ind,
		  int ncomps,
		  int idx)
{ 
  S.copy(utmp,idx,state_ind,1);
}


void
Godunov::ComputeAofs (const Box& grd, 
                      FArrayBox& areax,
                      FArrayBox& uedge,
                      FArrayBox& xflux,  
                      FArrayBox& areay,
                      FArrayBox& vedge,
                      FArrayBox& yflux,  
#if (BL_SPACEDIM == 3)                               
                      FArrayBox& areaz,
                      FArrayBox& wedge,
                      FArrayBox& zflux,
#endif
                      FArrayBox& vol,
                      FArrayBox& aofs,
                      int        aofs_ind,
                      int        iconserv )
{
  const int *lo = grd.loVect();
  const int *hi = grd.hiVect();

  FORT_ADV_FORCING(aofs.dataPtr(aofs_ind),ARLIM(aofs.loVect()), ARLIM(aofs.hiVect()),

		   xflux.dataPtr(), ARLIM(xflux.loVect()), ARLIM(xflux.hiVect()),
		   uedge.dataPtr(), ARLIM(uedge.loVect()), ARLIM(uedge.hiVect()),
		   areax.dataPtr(), ARLIM(areax.loVect()), ARLIM(areax.hiVect()),

		   yflux.dataPtr(), ARLIM(yflux.loVect()), ARLIM(yflux.hiVect()),
		   vedge.dataPtr(), ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
		   areay.dataPtr(), ARLIM(areay.loVect()), ARLIM(areay.hiVect()),
#if (BL_SPACEDIM == 3)                                                    
		   zflux.dataPtr(), ARLIM(zflux.loVect()), ARLIM(zflux.hiVect()),
		   wedge.dataPtr(), ARLIM(wedge.loVect()), ARLIM(wedge.hiVect()),
		   areaz.dataPtr(), ARLIM(areaz.loVect()), ARLIM(areaz.hiVect()),
#endif
		   vol.dataPtr(), ARLIM(vol.loVect()), ARLIM(vol.hiVect()),
		   lo, hi, &iconserv);
}

void
Godunov::ComputeAofsRmn (const Box& grd, 
			 FArrayBox& areax,
			 FArrayBox& uedge,
			 FArrayBox& xflux,  
			 FArrayBox& areay,
			 FArrayBox& vedge,
			 FArrayBox& yflux,  
#if (BL_SPACEDIM == 3)                               
			 FArrayBox& areaz,
			 FArrayBox& wedge,
			 FArrayBox& zflux,
#endif
			 FArrayBox& vol,
			 FArrayBox& aofs,
			 int        aofs_ind,
			 int        nscal)
{
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();

  FORT_ADV_RMN_FORCING(aofs.dataPtr(aofs_ind),ARLIM(aofs.loVect()), ARLIM(aofs.hiVect()),
		       xflux.dataPtr(), ARLIM(xflux.loVect()), ARLIM(xflux.hiVect()),
		       uedge.dataPtr(), ARLIM(uedge.loVect()), ARLIM(uedge.hiVect()),
		       areax.dataPtr(), ARLIM(areax.loVect()), ARLIM(areax.hiVect()),
			 
		       yflux.dataPtr(), ARLIM(yflux.loVect()), ARLIM(yflux.hiVect()),
		       vedge.dataPtr(), ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
		       areay.dataPtr(), ARLIM(areay.loVect()), ARLIM(areay.hiVect()),
#if (BL_SPACEDIM == 3)                                                    
		       zflux.dataPtr(), ARLIM(zflux.loVect()), ARLIM(zflux.hiVect()),
		       wedge.dataPtr(), ARLIM(wedge.loVect()), ARLIM(wedge.hiVect()),
		       areaz.dataPtr(), ARLIM(areaz.loVect()), ARLIM(areaz.hiVect()),
#endif
		       vol.dataPtr(), ARLIM(vol.loVect()), ARLIM(vol.hiVect()),
		       lo, hi, &nscal);
			   
}

//
// Sync advect a state component.
// This routine assumes uad,vad,wad have been precomputed.
//

void
Godunov::SyncAdvect (const Box&  grd,
                     const Real* dx,
                     Real        dt,
                     int         level,
                     const FArrayBox& areax,
                     FArrayBox& uedge,
                     FArrayBox& ucorr,
                     FArrayBox& xflux,
                     const FArrayBox& areay,
                     FArrayBox& vedge,
                     FArrayBox& vcorr,
                     FArrayBox& yflux,
#if (BL_SPACEDIM == 3)
                     const FArrayBox& areaz,
                     FArrayBox& wedge,
                     FArrayBox& wcorr,
                     FArrayBox& zflux,
#endif
                     FArrayBox& S,
                     FArrayBox& tforces,
                     FArrayBox& divu,
		     FArrayBox& rock_phi,
                     int        fab_ind,
                     FArrayBox& sync,
                     int        sync_ind,
                     int        iconserv,
                     int        state_ind,
                     const int* bc,
                     const FArrayBox& vol)
{
  int velpred = 0;
  //
  // Error block.
  //
  BL_ASSERT(S.box().contains(work_bx));

  BL_ASSERT(S.nComp()       >= fab_ind    );
  BL_ASSERT(tforces.nComp() >= fab_ind    );
  BL_ASSERT(sync.nComp()    >= sync_ind   );

  BL_ASSERT(ucorr.box()     == xflux_bx   );
  BL_ASSERT(ucorr.nComp()   >= 1          );

  BL_ASSERT(vcorr.box()     == yflux_bx   );
  BL_ASSERT(vcorr.nComp()   >= 1          );
#if (BL_SPACEDIM == 3)
  BL_ASSERT(wcorr.box()     == zflux_bx   );
  BL_ASSERT(wcorr.nComp()   >= 1          );
#endif    
  //
  // Compute the edge states.
  //
  edge_states(grd, dx, dt, velpred,
	      uedge, xflux,
	      vedge, yflux,
#if (BL_SPACEDIM == 3)     
	      wedge, zflux,
#endif
	      S, tforces, divu, fab_ind, state_ind, bc,
	      iconserv);

  //
  // Compute the advective tendency for the mac sync.
  //
  ComputeSyncAofs(grd,
		  areax, ucorr, xflux,  
		  areay, vcorr, yflux,  
#if (BL_SPACEDIM == 3)                             
		  areaz, wcorr, zflux,
#endif                     
		  vol, sync, sync_ind, iconserv);
}

//
// Compute the advective derivative of corrective fluxes for the mac sync.
//

void
Godunov::ComputeSyncAofs (const Box& grd,
                          const FArrayBox& areax,
                          FArrayBox& ucorr,
                          FArrayBox& xflux,  
                          const FArrayBox& areay,
                          FArrayBox& vcorr,
                          FArrayBox& yflux,  
#if (BL_SPACEDIM == 3)                             
                          const FArrayBox& areaz,
                          FArrayBox& wcorr,
                          FArrayBox& zflux,
#endif               
                          const FArrayBox& vol,
                          FArrayBox& sync,
                          int        sync_ind,
                          int        iconserv)
{
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();
  FORT_SYNC_ADV_FORCING(sync.dataPtr(sync_ind), ARLIM(sync.loVect()), ARLIM(sync.hiVect()),
                           
			xflux.dataPtr(),ARLIM(xflux.loVect()),ARLIM(xflux.hiVect()),
			ucorr.dataPtr(),ARLIM(ucorr.loVect()),ARLIM(ucorr.hiVect()),
			areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),

			yflux.dataPtr(),ARLIM(yflux.loVect()),ARLIM(yflux.hiVect()),
			vcorr.dataPtr(),ARLIM(vcorr.loVect()),ARLIM(vcorr.hiVect()),
			areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
#if (BL_SPACEDIM == 3)                                             
			zflux.dataPtr(),ARLIM(zflux.loVect()),ARLIM(zflux.hiVect()),
			wcorr.dataPtr(),ARLIM(wcorr.loVect()),ARLIM(wcorr.hiVect()),
			areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()),
#endif
			vol.dataPtr(), ARLIM(vol.loVect()), ARLIM(vol.hiVect()),
			lo, hi, &iconserv);
}  

void
Godunov::AdvectSyncRmn (const Box&  grd,
			const Real* dx,
			Real        dt, 
			FArrayBox&  areax,
			FArrayBox&  uedge,
			FArrayBox&  ucorr,
			FArrayBox&  xflux,
			FArrayBox&  kappax,
			FArrayBox&  areay,
			FArrayBox&  vedge,
			FArrayBox&  vcorr,
			FArrayBox&  yflux, 
			FArrayBox&  kappay,
#if (BL_SPACEDIM == 3)                               
			FArrayBox&  areaz,
			FArrayBox&  wedge,
			FArrayBox&  wcorr,
			FArrayBox&  zflux,
			FArrayBox&  kappaz,
#endif
			FArrayBox&  S,
			FArrayBox&  S_new,
			FArrayBox&  tforces,
			FArrayBox&  divu,
			int         fab_ind,
			FArrayBox&  aofs,
			int         aofs_ind,
			FArrayBox&  rock_phi,
			FArrayBox&  kappa,
			FArrayBox&  lbd_cc,
			FArrayBox&  dlbd_cc,
			FArrayBox&  kr_coef,
			const int   n_kr_coef,
			int         iconserv,
			int         state_ind,
			const int*  bc,
			FArrayBox&  vol,
			int         nscal)
{
  //
  // Compute edge states for an advected quantity.
  //
  edge_sync_rmn(grd, dx, dt,
		uedge, xflux, kappax,
		vedge, yflux, kappay,
#if (BL_SPACEDIM == 3)             
		wedge, zflux, kappaz,
#endif
		S, S_new, tforces, divu, rock_phi, kappa,
		lbd_cc, dlbd_cc,
		kr_coef, n_kr_coef, 
		fab_ind, state_ind, bc, 
		iconserv, nscal);
  //
  // Compute the advective tendency.  Do not divide by phi
  //
  ComputeSyncAofsRmn( grd,
		      areax, ucorr, xflux,  
		      areay, vcorr, yflux,  
#if (BL_SPACEDIM == 3)                             
		      areaz, wcorr, zflux,
#endif                     
		      vol, kr_coef, n_kr_coef, 
		      aofs, aofs_ind, nscal);
}

void
Godunov::ComputeSyncAofsRmn (const Box& grd,
			     const FArrayBox& areax,
			     FArrayBox& ucorr,
			     FArrayBox& xflux,  
			     const FArrayBox& areay,
			     FArrayBox& vcorr,
			     FArrayBox& yflux,  
#if (BL_SPACEDIM == 3)                             
			     const FArrayBox& areaz,
			     FArrayBox& wcorr,
			     FArrayBox& zflux,
#endif                     
			     const FArrayBox& vol,
			     const FArrayBox& kr_coef,
			     const int   n_kr_coef,
			     FArrayBox& sync,
			     int        sync_ind,
			     int        nscal)
{
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();

  FORT_SYNC_RMN_FORCING(sync.dataPtr(sync_ind), ARLIM(sync.loVect()), ARLIM(sync.hiVect()),
			
			xflux.dataPtr(),ARLIM(xflux.loVect()),ARLIM(xflux.hiVect()),
			ucorr.dataPtr(),ARLIM(ucorr.loVect()),ARLIM(ucorr.hiVect()),
			areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),
			
			yflux.dataPtr(),ARLIM(yflux.loVect()),ARLIM(yflux.hiVect()),
			vcorr.dataPtr(),ARLIM(vcorr.loVect()),ARLIM(vcorr.hiVect()),
			areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
#if (BL_SPACEDIM == 3)                                             
			zflux.dataPtr(),ARLIM(zflux.loVect()),ARLIM(zflux.hiVect()),
			wcorr.dataPtr(),ARLIM(wcorr.loVect()),ARLIM(wcorr.hiVect()),
			areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()),
#endif		
			kr_coef.dataPtr(), ARLIM(kr_coef.loVect()), ARLIM(kr_coef.hiVect()),&n_kr_coef,
			vol.dataPtr(), ARLIM(vol.loVect()), ARLIM(vol.hiVect()),
			lo, hi, &nscal);
}  

  

//
// Correct a conservatively-advected scalar for under-over shoots.
//
void
Godunov::ConservativeScalMinMax (FArrayBox& Sold,
                                 FArrayBox& Snew,
                                 int        ind_old_s, 
                                 int        ind_new_s, 
                                 const int* bc,
                                 const Box& grd)
{
  const int *slo        = Sold.loVect();
  const int *shi        = Sold.hiVect();
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();
  const Real *Sold_dat  = Sold.dataPtr(ind_old_s);
  const Real *Snew_dat  = Snew.dataPtr(ind_new_s);

#if (BL_SPACEDIM == 3)
  Box flatbox(grd);
  int zlen = flatbox.length(BL_SPACEDIM-1);
  flatbox.growHi(BL_SPACEDIM-1,3-zlen);
  FArrayBox smin(flatbox,1);
  FArrayBox smax(flatbox,1);
  const Real *smin_dat = smin.dataPtr();
  const Real *smax_dat = smax.dataPtr(); 
#endif

  FORT_CONSSCALMINMAX (Sold_dat, Snew_dat, 
		       ARLIM(slo), ARLIM(shi),
#if (BL_SPACEDIM == 3)
		       smin_dat, smax_dat,
		       ARLIM(lo), ARLIM(hi),
#endif
		       lo, hi, bc);
}

//
// Correct a convectively-advected scalar for under-over shoots.
//
void
Godunov::ConvectiveScalMinMax (FArrayBox& Sold,
                               FArrayBox& Snew,
                               int        ind_old, 
                               int        ind_new, 
                               const int* bc,
                               const Box& grd)
{
  const int *slo        = Sold.loVect();
  const int *shi        = Sold.hiVect();
  const int *snlo       = Snew.loVect();
  const int *snhi       = Snew.hiVect();
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();
  const Real *Sold_dat  = Sold.dataPtr(ind_old);
  const Real *Snew_dat  = Snew.dataPtr(ind_new);

#if (BL_SPACEDIM == 3)
  Box flatbox(grd);
  int zlen = flatbox.length(BL_SPACEDIM-1);
  flatbox.growHi(BL_SPACEDIM-1,3-zlen);
  FArrayBox smin(flatbox,1);
  FArrayBox smax(flatbox,1);
  const Real *smin_dat = smin.dataPtr();
  const Real *smax_dat = smax.dataPtr(); 
#endif

  FORT_CONVSCALMINMAX (Sold_dat, 
		       ARLIM(slo), ARLIM(shi),
		       Snew_dat,
		       ARLIM(snlo), ARLIM(snhi),
#if (BL_SPACEDIM == 3)
		       smin_dat, smax_dat,
		       ARLIM(lo), ARLIM(hi),
#endif
		       lo, hi, bc);
}

//
// Diagnostic functions follow
//

//
// Estimate the maximum eigenvalues.
//
void
Godunov::esteig (const Box&  grd,
		 const Real* dx,
		 FArrayBox&  uedge,
		 FArrayBox&  kappax,
		 FArrayBox&  vedge,
		 FArrayBox&  kappay,
#if (BL_SPACEDIM == 3)               
		 FArrayBox&  wedge,
		 FArrayBox&  kappaz,
#endif
		 FArrayBox&  S,
		 FArrayBox&  rock_phi,
		 FArrayBox&  kr_coef,		
		 const int   n_kr_coef,
		 const int*  bc,
		 Real*       eigmax)
{
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();

  const Real *s_dat     = S.dataPtr();
  const int *s_lo       = S.loVect();
  const int *s_hi       = S.hiVect();

  const Real *rphi_dat  = rock_phi.dataPtr();
  const int *rphi_lo    = rock_phi.loVect();
  const int *rphi_hi    = rock_phi.hiVect();

  const Real *kr_dat    = kr_coef.dataPtr();
  const int *kr_lo      = kr_coef.loVect();
  const int *kr_hi      = kr_coef.hiVect();

  const Box box = BoxLib::grow(grd,1);
  FArrayBox lbd(box,2);
  FArrayBox dlbd(box,3);
  FArrayBox d2lbd(box,3);
  
  const Real *l_dat   = lbd.dataPtr();
  const Real *dl_dat  = dlbd.dataPtr();
  const Real *d2l_dat = d2lbd.dataPtr();
  const int  *d_lo    = lbd.loVect();
  const int  *d_hi    = lbd.hiVect();

  FORT_EST_EIG(s_dat, ARLIM(s_lo), ARLIM(s_hi),
	       l_dat, dl_dat, d2l_dat, ARLIM(d_lo), ARLIM(d_hi),
	       rphi_dat, ARLIM(rphi_lo), ARLIM(rphi_hi),
	       kr_dat, ARLIM(kr_lo), ARLIM(kr_hi),&n_kr_coef,
	       uedge.dataPtr(),ARLIM(uedge.loVect()),ARLIM(uedge.hiVect()),
	       kappax.dataPtr(),ARLIM(kappax.loVect()),ARLIM(kappax.hiVect()),
	       vedge.dataPtr(),ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
	       kappay.dataPtr(),ARLIM(kappay.loVect()),ARLIM(kappay.hiVect()),
#if (BL_SPACEDIM == 3)
	       wedge.dataPtr(),ARLIM(wedge.loVect()),ARLIM(wedge.hiVect()),
	       kappaz.dataPtr(),ARLIM(kappaz.loVect()),ARLIM(kappaz.hiVect()),
#endif
	       lo, hi, dx, bc, eigmax);
}

//
// Estimate the maximum eigenvalues based on single phase flow
//
void
Godunov::esteig_lin (const Box&  grd,
		     FArrayBox&  uedge,
		     FArrayBox&  vedge,
#if (BL_SPACEDIM == 3)               
		     FArrayBox&  wedge,
#endif
		     FArrayBox&  rock_phi,
		     Real*       eigmax)
{
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();

  const Real *rphi_dat  = rock_phi.dataPtr();
  const int *rphi_lo    = rock_phi.loVect();
  const int *rphi_hi    = rock_phi.hiVect();
  
  FORT_EST_EIG_LIN(rphi_dat, ARLIM(rphi_lo), ARLIM(rphi_hi),
		   uedge.dataPtr(),ARLIM(uedge.loVect()),ARLIM(uedge.hiVect()),
		   vedge.dataPtr(),ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
#if (BL_SPACEDIM == 3)
		   wedge.dataPtr(),ARLIM(wedge.loVect()),ARLIM(wedge.hiVect()),
#endif
		   lo, hi, eigmax);
}

//
// Estimate the maximum eigenvalues based on tracer equation
//
void
Godunov::esteig_trc (const Box&  grd,
		     FArrayBox&  uedge,
		     FArrayBox&  vedge,
#if (BL_SPACEDIM == 3)               
		     FArrayBox&  wedge,
#endif
		     FArrayBox&  S,
		     const int   nc,
		     FArrayBox&  rock_phi,
		     Real*       eigmax)
{
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();

  const Real *s_dat     = S.dataPtr();
  const int *s_lo       = S.loVect();
  const int *s_hi       = S.hiVect();

  const Real *rphi_dat  = rock_phi.dataPtr();
  const int *rphi_lo    = rock_phi.loVect();
  const int *rphi_hi    = rock_phi.hiVect();
  
  FORT_EST_EIG_TRC(s_dat, ARLIM(s_lo), ARLIM(s_hi), &nc,
		   rphi_dat, ARLIM(rphi_lo), ARLIM(rphi_hi),
		   uedge.dataPtr(),ARLIM(uedge.loVect()),ARLIM(uedge.hiVect()),
		   vedge.dataPtr(),ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
#if (BL_SPACEDIM == 3)
		   wedge.dataPtr(),ARLIM(wedge.loVect()),ARLIM(wedge.hiVect()),
#endif
		   lo, hi, eigmax);
}

//
// Estimate the maximum eigenvalues with cpl term
//
void
Godunov::esteig_cpl (const Box&  grd,
		     const Real* dx,
		     FArrayBox&  uedge,
		     FArrayBox&  kappax,
		     FArrayBox&  vedge,
		     FArrayBox&  kappay,
#if (BL_SPACEDIM == 3)               
		     FArrayBox&  wedge,
		     FArrayBox&  kappaz,
#endif
		     FArrayBox&  S,
		     FArrayBox&  pc,
		     FArrayBox&  rock_phi,
		     FArrayBox&  kr_coef,		
		     const int   n_kr_coef,
		     const int*  bc,
		     Real*       eigmax)
{
  const int *lo         = grd.loVect();
  const int *hi         = grd.hiVect();

  const Real *s_dat     = S.dataPtr();
  const int *s_lo       = S.loVect();
  const int *s_hi       = S.hiVect();

  const Real *pc_dat    = pc.dataPtr();
  const int *pc_lo      = pc.loVect();
  const int *pc_hi      = pc.hiVect();

  const Real *rphi_dat  = rock_phi.dataPtr();
  const int *rphi_lo    = rock_phi.loVect();
  const int *rphi_hi    = rock_phi.hiVect();

  const Real *kr_dat    = kr_coef.dataPtr();
  const int *kr_lo      = kr_coef.loVect();
  const int *kr_hi      = kr_coef.hiVect();

  const Box box = BoxLib::grow(grd,1);
  FArrayBox lbd(box,2);
  FArrayBox dlbd(box,3);
  FArrayBox d2lbd(box,3);
  
  const Real *l_dat   = lbd.dataPtr();
  const Real *dl_dat  = dlbd.dataPtr();
  const Real *d2l_dat = d2lbd.dataPtr();
  const int  *d_lo    = lbd.loVect();
  const int  *d_hi    = lbd.hiVect();

  FORT_EST_EIG_CPL(s_dat, ARLIM(s_lo), ARLIM(s_hi),
	       l_dat, dl_dat, d2l_dat, ARLIM(d_lo), ARLIM(d_hi),
	       rphi_dat, ARLIM(rphi_lo), ARLIM(rphi_hi),
	       kr_dat, ARLIM(kr_lo), ARLIM(kr_hi),&n_kr_coef,
	       pc_dat, ARLIM(pc_lo), ARLIM(pc_hi),
	       uedge.dataPtr(),ARLIM(uedge.loVect()),ARLIM(uedge.hiVect()),
	       kappax.dataPtr(),ARLIM(kappax.loVect()),ARLIM(kappax.hiVect()),
	       vedge.dataPtr(),ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
	       kappay.dataPtr(),ARLIM(kappay.loVect()),ARLIM(kappay.hiVect()),
#if (BL_SPACEDIM == 3)
	       wedge.dataPtr(),ARLIM(wedge.loVect()),ARLIM(wedge.hiVect()),
	       kappaz.dataPtr(),ARLIM(kappaz.loVect()),ARLIM(kappaz.hiVect()),
#endif
	       lo, hi, dx, bc, eigmax);
}

//
// Source term functions follow
//

void
Godunov::Add_aofs_tf (FArrayBox& Sold,
                      FArrayBox& Snew,
                      int        start_ind,
                      int        num_comp,
                      FArrayBox& Aofs,
                      int        aofs_ind,
		      FArrayBox& tforces,
                      int        tf_ind,
		      FArrayBox& Rockphi,
                      const Box& grd,
                      Real       dt)
{
  BL_ASSERT(Snew.nComp()    >= start_ind + num_comp);
  BL_ASSERT(Sold.nComp()    >= start_ind + num_comp);
  BL_ASSERT(Aofs.nComp()    >= aofs_ind  + num_comp);

  const int *slo    = Sold.loVect();
  const int *shi    = Sold.hiVect();
  const int *alo    = Aofs.loVect();
  const int *ahi    = Aofs.hiVect();    
  const int *tlo    = tforces.loVect();
  const int *thi    = tforces.hiVect();
  const int *plo    = Rockphi.loVect();
  const int *phi    = Rockphi.hiVect();
  const int *lo     = grd.loVect();
  const int *hi     = grd.hiVect();
  const Real *SOdat = Sold.dataPtr(start_ind);
  const Real *SNdat = Snew.dataPtr(start_ind);
  const Real *AOdat = Aofs.dataPtr(aofs_ind);
  const Real *TFdat = tforces.dataPtr(tf_ind);
  const Real *pdat  = Rockphi.dataPtr();
    
  FORT_UPDATE_AOFS_TF(SOdat, ARLIM(slo), ARLIM(shi), 
		      SNdat, ARLIM(slo), ARLIM(shi),
		      AOdat, ARLIM(alo), ARLIM(ahi),
		      TFdat, ARLIM(tlo), ARLIM(thi),
		      pdat, ARLIM(plo), ARLIM(phi),
		      lo, hi, &dt, &num_comp);
}

void
Godunov::Add_aofs_tracer (FArrayBox& Sold,
			  FArrayBox& Snew,
			  int        start_ind,
			  int        num_comp,
			  FArrayBox& Aofs,
			  int        aofs_ind,
			  FArrayBox& tforces,
			  int        tf_ind,
			  FArrayBox& Rockphi,
			  const Box& grd,
			  Array<int> idx_total,
			  Real       dt)
{
  BL_ASSERT(Snew.nComp()    >= start_ind + num_comp);
  BL_ASSERT(Sold.nComp()    >= start_ind + num_comp);
  BL_ASSERT(Aofs.nComp()    >= aofs_ind  + num_comp);

  const int *slo    = Sold.loVect();
  const int *shi    = Sold.hiVect();
  const int *alo    = Aofs.loVect();
  const int *ahi    = Aofs.hiVect();    
  const int *tlo    = tforces.loVect();
  const int *thi    = tforces.hiVect();
  const int *plo    = Rockphi.loVect();
  const int *phi    = Rockphi.hiVect();
  const int *lo     = grd.loVect();
  const int *hi     = grd.hiVect();
  const int idx_n   = idx_total.size();
  const Real *SOdat = Sold.dataPtr(start_ind);
  const Real *SNdat = Snew.dataPtr(start_ind);
  const Real *AOdat = Aofs.dataPtr(aofs_ind);
  const Real *TFdat = tforces.dataPtr(tf_ind);
  const Real *pdat  = Rockphi.dataPtr();
    
    
  FORT_UPDATE_AOFS_TRACER(SOdat, ARLIM(slo), ARLIM(shi), 
			  SNdat, ARLIM(slo), ARLIM(shi),
			  AOdat, ARLIM(alo), ARLIM(ahi),
			  TFdat, ARLIM(tlo), ARLIM(thi),
			  pdat, ARLIM(plo), ARLIM(phi),
			  lo, hi, idx_total.dataPtr(), &idx_n,
			  &dt, &num_comp);
}

//
// Compute total source term for scalars.  Note for compatibility
// The switch iconserv, determines the form of the total source term
//
// iconserv==1   => tforces = tforces + visc - divU*S
//
// iconserv==0   => tforces = (tforces+ visc)
//

void
Godunov::Sum_tf_divu_visc (FArrayBox& S,
                           FArrayBox& tforces,
                           int        s_ind,
                           int        num_comp,
                           FArrayBox& visc,
                           int        v_ind,
                           FArrayBox& divu,
                           int        iconserv)
{
  BL_ASSERT(S.nComp()       >= s_ind+num_comp);
  BL_ASSERT(divu.nComp()    == 1             );
  BL_ASSERT(visc.nComp()    >= v_ind+num_comp);
    
  const int *slo    = S.loVect();
  const int *shi    = S.hiVect();
  const int *tlo    = tforces.loVect();
  const int *thi    = tforces.hiVect();
  const int *dlo    = divu.loVect();
  const int *dhi    = divu.hiVect();
  const int *vlo    = visc.loVect();
  const int *vhi    = visc.hiVect();
  const Real *Sdat  = S.dataPtr(s_ind);
  const Real *TFdat = tforces.dataPtr(s_ind);
  const Real *DUdat = divu.dataPtr();
  const Real *VIdat = visc.dataPtr(v_ind);
     
  FORT_SUM_TF_DIVU_VISC(Sdat,  ARLIM(slo), ARLIM(shi),
			TFdat, ARLIM(tlo), ARLIM(thi),
			DUdat, ARLIM(dlo), ARLIM(dhi),
			VIdat, ARLIM(vlo), ARLIM(vhi),
			tlo, thi, &num_comp, &iconserv);
}

bool
Godunov::are_any(const Array<AdvectionForm>& advectionType,
                 const AdvectionForm         testForm,
                 const int                   sComp,
                 const int                   nComp)
{
  for (int comp = sComp; comp < sComp + nComp; ++comp)
    {
      if (advectionType[comp] == testForm)
	return true;
    }

  return false;
}

int
Godunov::how_many(const Array<AdvectionForm>& advectionType,
                  const AdvectionForm         testForm,
                  const int                   sComp,
                  const int                   nComp)
{
  int counter = 0;

  for (int comp = sComp; comp < sComp + nComp; ++comp)
    {
      if (advectionType[comp] == testForm)
	++counter;
    }

  return counter;
}
