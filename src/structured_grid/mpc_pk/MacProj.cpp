//
// $Id: MacProj.cpp,v 1.77 2011-07-02 15:17:03 gpau Exp $
//
#include <winstd.H>

#include <LO_BCTYPES.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <MacProj.H>
#include <MacBndry.H>
#include <MacOpMacDrivers.H>
#include <PorousMedia.H>
#include <MACPROJ_F.H>
#include <MACOPERATOR_F.H>

#ifndef _PorousMedia_H_
enum StateType {State_Type=0, Press_Type, Vel_Type};
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)	\
  const int* fablo  = (fab).loVect();           \
  const int* fabhi  = (fab).hiVect();           \
  Real*      fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)	\
  const int* fablo   = (fab).loVect();		\
  const int* fabhi   = (fab).hiVect();		\
  const Real* fabdat = (fab).dataPtr();

#define DEF_BOX_LIMITS(box,boxlo,boxhi)		\
  const int* boxlo = (box).loVect();		\
  const int* boxhi = (box).hiVect();

#define GEOM_GROW 1
#define HYP_GROW 3

int  MacProj::verbose          = 0;
int  MacProj::do_any_diffuse   = 0;
bool MacProj::use_cg_solve     = false;
namespace
{
  bool        use_hypre_solve  = false;
  bool        use_fboxlib_mg   = false;
}
Real MacProj::mac_tol          = 1.0e-12;
Real MacProj::mac_abs_tol      = 1.0e-16;
Real MacProj::mac_sync_tol     = 1.0e-12;

//
// Setup functions follow
//

MacProj::MacProj (Amr*   _parent,
                  int    _finest_level,
                  BCRec* _phys_bc,
		  int    _do_any_diffuse)
  :
  parent(_parent),
  LevelData(_finest_level+1),
  phys_bc(_phys_bc), 
  phi_bcs(_finest_level+1),
  mac_reg(_finest_level+1, PArrayManage),
  mac_reg_rhoD(_finest_level+1, PArrayManage),
  volume(_finest_level+1),
  area(_finest_level+1),
  finest_level(_finest_level)
{
  read_params();
  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "Creating mac_projector with finest level " 
	      << finest_level << std::endl;

  do_any_diffuse = _do_any_diffuse;

  finest_level_allocated = finest_level;
}

MacProj::~MacProj () {}

void
MacProj::read_params ()
{
  ParmParse pp("mac");

  pp.query( "v",                verbose          );
  pp.query( "mac_tol",          mac_tol          );
  pp.query( "mac_sync_tol",     mac_sync_tol     );
  pp.query( "use_cg_solve",     use_cg_solve     );
#if MG_USE_HYPRE
  pp.query( "use_hypre_solve",  use_hypre_solve);
#endif
#if MG_USE_FBOXLIB
  pp.query( "use_fboxlib_mg",   use_fboxlib_mg);
#endif
  pp.query( "mac_abs_tol",      mac_abs_tol      );

  if ( use_cg_solve && use_hypre_solve )
    {
      BoxLib::Error("MacProj::read_params: cg_solve && .not. hypre_solve");
    }
  if ( use_cg_solve && use_fboxlib_mg )
    {
      BoxLib::Error("MacProj::read_params: cg_solve && .not. fboxlib_solve");
    }
}

void
MacProj::install_level (int                   level,
                        AmrLevel*             level_data,
                        MultiFab&             _volume,
                        MultiFab*             _area)
{
  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "Installing MacProj level " << level << '\n';

  if (parent->finestLevel() < finest_level)
    for (int lev = parent->finestLevel() + 1; lev <= finest_level; lev++)
      {
	mac_reg.clear(lev);
	mac_reg_rhoD.clear(lev);
      }

  finest_level = parent->finestLevel();

  if (level > finest_level_allocated)
    {
      finest_level_allocated = finest_level;
      LevelData.resize(finest_level+1);
      phi_bcs.resize(finest_level+1);
      mac_reg.resize(finest_level+1);
      mac_reg_rhoD.resize(finest_level+1);
      volume.resize(finest_level+1);
      area.resize(finest_level+1);
    }

  LevelData.clear(level);
  LevelData.set(level, level_data);
  volume.clear(level);
  volume.set(level, &_volume);
  area.set(level, _area);

  BuildPhiBC(level);

  if (level > 0)
    {
      mac_reg.clear(level);
      mac_reg.set(level,new FluxRegister(LevelData[level].boxArray(),
					 parent->refRatio(level-1),level,1));
      mac_reg_rhoD.clear(level);
      mac_reg_rhoD.set(level,new FluxRegister(LevelData[level].boxArray(),
					      parent->refRatio(level-1),level,1));
    }
}
void
MacProj::BuildPhiBC (int level)
{
  const BoxArray& grids   = LevelData[level].boxArray();
  const Geometry& geom    = parent->Geom(level);
  const int       ngrds   = grids.size();
  phi_bcs[level].resize(ngrds);
  const Box&      domain  = geom.Domain();
  const int*      domlo   = domain.loVect();
  const int*      domhi   = domain.hiVect();
  const int*      phys_lo = phys_bc->lo();
  const int*      phys_hi = phys_bc->hi();

  for (int i = 0; i < ngrds; i++)
    {
      BCRec&     bc = phi_bcs[level][i];
      const int* lo = grids[i].loVect();
      const int* hi = grids[i].hiVect();

      for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
	  if (lo[dir] == domlo[dir])
            {
	      //bc.setLo(dir,((phys_lo[dir]==Outflow) ||(phys_lo[dir]==Inflow)) ? LO_DIRICHLET : LO_NEUMANN);
	      bc.setLo(dir,((phys_lo[dir]==EXT_DIR)) ? LO_DIRICHLET : LO_NEUMANN);
            }
	  else
            {
	      bc.setLo(dir,LO_DIRICHLET);
            }
	  if (hi[dir] == domhi[dir])
            {
	      //bc.setHi(dir,((phys_hi[dir]==Outflow) ||(phys_hi[dir]==Inflow)) ? LO_DIRICHLET : LO_NEUMANN);
	      bc.setHi(dir,((phys_hi[dir]==EXT_DIR)) ? LO_DIRICHLET : LO_NEUMANN);
            }
	  else
            {
	      bc.setHi(dir,LO_DIRICHLET);
            }
        }
    }
}


//
// Projection functions follow ...
//
static
bool
grids_on_side_of_domain (const BoxArray&    grids,
                         const Box&         domain,
                         const Orientation& outFace)
{
  const int idir = outFace.coordDir();

  if (outFace.isLow())
    {
      for (int igrid = 0; igrid < grids.size(); igrid++)
        { 
	  if (grids[igrid].smallEnd(idir) == domain.smallEnd(idir))
            { 
	      return true;
            }
        }
    }
  
  if (outFace.isHigh())
    {
      for (int igrid = 0; igrid < grids.size(); igrid++)
        {
	  if (grids[igrid].bigEnd(idir) == domain.bigEnd(idir))
            {
	      return true;
            }
        }
    }

  return false;
}

//
// Compute the level advance mac projection.
//

void
MacProj::mac_project (int             level,
                      MultiFab*       u_mac,
                      MultiFab*       mac_coef,
		      MultiFab*       RhoD,
		      MultiFab*       forces,
                      MultiFab*       phi,
                      const MacBndry& mac_bndry,
		      const BCRec&    p_bc)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::mac_project()");

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... mac_project at level " << level << '\n';

  const BoxArray& grids      = LevelData[level].boxArray();
  const Geometry& geom       = parent->Geom(level);
  const Real*     dx         = geom.CellSize();
    
  //
  // Compute the nondivergent velocities, by creating the linop
  // and multigrid operator appropriate for the solved system.
  //

  MultiFab Rhs(grids,1,0);
  Rhs.setVal(0.);

  if (forces != 0) 
    {
      for (int kk=0; kk < forces->nComp(); kk++)
	MultiFab::Add(Rhs,*forces,kk,0,1,0);
    }

  int the_solver = 0;

  if (use_cg_solve)
    {
      the_solver = 1;
    }
  else if ( use_fboxlib_mg )
    {
      the_solver = 3;
    }

  const Box& domain = parent->Geom(level).Domain();
  const Real rhs_scale = 1.0;
  mac_level_driver(mac_bndry, p_bc, grids, the_solver, level,
		   dx, mac_tol, mac_abs_tol, rhs_scale, 
		   area[level], volume[level], mac_coef, Rhs, 
		   u_mac, RhoD, phi, domain, verbose);

  //
  // Test that u_mac is divergence free
  //
  if (verbose == 2)
    check_div_cond(level, u_mac, RhoD);

  print_min_max(u_mac);
}

void
MacProj::contribute_to_mac_reg (int       level,
                                MultiFab* u_mac)
{
  //
  // Store advection velocities in mac registers at crse/fine boundaries.
  //
  // Initialize advection velocity registers with coarse grid velocity.
  //
  if (level < finest_level)
    {
      FluxRegister& mr = mac_reg[level+1];

      mr.setVal(0.0);

      for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
	  mr.CrseInit(u_mac[dir],area[level][dir],dir,0,0,1,-1.0);
        }
    }
  //
  // Increment in fine grid velocity to velocity registers.
  //
  if (level > 0)
    {
      const Real mult = 1.0/parent->nCycle(level);

      for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
	  mac_reg[level].FineAdd(u_mac[dir],area[level][dir],dir,0,0,1,mult);
        }
    }
}

void
MacProj::contribute_to_mac_reg_rhoD (int       level,
				     MultiFab* rhod)
{
  //
  // Store advection velocities in mac registers at crse/fine boundaries.
  //
  // Initialize advection velocity registers with coarse grid velocity.
  //
  if (level < finest_level)
    {
      FluxRegister& mr = mac_reg_rhoD[level+1];

      mr.setVal(0.0);

      for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
	  mr.CrseInit(rhod[dir],area[level][dir],dir,0,0,1,-1.0);
        }
    }
  //
  // Increment in fine grid velocity to velocity registers.
  //
  if (level > 0)
    {
      const Real mult = 1.0/parent->nCycle(level);

      for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
	  mac_reg_rhoD[level].FineAdd(rhod[dir],area[level][dir],dir,0,0,1,mult);
        }
    }
}

void
MacProj::mac_sync_solve (int       level,
			 const BCRec& p_bc,
                         MultiFab* mac_coef,
                         MultiFab* mac_sync_phi,
			 MultiFab* mac_sync_u,
                         IntVect&  fine_ratio)
{

  BL_ASSERT(level < finest_level);

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... mac_sync_solve at level " << level << '\n';

  const BoxArray& grids      = LevelData[level].boxArray();
  const Geometry& geom       = parent->Geom(level);
  const Real*     dx         = geom.CellSize();
  const BoxArray& fine_boxes = LevelData[level+1].boxArray();
  IntVect         crse_ratio = level > 0 ? parent->refRatio(level-1)
    : IntVect::TheZeroVector();
  //
  // Alloc and define RHS by doing a reflux-like operation in coarse
  // grid cells adjacent to fine grids.  The values in these
  // cells should be SUM{MR/VOL} where the sum is taken over
  // all edges of a cell that adjoin fine grids, MR = value in
  // MAC register, VOL = cell volume.  All other cells have a
  // value of zero (including crse cells under fine grids).
  //
  MultiFab Rhs(grids,1,0);
  Rhs.setVal(0.0);
  //
  // Reflux subtracts values at hi edge of coarse cell and
  // adds values at lo edge.  We want the opposite here so
  // set scale to -1 & alloc space for Rhs.
  //
  FluxRegister& mr = mac_reg[level+1];
  const Real scale = 1.0;

  mr.Reflux(Rhs,volume[level],scale,0,0,1,geom);
  MultiFab* rhod_tmp;
  rhod_tmp = new MultiFab[BL_SPACEDIM];

  for (int dir=0; dir < BL_SPACEDIM; dir++)
    {
      BoxArray edge_grids(grids);
      edge_grids.surroundingNodes(dir);
      rhod_tmp[dir].define(edge_grids,1,0,Fab_allocate);
      rhod_tmp[dir].setVal(0.0);
    }

  if (do_any_diffuse)
    {  
      FluxRegister& mrd = mac_reg_rhoD[level+1];

      for (int dir=0; dir < BL_SPACEDIM; dir++)
	{
	  const Orientation lo_face(dir,Orientation::low);
	  const Orientation hi_face(dir,Orientation::high);
	
	  mrd[lo_face].copyTo(rhod_tmp[dir]);
	  mrd[hi_face].copyTo(rhod_tmp[dir]);
	}
    }

  for (int kf = 0, nfine = fine_boxes.size(); kf < nfine; kf++)
    {
      Box bf = BoxLib::coarsen(fine_boxes[kf],fine_ratio);

      for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
        {
	  BL_ASSERT(grids[Rhsmfi.index()] == Rhsmfi.validbox());

	  Box isect = Rhsmfi.validbox() & bf;

	  if (isect.ok())
            {
	      Rhs[Rhsmfi].setVal(0.0,isect,0);
            }
        }
    }
    
  mac_sync_phi->setVal(0.0);
  //
  // store the Dirichlet boundary condition for mac_sync_phi in mac_bndry
  //
  MacBndry mac_bndry(grids,1,geom);
  const int src_comp = 0;
  const int dest_comp = 0;
  const int num_comp = 1;
  if (level == 0)
    mac_bndry.setBndryValues(*mac_sync_phi,src_comp,dest_comp,num_comp,
			     p_bc);
  else
    {
      BoxArray crse_boxes(grids);
      crse_boxes.coarsen(crse_ratio);
      const int in_rad     = 0;
      const int out_rad    = 1;
      const int extent_rad = 2;
      BndryRegister crse_br(crse_boxes,in_rad,out_rad,extent_rad,num_comp);
      crse_br.setVal(0);
      mac_bndry.setBndryValues(crse_br,src_comp,*mac_sync_phi,src_comp,
			       dest_comp,num_comp,crse_ratio, p_bc);
    }
    
  //
  // Now define edge centered coefficients and adjust RHS for MAC solve.
  // 

  //
  // Solve the sync system.
  //
  int the_solver = 0;
  if (use_cg_solve)
    {
      the_solver = 1;
    }
  else if ( use_fboxlib_mg )
    {
      the_solver = 3;
    }
  const Box&  domain = parent->Geom(level).Domain();
  const Real rhs_scale = 1.0;
  mac_sync_driver(mac_bndry, p_bc, grids, the_solver, level, dx, 
		  mac_sync_tol, mac_abs_tol, rhs_scale, area[level],
		  volume[level], Rhs, mac_coef, mac_sync_phi, 
		  mac_sync_u, rhod_tmp, domain,verbose);
  delete [] rhod_tmp;

}

//
// After solving for mac_sync_phi in mac_sync_solve(), we
// can now do the sync advect step.  This consists of two steps
//
// 1. compute u_corr as the gradient of mac_sync_phi
// 2. compute advective tendency of u_corr and
//    add into Ssync
//
// If increment_sync is non-null, the (i-BL_SPACEDIM)-th component 
// of (*Ssync) is incremented only when increment[i]==1
// This is useful if that component gets incrmnted in a non-standard way.
//

void
MacProj::mac_sync_compute (int                   level,
                           MultiFab*             u_mac, 
			   MultiFab*             u_corr,
                           MultiFab*             Ssync,
                           MultiFab*             mac_coef,
			   MultiFab*             rock_phi,
			   MultiFab*             kappa,
			   MultiFab*             lambda_cc,
			   MultiFab*             dlambda_cc,
			   MultiFab*             kr_coef,
			   MultiFab*             kpedge,
                           MultiFab*             mac_sync_phi,
                           FluxRegister*         adv_flux_reg,
                           Array<AdvectionForm>& advectionType,
                           Real                  prev_time,
                           Real                  dt, 
			   int                   nscal,
			   Real                  be_cn_theta)
{
  //
  // Get parameters.
  //
  const BoxArray& grids               = LevelData[level].boxArray();
  const Geometry& geom                = parent->Geom(level);
  const Real*     dx                  = geom.CellSize();
  PorousMedia&    pm_level            = *(PorousMedia*) &(parent->getLevel(level));
  Godunov*        godunov             = pm_level.godunov;
  int model = pm_level.model;

  MultiFab* rhod_tmp;
  rhod_tmp = new MultiFab[BL_SPACEDIM];
  for (int dir=0; dir < BL_SPACEDIM; dir++)
    {

      BoxArray edge_grids(grids);
      edge_grids.surroundingNodes(dir);
      rhod_tmp[dir].define(edge_grids,1,0,Fab_allocate);
      rhod_tmp[dir].setVal(0.0);
    }

  if (do_any_diffuse)
    {
      FluxRegister& mrd = mac_reg_rhoD[level+1];
      for (int dir=0; dir < BL_SPACEDIM; dir++)
	{
	  const Orientation lo_face(dir,Orientation::low);
	  const Orientation hi_face(dir,Orientation::high);
	
	  mrd[lo_face].copyTo(rhod_tmp[dir]);
	  mrd[hi_face].copyTo(rhod_tmp[dir]);
	}
    }
    
  // 
  // Compute correction velocity
  // Note that u_corr = \lambda \nabla p_corr.  No minus sign.
  //
  Real vel_update_scale = 1.0;
  for (MFIter mfi(*mac_sync_phi); mfi.isValid(); ++mfi)
    {
      const int i     = mfi.index();
      mac_vel_sync(D_DECL(u_corr[0][i],u_corr[1][i],
			  u_corr[2][i]),
		   D_DECL(rhod_tmp[0][i],rhod_tmp[1][i],
			  rhod_tmp[2][i]),
		   D_DECL(area[level][0][i],area[level][1][i],
			  area[level][2][i]),
		   D_DECL(mac_coef[0][i],mac_coef[1][i],
			  mac_coef[2][i]),
		   (*mac_sync_phi)[mfi],grids[i], dx,vel_update_scale);
    }
  
  MultiFab* u_macG;
  u_macG = new MultiFab[BL_SPACEDIM];
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      BoxArray edge_grids(grids);
      edge_grids.surroundingNodes(dir).grow(1);
      u_macG[dir].define(edge_grids,1,0,Fab_allocate);
      u_macG[dir].copy(u_corr[dir]);

      // There are two options here: 
      //  I can either have u_mac + u_corr or just u_mac.
      //  The choice is not exactly clear.  I believe latter is a
      //  better option.
      if (model == 2)
	MultiFab::Copy(u_macG[dir],u_mac[dir],0,0,1,0);

    }

  //
  // Compute the component correction.
  //
  //
  // Get viscous forcing.
  //
    
  MultiFab scal_visc_terms(grids,nscal,1);
  scal_visc_terms.setVal(0.);
  if (be_cn_theta != 1.0) 
    pm_level.getViscTerms(scal_visc_terms,0,nscal,prev_time);

  Array<int> pm_level_bc;
  MultiFab* divu_fp = new MultiFab(grids,1,1);
  (*divu_fp).setVal(0.);

  FArrayBox tforces;
  FArrayBox flux[BL_SPACEDIM];
  MultiFab& S_new = pm_level.get_new_data(State_Type); 

  for (FillPatchIterator S_fpi(pm_level,scal_visc_terms,HYP_GROW,
			       prev_time,State_Type,0,nscal);
       S_fpi.isValid(); ++S_fpi)
    {
       pm_level.dirichletStateBC(S_fpi(),HYP_GROW,prev_time);
       const int i     = S_fpi.index();

      int state_ind     = 0;
      pm_level_bc       = pm_level.getBCArray(State_Type,i,state_ind,1);
      int use_conserv_diff = (advectionType[state_ind] == Conservative);

      // tforces is set to zero in PorousMedia::getForce, using prev_time, need to check
      pm_level.getForce(tforces,i,1,0,nscal,prev_time);
	
      // Compute total forcing terms.
      // iconserve set to 1.  
      godunov->Sum_tf_divu_visc(S_fpi(), tforces, state_ind, nscal,
				scal_visc_terms[i], state_ind, 
				(*divu_fp)[i], use_conserv_diff);

        
      // Set up the workspace for the godunov Box.
      
      godunov->Setup(grids[i], flux[0], flux[1], 
#if (BL_SPACEDIM == 3)
		     flux[2], 
#endif
		     nscal, model);

      if (model == 2)
	{
	  const int n_kr_coef = kr_coef->nComp();
	  godunov->AdvectSyncRmn(grids[i], dx, dt, 
				 area[level][0][i], u_macG[0][i], u_corr[0][i], 
				 flux[0], kpedge[0][i],
				 area[level][1][i], u_macG[1][i], u_corr[1][i], 
				 flux[1], kpedge[1][i],
#if (BL_SPACEDIM == 3)                        
				 area[level][2][i], u_macG[2][i], u_corr[2][i], 
				 flux[2], kpedge[2][i],
#endif
				 S_fpi(),S_new[i],tforces,
				 (*divu_fp)[i] , state_ind,
				 (*Ssync)[i]   , state_ind,
				 (*rock_phi)[i], (*kappa)[i],
				 (*lambda_cc)[i], (*dlambda_cc)[i],
				 (*kr_coef)[i],n_kr_coef,
				 use_conserv_diff,
				 state_ind,pm_level_bc.dataPtr(),volume[level][i],
				 nscal);
	}
      else
	{
	  godunov->AdvectStateLin(grids[i], dx, dt, 
				  area[level][0][i], u_macG[0][i], flux[0],
				  area[level][1][i], u_macG[1][i], flux[1], 
#if (BL_SPACEDIM == 3)                        
				  area[level][2][i], u_macG[2][i], flux[2],
#endif
				  S_fpi(),S_new[i],tforces, state_ind,
				  (*Ssync)[i]   , state_ind,
				  (*rock_phi)[i], state_ind,
				  pm_level_bc.dataPtr(),volume[level][i],nscal);	
	}
	
      if (level > 0)
	{
	  for (int dir = 0; dir < BL_SPACEDIM; dir++)
	    adv_flux_reg->FineAdd(flux[dir],dir,i,0,state_ind,nscal,-dt);

	  const Real mlt =  -1.0/( (double) parent->nCycle(level));
	  for (int dir = 0; dir < BL_SPACEDIM; dir++)
            mac_reg[level].FineAdd(u_corr[dir][i],area[level][dir][i],
				   dir,i,0,0,1,mlt);
        }
    }

  if (verbose > 2) 
    check_div_cond(level,u_corr,rhod_tmp);

  delete divu_fp;
  delete [] rhod_tmp;
  delete [] u_macG;
}

//
// Check the mac divergence.
//
void
MacProj::check_div_cond (int      level,
                         MultiFab U_edge[],
			 MultiFab RhoD[]) const
{
  const BoxArray& grids = LevelData[level].boxArray();

  Real sum = 0.0;

  FArrayBox dmac;
  Real npts = 0.0;
  for (MFIter U_edge0mfi(U_edge[0]); U_edge0mfi.isValid(); ++U_edge0mfi)
    {
      dmac.resize(grids[U_edge0mfi.index()],1);

      const FArrayBox& uxedge = U_edge[0][U_edge0mfi];
      const FArrayBox& uyedge = U_edge[1][U_edge0mfi];
      const FArrayBox& rhoDX  = RhoD[0][U_edge0mfi];
      const FArrayBox& rhoDY  = RhoD[1][U_edge0mfi];
      const FArrayBox& xarea  = area[level][0][U_edge0mfi];
      const FArrayBox& yarea  = area[level][1][U_edge0mfi];
      const FArrayBox& vol    = volume[level][U_edge0mfi];

      DEF_LIMITS(dmac,dmac_dat,dlo,dhi);
      DEF_CLIMITS(uxedge,ux_dat,uxlo,uxhi);
      DEF_CLIMITS(uyedge,uy_dat,uylo,uyhi);
      DEF_CLIMITS(rhoDX,rx_dat,rxlo,rxhi);
      DEF_CLIMITS(rhoDY,ry_dat,rylo,ryhi);
      DEF_CLIMITS(xarea,ax_dat,axlo,axhi);
      DEF_CLIMITS(yarea,ay_dat,aylo,ayhi);
      DEF_CLIMITS(vol,vol_dat,vlo,vhi);

#if (BL_SPACEDIM == 2)
      FORT_MACDIV(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
		  ux_dat,ARLIM(uxlo),ARLIM(uxhi),
		  uy_dat,ARLIM(uylo),ARLIM(uyhi),
		  rx_dat,ARLIM(rxlo),ARLIM(rxhi),
		  ry_dat,ARLIM(rylo),ARLIM(ryhi),
		  ax_dat,ARLIM(axlo),ARLIM(axhi), 
		  ay_dat,ARLIM(aylo),ARLIM(ayhi), 
		  vol_dat,ARLIM(vlo),ARLIM(vhi));
#else
      const FArrayBox& uzedge = U_edge[2][U_edge0mfi];
      DEF_CLIMITS(uzedge,uz_dat,uzlo,uzhi);
      const FArrayBox& rhoDZ = RhoD[2][U_edge0mfi];
      DEF_CLIMITS(rhoDZ,rz_dat,rzlo,rzhi);
      const FArrayBox& zarea = area[level][2][U_edge0mfi];
      DEF_CLIMITS(zarea,az_dat,azlo,azhi);

      FORT_MACDIV(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
		  ux_dat,ARLIM(uxlo),ARLIM(uxhi),
		  uy_dat,ARLIM(uylo),ARLIM(uyhi),
		  uz_dat,ARLIM(uzlo),ARLIM(uzhi),
		  rx_dat,ARLIM(rxlo),ARLIM(rxhi),
		  ry_dat,ARLIM(rylo),ARLIM(ryhi),
		  rz_dat,ARLIM(rzlo),ARLIM(rzhi),
		  ax_dat,ARLIM(axlo),ARLIM(axhi),
		  ay_dat,ARLIM(aylo),ARLIM(ayhi),
		  az_dat,ARLIM(azlo),ARLIM(azhi),
		  vol_dat,ARLIM(vlo),ARLIM(vhi));
#endif

      sum += dmac.sum(0);
      npts += U_edge0mfi.validbox().d_numPts();
    }
  ParallelDescriptor::ReduceRealSum(sum);
  ParallelDescriptor::ReduceRealSum(npts);
  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "integral of divergence = " << sum/npts << std::endl;
}

//
// Check discrete divergence of each cell
//
void
MacProj::check_div_cond (int       level,
			 MultiFab& divu,
                         MultiFab  uedge[],
			 MultiFab  RhoD[]) 
{
  for (MFIter mfi(uedge[0]); mfi.isValid(); ++mfi)
    {
      FArrayBox& dmac = divu[mfi];

      const FArrayBox& uxedge = uedge[0][mfi];
      const FArrayBox& uyedge = uedge[1][mfi];
      const FArrayBox& rhoDX  = RhoD[0][mfi];
      const FArrayBox& rhoDY  = RhoD[1][mfi];
      const FArrayBox& xarea  = area[level][0][mfi];
      const FArrayBox& yarea  = area[level][1][mfi];
      const FArrayBox& vol    = volume[level][mfi];

      DEF_LIMITS(dmac,dmac_dat,dlo,dhi);
      DEF_CLIMITS(uxedge,ux_dat,uxlo,uxhi);
      DEF_CLIMITS(uyedge,uy_dat,uylo,uyhi);
      DEF_CLIMITS(rhoDX,rx_dat,rxlo,rxhi);
      DEF_CLIMITS(rhoDY,ry_dat,rylo,ryhi);
      DEF_CLIMITS(xarea,ax_dat,axlo,axhi);
      DEF_CLIMITS(yarea,ay_dat,aylo,ayhi);
      DEF_CLIMITS(vol,vol_dat,vlo,vhi);

#if (BL_SPACEDIM == 2)
      FORT_MACDIV(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
		  ux_dat,ARLIM(uxlo),ARLIM(uxhi),
		  uy_dat,ARLIM(uylo),ARLIM(uyhi),
		  rx_dat,ARLIM(rxlo),ARLIM(rxhi),
		  ry_dat,ARLIM(rylo),ARLIM(ryhi),
		  ax_dat,ARLIM(axlo),ARLIM(axhi), 
		  ay_dat,ARLIM(aylo),ARLIM(ayhi), 
		  vol_dat,ARLIM(vlo),ARLIM(vhi));
#else
      const FArrayBox& uzedge = uedge[2][mfi];
      DEF_CLIMITS(uzedge,uz_dat,uzlo,uzhi);
      const FArrayBox& rhoDZ = RhoD[2][mfi];
      DEF_CLIMITS(rhoDZ,rz_dat,rzlo,rzhi);
      const FArrayBox& zarea = area[level][2][mfi];
      DEF_CLIMITS(zarea,az_dat,azlo,azhi);

      FORT_MACDIV(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
		  ux_dat,ARLIM(uxlo),ARLIM(uxhi),
		  uy_dat,ARLIM(uylo),ARLIM(uyhi),
		  uz_dat,ARLIM(uzlo),ARLIM(uzhi),
		  rx_dat,ARLIM(rxlo),ARLIM(rxhi),
		  ry_dat,ARLIM(rylo),ARLIM(ryhi),
		  rz_dat,ARLIM(rzlo),ARLIM(rzhi),
		  ax_dat,ARLIM(axlo),ARLIM(axhi),
		  ay_dat,ARLIM(aylo),ARLIM(ayhi),
		  az_dat,ARLIM(azlo),ARLIM(azhi),
		  vol_dat,ARLIM(vlo),ARLIM(vhi));
#endif
    }
}

void
MacProj::set_dirichlet_bcs (int             level,
			    MultiFab*       mac_phi,
			    const MultiFab& RhoG,
			    const BCRec&    p_bc,
			    const Array<Real>&  press_lo,
			    const Array<Real>&  press_hi)
{

  //
  // Set dirichlet boundary condition with hydrostatic condition
  //

  bool hasDirichlet;
  Orientation Faces[2*BL_SPACEDIM];
  int numDirichlet;

  GetDirichletFaces(hasDirichlet,Faces,p_bc,numDirichlet);

  const BoxArray&   grids  = LevelData[level].boxArray();
  const Geometry&   geom   = parent->Geom(level);
  const Box&        domain = parent->Geom(level).Domain();

  //
  // Create 1-wide cc box just outside boundary to hold phi.
  // Create 2-wide cc box just inside  boundary to hold rho
  //

  BoxList  ccBoxList;
  BoxList phiBoxList;

  Array<IntVect> idx(numDirichlet);

  PorousMedia* pm_level = dynamic_cast<PorousMedia*>(&parent->getLevel(level));
  Real gravity = pm_level->getGravity();

  for (int iface = 0; iface < numDirichlet; iface++)
    {
      idx[iface] = IntVect::TheZeroVector();

      if (grids_on_side_of_domain(grids,geom.Domain(),Faces[iface])) 
	{
	  const int outDir    = Faces[iface].coordDir();

	  Box ccBndBox;
	  if (Faces[iface].faceDir() == Orientation::high)
	    {
	      ccBndBox = Box(BoxLib::adjCellHi(domain,outDir,2));
	      ccBndBox.shift(outDir,-2);
	      idx[iface].setVal(outDir,0);
	    } 
	  else 
	    {
	      ccBndBox = Box(BoxLib::adjCellLo(domain,outDir,2));
	      ccBndBox.shift(outDir,2);
	      idx[iface].setVal(outDir,0);
	    }
	  ccBndBox.grow(1);
	  ccBoxList.push_back(ccBndBox);

	  Box phiBox  = BoxLib::adjCell(domain,Faces[iface],1);
	  phiBoxList.push_back(phiBox);

	  const Box      valid_ccBndBox       = ccBndBox & domain;
	  const BoxArray uncovered_ba = BoxLib::complementIn(valid_ccBndBox,grids);
	  if (uncovered_ba.size() && 
	      BoxLib::intersect(grids,valid_ccBndBox).size() )
	    BoxLib::Error("MacProj: Cannot yet handle partially refined outflow");
	}
    }

  if ( !ccBoxList.isEmpty() ) 
    {
      BoxArray  ccBoxArray( ccBoxList);
      BoxArray phiBoxArray(phiBoxList);
 
      FArrayBox rhodat[2*BL_SPACEDIM];
      FArrayBox phidat[2*BL_SPACEDIM];

      BoxArray batmp = grids;
      batmp.grow(1);
      MultiFab rhotmp(batmp,1,0);
      for (MFIter mfi(rhotmp); mfi.isValid(); ++mfi)
	rhotmp[mfi].copy(RhoG[mfi]);

      for ( int iface = 0; iface < numDirichlet; ++iface) 
	{
	  rhodat[iface].resize(ccBoxArray[iface], 1);
	  phidat[iface].resize(phiBoxArray[iface], 1);
	
	  phidat[iface].setVal(0.0);
	  rhodat[iface].setVal(0.0);

	  rhotmp.shift(idx[iface]);
	  rhotmp.copy(rhodat[iface], 0, 0, 1);
	  rhotmp.shift(-idx[iface]);
	
	}
      
      const int* lo_bc = phys_bc->lo();
      const int* hi_bc = phys_bc->hi();

      computeRhoG (rhodat, phidat, geom, Faces, numDirichlet, gravity,
		   lo_bc, hi_bc, press_lo.dataPtr(), press_hi.dataPtr());

      // Must do this kind of copy instead of mac_phi->copy(phidat);
      //   because we're copying onto the ghost cells of the FABs,
      //   not the valid regions.
      for ( int iface =0; iface < numDirichlet; ++iface )
	{
	  for (MFIter mfi(*mac_phi); mfi.isValid(); ++mfi)
	    {
	      Box ovlp = (*mac_phi)[mfi].box() & phidat[iface].box();
	      if (ovlp.ok()) 
		(*mac_phi)[mfi].copy(phidat[iface],ovlp,0,ovlp,0,1);
	    }
	} 
    }
}

void
MacProj::computeRhoG (FArrayBox*         rhoG,
                      FArrayBox*         phiMF,
                      const Geometry&    geom,
                      Orientation*       outFaces,
                      int                numOutFlowFaces,
                      Real               gravity,
                      const int*         lo_bc,
                      const int*         hi_bc,
		      const Real*        press_in,
		      const Real*        press_out)

{
  const Real* dx    = geom.CellSize();
  const Box& domain = geom.Domain();
  const int* domlo  = domain.loVect();
  const int* domhi  = domain.hiVect();

  for (int iface = 0; iface < numOutFlowFaces; iface++) 
    {

      int face          = int(outFaces[iface]);

      DEF_LIMITS(phiMF[iface], phiPtr,philo,phihi);
      DEF_LIMITS(rhoG[iface], rhoPtr,rholo,rhohi);
      
      FORT_RHOGBC(rhoPtr,ARLIM(rholo),ARLIM(rhohi),
		  phiPtr,ARLIM(philo),ARLIM(phihi),
		  &face,&gravity,dx,domlo,domhi,
		  lo_bc,hi_bc,press_in,press_out);
      
    }
}

void
MacProj::GetOutFlowFaces (bool&        haveOutFlow,
                          Orientation* outFaces,
                          BCRec*       _phys_bc,
			  const BCRec& p_bc,
                          int&         numOutFlowBC)
{
  haveOutFlow = false;

  numOutFlowBC = 0;

  for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
      if ((_phys_bc->lo(idir) == Outflow || _phys_bc->lo(idir) == SEEPAGE)
	  &&  p_bc.lo(idir) == EXT_DIR)
        {
	  haveOutFlow = true;
	  outFaces[numOutFlowBC] = Orientation(idir,Orientation::low);
	  numOutFlowBC++;
        }

      if ((_phys_bc->hi(idir) == Outflow || _phys_bc->hi(idir) == SEEPAGE)
	  &&  p_bc.hi(idir) == EXT_DIR)
        {
	  haveOutFlow = true;
	  outFaces[numOutFlowBC] = Orientation(idir,Orientation::high);
	  numOutFlowBC++;
        }

    }
}

void
MacProj::GetInFlowFaces  (bool&        haveInFlow,
                          Orientation* inFaces,
                          BCRec*       _phys_bc,
			  const BCRec& p_bc,
                          int&         numInFlowBC)
{
  haveInFlow = false;

  numInFlowBC = 0;

  for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
      if (_phys_bc->lo(idir) == Inflow &&  p_bc.lo(idir) == EXT_DIR)
        {
	  haveInFlow = true;
	  inFaces[numInFlowBC] = Orientation(idir,Orientation::low);
	  numInFlowBC++;
        }

      if (_phys_bc->hi(idir) == Inflow &&  p_bc.hi(idir) == EXT_DIR)
        {
	  haveInFlow = true;
	  inFaces[numInFlowBC] = Orientation(idir,Orientation::high);
	  numInFlowBC++;
        }

    }
}

void
MacProj::GetDirichletFaces  (bool&        haveDirichlet,
			     Orientation* Faces,
			     const BCRec& p_bc,
			     int&         numDirichlet)
{
  // 
  //  find the dirichlet faces.
  //

  haveDirichlet = false;

  numDirichlet = 0;

  for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
      if (p_bc.lo(idir) == EXT_DIR)
        {
	  haveDirichlet = true;
	  Faces[numDirichlet] = Orientation(idir,Orientation::low);
	  numDirichlet++;
        }

      if (p_bc.hi(idir) == EXT_DIR)
        {
	  haveDirichlet = true;
	  Faces[numDirichlet] = Orientation(idir,Orientation::high);
	  numDirichlet++;
        }
    }
}

void
MacProj::print_min_max (MultiFab* u_mac)
{

  // Write out the min and max of the MAC velocities
  Real umax =  -1.e20;
  Real vmax =  -1.e20;
  Real umin =   1.e20;
  Real vmin =   1.e20;
#if(BL_SPACEDIM == 3)
  Real wmax =  -1.e20;
  Real wmin =   1.e20;
#endif
  if (verbose) {
    for (MFIter mfi(u_mac[0]); mfi.isValid(); ++mfi)
      {
	int i = mfi.index();
	umax = std::max(umax,u_mac[0][i].max(u_mac[0].boxArray()[i]));
	umin = std::min(umin,u_mac[0][i].min(u_mac[0].boxArray()[i]));
	vmax = std::max(vmax,u_mac[1][i].max(u_mac[1].boxArray()[i]));
	vmin = std::min(vmin,u_mac[1][i].min(u_mac[1].boxArray()[i]));
#if(BL_SPACEDIM == 3)
	wmax = std::max(wmax,u_mac[2][i].max(u_mac[2].boxArray()[i]));
	wmin = std::min(wmin,u_mac[2][i].min(u_mac[2].boxArray()[i]));
#endif
      }
    ParallelDescriptor::ReduceRealMax(umax);
    ParallelDescriptor::ReduceRealMin(umin);
    ParallelDescriptor::ReduceRealMax(vmax);
    ParallelDescriptor::ReduceRealMin(vmin);
#if(BL_SPACEDIM == 3)
    ParallelDescriptor::ReduceRealMax(wmax);
    ParallelDescriptor::ReduceRealMin(wmin);
#endif
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "   UMAC MAX/MIN  " << umax << "  " << umin << std::endl;
      std::cout << "   VMAC MAX/MIN  " << vmax << "  " << vmin << std::endl;
#if(BL_SPACEDIM == 3)
      std::cout << "   WMAC MAX/MIN  " << wmax << "  " << wmin << std::endl;
#endif
    }
  }
}

