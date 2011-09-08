//
// Comment out this line to use diffusion class outside
// the context of PorousMedia and classes derived from it.
//
#include <winstd.H>

#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <ErrorList.H>

#include <PorousMedia.H>
#include <Diffusion.H>
#include <MultiGrid.H>
#include <CGSolver.H>

#include <DIFFUSION_F.H>
#include <VISCOPERATOR_F.H>

#include <algorithm>
#include <cfloat>

#if defined(BL_OSF1)
#if defined(BL_USE_DOUBLE)
const Real BL_BOGUS      = DBL_QNAN;
#else
const Real BL_BOGUS      = FLT_QNAN;
#endif
#else
const Real BL_BOGUS      = 1.e200;
#endif

const Real BL_SAFE_BOGUS = -666.e200;

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)	\
  const int* fablo = (fab).loVect();		\
  const int* fabhi = (fab).hiVect();		\
  Real* fabdat     = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)	\
  const int* fablo   = (fab).loVect();		\
  const int* fabhi   = (fab).hiVect();		\
  const Real* fabdat = (fab).dataPtr();

#define GEOM_GROW 1
#define HYP_GROW 3

namespace
{
    bool initialized = false;
}
//
// Set defaults in Initialize() !!!
//
int  Diffusion::verbose;
Real Diffusion::visc_tol;
int  Diffusion::do_reflux;
int  Diffusion::max_order;
int  Diffusion::scale_abec;
Real Diffusion::visc_abs_tol;
int  Diffusion::use_cg_solve;
bool Diffusion::use_mg_precond_flag;

#ifdef MG_USE_FBOXLIB
namespace
{
    bool use_fboxlib_mg;
}
#endif

Array<Real> Diffusion::visc_coef;
Array<int>  Diffusion::is_diffusive;

void
Diffusion::Finalize ()
{
    visc_coef.clear();
    is_diffusive.clear();

    initialized = false;
}

void
Diffusion::Initialize ()
{
    if (initialized) return;

    Diffusion::verbose             = 0;
    Diffusion::visc_tol            = 1.0e-16;      
    Diffusion::do_reflux           = 1;
    Diffusion::max_order           = 4;
    Diffusion::scale_abec          = 0;
    Diffusion::visc_abs_tol        = 1.0e-16;
    Diffusion::use_cg_solve        = 0;
    Diffusion::use_mg_precond_flag = false;

#ifdef MG_USE_FBOXLIB
    use_fboxlib_mg = false;
#endif

    BoxLib::ExecOnFinalize(Diffusion::Finalize);

    initialized = true;
}

static
Real
dotxy (const MultiFab& r,
       const MultiFab& z,
       int             idxr = 0,
       int             idxz = 0,
       bool            local = false)
{
  BL_PROFILE("Diffusion::dotxy");
  BL_ASSERT(idxr < r.nComp() && idxz < z.nComp());
  int  nz  = z.nComp();
  int  nr  = r.nComp();
  Real rho = 0.0;
  for (MFIter mfi(r); mfi.isValid(); ++mfi)
    {
      const int k = mfi.index();
      Real trho;
      int  idxr_f = idxr + 1;
      int  idxz_f = idxz + 1;
      
      FORT_DFXDOTY(&trho,
		   z[k].dataPtr(),
		   ARLIM(z[k].loVect()),ARLIM(z[k].hiVect()),
		   &nz, &idxz_f,
		   r[k].dataPtr(),
		   ARLIM(r[k].loVect()),ARLIM(r[k].hiVect()),
		   &nr, &idxr_f,
		   r.box(k).loVect(),r.box(k).hiVect());
        rho += trho;
    }

    if (!local)
        ParallelDescriptor::ReduceRealSum(rho);

    return rho;
}


Diffusion::Diffusion (Amr*               Parent,
                      AmrLevel*          Caller,
                      Diffusion*         Coarser,
                      int                num_state,
                      FluxRegister*      Viscflux_reg,
                      MultiFab&          Volume,
                      MultiFab*          Area,
                      const Array<int>&  _is_diffusive,
                      const Array<Real>& _visc_coef)
  :
  parent(Parent),
  caller(Caller),
  grids(caller->boxArray()),
  level(caller->Level()),
  volume(Volume),
  area(Area),
  coarser(Coarser),
  finer(0),
  NUM_SCALARS(num_state),
  viscflux_reg(Viscflux_reg)
{
  if (!initialized)
    {
      Initialize();

      ParmParse ppdiff("diffuse");

      ppdiff.query("v",            verbose);
      ppdiff.query("use_cg_solve", use_cg_solve);

#ifdef MG_USE_FBOXLIB
      ppdiff.query("use_fboxlib_mg", use_fboxlib_mg);
      if ( use_cg_solve && use_fboxlib_mg )
	{
	  BoxLib::Error("Diffusion::read_params: cg_solve && .not. fboxlib_solve");
	}
#endif
      int use_mg_precond = 0;

      ppdiff.query("max_order",      max_order);
      ppdiff.query("scale_abec",     scale_abec);
      ppdiff.query("use_mg_precond", use_mg_precond);

      use_mg_precond_flag = (use_mg_precond ? true : false);

      ParmParse pp("prob");

      pp.query("visc_tol",     visc_tol);
      pp.query("do_reflux",    do_reflux);
      pp.query("visc_abs_tol", visc_abs_tol);

      do_reflux = (do_reflux ? 1 : 0);

      const int n_visc = _visc_coef.size();
      const int n_diff = _is_diffusive.size();
      //
      // Check whether number of diffusion coefficients is sufficient.  
      // Coefficients defined for i > NUM_SCALARS are ignored.
      //
      if (n_diff < NUM_SCALARS){
	BoxLib::Abort("Diffusion::Diffusion(): is_diffusive array is not long enough");
      }

      if (n_visc < NUM_SCALARS){
	std::cout << n_visc << std::endl;
	std::cout << NUM_SCALARS << std::endl ;
	BoxLib::Abort("Diffusion::Diffusion(): visc_coef array is not long enough");
      }
      visc_coef.resize(NUM_SCALARS);
      is_diffusive.resize(NUM_SCALARS);
      for (int i = 0; i < NUM_SCALARS; i++)
        {
	  is_diffusive[i] = _is_diffusive[i];
	  visc_coef[i] = _visc_coef[i];
        }

      echo_settings();
    }

  if (level > 0)
    {
      crse_ratio = parent->refRatio(level-1);
      coarser->finer = this;
    }

  rho = new MultiFab(grids,1,1);  
  rho->setVal(1.);

}

Diffusion::Diffusion (Amr*               Parent,
                      AmrLevel*          Caller,
                      Diffusion*         Coarser,
                      int                num_state,
                      FluxRegister*      Viscflux_reg,
                      MultiFab&          Volume,
                      MultiFab*          Area)
  :
  parent(Parent),
  caller(Caller),
  grids(caller->boxArray()),
  level(caller->Level()),
  volume(Volume),
  area(Area),
  coarser(Coarser),
  finer(0),
  NUM_SCALARS(num_state),
  viscflux_reg(Viscflux_reg)
{
  if (!initialized)
    {
      Initialize();

      ParmParse ppdiff("diffuse");

      ppdiff.query("v",            verbose);
      ppdiff.query("use_cg_solve", use_cg_solve);

#ifdef MG_USE_FBOXLIB
      ppdiff.query("use_fboxlib_mg", use_fboxlib_mg);
      if ( use_cg_solve && use_fboxlib_mg )
	{
	  BoxLib::Error("Diffusion::read_params: cg_solve && .not. fboxlib_solve");
	}
#endif
      int use_mg_precond = 0;

      ppdiff.query("max_order",      max_order);
      ppdiff.query("scale_abec",     scale_abec);
      ppdiff.query("use_mg_precond", use_mg_precond);

      use_mg_precond_flag = (use_mg_precond ? true : false);

      ParmParse pp("ns");

      pp.query("visc_tol",  visc_tol);
      pp.query("do_reflux", do_reflux);

      do_reflux = (do_reflux ? 1 : 0);

      echo_settings();
    }

  if (level > 0)
    {
      crse_ratio = parent->refRatio(level-1);
      coarser->finer = this;
    }

  rho = new MultiFab(grids,1,1);   
  (*rho).setVal(1.0);
}

Diffusion::~Diffusion () 
{
  delete rho;
}

FluxRegister*
Diffusion::viscFluxReg ()
{
  return viscflux_reg;
}

int
Diffusion::maxOrder() const
{
  return max_order;
}

void
Diffusion::echo_settings () const
{
  //
  // Print out my settings.
  //
  if (verbose && ParallelDescriptor::IOProcessor())
    {
      std::cout << "Diffusion settings...\n";
      std::cout << "  From diffuse:\n";
      std::cout << "   use_cg_solve =        " << use_cg_solve << '\n';
      std::cout << "   use_mg_precond_flag = " << use_mg_precond_flag << '\n';
      std::cout << "   max_order =           " << max_order << '\n';
      std::cout << "   scale_abec =          " << scale_abec << '\n';

      std::cout << "\n\n  From ns:\n";
      std::cout << "   do_reflux =           " << do_reflux << '\n';
      std::cout << "   visc_tol =            " << visc_tol << '\n';
    
      std::cout << "   is_diffusive =";
      for (int i =0; i < NUM_SCALARS; i++)
	std::cout << "  " << is_diffusive[i];
    
      std::cout << "\n   visc_coef =";
      for (int i = 0; i < NUM_SCALARS; i++)
	std::cout << "  " << visc_coef[i];

      std::cout << '\n';
    }
}

Real
Diffusion::get_scaled_abs_tol (const MultiFab& rhs,
                               Real            reduction) const
{
  Real norm_est = 0;
  for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    norm_est = std::max(norm_est, rhs[mfi].norm(0));
  ParallelDescriptor::ReduceRealMax(norm_est);
  return norm_est * reduction;
}

void
Diffusion::set_rho(const MultiFab* rho_in)
{
  MultiFab::Copy(*rho,*rho_in,0,0,1,1);
}

void
Diffusion::set_rho(Real rho_in)
{
  rho->setVal(rho_in);
}


void
Diffusion::diffuse_scalar (Real                   dt,
			   int                    sigma,
			   Real                   be_cn_theta,
			   MultiFab* const*       fluxn,
			   MultiFab* const*       fluxnp1,
			   int                    dataComp,
			   MultiFab*              delta_rhs, 
			   const MultiFab*        alpha, 
			   const MultiFab* const* betan, 
			   const MultiFab* const* betanp1,
			   const SolveMode&       solve_mode)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::diffuse_scalar()");
  //
  // This routine expects that physical BC's have been loaded into
  // the grow cells of the old and new state at this level.  If rho_flag==2,
  // the values there are rho.phi, where phi is the quantity being diffused.
  // Values in these cells will be preserved.  Also, if there are any
  // explicit update terms, these have already incremented the new state
  // on the valid region (i.e., on the valid region the new state is the old
  // state + dt*Div(explicit_fluxes), e.g.)
  //
    
  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... diffusing scalar " 
	      << caller->get_desc_lst()[State_Type].name(sigma) << '\n';

  int allnull, allthere;
  checkBetas(betan, betanp1, allthere, allnull);

  BL_ASSERT(solve_mode==ONEPASS || (delta_rhs && delta_rhs->boxArray()==grids));

  //
  // At this point, S_old has bndry at time N, S_new has bndry at time N+1
  //
  MultiFab& S_old = caller->get_old_data(State_Type);
  MultiFab& S_new = caller->get_new_data(State_Type);
    
  MultiFab Rhs(grids,1,0);
  MultiFab Soln(grids,1,1);
  //
  // Set up Rhs.
  //
  Real a = 0.0;
  Real b = -(1.0-be_cn_theta)*dt;
  if (allnull)
    b *= visc_coef[sigma];
  ViscBndry visc_bndry_0;
  const Real prev_time   = caller->get_state_data(State_Type).prevTime();
  ABecLaplacian* visc_op = getViscOp(sigma,a,b,prev_time,visc_bndry_0,
				     0,dataComp,betan);
  visc_op->maxOrder(max_order);
  //
  // Copy to single-component multifab, then apply op to rho-scaled state
  //
  MultiFab::Copy(Soln,S_old,sigma,0,1,1);
  Soln.divide(*rho,0,1,1);
    
  visc_op->apply(Rhs,Soln);


  visc_op->compFlux(D_DECL(*fluxn[0],*fluxn[1],*fluxn[2]),Soln);
  for (int i = 0; i < BL_SPACEDIM; ++i)
    (*fluxn[i]).mult(-b/(dt*caller->Geom().CellSize()[i]));
  delete visc_op;

  //
  // If this is a predictor step, put "explicit" updates passed via S_new
  // into delta_rhs, so they dont get lost,
  // pull it off S_new to avoid double counting
  //   (for rho_flag == 1:
  //       S_new = S_old - dt.(U.Grad(phi)); want Rhs -= (U.Grad(phi)),
  //    else
  //       S_new = S_old - dt.Div(U.Phi),   want Rhs -= Div(U.Phi) )
  //

  FArrayBox tmpfab;
    
  if (solve_mode == PREDICTOR)
    {
      for (MFIter Smfi(S_new); Smfi.isValid(); ++Smfi)
        {
	  const Box& box = Smfi.validbox();
	  tmpfab.resize(box,1);
	  tmpfab.copy(S_new[Smfi],box,sigma,box,0,1);
	  tmpfab.minus(S_old[Smfi],box,sigma,0,1);
	  S_new[Smfi].minus(tmpfab,box,0,sigma,1); // Remove this term from S_new
	  tmpfab.mult(1.0/dt,box,0,1);
	  if (alpha!=0)
	    tmpfab.mult((*alpha)[Smfi],box,dataComp,0,1);
	  (*delta_rhs)[Smfi].plus(tmpfab,box,0,dataComp,1);
        }
    }

  //
  // Add body sources
  //
  if (delta_rhs != 0)
    {
      for (MFIter mfi(*delta_rhs); mfi.isValid(); ++mfi)
        {
	  const Box& box = mfi.validbox();
	  tmpfab.resize(box,1);
	  tmpfab.copy((*delta_rhs)[mfi],box,dataComp,box,0,1);
	  tmpfab.mult(dt,box,0,1);
	  tmpfab.mult(volume[mfi],box,0,0,1);
	  Rhs[mfi].plus(tmpfab,box,0,0,1);
        }
    }
   
  MultiFab::Copy(Soln,S_new,sigma,0,1,1);
  //Soln.divide(*rho,0,1,1);
  for (MFIter mfi(Soln); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.validbox();
      Soln[mfi].mult(volume[mfi],box,0,0,1);
      if (alpha!=0)
	Soln[mfi].mult((*alpha)[mfi],box,dataComp,0,1);
      Rhs[mfi].plus(Soln[mfi],box,0,0,1);
    }
  //
  // Make a good guess for Soln
  //
  MultiFab::Copy(Soln,S_new,sigma,0,1,0);
  Soln.divide(*rho,0,1,1);

  //
  // Construct viscous operator with bndry data at time N+1.
  //
  a = 1.0;
  b = be_cn_theta*dt;
  if (allnull) {
    b *= visc_coef[sigma];
  }
  ViscBndry  visc_bndry;
  const Real cur_time = caller->get_state_data(State_Type).curTime();
  Real       rhsscale = 1.0;
  visc_op  = getViscOp(sigma,a,b,cur_time,visc_bndry,
		       &rhsscale,dataComp,betanp1,alpha);
  Rhs.mult(rhsscale,0,1);
  visc_op->maxOrder(max_order);
  //
  // Construct solver and call it.
  //
  const Real S_tol     = visc_tol;
  const Real S_tol_abs = get_scaled_abs_tol(Rhs, visc_tol);
  if (use_cg_solve)
    {
      CGSolver cg(*visc_op,use_mg_precond_flag);
      cg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }
#ifdef MG_USE_FBOXLIB
  else if ( use_fboxlib_mg )
    {
      std::vector<BoxArray> bav(1);
      bav[0] = S_new.boxArray();
      std::vector<DistributionMapping> dmv(1);
      dmv[0] = Rhs.DistributionMap();
      bool nodal = false;
      std::vector<Geometry> geom(1);
      geom[0] = visc_bndry.getGeom();
      const BCRec& scal_bc = caller->get_desc_lst()[State_Type].getBC(sigma);
 
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
	      mg_bc[i*2 + 0] = scal_bc.lo(i)==EXT_DIR? MGT_BC_DIR : MGT_BC_NEU;
	      mg_bc[i*2 + 1] = scal_bc.hi(i)==EXT_DIR? MGT_BC_DIR : MGT_BC_NEU;
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
	const Real* dx_crse   = parent->Geom(level-1).CellSize();
	for ( int i = 0; i < BL_SPACEDIM; ++i ) {
	  xa[0][i] = 0.5 * dx_crse[i];
	  xb[0][i] = 0.5 * dx_crse[i];
	}
      }

      MGT_Solver mgt_solver(geom, mg_bc, bav, dmv, nodal);
      mgt_solver.set_maxorder(max_order);

      const MultiFab* aa_p[1];
      aa_p[0] = &(visc_op->aCoefficients());
      const MultiFab* bb_p[1][BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i )
        {
	  bb_p[0][i] = &(visc_op->bCoefficients(i));
        }
      mgt_solver.set_visc_coefficients(aa_p, bb_p, b, xa, xb);

      MultiFab* phi_p[1];
      MultiFab* Rhs_p[1];
      phi_p[0] = &Soln;
      Rhs_p[0] = &Rhs;
      Real final_resnorm;
      mgt_solver.solve(phi_p, Rhs_p, S_tol, S_tol_abs, visc_bndry,final_resnorm);
    }
#endif
  else
    {
      MultiGrid mg(*visc_op);
      mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }
  //
  // Get extensivefluxes from new-time op
  //
  visc_op->compFlux(D_DECL(*fluxnp1[0],*fluxnp1[1],*fluxnp1[2]),Soln);
  for (int i = 0; i < BL_SPACEDIM; ++i)
    (*fluxnp1[i]).mult(b/(dt*caller->Geom().CellSize()[i]));
  delete visc_op;
  //
  // Copy into state variable at new time, without bc's
  //
  MultiFab::Multiply(Soln,*rho,0,0,1,1);
  MultiFab::Copy(S_new,Soln,0,sigma,1,0);
}


void Diffusion::diffuse_init_CPL(Real                   dt,
				 int                    nc,
				 Real                   be_cn_theta,
				 MultiFab* const*       fluxn,
				 int                    dataComp,
				 MultiFab*              delta_rhs, 
				 const MultiFab*        alpha, 
				 const MultiFab* const* betan, 
				 const MultiFab*        pcn,
				 MultiFab*              Rhs)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::diffuse_init_CPL()");

  //
  // Set up initial RHS of the Newton-solve that depends only on
  // solutions at time n and not n+1.
  // RHS = s^{n+1,0} - \nabla \rho D \nabla pc^n + delta_rhs
  //

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "construct initial RHS for capillary solve : " 
	      << caller->get_desc_lst()[State_Type].name(nc) << '\n';
  
  MultiFab Tmpn(grids,1,0);
  MultiFab Soln(grids,1,1);

  Tmpn.setVal(0.);
  Soln.setVal(0.);

  // Rhs is initially s^{n+1,0}
  // Multiply by phi and volume to get phi*s^{n} + dt*(div F) 
  for (MFIter mfi(*Rhs); mfi.isValid(); ++mfi)   
    {
      const Box& box = mfi.validbox();
      (*Rhs)[mfi].mult(volume[mfi],box,0,0,1);
      if (alpha!=0)
	(*Rhs)[mfi].mult((*alpha)[mfi],box,0,0,1);
    }
  
  // Determine - (dt/2) \nabla H^n \nabla P_c^n
  Real a = 0.0;
  Real b = (1.0-be_cn_theta)*dt;
  ViscBndry visc_bndry_0;
  const BCRec& bc      = caller->get_desc_lst()[State_Type].getBC(nc);
  const Real prev_time = caller->get_state_data(State_Type).prevTime();
  IntVect ref_ratio    = level > 0 ? 
    parent->refRatio(level-1) : IntVect::TheUnitVector();
  getBndryData(visc_bndry_0,nc,1,prev_time);
  visc_bndry_0.setScalarValues(bc,ref_ratio,pcn);
  ABecLaplacian* visc_op = getViscOp(nc,a,b,prev_time,visc_bndry_0,
				     0,dataComp,betan,alpha,true);
  MultiFab::Copy(Soln,*pcn,0,0,1,1);  
  visc_op->apply(Tmpn,Soln);
  visc_op->compFlux(D_DECL(*fluxn[0],*fluxn[1],*fluxn[2]),Soln);
  for (int i = 0; i < BL_SPACEDIM; ++i)
    (*fluxn[i]).mult(b/dt/caller->Geom().CellSize()[i]);
  delete visc_op;

  //RHS = RHS - (dt/2)*(\nabla H^n \nabla P_c^n) 
  Rhs->plus(Tmpn,0,1,0);

  // RHS = RHS + delta_rhs
  if (delta_rhs != 0)
    {
      MultiFab::Copy(Soln,*delta_rhs,0,0,1,1);
      for (MFIter mfi(Soln); mfi.isValid(); ++mfi)
	{
	  const Box& box = mfi.validbox();
	  Soln[mfi].mult(volume[mfi],box,0,0,1);
	  (*Rhs)[mfi].plus(Soln[mfi],box,0,0,1);
	}
    }
}

void
Diffusion::diffuse_iter_CPL (Real                   dt,
			     int                    nc,
			     int                    ncomps,
			     Real                   be_cn_theta,
			     int                    dataComp,
			     const MultiFab*        alpha, 
			     const MultiFab* const* betanp1,
			     const MultiFab* const* betanp1_dp,
			     MultiFab*              pcnp1,
			     MultiFab*              S_nwt,
			     Real*                  err_nwt)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::diffuse_iter_CPL()");
  //
  // - This routine expects that physical BC's have been loaded into
  // the grow cells of the old and new state at this level.  
  // - If there are any explicit update terms, these are reflected in 
  // the new state (i.e., on the valid region the new state is the old
  // state + dt*Div(explicit_fluxes), e.g.)
  //
  MultiFab& S_new = caller->get_new_data(State_Type);
  // setup multifabs for solvers
  MultiFab Rhs(grids,1,0);
  MultiFab Soln(grids,1,1);

  Rhs.setVal(0.);
  Soln.setVal(0.);

  // Compute residual
  ViscBndry visc_bndry_0;
  const BCRec& bc     = caller->get_desc_lst()[State_Type].getBC(nc);
  const Real cur_time = caller->get_state_data(State_Type).curTime();
  IntVect ref_ratio   = level > 0 ? 
    parent->refRatio(level-1) : IntVect::TheUnitVector();
  getBndryData(visc_bndry_0,nc,1,cur_time);
  residual_CPL(visc_bndry_0,
	       nc,dt,be_cn_theta,Rhs,S_nwt,pcnp1,alpha,betanp1);  

  // Make a good guess for Soln
  Soln.setVal(0);
    
  // Construct viscous operator with bndry data at time N+1.
  int a = 1.0;
  int b = -be_cn_theta*dt;  
  Real rhsscale = 1.0;
  visc_bndry_0.setdeltaSValues(bc,ref_ratio);
  ABecLaplacian* visc_op_dp = getViscOp(nc,a,b,cur_time,visc_bndry_0,
       &rhsscale,dataComp,betanp1_dp,alpha,true);
  //
  // Construct solver and call it.
  //
  const Real S_tol     = visc_tol;
  const Real S_tol_abs = get_scaled_abs_tol(Rhs, visc_tol);
  if (use_cg_solve)
    {
      CGSolver cg(*visc_op_dp,use_mg_precond_flag);
      cg.solve(Soln,Rhs,S_tol,S_tol_abs);	
    }

#ifdef MG_USE_FBOXLIB
  else if ( use_fboxlib_mg )
    {
      std::vector<BoxArray> bav(1);
      bav[0] = S_new.boxArray();
      std::vector<DistributionMapping> dmv(1);
      dmv[0] = Rhs.DistributionMap();
      bool nodal = false;
      std::vector<Geometry> geom(1);
      geom[0] = visc_bndry_0.getGeom();
      const BCRec& scal_bc = caller->get_desc_lst()[State_Type].getBC(nc);
 
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
	      mg_bc[i*2 + 0] = scal_bc.lo(i)==EXT_DIR? MGT_BC_DIR : MGT_BC_NEU;
	      mg_bc[i*2 + 1] = scal_bc.hi(i)==EXT_DIR? MGT_BC_DIR : MGT_BC_NEU;
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
	const Real* dx_crse   = parent->Geom(level-1).CellSize();
	for ( int i = 0; i < BL_SPACEDIM; ++i ) {
	  xa[0][i] = 0.5 * dx_crse[i];
	  xb[0][i] = 0.5 * dx_crse[i];
	}
      }

      MGT_Solver mgt_solver(geom, mg_bc, bav, dmv, nodal);
      mgt_solver.set_maxorder(max_order);

      const MultiFab* aa_p[1];
      aa_p[0] = &(visc_op_dp->aCoefficients());
      const MultiFab* bb_p[1][BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i )
        {
	  bb_p[0][i] = &(visc_op_dp->bCoefficients(i));
        }
      mgt_solver.set_visc_coefficients(aa_p, bb_p, b, xa, xb);

      MultiFab* phi_p[1];
      MultiFab* Rhs_p[1];
      phi_p[0] = &Soln;
      Rhs_p[0] = &Rhs;
      Real final_resnorm;
      mgt_solver.solve(phi_p, Rhs_p, S_tol, S_tol_abs, 
		       visc_bndry_0,final_resnorm);
    }
#endif
  else
    {
      MultiGrid mg(*visc_op_dp);
      mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }

  delete visc_op_dp;

  // Multiply by rho to get \delta n
  MultiFab::Multiply(Soln,*rho,0,0,1,0);

  // Do line search
  PorousMedia* pm_level = dynamic_cast<PorousMedia*>(&parent->getLevel(level));
  MultiFab** betatmp  = 0;
  allocFluxBoxesLevel(betatmp,0,1);
  // Store old S_new 
  MultiFab Stmp(grids,1,1);
  MultiFab Rhsp1(grids,1,0);
  MultiFab::Copy(Stmp,S_new,nc,0,1,1);
  MultiFab::Add(S_new,Soln,0,nc,1,0);
  pm_level->scalar_adjust_constraint(0,ncomps-1);
  pm_level->FillStateBndry(cur_time,State_Type,0,ncomps);
  pm_level->calcCapillary(cur_time);
  pm_level->calcLambda(cur_time);
  pm_level->calcDiffusivity_CPL(betatmp,pm_level->lambdap1_cc);
  residual_CPL(visc_bndry_0,
	       nc,dt,be_cn_theta,Rhsp1,S_nwt,pm_level->pcnp1_cc,alpha,betatmp);

  Real rhs_norm = Rhs.norm2(0);
  Real rhsp1_norm = Rhsp1.norm2(0);
  Real alphak = 1.0;
  int iter = 0;
  while (iter < 10 && rhsp1_norm > rhs_norm && alphak > 1.e-3) 
    {
      
      alphak = 0.33*alphak;
      MultiFab::Copy(S_new,Stmp,0,nc,1,1);
      MultiFab::Add(S_new,Soln,0,nc,1,0);
      Soln.mult(0.33);
      pm_level->scalar_adjust_constraint(0,ncomps-1);
      pm_level->FillStateBndry(cur_time,State_Type,0,ncomps);
      pm_level->calcCapillary(cur_time);
      pm_level->calcLambda(cur_time);
      pm_level->calcDiffusivity_CPL(betatmp,pm_level->lambdap1_cc);
      residual_CPL(visc_bndry_0,
		   nc,dt,be_cn_theta,Rhsp1,S_nwt,
		   pm_level->pcnp1_cc,alpha,betatmp);
      rhsp1_norm = Rhsp1.norm2(0);
      iter += 1;
    }
    
  // Compute the err_nwt
  //*err_nwt = Soln.norm2(0)/S_new.norm2(0);
  *err_nwt = Rhs.norm2(0)/S_new.norm2(0);

  removeFluxBoxesLevel(betatmp);
}

void 
Diffusion::residual_CPL (ViscBndry              visc_bndry,
			 int                    nc,
			 Real                   dt,
			 Real                   be_cn_theta,	
			 MultiFab&              Rhs,
			 MultiFab*              S_nwt,
			 const MultiFab*        pcnp1,
			 const MultiFab*        alpha,
			 const MultiFab* const* betanp1,
			 MultiFab* const*       fluxnp1)
{ 
  MultiFab Soln(grids,1,1);
  Rhs.setVal(0.);
  Soln.setVal(0.);

  // Set up Rhs.
  Real a = 0.0;
  Real b = (1.0-be_cn_theta)*dt;

  // Rhs = - (dt/2) * \nabla H^{n+1} \nabla P_c^{n+1}
  const BCRec& bc     = caller->get_desc_lst()[State_Type].getBC(nc);
  const Real cur_time = caller->get_state_data(State_Type).curTime();
  IntVect ref_ratio   = level > 0 ? 
    parent->refRatio(level-1) : IntVect::TheUnitVector();
  visc_bndry.setScalarValues(bc,ref_ratio,pcnp1);
  ABecLaplacian* visc_op = getViscOp(nc,a,b,cur_time,visc_bndry,
		0,0,betanp1,alpha,true);
  MultiFab::Copy(Soln,*pcnp1,0,0,1,1);
  visc_op->apply(Rhs,Soln);

  if (fluxnp1 !=0)
    {
      visc_op->compFlux(D_DECL(*fluxnp1[0],*fluxnp1[1],*fluxnp1[2]),Soln);
      for (int i = 0; i < BL_SPACEDIM; ++i)
	(*fluxnp1[i]).mult(b/dt/caller->Geom().CellSize()[i]);
    }

  delete visc_op;

  // RHS = RHS + time n data - \phi s^{n+1,k}
  MultiFab& S_new = caller->get_new_data(State_Type);
  MultiFab::Copy(Soln,S_new,nc,0,1,1);
  for (MFIter mfi(Soln); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.validbox();
      Rhs[mfi].plus((*S_nwt)[mfi],box,0,0,1);
      Soln[mfi].mult(volume[mfi],box,0,0,1);
      if (alpha!=0)
	Soln[mfi].mult((*alpha)[mfi],box,0,0,1);
      Rhs[mfi].minus(Soln[mfi],box,0,0,1);
    }
}

#ifdef MG_USE_FBOXLIB
void 
Diffusion::residual_richard (ABecLaplacian*         visc_op,
			     Real                   dt,
			     Real                   gravity,
			     Array<Real>            density,
			     MultiFab&              res,
			     const MultiFab*        pc,
			     const MultiFab* const* beta,
			     const MultiFab*        alpha,
			     const MultiFab*        Stmp)
{ 
  MultiFab tmp(grids,1,0);
  MultiFab Soln(grids,1,1);
  tmp.setVal(0.);
  Soln.setVal(0.);

  // res =  res_fix - \nabla H \nabla P_c
  MultiFab::Copy(Soln,*pc,0,0,1,1);
  visc_op->apply(tmp,Soln);
  MultiFab::Add(res,tmp,0,0,1,0);
  // res = res + \nabla H \rho |g| 
  const Real* dx   = caller->Geom().CellSize();
  for (MFIter mfi(res); mfi.isValid(); ++mfi)
    {
      const int* lo   = mfi.validbox().loVect();
      const int* hi   = mfi.validbox().hiVect();
	
      FArrayBox& rg       = res[mfi];     
      const FArrayBox& bx = (*beta[0])[mfi];
      const FArrayBox& by = (*beta[1])[mfi];
	
      DEF_LIMITS (rg,rg_dat,rglo,rghi);
      DEF_CLIMITS(bx,bx_dat,bxlo,bxhi);
      DEF_CLIMITS(by,by_dat,bylo,byhi);
	
#if (BL_SPACEDIM==3)
	const FArrayBox& bz = (*beta[2])[mfi];
	DEF_CLIMITS(bz,bz_dat,bzlo,bzhi);
#endif
	FORT_DRHOG_RICHARD (rg_dat, ARLIM(rglo), ARLIM(rghi),
			    bx_dat, ARLIM(bxlo), ARLIM(bxhi),
			    by_dat, ARLIM(bylo), ARLIM(byhi),
#if (BL_SPACEDIM==3)
			    bz_dat, ARLIM(bzlo), ARLIM(bzhi),
#endif
			    &density[0],lo,hi,dx,&gravity,&dt);
    }

  // res  = res + \phi (n(t+dt))
  if (alpha)
    {
      if (Stmp)
	{
	  MultiFab::Copy(Soln,*Stmp,0,0,1,0);
	}
      else
	{
	  MultiFab& S_new = caller->get_new_data(State_Type);
	  MultiFab::Copy(Soln,S_new,0,0,1,0);
	}
      for (MFIter mfi(Soln); mfi.isValid(); ++mfi)
	{
	  const Box& box = mfi.validbox();
	  Soln[mfi].mult((*alpha)[mfi],box,0,0,1);
	  res[mfi].plus(Soln[mfi],box,0,0,1);
	}
    }
}

void 
Diffusion::residual_richard (MGT_Solver&               mgt_solver,
			     Real                      dt,
			     Real                      gravity,
			     Array<Real>               density,
			     MultiFab**                Rhs,
			     PArray<MultiFab>&         pc,
			     Array<PArray<MultiFab> >& beta,
			     PArray<MultiFab>&         alpha,
			     PArray<MultiFab>&         res_fix,
			     MultiFab**                S,
			     ViscBndry                 visc_bndry)
{ 
  int nlevs = alpha.size();
  MultiFab* soln[alpha.size()];
  
  for (int lev=0; lev<nlevs;lev++)
    soln[lev] = &pc[lev];

  mgt_solver.applyop(soln,Rhs,visc_bndry);

  for (int lev = 0; lev < nlevs; lev++) 
    { 
      MultiFab::Add(*Rhs[lev],res_fix[lev],0,0,1,0);
      // gravity term
      const Real* dx   = parent->Geom(lev).CellSize();
      for (MFIter mfi(*Rhs[lev]); mfi.isValid(); ++mfi)
	{
	  const int* lo   = mfi.validbox().loVect();
	  const int* hi   = mfi.validbox().hiVect();
	  
	  FArrayBox& rg       = (*Rhs[lev])[mfi];  
	  const FArrayBox& bx = beta[0][lev][mfi];
	  const FArrayBox& by = beta[1][lev][mfi];
	  
	  DEF_LIMITS (rg,rg_dat,rglo,rghi);
	  DEF_CLIMITS(bx,bx_dat,bxlo,bxhi);
	  DEF_CLIMITS(by,by_dat,bylo,byhi);
#if (BL_SPACEDIM==3)
	  const FArrayBox& bz = beta[2][lev][mfi];
	  DEF_CLIMITS(bz,bz_dat,bzlo,bzhi);
#endif
	  FORT_DRHOG_RICHARD (rg_dat, ARLIM(rglo), ARLIM(rghi),
			      bx_dat, ARLIM(bxlo), ARLIM(bxhi),
			      by_dat, ARLIM(bylo), ARLIM(byhi),
#if (BL_SPACEDIM==3)
			      bz_dat, ARLIM(bzlo), ARLIM(bzhi),
#endif
			      &density[0],lo,hi,dx,&gravity,&dt);
	}

      // res  = res + \phi (n(t+dt))
      const BoxArray& grd = alpha[lev].boxArray();
      MultiFab Tmp(grd,1,0);
      MultiFab::Copy(Tmp,*S[lev],0,0,1,0);
      for (MFIter mfi(Tmp); mfi.isValid(); ++mfi)
	{
	  const Box& box = mfi.validbox();
	  Tmp[mfi].mult(alpha[lev][mfi],box,0,0,1);
	  (*Rhs[lev])[mfi].plus(Tmp[mfi],box,0,0,1);
	}
    }

    for (int lev = nlevs-2; lev >= 0; lev--)
    {
      PorousMedia* pm = dynamic_cast<PorousMedia*>(&parent->getLevel(lev));
      pm->avgDown(Rhs[lev],lev,Rhs[lev+1],lev+1);
    }
}


#endif

void 
Diffusion::compute_flux( int                    nc,
			 Real                   dt,
			 Real                   be_cn_theta,
			 MultiFab* const*       flux,
			 const MultiFab*        pc,
			 const MultiFab* const* beta)
{
  // 
  // Compute flux 
  //
  MultiFab Soln(grids,1,1);
  Soln.setVal(0.);
  Real a = 0.0;
  Real b = (1.0-be_cn_theta)*dt;
  ViscBndry visc_bndry;
  const BCRec& bc     = caller->get_desc_lst()[State_Type].getBC(nc);
  const Real cur_time = caller->get_state_data(State_Type).curTime();
  IntVect ref_ratio   = level > 0 ? 
    parent->refRatio(level-1) : IntVect::TheUnitVector();
  getBndryData(visc_bndry,nc,1,cur_time);
  visc_bndry.setScalarValues(bc,ref_ratio,pc);
  MultiFab* alpha = 0;
  ABecLaplacian* visc_op = getViscOp(nc,a,b,cur_time,visc_bndry,
				     0,0,beta,alpha,true);
  MultiFab::Copy(Soln,*pc,0,0,1,1);
  visc_op->compFlux(D_DECL(*flux[0],*flux[1],*flux[2]),Soln);
  for (int i = 0; i < BL_SPACEDIM; ++i)
    (*flux[i]).mult(b/dt/caller->Geom().CellSize()[i]);

  delete visc_op;
}

#ifdef MG_USE_FBOXLIB
void
Diffusion::richard_iter (Real                   dt,
			 int                    nc,
			 Real                   gravity,
			 Array<Real>            density,
			 MultiFab&              res_fix,
			 const MultiFab*        alpha, 
			 const MultiFab* const* beta,
			 const MultiFab* const* beta_dp,
			 MultiFab*              pc,
			 MultiFab*              umac,
			 const bool             do_upwind,
			 Real*                  err_nwt)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::richard_iter()");
  //
  // This routine solves the time-dependent richards equation
  //
  MultiFab& S_new = caller->get_new_data(State_Type);
  Real snorm = S_new.norm2(0);
  // setup multifabs for solvers
  MultiFab Rhs(grids,1,0);
  MultiFab Soln(grids,1,1);
  Rhs.setVal(0.);
  Soln.setVal(0.);
  PorousMedia* pm_level = dynamic_cast<PorousMedia*>(&parent->getLevel(level));
  // Compute residual
  Real a = 0;
  Real b = -dt*density[0];
  ViscBndry visc_bndry;

  const BCRec& bc     = caller->get_desc_lst()[Press_Type].getBC(0);
  const Real cur_time = caller->get_state_data(State_Type).curTime();
  IntVect ref_ratio   = level > 0 ? 
    parent->refRatio(level-1) : IntVect::TheUnitVector();
  getBndryData(visc_bndry,nc,1,cur_time);
  visc_bndry.setScalarValues(bc,ref_ratio,pc);
  ABecLaplacian* visc_op = getViscOp(visc_bndry,a,b,0,0,beta,alpha,false);
  MultiFab::Copy(Rhs,res_fix,0,0,1,0);
  residual_richard(visc_op,dt*density[0],gravity,density,Rhs,pc,beta,alpha);
  Real prev_res_norm = Rhs.norm2(0);

  // preconditioning the residual.
  // If beta_dp is the exact Jacobian, then this is the Newton step.
  if (beta_dp != 0)
    {
      Rhs.mult(-1.0);
      MultiFab* phi_p[1];
      MultiFab* Rhs_p[1];
      phi_p[0] = &Soln;
      Rhs_p[0] = &Rhs;
      int nc_opt = 2;
      Real a_dp = 0.0;
      if (alpha) a_dp = 1.0;
      Real b_dp = dt*density[0];  
      visc_bndry.setdeltaSValues(bc,ref_ratio);

      Array< Array<Real> > xa(1);
      Array< Array<Real> > xb(1);
      PArray<MultiFab> a2_p;
      PArray<MultiFab> a1_p;
      Array<PArray<MultiFab> > bb_p(BL_SPACEDIM);    
      const Real S_tol     = visc_tol;
      const Real S_tol_abs = get_scaled_abs_tol(Rhs, visc_abs_tol);
      Real final_resnorm;
      int  fill_bcs_for_gradient = 1;
      MGT_Solver mgt_solver = getOp(0,1,xa,xb,cur_time,visc_bndry,bc,true);
      coefs_fboxlib_mg (a1_p, bb_p, alpha, beta_dp, a_dp, false);
      mgt_solver.set_porous_coefficients(a1_p, a2_p, bb_p, b_dp, xa, xb, nc_opt);    
      mgt_solver.solve(phi_p, Rhs_p, S_tol, S_tol_abs, 
		       visc_bndry, fill_bcs_for_gradient, final_resnorm);
    }
  
//   // Construct krylov subspace
//   bool do_newton_krylov = false;
//   if (do_newton_krylov)
//     {
//       MultiFab::Copy(Rhs,Soln,0,0,1,0);

//       int krylov_dim = 15;
//       int H_dim = krylov_dim - 1;
//       MultiFab krylov(grids,krylov_dim,0);
//       MultiFab krylovn(grids,1,0);
//       MultiFab krylovnp1(grids,1,0);
//       Real H[H_dim*H_dim];
//       for (int i = 0; i<H_dim*H_dim; i++)
// 	H[i] = 0.0;
//       Real rhs_norm = sqrt(dotxy(Rhs,Rhs));
//       MultiFab::Copy(krylov,Rhs,0,0,1,0);
//       for (MFIter mfi(krylov); mfi.isValid(); ++mfi)
// 	krylov[mfi].mult(1.0/rhs_norm,0);
//       for (int j=0; j<H_dim;j++)
// 	{
// 	  Soln.setVal(0.);
// 	  MultiFab::Copy(Soln,krylov,j,0,1,0);     
// 	  jac_richard(nc,gravity,density,res_fix,Soln,Rhs,visc_op,visc_op_dp);
// 	  MultiFab::Copy(krylov,Soln,0,j+1,1,0);
// 	  MultiFab::Copy(krylovnp1,Soln,0,0,1,0);
	  
// 	  // note that H is column first here
// 	  for (int i=0;i<=j;i++)
// 	    {
// 	      MultiFab::Copy(krylovn,krylov,i,0,1,0);
// 	      H[i+j*H_dim] = dotxy(krylovnp1,krylovn);
// 	      krylovn.mult(-H[i+j*H_dim]);
// 	      MultiFab::Add(krylov,krylovn,0,j+1,1,0);
// 	    }
// 	  Real k_norm = sqrt(dotxy(krylov,krylov,j+1,j+1));
// 	  for (MFIter mfi(krylov); mfi.isValid(); ++mfi)
// 	    krylov[mfi].mult(1.0/k_norm,j+1);
// 	  if (j+1 < H_dim)
// 	    H[j+1+j*H_dim] = k_norm;
// 	}
//       // Arnoldi step
//       Real ak[H_dim];
//       FORT_ARNOLDI(ak,H,&rhs_norm,&H_dim);
//       Soln.setVal(0.);
//       for (MFIter mfi(krylov); mfi.isValid(); ++mfi)
// 	{
// 	  for  (int i=0;i<H_dim;i++)
// 	    {
// 	      krylov[mfi].mult(ak[i],i);
// 	      Soln[mfi].plus(krylov[mfi],mfi.validbox(),i,0,1);
// 	    }
// 	}
//     }
  // line search 
  MultiFab Stmp(grids,1,1);
  MultiFab Rhsp1(grids,1,0);
  Stmp.setVal(0.);
  Rhsp1.setVal(0.);
  MultiFab::Copy(Stmp,S_new,nc,0,1,1);
  MultiFab::Add(Stmp,Soln,0,0,1,0);

  // determie new capillary pressure
  MultiFab pctmp(grids,1,1);
  pm_level->calcCapillary(&pctmp,Stmp);
  // determine new lambda_1 and beta
  MultiFab lambda(grids,pm_level->ncomps,1);
  pm_level->calcLambda(&lambda,Stmp);
  MultiFab** betatmp;
  allocFluxBoxesLevel(betatmp,0,1);
  pm_level->calc_richard_coef(betatmp,&lambda,umac,0,do_upwind);
  // Compute residual
  setBeta (visc_op,0,betatmp,false);
  MultiFab::Copy(Rhsp1,res_fix,0,0,1,0);
  residual_richard(visc_op,dt*density[0],gravity,density,Rhsp1,&pctmp,betatmp,alpha,&Stmp);
  Real rhsp1_norm = Rhsp1.norm2(0);
  Real alphak = 1.0;
  int iter = 0;
  while (iter < 10 && rhsp1_norm > prev_res_norm && alphak > 1.e-3) 
    {
      Rhsp1.setVal(0.);
      alphak = 0.1*alphak;
      Soln.mult(0.1);
      MultiFab::Copy(Stmp,S_new,nc,0,1,1);
      MultiFab::Add(Stmp,Soln,0,0,1,0);
      pm_level->calcCapillary(&pctmp,Stmp);
      pm_level->calcLambda(&lambda,Stmp);
      pm_level->calc_richard_coef(betatmp,&lambda,umac,0,do_upwind);
      setBeta (visc_op,0,betatmp,false);
      MultiFab::Copy(Rhsp1,res_fix,0,0,1,0);
      residual_richard(visc_op,dt*density[0],gravity,density,Rhsp1,&pctmp,betatmp,alpha,&Stmp);
      rhsp1_norm = Rhsp1.norm2(0);
      iter += 1;
    }

  MultiFab::Copy(S_new,Stmp,0,nc,1,0);
  MultiFab::Copy(S_new,Rhsp1,0,pm_level->ncomps+pm_level->ntracers,1,0);

  // Compute the err_nwt
  *err_nwt = rhsp1_norm/snorm;

  removeFluxBoxesLevel(betatmp);
  delete visc_op;
}

void
Diffusion::richard_composite_iter (Real                      dt,
				   int                       nlevs,
				   int                       nc,
				   Real                      gravity,
				   Array<Real>               density,
				   PArray<MultiFab>&         res_fix,
				   PArray<MultiFab>&         alpha, 
				   Array<PArray<MultiFab> >& beta,
				   Array<PArray<MultiFab> >& beta_dp,
				   PArray<MultiFab>&         pc,
				   const bool                do_upwind,
				   Real*                     err_nwt)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::richard_composite_iter()");
  //
  // This routine solves the time-dependent richards equation based on 
  // composite solve.
  //

  // setup multilevel objects
  Real b = -dt*density[0];
  MultiFab* Rhs[nlevs];
  MultiFab* Soln[nlevs];
  std::vector<BoxArray> bav(nlevs);
  std::vector<DistributionMapping> dmv(nlevs);
  std::vector<Geometry> geom(nlevs);    
  for (int lev = 0; lev < nlevs; lev++) 
    {
      MultiFab& S = parent->getLevel(lev).get_new_data(State_Type);
      bav[lev]  = S.boxArray();
      dmv[lev]  = S.DistributionMap();
      geom[lev] = parent->Geom(lev+level);

      Rhs[lev]  = new MultiFab(bav[lev],1,0);
      Rhs[lev]->setVal(0.);
      Soln[lev] = new MultiFab(bav[lev],1,1);
      MultiFab::Copy(*Soln[lev],S,nc,0,1,1);
    }
  Real snorm = (*Soln[0]).norm2(0);
  
  ViscBndry visc_bndry;     
  const BCRec& bc     = caller->get_desc_lst()[Press_Type].getBC(0);
  const Real cur_time = caller->get_state_data(State_Type).curTime();
  IntVect ref_ratio   = level > 0 ? 
    parent->refRatio(level-1) : IntVect::TheUnitVector();
  getBndryData(visc_bndry,nc,1,cur_time);
  visc_bndry.setScalarValues(bc,ref_ratio,&pc[0]);

  int nc_opt = 2;
  Array< Array<Real> > xa(nlevs);
  Array< Array<Real> > xb(nlevs);
  PArray<MultiFab> a1_p(nlevs,PArrayManage);
  PArray<MultiFab> a2_p;

  for (int lev=0; lev<nlevs; lev++)
    {
      a1_p.set(lev,new MultiFab(alpha[lev].boxArray(),alpha[lev].nComp(),
				alpha[lev].nGrow()));
      a1_p[lev].setVal(0.);
    }

  MGT_Solver mgt_solver = getOp(0,nlevs,geom,bav,dmv, xa, xb,cur_time,
				visc_bndry,bc,true);
  
  // Setup RHS
  mgt_solver.set_maxorder(2);
  mgt_solver.set_porous_coefficients(a1_p, a2_p, beta, b, xa, xb, nc_opt); 
  residual_richard(mgt_solver,dt*density[0],gravity,density,Rhs,
		   pc,beta,alpha,res_fix,Soln,visc_bndry);
  Real prev_res_norm = (*Rhs[0]).norm2(0);

  // preconditioning the residual.
  // If beta_dp is the exact Jacobian, then this is the Newton step.
  
  if (beta_dp.size() > 0)
    {
      for (int lev=0; lev<nlevs; lev++)
	{
	  (*Rhs[lev]).mult(-1.0);
	  (*Soln[lev]).setVal(0.);
	}

      Real b_dp = dt*density[0];   
      const Real S_tol     = visc_tol;
      const Real S_tol_abs = visc_abs_tol;
      Real final_resnorm;
      int  fill_bcs_for_gradient = 1; 
      mgt_solver.set_maxorder(2);
      mgt_solver.set_porous_coefficients(alpha, a2_p, beta_dp, b_dp, xa, xb, nc_opt);    
      mgt_solver.solve(Soln, Rhs, S_tol, S_tol_abs, 
		       fill_bcs_for_gradient, final_resnorm);

      for (int lev = nlevs-2; lev >= 0; lev--)
	{
	  PorousMedia* pm = dynamic_cast<PorousMedia*>(&parent->getLevel(lev));
	  pm->avgDown(Soln[lev],lev,Soln[lev+1],lev+1);
	}
    }
  
  // line search 
  MultiFab* Stmp[nlevs];
  MultiFab* Rhsp1[nlevs];
  PArray<MultiFab> pctmp(nlevs, PArrayManage);
  Array<PArray<MultiFab> > betatmp(BL_SPACEDIM);

  for (int dir=0; dir<BL_SPACEDIM; dir++)
    betatmp[dir].resize(nlevs,PArrayNoManage);
  
  for (int lev=0; lev<nlevs; lev++)
    {  
      PorousMedia* pm = dynamic_cast<PorousMedia*>(&parent->getLevel(lev));

      pctmp.set(lev,new MultiFab(bav[lev],1,1));

      for (int dir=0; dir<BL_SPACEDIM; dir++)
	{
	  BoxArray ba = bav[lev];
	  ba.surroundingNodes(dir);
	  betatmp[dir].set(lev, new MultiFab(ba,beta[dir][lev].nComp(),0));
	}

      Stmp[lev]  = new MultiFab(bav[lev],1,1);
      MultiFab& S = pm->get_new_data(State_Type);
      MultiFab::Copy(*Stmp[lev],S,nc,0,1,1);
      MultiFab::Add(*Stmp[lev],*Soln[lev],0,0,1,0);
      
      MultiFab* tmp_betatmp[BL_SPACEDIM];
      MultiFab* tmp_umac = pm->u_mac_curr;
      for (int dir=0;dir<BL_SPACEDIM;dir++)
	tmp_betatmp[dir] = &betatmp[dir][lev];
      pm->calcCapillary(&pctmp[lev],*Stmp[lev]);
      MultiFab lambda(bav[lev],pm->ncomps,1);
      pm->calcLambda(&lambda,*Stmp[lev]);
      pm->calc_richard_coef(tmp_betatmp,&lambda,tmp_umac,0,do_upwind);
      Rhsp1[lev] = new MultiFab(bav[lev],1,0);
    }
  mgt_solver.set_maxorder(2);
  mgt_solver.set_porous_coefficients(a1_p, a2_p, betatmp, b, xa, xb, nc_opt); 
  residual_richard(mgt_solver,dt*density[0],gravity,density,Rhsp1,
		   pctmp,betatmp,alpha,res_fix,Stmp,visc_bndry);

  Real rhsp1_norm = (*Rhsp1[0]).norm2(0);
  Real alphak = 1.0;
  int iter = 0;
  while (iter < 10 && rhsp1_norm > prev_res_norm && alphak > 1.e-3) 
    {
        for (int lev=0; lev<nlevs; lev++)
	  {
	    (*Rhsp1[lev]).setVal(0.);
	    alphak = 0.1*alphak;
	    (*Soln[lev]).mult(0.1);
	    MultiFab& S = parent->getLevel(lev).get_new_data(State_Type);
	    MultiFab::Copy(*Stmp[lev],S,nc,0,1,1);
	    MultiFab::Add(*Stmp[lev],*Soln[lev],0,0,1,0);
      
	    PorousMedia* pm = dynamic_cast<PorousMedia*>(&parent->getLevel(lev));

	    MultiFab* tmp_betatmp[BL_SPACEDIM];
	    MultiFab* tmp_umac = pm->u_mac_curr;
	    for (int dir=0;dir<BL_SPACEDIM;dir++)
	      tmp_betatmp[dir] = &betatmp[dir][lev];

	    MultiFab lambda(bav[lev],pm->ncomps,1);
	    pm->calcLambda(&lambda,*Stmp[lev]);
	    pm->calcCapillary(&pctmp[lev],*Stmp[lev]);
	    pm->calc_richard_coef(tmp_betatmp,&lambda,tmp_umac,0,do_upwind);
	  }
	mgt_solver.set_porous_coefficients(a1_p, a2_p, betatmp, b, xa, xb, nc_opt); 
	residual_richard(mgt_solver,dt*density[0],gravity,density,
			 Rhsp1,pctmp,betatmp,alpha,res_fix,Stmp,visc_bndry);

	rhsp1_norm = (*Rhsp1[0]).norm2(0);
    }

  for (int lev=0; lev<nlevs;lev++)
    {
      MultiFab& S = parent->getLevel(lev).get_new_data(State_Type);
      MultiFab::Copy(S,*Stmp[lev],0,nc,1,1);
    }

  // Compute the err_nwt
  *err_nwt = rhsp1_norm/snorm;

  for (int lev = 0; lev < nlevs; lev++)
    {
      delete Soln[lev];
      delete Rhs[lev];
      delete Stmp[lev];
      delete Rhsp1[lev];
    }

}

void
Diffusion::richard_iter_p (Real                   dt,
			   int                    nc,
			   Real                   gravity,
			   Array<Real>            density,
			   MultiFab&              res_fix,
			   const MultiFab*        alpha, 
			   const MultiFab*        dalpha,
			   const MultiFab* const* beta,
			   const MultiFab* const* beta_dp,
			   MultiFab*              pc,
			   MultiFab*              umac,
			   const bool             do_upwind,
			   Real*                  err_nwt)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::richard_iter_p()");
 
 //
  // This routine solves the time-dependent richards equation
  //
  MultiFab& S_new = caller->get_new_data(State_Type);
  Real snorm = S_new.norm2(0);

  // setup multifabs for solvers
  MultiFab Rhs(grids,1,0);
  MultiFab Soln(grids,1,1);
  Rhs.setVal(0.);
  Soln.setVal(0.);

  PorousMedia* pm_level = dynamic_cast<PorousMedia*>(&parent->getLevel(level));
  // Compute residual
  Real a = 0;
  Real b = -dt*density[0];
  ViscBndry visc_bndry;

  const BCRec& bc     = caller->get_desc_lst()[Press_Type].getBC(0);
  const Real cur_time = caller->get_state_data(State_Type).curTime();
  IntVect ref_ratio   = level > 0 ? 
    parent->refRatio(level-1) : IntVect::TheUnitVector();
  getBndryData(visc_bndry,nc,1,cur_time);
  visc_bndry.setScalarValues(bc,ref_ratio,pc);
  ABecLaplacian* visc_op = getViscOp(visc_bndry,a,b,0,0,beta,alpha,false);
  MultiFab::Copy(Rhs,res_fix,0,0,1,0);
  residual_richard(visc_op,dt*density[0],gravity,density,Rhs,pc,beta,alpha);
  Real prev_res_norm = Rhs.norm2(0);

  // preconditioning the residual.
  // If beta_dp is the exact Jacobian, then this is the Newton step.
  if (beta_dp != 0)
    {
      Rhs.mult(-1.0);
      MultiFab* phi_p[1];
      MultiFab* Rhs_p[1];
      phi_p[0] = &Soln;
      Rhs_p[0] = &Rhs;
      int nc_opt = 2;
      Real a_dp = 0.0;
      if (dalpha) a_dp = 1.0;
      Real b_dp = dt*density[0];  
      visc_bndry.setdeltaSValues(bc,ref_ratio);

      Array< Array<Real> > xa(1);
      Array< Array<Real> > xb(1);
      PArray<MultiFab> a2_p;
      PArray<MultiFab> a1_p;
      Array<PArray<MultiFab> > bb_p(BL_SPACEDIM);    
      const Real S_tol     = visc_tol;
      const Real S_tol_abs = get_scaled_abs_tol(Rhs, visc_abs_tol);
      Real final_resnorm;
      int  fill_bcs_for_gradient = 1;
      MGT_Solver mgt_solver = getOp(0,1,xa,xb,cur_time,visc_bndry,bc,true);
      coefs_fboxlib_mg (a1_p, bb_p, dalpha, beta_dp, a_dp, false);
      mgt_solver.set_porous_coefficients(a1_p, a2_p, bb_p, b_dp, xa, xb, nc_opt);    
      mgt_solver.solve(phi_p, Rhs_p, S_tol, S_tol_abs, 
		       visc_bndry, fill_bcs_for_gradient, final_resnorm);
    }

  // line search 
  MultiFab Stmp(grids,1,1);
  MultiFab pctmp(grids,1,1);
  MultiFab Rhsp1(grids,1,0);
  Stmp.setVal(0.);
  Rhsp1.setVal(0.);
  MultiFab::Copy(pctmp,*pc,nc,0,1,1);
  MultiFab::Add(pctmp,Soln,0,0,1,0);
  int iter = 0;
  // Real pcmin = 1.e20;
//   for (MFIter mfi(pctmp); mfi.isValid(); ++mfi)
//     pcmin = std::min(pcmin,pctmp[mfi].min(mfi.validbox(),0));
//   const int IOProc = ParallelDescriptor::IOProcessorNumber();
//   ParallelDescriptor::ReduceRealMin(&pcmin, 1, IOProc);
//   while (pcmin < 0. && iter < 10)
//     {
//       Soln.mult(0.3);
//       MultiFab::Copy(pctmp,*pc,nc,0,1,1);
//       MultiFab::Add(pctmp,Soln,0,0,1,0);

//       pcmin = 1.e20;
//       for (MFIter mfi(pctmp); mfi.isValid(); ++mfi)
// 	pcmin = std::min(pcmin,pctmp[mfi].min(mfi.validbox(),0));
//       ParallelDescriptor::ReduceRealMin(&pcmin, 1, IOProc);
//       iter++;
//     }

  pm_level->calcInvCapillary(Stmp,pctmp);
  // determine new lambda_1 and beta
  MultiFab lambda(grids,pm_level->ncomps,1);
  pm_level->calcLambda(&lambda,Stmp);
  MultiFab** betatmp;
  allocFluxBoxesLevel(betatmp,0,1);
  pm_level->calc_richard_coef(betatmp,&lambda,umac,0,do_upwind);
  // Compute residual
  setBeta (visc_op,0,betatmp,false);
  MultiFab::Copy(Rhsp1,res_fix,0,0,1,0);
  residual_richard(visc_op,dt*density[0],gravity,density,Rhsp1,&pctmp,betatmp,alpha,&Stmp);
  Real rhsp1_norm = Rhsp1.norm2(0);
  Real alphak = 1.0;
  iter = 0;
  while (iter < 10 && rhsp1_norm > prev_res_norm && alphak > 1.e-3) 
    {
      Rhsp1.setVal(0.);
      alphak = 0.1*alphak;
      Soln.mult(0.1);
      MultiFab::Copy(pctmp,*pc,nc,0,1,1);
      MultiFab::Add(pctmp,Soln,0,0,1,0);
      pm_level->calcInvCapillary(Stmp,pctmp);
      pm_level->calcLambda(&lambda,Stmp);
      pm_level->calc_richard_coef(betatmp,&lambda,umac,0,do_upwind);
      setBeta (visc_op,0,betatmp,false);
      MultiFab::Copy(Rhsp1,res_fix,0,0,1,0);
      residual_richard(visc_op,dt*density[0],gravity,density,Rhsp1,&pctmp,betatmp,alpha,&Stmp);
      rhsp1_norm = Rhsp1.norm2(0);
      iter += 1;
    }

  MultiFab::Copy(S_new,Stmp,0,nc,1,0);
  MultiFab::Copy(S_new,Rhsp1,0,pm_level->ncomps+pm_level->ntracers,1,0);

  // Compute the err_nwt
  *err_nwt = rhsp1_norm/snorm;

  removeFluxBoxesLevel(betatmp);
  delete visc_op;
}

void
Diffusion::richard_iter_eqb (int                    nc,
			     Real                   gravity,
			     Array<Real>            density,
			     MultiFab&              res_fix,
			     const MultiFab* const* beta,
			     const MultiFab* const* beta_dp,
			     MultiFab*              pc,
			     MultiFab*              umac,
			     const bool             do_upwind,
			     Real*                  err_nwt)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::richard_iter_eqb()");
  //
  // This routine solves the equilibrium richard's equation
  //

  MultiFab* alpha = 0;
  Real dt = 1.0;
  richard_iter (dt,nc,gravity,density,res_fix,alpha, 
		beta,beta_dp,pc,umac,do_upwind,err_nwt);
  
}

void
Diffusion::jac_richard (int                    nc,
			Real                   gravity,
			Array<Real>            density,
			MultiFab&              res_fix,
			MultiFab&              soln,
			MultiFab&              soln_old,
			ABecLaplacian*         visc_op,
			ABecLaplacian*         visc_op_dp)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::jac_richard()");
  //
  // This routine computes jacobian*soln and preconditioned it with beta_dp
  //
  MultiFab& S_new = caller->get_new_data(State_Type);
  MultiFab res_old(grids,1,0);
  MultiFab::Copy(res_old,soln_old,0,0,1,0);
  res_old.mult(-1.0);

  // determine perturbed solution
  Real deps = 1.e-8;
  soln.mult(deps);
  MultiFab::Add(soln,S_new,nc,0,1,1);
  // determie new capillary pressure
  PorousMedia* pm_level = dynamic_cast<PorousMedia*>(&parent->getLevel(level));
  MultiFab pc(grids,1,1);
  pm_level->calcCapillary(&pc,soln);
  // determine new lambda_1 and beta
  MultiFab lambda(grids,pm_level->ncomps,1);
  pm_level->calcLambda(&lambda,soln);
  MultiFab** betatmp;
  allocFluxBoxesLevel(betatmp,0,1);
  pm_level->calc_richard_coef(betatmp,&lambda,0);
  // Compute residual
  setBeta (visc_op,0,betatmp);
  MultiFab::Copy(soln,res_fix,0,0,1,0);
  residual_richard(visc_op,1.0,gravity,density,soln,&pc,betatmp);
  removeFluxBoxesLevel(betatmp);
  if (visc_op_dp != 0)
    {
      // preconditioning the residual
      MultiFab Soln(grids,1,0);
      Soln.setVal(0.); 
      const BCRec& bc      = caller->get_desc_lst()[State_Type].getBC(nc);
      const Real S_tol     = visc_tol;
      const Real S_tol_abs = get_scaled_abs_tol(soln, visc_tol);
      solve(visc_op_dp,Soln,soln,bc,S_tol,S_tol_abs);
      MultiFab::Copy(soln,Soln,0,0,1,0);
    }

  // Jr \approx (F(s+\epsilon*r) - F(s))/\epsilon
  MultiFab::Add(soln,res_old,0,0,1,0);
  soln.mult(1.0/deps);
}

void 
Diffusion::richard_flux( int                    nc,
			 Real                   be_cn_theta,
			 Real                   gravity,
			 Array<Real>&           density,
			 MultiFab* const*       flux,
			 const MultiFab*        pc,
			 const MultiFab* const* beta)
{
  // 
  // Compute flux 
  //
  MultiFab Soln(grids,1,1);
  Soln.setVal(0.);
  Real a = 0.0;
  Real b = be_cn_theta;
  ViscBndry visc_bndry;

  const BCRec& bc     = caller->get_desc_lst()[Press_Type].getBC(0);
  const Real cur_time = caller->get_state_data(State_Type).curTime();
  IntVect ref_ratio   = level > 0 ? 
    parent->refRatio(level-1) : IntVect::TheUnitVector();
  getBndryData(visc_bndry,nc,1,cur_time);
  visc_bndry.setScalarValues(bc,ref_ratio,pc);
  MultiFab* alpha = 0;
  ABecLaplacian* visc_op = getViscOp(visc_bndry,a,b,0,0,beta,alpha,false);
  MultiFab::Copy(Soln,*pc,0,0,1,1);
  visc_op->compFlux(D_DECL(*flux[0],*flux[1],*flux[2]),Soln);
  for (int i = 0; i < BL_SPACEDIM; ++i)
    (*flux[i]).mult(b);
  delete visc_op;

  const Real* dx   = caller->Geom().CellSize();
  for (MFIter mfi(Soln); mfi.isValid(); ++mfi)
    {
      const int* lo   = mfi.validbox().loVect();
      const int* hi   = mfi.validbox().hiVect();
	
      FArrayBox& fx = (*flux[0])[mfi];
      FArrayBox& fy = (*flux[1])[mfi];
      const FArrayBox& bx = (*beta[0])[mfi];
      const FArrayBox& by = (*beta[1])[mfi];
	
      DEF_LIMITS(fx,fx_dat,fxlo,fxhi);
      DEF_LIMITS(fy,fy_dat,fylo,fyhi);
      DEF_CLIMITS(bx,bx_dat,bxlo,bxhi);
      DEF_CLIMITS(by,by_dat,bylo,byhi);
	
#if (BL_SPACEDIM==3)
      FArrayBox& fz = (*flux[2])[mfi];
      const FArrayBox& bz = (*beta[2])[mfi];
      DEF_LIMITS(fz,fz_dat,fzlo,fzhi);
      DEF_CLIMITS(bz,bz_dat,bzlo,bzhi);
#endif
	FORT_FRHOG_RICHARD (fx_dat, ARLIM(fxlo), ARLIM(fxhi),
			    bx_dat, ARLIM(bxlo), ARLIM(bxhi),
			    fy_dat, ARLIM(fylo), ARLIM(fyhi),
			    by_dat, ARLIM(bylo), ARLIM(byhi),
#if (BL_SPACEDIM==3)
			    fz_dat, ARLIM(fzlo), ARLIM(fzhi),
			    bz_dat, ARLIM(bzlo), ARLIM(bzhi),
#endif
			    &density[0],lo,hi,dx,&gravity,&b);
    }

  for (int i = 0; i < BL_SPACEDIM; ++i)
    {
      for (MFIter fmfi(*flux[i]); fmfi.isValid(); ++fmfi)
	{
	  const int idx = fmfi.index();
	  (*flux[i])[idx].mult((area[i])[idx],0,0,1);
	}
    }
}

void
Diffusion::check_consistency  (Real                   dt,
			       int                    nc,
			       Real                   pm,
			       Real                   be_cn_theta,
			       int                    dataComp,
			       const MultiFab*        alpha, 
			       const MultiFab* const* betan, 
			       const MultiFab* const* betanp1,
			       const MultiFab*        pcn, 
			       const MultiFab*        pcnp1,
			       MultiFab*              S_nwt,
			       Real*                  err_nwt)
{

  // At this point, S_old has bndry at time N, S_new has bndry at time N+1
  MultiFab& S_new = caller->get_new_data(State_Type);

  MultiFab Rhs(grids,1,0);
  MultiFab Tmpn(grids,1,0);
  MultiFab Tmpnp1(grids,1,0);
  MultiFab Soln(grids,1,1);

  Rhs.setVal(0.);
  Tmpn.setVal(0.);
  Tmpnp1.setVal(0.);
  Soln.setVal(0.);

  // Set up Rhs.
  Real a = 0.0;
  Real b = pm*(1.0-be_cn_theta)*dt;

  // (dt/2) * \nabla H^n \nabla P_c^n
  ViscBndry visc_bndry_0;
  const BCRec& bc      = caller->get_desc_lst()[State_Type].getBC(nc);
  const Real prev_time = caller->get_state_data(State_Type).prevTime();
  IntVect ref_ratio    = level > 0 ? parent->refRatio(level-1) : IntVect::TheUnitVector();
  getBndryData(visc_bndry_0,nc,1,prev_time);
    
  visc_bndry_0.setScalarValues(bc,ref_ratio,pcn);
  ABecLaplacian* visc_op = getViscOp(nc,a,b,prev_time,visc_bndry_0,
				     0,dataComp,betan,alpha,true);
  visc_op->maxOrder(max_order);
  MultiFab::Copy(Soln,*pcn,0,0,1,1);
  visc_op->apply(Tmpn,Soln);
  delete visc_op;
  // (dt/2) * \nabla H^{n+1} \nabla P_c^{n+1}
  const Real cur_time = caller->get_state_data(State_Type).curTime();
  visc_bndry_0.setScalarValues(bc,ref_ratio,pcnp1);
  visc_op = getViscOp(nc,a,b,cur_time,visc_bndry_0,
		      0,dataComp,betanp1,alpha,true);
  visc_op->maxOrder(max_order);
  MultiFab::Copy(Soln,*pcnp1,0,0,1,1);
  visc_op->apply(Tmpnp1,Soln);
    
  delete visc_op;

  // RHS = (dt/2) * (\nabla H^n \nabla P_c^n + \nabla H^{n+1} \nabla P_c^{n+1}).
  for (MFIter mfi(Rhs); mfi.isValid(); ++mfi)   {
    const Box& box = mfi.validbox();
    Rhs[mfi].plus(Tmpn[mfi],box,0,0,1);
    Rhs[mfi].plus(Tmpnp1[mfi],box,0,0,1);
  }

  // RHS = RHS - s^{n+1,k}
  MultiFab::Copy(Soln,S_new,nc,0,1,1);
  for (MFIter mfi(Soln); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.validbox();
      Soln[mfi].mult(volume[mfi],box,0,0,1);
      if (alpha!=0)
	Soln[mfi].mult((*alpha)[mfi],box,0,0,1);
      Rhs[mfi].minus(Soln[mfi],box,0,0,1);
    }
    
    
  // RHS = RHS + phi*s^{n} + dt*(div F)
  MultiFab::Copy(Soln,*S_nwt,0,0,1,1);
  for (MFIter mfi(Soln); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.validbox();
      Soln[mfi].mult(volume[mfi],box,0,0,1);
      if (alpha!=0)
	Soln[mfi].mult((*alpha)[mfi],box,0,0,1);
      Rhs[mfi].plus(Soln[mfi],box,0,0,1);
    }

  // Compute the err_nwt
  *err_nwt = Rhs.norm2(0);   
}
#endif

void
Diffusion::diffuse_Ssync (MultiFab*              Ssync,
                          int                    nc,
                          Real                   dt,
                          Real                   be_cn_theta,
                          MultiFab* const*       flux,
                          int                    dataComp,
                          const MultiFab* const* beta,
                          const MultiFab*        alpha)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::diffuse_Ssync()");

  const int IOProc    = ParallelDescriptor::IOProcessorNumber();

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "Diffusion::diffuse_Ssync: "
	      << caller->get_desc_lst()[State_Type].name(nc) << '\n';

  int allnull, allthere;
  checkBeta(beta, allthere, allnull);

  MultiFab Soln(grids,1,1);
  MultiFab Rhs(grids,1,0);

  Soln.setVal(0);
  MultiFab::Copy(Rhs,*Ssync,nc,0,1,0);

  if (verbose > 1)
    {
      Real r_norm = Rhs.norm0(0);
      ParallelDescriptor::ReduceRealMax(r_norm,IOProc);

      if (ParallelDescriptor::IOProcessor())
	std::cout << "Original max of Ssync " << r_norm << '\n';
    }

  //
  // SET UP COEFFICIENTS FOR VISCOUS SOLVER.
  //
  const Real a = 1.0;
  Real       b = be_cn_theta*dt;
  if (allnull)
    b *= visc_coef[nc];
  Real           rhsscale = 1.0;
  ABecLaplacian* visc_op  = getViscOp(nc,a,b,&rhsscale,dataComp,beta,alpha);
  //
  // Compute RHS.
  //
  for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
      Rhs[Rhsmfi].mult(volume[Rhsmfi]); 
    }
  Rhs.mult(rhsscale,0,1);
  //
  // Construct solver and call it.
  //
  const Real S_tol     = visc_tol;
  const Real S_tol_abs = get_scaled_abs_tol(Rhs, visc_tol);
  if (use_cg_solve)
    {
      CGSolver cg(*visc_op,use_mg_precond_flag);
      cg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }
  else
    {
      MultiGrid mg(*visc_op);
      mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }
  int flux_allthere, flux_allnull;
  checkBeta(flux, flux_allthere, flux_allnull);
  if (flux_allthere)
    {
      visc_op->compFlux(D_DECL(*flux[0],*flux[1],*flux[2]),Soln);
      for (int i = 0; i < BL_SPACEDIM; ++i)
	(*flux[i]).mult(b/(dt*caller->Geom().CellSize()[i]),0);
    }

  MultiFab::Multiply(Soln,*rho,0,0,1,0);
  MultiFab::Copy(*Ssync,Soln,0,nc,1,0);
  if (verbose > 1)
    {
      Real s_norm = Soln.norm0(0);
      ParallelDescriptor::ReduceRealMax(s_norm,IOProc);
      if (ParallelDescriptor::IOProcessor())
	std::cout << "Final max of Ssync " << s_norm << '\n';
    }
    
  delete visc_op;
}

void
Diffusion::diffuse_Ssync_CPL (MultiFab*              Ssync,
			      int                    nc,
			      Real                   dt,
			      Real                   be_cn_theta,
			      MultiFab* const*       flux,
			      int                    dataComp,
			      const MultiFab*        alpha,
			      const MultiFab* const* betan_dp)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::diffuse_Ssync_CPL()");

  const int IOProc    = ParallelDescriptor::IOProcessorNumber();

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "Diffusion::diffuse_Ssync_CPL: "
	      << caller->get_desc_lst()[State_Type].name(nc) << '\n';

  MultiFab Soln(grids,1,1);
  MultiFab Rhs(grids,1,0);
  Soln.setVal(0);
  MultiFab::Copy(Rhs,*Ssync,nc,0,1,0);
  if (verbose > 1)
    {
      Real r_norm = Rhs.norm0(0);  
      ParallelDescriptor::ReduceRealMax(r_norm,IOProc);
      if (ParallelDescriptor::IOProcessor())
	std::cout << "Original max of Ssync " << r_norm << '\n';
    }

  //
  // SET UP COEFFICIENTS FOR VISCOUS SOLVER.
  //
  const Real a = 1.0;
  Real       b = -be_cn_theta*dt;
  ViscBndry visc_bndry_0;
  const Real prev_time = caller->get_state_data(State_Type).prevTime();
  getBndryData(visc_bndry_0,nc,1,prev_time);
  ABecLaplacian* visc_op = getViscOp(nc,a,b,prev_time,visc_bndry_0,
				     0,dataComp,betan_dp,alpha,true);
    
  //
  // Compute RHS.
  //
  for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
      Rhs[Rhsmfi].mult(volume[Rhsmfi]); 
    }

  //
  // Construct solver and call it.
  //
  const Real S_tol     = visc_tol;
  const Real S_tol_abs = get_scaled_abs_tol(Rhs, visc_tol);
  if (use_cg_solve)
    {
      CGSolver cg(*visc_op,use_mg_precond_flag);
      cg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }
  else
    {
      MultiGrid mg(*visc_op);
      mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }

  //delete visc_op;
  //visc_bndry_0.setScalarValues(bc,ref_ratio,pcn);
  //visc_op = getViscOp(nc,a,b,prev_time,visc_bndry_0,
  //			rho_flag,0,dataComp,betan,alpha,true);
    
  visc_op->compFlux(D_DECL(*flux[0],*flux[1],*flux[2]),Soln);
  for (int i = 0; i < BL_SPACEDIM; ++i)
    (*flux[i]).mult(-b/(dt*caller->Geom().CellSize()[i]),0);
    
  MultiFab::Copy(*Ssync,Soln,0,nc,1,0);
    
  if (verbose > 1)
    {
      Real s_norm = Soln.norm0(0);
      ParallelDescriptor::ReduceRealMax(s_norm,IOProc);
      if (ParallelDescriptor::IOProcessor())
	std::cout << "Final max of Ssync " << s_norm << '\n';
    }
    
  delete visc_op;
}

void
Diffusion::solve(ABecLaplacian*    visc_op,
		 MultiFab&         Soln,
		 MultiFab&         Rhs,
		 const BCRec&      bc, 
		 Real              S_tol,
		 Real              S_tol_abs)
{
  //
  // Construct solver and call it.
  //
  if (use_cg_solve)
    {
      CGSolver cg(*visc_op,use_mg_precond_flag);
      cg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }

#ifdef MG_USE_FBOXLIB
  else if ( use_fboxlib_mg )
    {
      std::vector<BoxArray> bav(1);
      bav[0] = Rhs.boxArray();
      std::vector<DistributionMapping> dmv(1);
      dmv[0] = Rhs.DistributionMap();
      bool nodal = false;
      std::vector<Geometry> geom(1);
      const BndryData visc_bndry = (*visc_op).bndryData();
      geom[0] = visc_bndry.getGeom();
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
	      mg_bc[i*2 + 0] = bc.lo(i)==EXT_DIR? MGT_BC_DIR : MGT_BC_NEU;
	      mg_bc[i*2 + 1] = bc.hi(i)==EXT_DIR? MGT_BC_DIR : MGT_BC_NEU;
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
	const Real* dx_crse   = parent->Geom(level-1).CellSize();
	for ( int i = 0; i < BL_SPACEDIM; ++i ) {
	  xa[0][i] = 0.5 * dx_crse[i];
	  xb[0][i] = 0.5 * dx_crse[i];
	}
      }

      MGT_Solver mgt_solver(geom, mg_bc, bav, dmv, nodal);
      mgt_solver.set_maxorder(max_order);

      const MultiFab* aa_p[1];
      aa_p[0] = &(visc_op->aCoefficients());
      const MultiFab* bb_p[1][BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i )
        {
	  bb_p[0][i] = &(visc_op->bCoefficients(i));
        }
      Real b = visc_op->get_beta();
      mgt_solver.set_visc_coefficients(aa_p, bb_p, b, xa, xb);

      MultiFab* phi_p[1];
      MultiFab* Rhs_p[1];
      phi_p[0] = &Soln;
      Rhs_p[0] = &Rhs;
      Real final_resnorm;
      mgt_solver.solve(phi_p, Rhs_p, S_tol, S_tol_abs, 
		       visc_bndry,final_resnorm);
    }
#endif
  else
    {
      MultiGrid mg(*visc_op);
      mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }
}

ABecLaplacian*
Diffusion::getViscOp (int                    comp,
                      Real                   a,
                      Real                   b,
                      Real                   time,
                      ViscBndry&             visc_bndry,
                      Real*                  rhsscale,
                      int                    dataComp,
                      const MultiFab* const* beta,
                      const MultiFab*        alpha_in,
                      bool		     bndry_already_filled)
{
  int allnull, allthere;
  checkBeta(beta, allthere, allnull);

  const Real* dx = caller->Geom().CellSize();

  if (!bndry_already_filled)
    getBndryData(visc_bndry,comp,1,time);
  ABecLaplacian* visc_op = new ABecLaplacian(visc_bndry,dx);
  visc_op->maxOrder(max_order);

  //
  // alpha should be the same size as volume.
  //
  MultiFab alpha(grids,1,GEOM_GROW);
  if (a != 0) 
    {
      for (MFIter alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
	{
	  const int  i  = alphamfi.index();
	  const Box& bx = alphamfi.validbox();

	  const int*  lo      = bx.loVect();
	  const int*  hi      = bx.hiVect();
	  Real*       dat     = alpha[alphamfi].dataPtr();
	  Box         abx     = BoxLib::grow(bx,alpha.nGrow());
	  const int*  alo     = abx.loVect();
	  const int*  ahi     = abx.hiVect();
	  const Real* voli    = volume[alphamfi].dataPtr();
	  Box         vbox    = BoxLib::grow(volume.box(i),volume.nGrow());
	  const int*  vlo     = vbox.loVect();
	  const int*  vhi     = vbox.hiVect();
	  
	  Real*       rhodat  = (*rho)[alphamfi].dataPtr();
	  Box         rbx     = BoxLib::grow(bx,rho->nGrow());
	  const int*  rlo     = rbx.loVect();
	  const int*  rhi     = rbx.hiVect();

	  FORT_SETALPHA(dat, ARLIM(alo), ARLIM(ahi),
			rhodat, ARLIM(rlo), ARLIM(rhi),
			lo, hi, voli, ARLIM(vlo), ARLIM(vhi));

	  if (alpha_in !=0)
	    alpha[i].mult((*alpha_in)[i],bx,dataComp,0,1);
	}
    }
  else
    alpha.setVal(0.);

  if (rhsscale != 0)
    {
      *rhsscale = scale_abec ? 1.0/alpha.max(0) : 1.0;
      visc_op->setScalars(a*(*rhsscale),b*(*rhsscale));
    }
  else
    {
      visc_op->setScalars(a,b);
    }
  if (allnull)
    {
      visc_op->aCoefficients(alpha);
      for (int n = 0; n < BL_SPACEDIM; n++)
        {
	  MultiFab bcoeffs(area[n].boxArray(),1,0);
	  bcoeffs.copy(area[n]);
	  bcoeffs.mult(dx[n]);
	  visc_op->bCoefficients(bcoeffs,n);
	    
        }
    }
  else
    {
      visc_op->aCoefficients(alpha);
      for (int n = 0; n < BL_SPACEDIM; n++)
        {
	  MultiFab bcoeffs(area[n].boxArray(),1,0);
	  bcoeffs.copy(area[n]);
	  for (MFIter bcoeffsmfi(bcoeffs); bcoeffsmfi.isValid(); ++bcoeffsmfi)
            {
	      const int i = bcoeffsmfi.index();
	      bcoeffs[i].mult((*beta[n])[i],dataComp,0,1);
	      bcoeffs[i].mult(dx[n]);
            }	    
	  visc_op->bCoefficients(bcoeffs,n);	    
        }
    }
  return visc_op;
}

ABecLaplacian*
Diffusion::getViscOp (ViscBndry&             bndry,
                      Real                   a,
                      Real                   b,
                      Real*                  rhsscale,
                      int                    dataComp,
                      const MultiFab* const* beta,
                      const MultiFab*        alpha_in,
		      bool                   do_volume)
{
  //
  // Note: Construct the ABecLaplacian based on given ViscBndry.
  //
  int allnull, allthere;
  checkBeta(beta, allthere, allnull);

  const Geometry& geom = caller->Geom();
  const Real*  dx      = geom.CellSize();
  ABecLaplacian* visc_op = new ABecLaplacian(bndry,dx);
  visc_op->maxOrder(2);
  //
  // alpha should be the same size as volume.
  //
  MultiFab alpha(grids,1,GEOM_GROW);
  if (a != 0) 
    {
      for (MFIter alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
	{
	  const int   i       = alphamfi.index();
	  const Box&  bx      = alphamfi.validbox();

	  if (do_volume)
	    {
	      const int*  lo      = bx.loVect();
	      const int*  hi      = bx.hiVect();
	      Real*       dat     = alpha[alphamfi].dataPtr();
	      Box         abx     = BoxLib::grow(bx,alpha.nGrow());
	      const int*  alo     = abx.loVect();
	      const int*  ahi     = abx.hiVect();
	      const Real* voli    = volume[alphamfi].dataPtr();
	      Box         vbox    = BoxLib::grow(volume.box(i),volume.nGrow());
	      const int*  vlo     = vbox.loVect();
	      const int*  vhi     = vbox.hiVect();
	      
	      Real*       rhodat  = (*rho)[alphamfi].dataPtr();
	      Box         rbx     = BoxLib::grow(bx,rho->nGrow());
	      const int*  rlo     = rbx.loVect();
	      const int*  rhi     = rbx.hiVect();
	      
	      FORT_SETALPHA(dat, ARLIM(alo), ARLIM(ahi),
			    rhodat, ARLIM(rlo), ARLIM(rhi),
			    lo, hi, voli, ARLIM(vlo), ARLIM(vhi));
	    }
	  if (alpha_in !=0)
	    alpha[i].mult((*alpha_in)[i],bx,dataComp,0,1);
	}  
    }
  else
    alpha.setVal(0.);


  if (rhsscale != 0)
    {
      *rhsscale = scale_abec ? 1.0/alpha.max(0) : 1.0;

      visc_op->setScalars(a*(*rhsscale),b*(*rhsscale));
    }
  else
    {
      visc_op->setScalars(a,b);
    }
  visc_op->aCoefficients(alpha);

  if (allnull)
    {
      for (int n = 0; n < BL_SPACEDIM; n++)
        {
	  MultiFab bcoeffs(area[n].boxArray(),1,0);
	  bcoeffs.copy(area[n]);
	  bcoeffs.mult(dx[n],0,1,0);
	  visc_op->bCoefficients(bcoeffs,n);
        }
    }
  else
    {
      for (int n = 0; n < BL_SPACEDIM; n++)
        {
	  MultiFab bcoeffs(area[n].boxArray(),1,0);
	  bcoeffs.setVal(1.);
	  if (do_volume) bcoeffs.copy(area[n]);
	  for (MFIter bcoeffsmfi(bcoeffs); bcoeffsmfi.isValid(); ++bcoeffsmfi)
            {
	      const int i = bcoeffsmfi.index();
	      bcoeffs[i].mult((*beta[n])[i],dataComp,0,1);
            }
	  if (do_volume) bcoeffs.mult(dx[n],0,1,0);
	  visc_op->bCoefficients(bcoeffs,n);
        }
    }

  return visc_op;
}

ABecLaplacian*
Diffusion::getViscOp (int                    comp,
                      Real                   a,
                      Real                   b,
                      Real*                  rhsscale,
                      int                    dataComp,
                      const MultiFab* const* beta,
                      const MultiFab*        alpha_in)
{
  //
  // Note: This assumes that the "NEW" density is to be used, if rho_flag==2
  //
  int allnull, allthere;
  checkBeta(beta, allthere, allnull);

  const Geometry& geom = caller->Geom();
  const Real*  dx      = geom.CellSize();
  const BCRec& bc      = caller->get_desc_lst()[State_Type].getBC(comp);

  IntVect ref_ratio = level > 0 ? parent->refRatio(level-1) : IntVect::TheUnitVector();

  ViscBndry bndry(grids,1,geom);
  bndry.setHomogValues(bc, ref_ratio);

  ABecLaplacian* visc_op = new ABecLaplacian(bndry,dx);
  visc_op->maxOrder(max_order);
  //
  // alpha should be the same size as volume.
  //
  MultiFab alpha(grids,1,GEOM_GROW);
  if (a != 0) 
    {
      for (MFIter alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
	{
	  const int   i       = alphamfi.index();
	  const Box&  bx      = alphamfi.validbox();
	  const int*  lo      = bx.loVect();
	  const int*  hi      = bx.hiVect();
	  Real*       dat     = alpha[alphamfi].dataPtr();
	  Box         abx     = BoxLib::grow(bx,alpha.nGrow());
	  const int*  alo     = abx.loVect();
	  const int*  ahi     = abx.hiVect();
	  const Real* voli    = volume[alphamfi].dataPtr();
	  Box         vbox    = BoxLib::grow(volume.box(i),volume.nGrow());
	  const int*  vlo     = vbox.loVect();
	  const int*  vhi     = vbox.hiVect();
	  
	  Real*       rhodat  = (*rho)[alphamfi].dataPtr();
	  Box         rbx     = BoxLib::grow(bx,rho->nGrow());
	  const int*  rlo     = rbx.loVect();
	  const int*  rhi     = rbx.hiVect();
	  
	  FORT_SETALPHA(dat, ARLIM(alo), ARLIM(ahi),
			rhodat, ARLIM(rlo), ARLIM(rhi),
			lo, hi, voli, ARLIM(vlo), ARLIM(vhi));

	  if (alpha_in !=0)
	    alpha[i].mult((*alpha_in)[i],bx,dataComp,0,1);
	}  
    }
  else
    alpha.setVal(0.);


  if (rhsscale != 0)
    {
      *rhsscale = scale_abec ? 1.0/alpha.max(0) : 1.0;

      visc_op->setScalars(a*(*rhsscale),b*(*rhsscale));
    }
  else
    {
      visc_op->setScalars(a,b);
    }
  visc_op->aCoefficients(alpha);

  if (allnull)
    {
      for (int n = 0; n < BL_SPACEDIM; n++)
        {
	  MultiFab bcoeffs(area[n].boxArray(),1,0);
	  bcoeffs.copy(area[n]);
	  bcoeffs.mult(dx[n],0,1,0);
	  visc_op->bCoefficients(bcoeffs,n);
        }
    }
  else
    {
      for (int n = 0; n < BL_SPACEDIM; n++)
        {
	  MultiFab bcoeffs(area[n].boxArray(),1,0);
	  bcoeffs.copy(area[n]);
	  for (MFIter bcoeffsmfi(bcoeffs); bcoeffsmfi.isValid(); ++bcoeffsmfi)
            {
	      const int i = bcoeffsmfi.index();
	      bcoeffs[i].mult((*beta[n])[i],dataComp,0,1);
            }
	  bcoeffs.mult(dx[n],0,1,0);
	  visc_op->bCoefficients(bcoeffs,n);
        }
    }

  return visc_op;
}

void
Diffusion::setBeta (ABecLaplacian*         visc_op,
		    int                    dataComp,
		    const MultiFab* const* beta,
		    bool                   do_volume)
{
  //
  // Note: set Beta of visc_op
  //
  const Geometry& geom = caller->Geom();
  const Real*  dx      = geom.CellSize();
  for (int n = 0; n < BL_SPACEDIM; n++)
    {
      MultiFab bcoeffs(area[n].boxArray(),1,0);
      bcoeffs.setVal(1.);
      if (do_volume) bcoeffs.copy(area[n]);
      for (MFIter mfi(bcoeffs); mfi.isValid(); ++mfi)
	{
	  const int i = mfi.index();
	  bcoeffs[i].mult((*beta[n])[i],dataComp,0,1);
	}
      if (do_volume) bcoeffs.mult(dx[n],0,1,0);
      visc_op->bCoefficients(bcoeffs,n);
    }
}

#ifdef MG_USE_FBOXLIB
MGT_Solver
Diffusion::getOp (int                   comp,
		  int                   nlevs,
		  Array< Array<Real> >& xa,
		  Array< Array<Real> >& xb,
		  Real                  time,
		  ViscBndry&            visc_bndry,
		  const BCRec&          bc, 
		  bool	                bndry_already_filled)
{

  bool nodal = false;

  std::vector<Geometry> geom(nlevs);
  std::vector<BoxArray> bav(nlevs);
  std::vector<DistributionMapping> dmv(nlevs);

  if (!bndry_already_filled)
    getBndryData(visc_bndry,comp,1,time);

  for (int lev = 0; lev < nlevs; lev++) 
  {
    MultiFab& P = parent->getLevel(lev+level).get_data(State_Type,time);
    bav[lev]  = P.boxArray();
    dmv[lev]  = P.DistributionMap();
    geom[lev] = parent->Geom(lev+level);
  }
  
  int mg_bc[2*BL_SPACEDIM];
  for ( int dir = 0; dir < BL_SPACEDIM; ++dir )
  {
    if ( geom[0].isPeriodic(dir))
    {
      mg_bc[dir*2 + 0] = 0;
      mg_bc[dir*2 + 1] = 0;
    }
    else
    {
      mg_bc[dir*2 + 0] = bc.lo(dir)==EXT_DIR? MGT_BC_DIR : MGT_BC_NEU;
      mg_bc[dir*2 + 1] = bc.hi(dir)==EXT_DIR? MGT_BC_DIR : MGT_BC_NEU;
    }
  }

  for (int lev = 0; lev < nlevs; lev++) 
  {
    xa[lev].resize(BL_SPACEDIM);
    xb[lev].resize(BL_SPACEDIM);
    if (lev+level == 0) 
    { 
      for ( int dir = 0; dir < BL_SPACEDIM; ++dir ) 
      {
	xa[lev][dir] = 0.;
	xb[lev][dir] = 0.; 
      }    
    } 
    else 
    {
      const Real* dx_crse = parent->Geom(lev+level-1).CellSize();
      for ( int dir = 0; dir < BL_SPACEDIM; ++dir ) 
      {
	xa[lev][dir] = 0.5 * dx_crse[dir];
	xb[lev][dir] = 0.5 * dx_crse[dir];
      }
    }
  }

  MGT_Solver mgt_solver(geom, mg_bc, bav, dmv, nodal);
  mgt_solver.set_maxorder(2);
    
  return mgt_solver;
}

MGT_Solver
Diffusion::getOp (int                              comp,
		  int                              nlevs,
		  std::vector<Geometry>            geom,
		  std::vector<BoxArray>            bav,
		  std::vector<DistributionMapping> dmv,
		  Array< Array<Real> >&            xa,
		  Array< Array<Real> >&            xb,
		  Real                             time,
		  ViscBndry&                       visc_bndry,
		  const BCRec&                     bc, 
		  bool	                           bndry_already_filled)
{
  bool nodal = false;

  if (!bndry_already_filled)
    getBndryData(visc_bndry,comp,1,time);

  int mg_bc[2*BL_SPACEDIM];
  for ( int dir = 0; dir < BL_SPACEDIM; ++dir )
  {
    if (geom[0].isPeriodic(dir))
    {
      mg_bc[dir*2 + 0] = 0;
      mg_bc[dir*2 + 1] = 0;
    }
    else
    {
      mg_bc[dir*2 + 0] = bc.lo(dir)==EXT_DIR? MGT_BC_DIR : MGT_BC_NEU;
      mg_bc[dir*2 + 1] = bc.hi(dir)==EXT_DIR? MGT_BC_DIR : MGT_BC_NEU;
    }
  }

  for (int lev = 0; lev < nlevs; lev++) 
  {
    xa[lev].resize(BL_SPACEDIM);
    xb[lev].resize(BL_SPACEDIM);
    if (lev == 0) 
    { 
      for ( int dir = 0; dir < BL_SPACEDIM; ++dir ) 
      {
	xa[lev][dir] = 0.;
	xb[lev][dir] = 0.; 
      }    
    } 
    else 
    {
      const Real* dx_crse   = geom[lev-1].CellSize();
      for ( int dir = 0; dir < BL_SPACEDIM; ++dir ) 
      {
	xa[lev][dir] = 0.5 * dx_crse[dir];
	xb[lev][dir] = 0.5 * dx_crse[dir];
      }
    }
  }

  MGT_Solver mgt_solver(geom, mg_bc, bav, dmv, nodal);
  mgt_solver.set_maxorder(max_order);
    
  return mgt_solver;
}

void
Diffusion::coefs_fboxlib_mg (PArray<MultiFab>&          acoefs,
			     Array< PArray<MultiFab> >& bcoefs,
			     const MultiFab*            alpha_in,
			     const MultiFab* const*     beta,
			     const Real                 a,
			     bool                       do_volume)
{		
  //
  // Compute the alpha and beta for fboxlib.
  // - assume the alpha_in and beta has not been multiplied by volume
  //
  if (a > 0) 
    {
      if (acoefs.size() == 0)
	{
	  acoefs.resize(1,PArrayManage);
	  acoefs.set(0,new MultiFab((*alpha_in).boxArray(),1,0,Fab_allocate));
	}    
      acoefs[0].setVal(1.0);
      for (MFIter mfi(acoefs[0]); mfi.isValid(); ++mfi)
	{
	  const int i = mfi.index();
	  BL_ASSERT(grids[i] == mfi.validbox());
	  if (do_volume)
	    acoefs[0][i].copy(volume[i],mfi.validbox(),0,mfi.validbox(),0,1);
	  acoefs[0][i].mult((*alpha_in)[i],mfi.validbox(),0,0,1);
	}
    }
  else
    {
      if (acoefs.size() == 0)
	{
	  acoefs.resize(1,PArrayManage);
	  acoefs.set(0,new MultiFab(grids,1,0,Fab_allocate));
	} 
      acoefs[0].setVal(0.);
    }


  int n_comps = beta[0][0].nComp();
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      if (bcoefs[dir].size()==0)
	bcoefs[dir].resize(1,PArrayManage);
      bcoefs[dir].set(0,new MultiFab((*beta[dir]).boxArray(),
				     n_comps,0,Fab_allocate));
    }


  // Do multicomponent-beta.
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {  
    bcoefs[dir][0].setVal(1.0);
    // Initialize coefficients with area    
    if (do_volume)
      {
	for (int n = 0; n < n_comps; n++)
	  bcoefs[dir][0].copy(area[dir],0,n,1);
      }

    // Multiply coefficients by dx
    const Real* dx = caller->Geom().CellSize();

    // Multiply each coefficient by that component's lambda
    for (MFIter mfi(bcoefs[dir][0]); mfi.isValid(); ++mfi)
    {
      const int idx = mfi.index();
      for (int n = 0; n < n_comps; n++) 
	bcoefs[dir][0][idx].mult((*beta[dir])[idx],n,n,1);
      if (do_volume)
	bcoefs[dir][0][idx].mult(dx[dir]);
    }
  }
}

void
Diffusion::coefs_fboxlib_mg (PArray<MultiFab>&          acoefs,
			     Array< PArray<MultiFab> >& bcoefs,
			     PArray<MultiFab>&          alpha,
			     Array< PArray<MultiFab> >& beta,
			     const Real                 a,
			     int                        nlevs)
{		
  //
  // Compute the alpha and beta for fboxlib.
  // Quantities are not multiplied by volume.
  //
  if (acoefs.size() == 0)
    {
      acoefs.resize(nlevs,PArrayManage);
      for (int lev=0; lev<nlevs; lev++)
	{
	  acoefs.set(lev,new MultiFab(alpha[lev].boxArray(),1,0,Fab_allocate));
	  acoefs[lev].setVal(0.);
	}
    }  
 
  if (a > 0) 
    {
      for (int lev=0; lev<nlevs; lev++)
	MultiFab::Copy(acoefs[lev],alpha[lev],0,0,1,0);
    }  


  int n_comps = beta[0][0].nComp();
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      if (bcoefs[dir].size()==0)
	{
	  bcoefs[dir].resize(nlevs,PArrayManage);
	  for (int lev=0; lev<nlevs; lev++)
	    {
	      bcoefs[dir].set(lev,new MultiFab(beta[dir][lev].boxArray(),
					       n_comps,0,Fab_allocate));
	      MultiFab::Copy(bcoefs[dir][lev],beta[dir][lev],0,0,n_comps,0);
	    }
	}
    }

}
#endif

void
Diffusion::getViscTerms (MultiFab&              visc_terms,
                         int                    src_comp,
                         int                    comp, 
                         Real                   time,
                         int                    dataComp,
                         const MultiFab* const* beta)
{
  int allnull, allthere;
  checkBeta(beta, allthere, allnull);
  //
  // Before computing the godunov predictors we may have to
  // precompute the viscous source terms.  To do this we must
  // construct a Laplacian operator, set the coeficients and apply
  // it to the time N data.  First, however, we must precompute the
  // fine N bndry values.  We will do this for each scalar that diffuses.
  //
  // Note: This routine DOES NOT fill grow cells
  //
  const Real* dx = caller->Geom().CellSize();
  MultiFab&   S  = caller->get_data(State_Type,time);
  int ngrow = visc_terms.nGrow();
  visc_terms.setVal(0.0,comp-src_comp,1,ngrow);
  //
  // FIXME
  // LinOp classes cannot handle multcomponent MultiFabs yet,
  // construct the components one at a time and copy to visc_terms.
  //
  MultiFab visc_tmp(grids,1,1);
  MultiFab s_tmp(grids,1,1);

  if (is_diffusive[comp])
    {
      ViscBndry visc_bndry;
      getBndryData(visc_bndry,comp,1,time);
      //
      // Set up operator and apply to compute viscous terms.
      //
      const Real a = 0.0;
      const Real b = allnull ? -visc_coef[comp] : -1.0;

      ABecLaplacian visc_op(visc_bndry,dx);

      visc_op.setScalars(a,b);
      visc_op.maxOrder(max_order);

      if (allnull)
        {
	  for (int n = 0; n < BL_SPACEDIM; n++)
            {
	      MultiFab bcoeffs(area[n].boxArray(),1,0);
	      bcoeffs.copy(area[n]);
	      bcoeffs.mult(dx[n]);
	      visc_op.bCoefficients(bcoeffs,n);
	    }
        }
      else
        {
	  for (int n = 0; n < BL_SPACEDIM; n++)
            {
	      MultiFab bcoeffs(area[n].boxArray(),1,0);
	      bcoeffs.copy(area[n]);
 
	      for (MFIter bcoeffsmfi(bcoeffs); bcoeffsmfi.isValid(); ++bcoeffsmfi)
                {
		  const int i = bcoeffsmfi.index();
		  bcoeffs[i].mult((*beta[n])[i],dataComp,0,1);
		  bcoeffs[i].mult(dx[n]);
                }
	      visc_op.bCoefficients(bcoeffs,n);
            }
        }
      //
      // Copy to single component multifab for operator classes.
      //

      MultiFab::Copy(s_tmp,S,comp,0,1,1);
      visc_op.apply(visc_tmp,s_tmp);
      //
      // Must divide by volume.
      //
      for (MFIter visc_tmpmfi(visc_tmp); visc_tmpmfi.isValid(); ++visc_tmpmfi)
        {
	  const int i = visc_tmpmfi.index();
	  BL_ASSERT(grids[i] == visc_tmpmfi.validbox());
	  visc_tmp[i].divide(volume[i],visc_tmpmfi.validbox(),0,0,1);
        }
      MultiFab::Copy(visc_terms,visc_tmp,0,comp-src_comp,1,0);
    }
}

void
Diffusion::getCplViscTerms (MultiFab&              visc_terms,
			    int                    src_comp,
			    Real                   time,
			    const Real*            density,
			    int                    dataComp,
			    const MultiFab* const* beta,
			    const MultiFab*        pc)
{ 
  //
  // Before computing the godunov predictors we may have to
  // precompute the viscous source terms due to capillary pressure.  
  // We construct a Laplacian operator, set the coeficients and apply
  // it to the time N data.  First, however, we must precompute the
  // fine N bndry values.  We will do this for each scalar that diffuses.
  //
  // Note: This routine DOES NOT fill grow cells
  //

  
  MultiFab* alpha  = 0;
  const BCRec& bc   = caller->get_desc_lst()[State_Type].getBC(src_comp);
  IntVect ref_ratio = level > 0 ? 
    parent->refRatio(level-1) : IntVect::TheUnitVector();

  const Real a =  0.0;
  const Real b =  1.0;
  
  ViscBndry visc_bndry;
  getBndryData(visc_bndry,src_comp,1,time);
  visc_bndry.setScalarValues(bc,ref_ratio,pc);
  ABecLaplacian* visc_op = getViscOp(src_comp,a,b,time,visc_bndry,
				     0,0,beta,alpha,true);
  MultiFab Tmpn(grids,1,0);
  MultiFab Soln(grids,1,1);
  MultiFab::Copy(Soln,*pc,0,0,1,1);
  visc_op->apply(Tmpn,Soln);

  delete visc_op;

  for (MFIter mfi(Tmpn); mfi.isValid(); ++mfi)
    {
      const int i = mfi.index();
      BL_ASSERT(grids[i] == mfi.validbox());
      Tmpn[i].divide(volume[i],mfi.validbox(),0,0,1);
    }

  MultiFab::Add(visc_terms,Tmpn,0,dataComp+src_comp,1,0);

  int nd = 1;
  if (src_comp == 1)
    nd = 0;
  
  Tmpn.mult(-density[nd]/density[src_comp]);
  MultiFab::Add(visc_terms,Tmpn,0,dataComp+nd,1,0);

}

void
Diffusion::getBndryData (ViscBndry& bndry,
                         int        src_comp,
                         int        num_comp,
                         Real       time)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getBndryData()");

  BL_ASSERT(num_comp == 1);
  //
  // Fill phys bndry vals of grow cells of (tmp) multifab passed to bndry.
  //
    
  const int     nGrow = 1;
  const BCRec&  bc    = caller->get_desc_lst()[State_Type].getBC(src_comp);
  MultiFab S(grids, num_comp, nGrow);
  S.setVal(BL_SAFE_BOGUS);
  bndry.define(grids,num_comp,caller->Geom());
  PorousMedia* pm_level = dynamic_cast<PorousMedia*>(&parent->getLevel(level));
  MFIter Rho_mfi(*rho);

  for (FillPatchIterator Phi_fpi(*caller,S,nGrow,time,State_Type,0,NUM_SCALARS);
       Rho_mfi.isValid() && Phi_fpi.isValid(); 
       ++Rho_mfi, ++Phi_fpi)
    {
      pm_level->dirichletStateBC(Phi_fpi(),nGrow,time);
      const BoxList gCells = BoxLib::boxDiff(Phi_fpi().box(),Phi_fpi.validbox());
      for (BoxList::const_iterator bli = gCells.begin();
	   bli != gCells.end();
	   ++bli)
        {
	  S[Phi_fpi].copy(Phi_fpi(),*bli,src_comp,*bli,0,num_comp);
	  for (int n=0; n<num_comp; ++n)
	    S[Phi_fpi].divide((*rho)[Rho_mfi],*bli,0,n,1);
        }
    }
  if (level == 0)
    {
      bndry.setBndryValues(S,0,0,num_comp,bc);
    }
  else
    {
      BoxArray cgrids = grids;
      cgrids.coarsen(crse_ratio);
      BndryRegister crse_br(cgrids,0,1,2,num_comp);
      //
      // interp for solvers over ALL c-f brs, need safe data.
      //
      crse_br.setVal(BL_BOGUS);
      coarser->FillBoundary(crse_br,src_comp,0,num_comp,time);
      bndry.setBndryValues(crse_br,0,S,0,0,num_comp,crse_ratio,bc);
    }
}

void
Diffusion::FillBoundary (BndryRegister& bdry,
                         int            state_ind,
                         int            dest_comp,
                         int            num_comp,
                         Real           time)
{
  //
  // Need one grow cell filled for linear solvers.
  // We assume filPatch gets this right, where possible.
  //
  const int     nGrow = 1;
    
  MultiFab S(caller->boxArray(),num_comp,nGrow);
  S.setVal(BL_SAFE_BOGUS);
  MFIter Rho_mfi(*rho);
  PorousMedia* pm_level = dynamic_cast<PorousMedia*>(&parent->getLevel(level));
  for (FillPatchIterator S_fpi(*caller,S,nGrow,time,State_Type,0,NUM_SCALARS);
       Rho_mfi.isValid() && S_fpi.isValid(); ++Rho_mfi, ++S_fpi)
    {
      pm_level->dirichletStateBC(S_fpi(),nGrow,time);

      S[S_fpi.index()].copy(S_fpi(),state_ind,0,num_comp);
      for (int n=0; n<num_comp; ++n)
	S[S_fpi.index()].divide((*rho)[Rho_mfi],0,n,1);	
    }
  //
  // Copy into boundary register.
  //
  bdry.copyFrom(S,nGrow,0,dest_comp,num_comp);    
}


void
Diffusion::checkBetas (const MultiFab* const* beta1, 
                       const MultiFab* const* beta2,
                       int&                   allthere,
                       int&                   allnull) const
{
  int allnull1, allnull2, allthere1, allthere2;

  checkBeta(beta1,allthere1,allnull1);
  checkBeta(beta2,allthere2,allnull2);
  allnull  = allnull1 && allnull2;
  allthere = allthere1 && allthere2;

  if (!(allthere || allnull))
    BoxLib::Abort("Diffusion::checkBetas(): betas must either be all 0 or all non-0");
}

void
Diffusion::checkBeta (const MultiFab* const* beta,
                      int&                   allthere,
                      int&                   allnull) const
{
  allnull  = 1;
  allthere = beta != 0;

  if (allthere)
    {
      for (int d = 0; d < BL_SPACEDIM; d++)
        {
	  allnull = allnull && beta[d] == 0;
	  allthere = allthere && beta[d] != 0;
        }
    }

  if (!(allthere || allnull))
    BoxLib::Abort("Diffusion::checkBeta(): betas must be all 0 or all non-0");
}

void
Diffusion::checkBeta (const MultiFab* const* beta,
                      int&                   allthere) const
{
  allthere = beta != 0;

  if (allthere)
    {
      for (int d = 0; d < BL_SPACEDIM; d++)
	allthere = allthere && beta[d] != 0;
    }

  if (!allthere)
    BoxLib::Abort("Diffusion::checkBeta(): betas must be all non-0");
}

void
Diffusion::allocFluxBoxesLevel (MultiFab**& fluxbox, 
                                int         nghost,
                                int         nvar)
{
  fluxbox = new MultiFab*[BL_SPACEDIM];
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      BoxArray edge_boxes(grids);
      edge_boxes.surroundingNodes(dir);
      fluxbox[dir] = new MultiFab(edge_boxes,nvar,nghost);
    }
}

void
Diffusion::removeFluxBoxesLevel (MultiFab**& fluxbox) 
{
  if (fluxbox != 0)
    {
      for (int i = 0; i<BL_SPACEDIM; i++)
	delete fluxbox[i];
      delete [] fluxbox;
      fluxbox = 0;
    }
}

bool
Diffusion::are_any(const Array<DiffusionForm>& diffusionType,
                   const DiffusionForm         testForm,
                   const int                   sComp,
                   const int                   nComp)
{
  for (int comp = sComp; comp < sComp + nComp; ++comp)
    {
      if (diffusionType[comp] == testForm)
	return true;
    }

  return false;
}

int
Diffusion::how_many(const Array<DiffusionForm>& diffusionType,
                    const DiffusionForm         testForm,
                    const int                   sComp,
                    const int                   nComp)
{
  int counter = 0;

  for (int comp = sComp; comp < sComp + nComp; ++comp)
    {
      if (diffusionType[comp] == testForm) 
	++counter;
    }

  return counter;
}

