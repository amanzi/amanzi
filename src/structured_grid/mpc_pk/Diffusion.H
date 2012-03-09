
#ifndef _Diffusion_H_
#define _Diffusion_H_

//
// Comment out this line to use diffusion class outside
// the context of PorousMedia and classes derived from it.
//
#define USE_POROUSMEDIA 1

#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <LevelBld.H>
#include <BC_TYPES_PM.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>
#include <ABecLaplacian.H>
#include <ViscBndry.H>

#ifdef MG_USE_FBOXLIB
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#endif

#ifndef _PorousMedia_H_
enum StateType {State_Type=0}; 
enum StateNames  { N_w_g=0, N_w_l=0, N_o_g=0, N_o_l,  Enthalpy};
//enum StateNames  { N_w, N_o, Enthalpy};
#endif
//
// Useful enumeration of the different forms of the diffusion terms
//
enum DiffusionForm { RhoInverse_Laplacian_S, Laplacian_SoverRho, Laplacian_S };

class Diffusion 
{
public:

    enum SolveMode {PREDICTOR, CORRECTOR, ONEPASS};

    Diffusion ();

    Diffusion (Amr*               Parent,
               AmrLevel*          Caller,
               Diffusion*         coarser,
               int                num_state,
               FluxRegister*      Viscflux_reg,
               MultiFab&          Volume,
               MultiFab*          Area,
               const Array<int>&  _is_diffusive,
               const Array<Real>& _visc_coef);


    Diffusion (Amr*               Parent,
               AmrLevel*          Caller,
               Diffusion*         coarser,
               int                num_state,
               FluxRegister*      Viscflux_reg,
               MultiFab&          Volume,
               MultiFab*          Area);

    virtual ~Diffusion ();

    void echo_settings () const;

    FluxRegister* viscFluxReg ();

    Real get_scaled_abs_tol (const MultiFab& rhs,
                             Real            reduction) const;

    void set_rho(const MultiFab* rho_in);

    void set_rho(Real rho_in);

    void diffuse_scalar (Real                   dt,
			 int                    sigma,
			 Real                   be_cn_theta,
			 MultiFab* const*       fluxn,
			 MultiFab* const*       fluxnp1,
			 int                    dataComp = 0,
			 MultiFab*              delta_rhs = 0, 
			 const MultiFab*        alpha = 0, 
			 const MultiFab* const* betan = 0, 
			 const MultiFab* const* betanp1 = 0,
			 const SolveMode&       solve_mode = ONEPASS);

    void diffuse_init_CPL(Real                   dt,
			  int                    nc,
			  Real                   be_cn_theta,
			  MultiFab* const*       fluxn,
			  int                    dataComp = 0,
			  MultiFab*              delta_rhs = 0, 
			  const MultiFab*        alpha = 0, 
			  const MultiFab* const* betan = 0, 
			  const MultiFab*        pcn = 0,
			  MultiFab*              Rhs = 0);

    void diffuse_iter_CPL (Real                   dt,
			   int                    nc,
			   int                    ncomps,
			   Real                   be_cn_theta,
			   int                    dataComp = 0,
			   const MultiFab*        alpha = 0, 
			   const MultiFab* const* betanp1 = 0,
			   const MultiFab* const* betanp1_dp = 0,
			   MultiFab*              pcnp1 = 0,
			   MultiFab*              S_nwt = 0,
			   Real*                  err_nwt = 0);

#ifdef MG_USE_FBOXLIB
    void richard_iter (Real                   dt,
		       int                    nc,
		       Real                   gravity,
		       Array<Real>            density,
		       MultiFab&              res_fix,
		       const MultiFab*        alpha,
		       const MultiFab* const* beta,
		       const MultiFab* const* beta_dp,
		       MultiFab*              umac,
		       const bool             do_upwind,
		       Real*                  err_nwt); 

    void richard_composite_iter (Real                      dt,
				 int                       nlevs,
				 int                       nc,
				 Real                      gravity,
				 Array<Real>               density,
				 PArray<MultiFab>&         res_fix,
				 PArray<MultiFab>&         alpha, 
				 Array<PArray<MultiFab> >& beta,
				 Array<PArray<MultiFab> >& beta_dp,
				 const bool                do_upwind,
				 Real*                     err_nwt);

    void richard_iter_p(Real                   dt,
			int                    nc,
			Real                   gravity,
			Array<Real>            density,
			MultiFab&              res_fix,
			const MultiFab*        alpha,
			const MultiFab*        dalpha,
			const MultiFab* const* beta,
			const MultiFab* const* beta_dp,
			MultiFab*              umac,
			const bool             do_upwind,
			Real*                  err_nwt); 

    void richard_composite_iter_p (Real                      dt,
				   int                       nlevs,
				   int                       nc,
				   Real                      gravity,
				   Array<Real>               density,
				   PArray<MultiFab>&         res_fix,
				   PArray<MultiFab>&         alpha, 
				   PArray<MultiFab>&         dalpha, 
				   Array<PArray<MultiFab> >& beta,
				   Array<PArray<MultiFab> >& beta_dp,
				   const bool                do_upwind,
				   Real*                     err_nwt);
 
    void richard_iter_eqb (int                    nc,
			   Real                   gravity,
			   Array<Real>            density,
			   MultiFab&              res_fix,
			   const MultiFab* const* beta,
			   const MultiFab* const* beta_dp,
			   MultiFab*              umac,
			   const bool             do_upwind,
			   Real*                  err_nwt);    
  
    void jac_richard (int                    nc,
		      Real                   gravity,
		      Array<Real>            density,
		      MultiFab&              res_fix,
		      MultiFab&              soln,
		      MultiFab&              soln_old, 
		      ABecLaplacian*         visc_op,
		      ABecLaplacian*         visc_op_dp = 0);

    void richard_flux( int                    nc,
		       Real                   be_cn_theta,
		       Real                   gravity,
		       Array<Real>&           density,
		       MultiFab* const*       flux,
		       const MultiFab*        pc,
		       const MultiFab* const* beta);
  
    void residual_richard (ABecLaplacian*         visc_op,
			   Real                   dt,
			   Real                   gravity,
			   Array<Real>            density,
			   MultiFab&              res,
			   const MultiFab*        pc,
			   const MultiFab* const* beta,
			   const MultiFab*        alpha = 0,
			   const MultiFab*        Stmp = 0 );

    void residual_richard (MGT_Solver&               mgt_solver,
			   Real                      dt,
			   Real                      gravity,
			   Array<Real>               density,
			   MultiFab**                Rhs,
			   PArray<MultiFab>&         pc,
			   Array<PArray<MultiFab> >& beta,
			   PArray<MultiFab>&         alpha,
			   PArray<MultiFab>&         res_fix,
			   MultiFab**                S,
			   ViscBndry                 visc_bndry);

#endif  

    void residual_CPL (ViscBndry             visc_bndry,
		       int                   nc,
		       Real                  dt,
		       Real                  be_cn_theta,	
		       MultiFab&             Rhs,
		       MultiFab*             S_nwt,
		       const MultiFab*       pcnp1,
		       const MultiFab*       alpha,		      
		       const MultiFab* const* betanp1,
		       MultiFab* const* fluxnp1 = 0);

    void diffuse_scalar_eqb (Real                   dt,
			     int                    sigma,
			     Real                   be_cn_theta,
			     MultiFab* const*       fluxn,
			     MultiFab* const*       fluxnp1,
			     int                    dataComp = 0,
			     MultiFab*              delta_rhs = 0, 
			     const MultiFab*        alpha = 0, 
			     const MultiFab* const* betan = 0, 
			     const MultiFab* const* betanp1 = 0,
			     const MultiFab* const* betanp1_dp =0,
			     MultiFab*              S_nwt = 0,
			     Real*                  err_nwt=0);

    void check_consistency  (Real                   dt,
			     int                    nc,
			     Real                   pm,
			     Real                   be_cn_theta,
			     int                    dataComp=0,
			     const MultiFab*        alpha=0, 
			     const MultiFab* const* betan=0, 
			     const MultiFab* const* betanp1=0,
			     const MultiFab*        pcn=0, 
			     const MultiFab*        pcnp1=0,
			     MultiFab*              S_nwt=0,
			     Real*                  err_nwt=0);
    
    void diffuse_Ssync (MultiFab*              Ssync,
                        int                    sigma,
                        Real                   dt,
                        Real                   be_cn_theta,
			MultiFab* const*       flux,
			int                    dataComp = 0,
                        const MultiFab* const* beta = 0,
			const MultiFab*        alpha = 0);

    void diffuse_Ssync_CPL (MultiFab*              Ssync,
			    int                    nc,
			    Real                   dt,
			    Real                   be_cn_theta,
			    MultiFab* const*       flux,
			    int                    dataComp = 0,
			    const MultiFab*        alpha = 0,
			    const MultiFab* const* betanp1=0);

    ABecLaplacian* getViscOp (int                    src_comp,
                              Real                   a,
                              Real                   b, 
                              Real                   time,
                              ViscBndry&             visc_bndry,
                              Real*                  rhsscale = 0,
			      int                    dataComp = 0,
                              const MultiFab* const* beta = 0,
                              const MultiFab*        alpha_in = 0,
                              bool		     bndry_already_filled = false);  

    ABecLaplacian* getViscOp (ViscBndry&             visc_bndry,
                              Real                   a,
                              Real                   b,
                              Real*                  rhsscale = 0,
			      int                    dataComp = 0,
                              const MultiFab* const* beta = 0,
                              const MultiFab*        alpha_in = 0,
			      bool                   do_volume = true);

    ABecLaplacian* getViscOp (int                    src_comp,
                              Real                   a,
                              Real                   b,
                              Real*                  rhsscale = 0,
			      int                    dataComp = 0,
                              const MultiFab* const* beta = 0,
                              const MultiFab*        alpha_in = 0);

#ifdef MG_USE_FBOXLIB
    MGT_Solver getOp (int                   comp,
		      int                   nlevs,
		      Array< Array<Real> >& xa,
		      Array< Array<Real> >& xb,
		      Real                  time,
		      ViscBndry&            visc_bndry,
		      const BCRec&          bc,
		      bool bndry_already_filled = true);

    MGT_Solver getOp (int                   comp,
		      int                   nlevs,
		      std::vector<Geometry> geom,
		      std::vector<BoxArray> bav,
		      std::vector<DistributionMapping> dmv,
		      Array< Array<Real> >& xa,
		      Array< Array<Real> >& xb,
		      Real                  time,
		      ViscBndry&            visc_bndry,
		      const BCRec&          bc,
		      bool	 bndry_already_filled = true);

    void coefs_fboxlib_mg (PArray<MultiFab>& acoefs,
			   Array< PArray<MultiFab> >& bcoefs,
			   const MultiFab* alpha_in, 
			   const MultiFab* const* beta,
			   const Real      a,
			   bool  do_volume = true);

    void coefs_fboxlib_mg (PArray<MultiFab>&          acoefs,
			   Array< PArray<MultiFab> >& bcoefs,
			   PArray<MultiFab>&          alpha,
			   Array< PArray<MultiFab> >& beta,
			   const Real                 a,
			   int                        nlevs);
#endif

    void setBeta (ABecLaplacian*         visc_op,
		  int                    dataComp,
		  const MultiFab* const* beta,
		  bool                   do_volume = true);
  
    void solve (ABecLaplacian*    visc_op,
		MultiFab&         Soln,
	        MultiFab&         Rhs,
		const BCRec&      bc, 
		Real              S_tol,
		Real              S_tol_abs);

    void compute_flux( int                    nc,
		       Real                   dt,
		       Real                   be_cn_theta,
		       MultiFab* const*       flux,
		       const MultiFab*        pc,
		       const MultiFab* const* beta);

    void getViscTerms (MultiFab&              visc_terms,
                       int                    src_comp,
                       int                    comp,
                       Real                   time,
		       int                    dataComp = 0,
                       const MultiFab* const* beta = 0);

    void getCplViscTerms (MultiFab&              visc_terms,
			  int                    src_comp,
			  Real                   time,
			  const Real*            density,
			  int                    dataComp = 0,
			  const MultiFab* const* beta = 0,
			  const MultiFab*        pc = 0);

    void getBndryData (ViscBndry&         bndry,
                       int                state_ind,
                       int                num_comp,
                       Real               time);

    void getBndryData (ViscBndry&          bndry,
                       int                 state_ind,
                       int                 num_comp,
                       AmrLevel::TimeLevel the_time,
                       int                 rho_flag);

    void FillBoundary (BndryRegister&     bdry,
                       int                src_comp,
                       int                dest_comp,
                       int                num_comp,
                       Real               time);
  
    void checkBetas (const MultiFab* const* beta1,
                     const MultiFab* const* beta2,
                     int&                   allthere,
                     int&                   allnull) const;

    void checkBeta (const MultiFab* const* beta,
                    int&                   allthere,
                    int&                   allnull) const;

    void checkBeta (const MultiFab* const* beta,
                    int&                   allthere) const;

    void allocFluxBoxesLevel (MultiFab**& fluxbox, 
                              int         nghost = 0,
                              int         nvar = 1);

    void removeFluxBoxesLevel (MultiFab**& fluxbox);

    int maxOrder () const;

    static bool are_any (const Array<DiffusionForm>& diffusionType,
                         const DiffusionForm         testForm,
                         const int                   sComp,
                         const int                   nComp);

    static int how_many (const Array<DiffusionForm>& diffusionType,
                         const DiffusionForm         testForm,
                         const int                   sComp,
                         const int                   nComp);
protected:

    static void Initialize ();
    static void Finalize ();
    //
    // Data Required by Derived Classes
    //
    Amr*            parent;
    AmrLevel*       caller;
    const BoxArray& grids;
    const int       level;
    //
    // Volume and area fractions.
    //
    const MultiFab& volume;
    const MultiFab* area;
    //
    // Static data.
    //
    static bool        use_mg_precond_flag;
    static int         use_cg_solve;
    static int         scale_abec;
    static Array<int>  is_diffusive;    // Does variable diffuse?
    static Array<Real> visc_coef;       // Const coef viscosity terms
    static int         verbose;
    static Real        visc_tol;
    static Real        visc_abs_tol;
private:
    //
    // The data.
    //
    Diffusion*    coarser;
    Diffusion*    finer;
    int           NUM_SCALARS;
    IntVect       crse_ratio;
    FluxRegister* viscflux_reg;

    //
    // Density at diffusion solve.
    //
    MultiFab* rho;
    //
    // Static data.
    //
    static int do_reflux;
    static int max_order;
};

#endif

