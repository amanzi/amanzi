#ifndef _GODUNOV_CLASS_H_
#define _GODUNOV_CLASS_H_ 

//
// A useful enumeration of the forms for the advective terms
//
enum AdvectionForm {Conservative = 0, NonConservative};

//===========================================================
// This header file defines the Godunov driver class.  Its
// interface simplifies the calling sequence for the Godunov
// box
//===========================================================

class Godunov {
public:      

  // =============================================================
  // Setup functions follow
  // =============================================================

  // construct the advection object
  Godunov();
    
  // size 1d work space explicitly
  Godunov( int max_size );

  static void Initialize();
  static void Finalize();
    
  // destruct the advection object
  ~Godunov();
 
  // set up the Godunov Box.  Assign work space.
  virtual void Setup( const Box &grd, 
		      FArrayBox &xflux, 
		      FArrayBox &yflux,
#if (BL_SPACEDIM == 3 )
		      FArrayBox &zflux, 
#endif
		      const int nscal,
		      const int model);    

  void Setup_tracer(  const Box &grd,
		      FArrayBox &xflux,
		      FArrayBox &yflux,
#if (BL_SPACEDIM == 3 )
		      FArrayBox &zflux,
#endif
		      const int nscal); 

  // =============================================================
  // Data access functions
  // ==============================================================

  int useForcesInTrans() {return use_forces_in_trans;}
    
  // =============================================================
  // Advection functions follow
  // =============================================================

  // compute the edge states
  void edge_states( const Box &grd, const Real *dx, Real dt, int velpred,
		    FArrayBox &uedge, FArrayBox &stx,
		    FArrayBox &vedge, FArrayBox &sty,
#if (BL_SPACEDIM == 3 )
		    FArrayBox &wedge, FArrayBox &stz,
#endif
		    FArrayBox &S, 
		    FArrayBox &tforces, FArrayBox& divu,
		    int fab_ind, int state_ind, const int *bc,
		    int use_conserv);

  // compute the edge states for single-phase flow (linear advection)
  void edge_states_lin(const Box &grd, const Real *dx, Real dt, 
		       FArrayBox &uedge, FArrayBox &stx, 
		       FArrayBox &vedge, FArrayBox &sty, 
#if (BL_SPACEDIM == 3 )
		       FArrayBox &wedge, FArrayBox &stz, 
#endif
		       FArrayBox &S,  FArrayBox &S_new,
		       FArrayBox &tforces, FArrayBox& rockphi,
		       int fab_ind, int state_ind, const int *bc, int nscal);

  // compute the edge states involving riemann solver
  void edge_states_rmn(const Box &grd, const Real *dx, Real dt, 
		       FArrayBox &uedge, FArrayBox &stx, FArrayBox &kappax,
		       FArrayBox &vedge, FArrayBox &sty, FArrayBox &kappay,
#if (BL_SPACEDIM == 3 )
		       FArrayBox &wedge, FArrayBox &stz, FArrayBox &kappaz,
#endif
		       FArrayBox& S,  FArrayBox& S_new,
		       FArrayBox& tforces, FArrayBox& divu,
		       FArrayBox& rockphi, FArrayBox& kappa,
		       FArrayBox& lambda_cc, FArrayBox& dlambda_cc,
		       FArrayBox&  kr_coef, const int n_kr_coef,
		       int fab_ind, int state_ind, const int *bc,
		       int use_conserv, int nscal);

  void edge_states_cpl(const Box &grd, const Real *dx, Real dt, 
		       FArrayBox &uedge, FArrayBox &stx, FArrayBox &kappax, FArrayBox &lambdax,
		       FArrayBox &vedge, FArrayBox &sty, FArrayBox &kappay, FArrayBox &lambday,
#if (BL_SPACEDIM == 3 )
		       FArrayBox &wedge, FArrayBox &stz, FArrayBox &kappaz, FArrayBox &lambdaz,
#endif
		       FArrayBox& S,  FArrayBox& S_new,
		       FArrayBox& tforces, FArrayBox& divu,
		       FArrayBox& rockphi, FArrayBox& kappa, FArrayBox& pc,
		       FArrayBox& lambda_cc, FArrayBox& dlambda_cc,
		       FArrayBox&  kr_coef, const int n_kr_coef,
		       int fab_ind, int state_ind, const int *bc,
		       int use_conserv, int nscal);

  void edge_states_pmr(const Box &grd, const Real *dx, Real dt, 
		       FArrayBox &uedge, FArrayBox &stx, FArrayBox &kappax,
		       FArrayBox &vedge, FArrayBox &sty, FArrayBox &kappay,
#if (BL_SPACEDIM == 3 )
		       FArrayBox &wedge, FArrayBox &stz, FArrayBox &kappaz,
#endif
		       FArrayBox &S,  FArrayBox &S_new,
		       FArrayBox &tforces, 
		       FArrayBox& divu,
		       FArrayBox& rockphi,
		       FArrayBox& kappa,
		       int fab_ind, int state_ind, const int *bc,
		       int use_conserv, int nscal, Real gravity, Real* eigmax);

  void edge_states_tracer (const Box& grd, const Real* dx, Real dt,
			   FArrayBox&  uedge, FArrayBox&  stx,
			   FArrayBox&  xlo, FArrayBox&  xhi,
			   FArrayBox&  vedge, FArrayBox&  sty,
			   FArrayBox&  ylo, FArrayBox&  yhi,
#if (BL_SPACEDIM == 3)               
			   FArrayBox&  wedge, FArrayBox&  stz,
			   FArrayBox&  zlo, FArrayBox&  zhi,
#endif
			   FArrayBox&  S, FArrayBox&  S_new,
			   FArrayBox&  St, FArrayBox&  St_new,
			   FArrayBox&  rock_phi, const int*  bc, int nscal);

  void edge_sync_rmn(const Box &grd, const Real *dx, Real dt, 
		     FArrayBox &uedge, FArrayBox &stx, FArrayBox &kappax,
		     FArrayBox &vedge, FArrayBox &sty, FArrayBox &kappay,
#if (BL_SPACEDIM == 3 )
		     FArrayBox &wedge, FArrayBox &stz, FArrayBox &kappaz,
#endif
		     FArrayBox &S,  FArrayBox &S_new,
		     FArrayBox &tforces, 
		     FArrayBox& divu,
		     FArrayBox& rockphi,
		     FArrayBox& kappa,
		     FArrayBox& lambda_cc, FArrayBox& dlambda_cc,
		     FArrayBox&  kr_coef, const int n_kr_coef,  
		     int fab_ind, int state_ind, const int *bc,
		     int use_conserv, int nscal);
  
  // advect a state component
  virtual void AdvectState( const Box &grd, const Real *dx, Real dt, 
			    FArrayBox &areax, FArrayBox &uedge, FArrayBox &xflux,  
			    FArrayBox &areay, FArrayBox &vedge, FArrayBox &yflux,  
#if (BL_SPACEDIM == 3 )                               
			    FArrayBox &areaz, FArrayBox &wedge, FArrayBox &zflux,
#endif
			    FArrayBox &S, 
			    FArrayBox &tforces, FArrayBox& divu,
			    int fab_ind,
			    FArrayBox &aofs,                  
			    int aofs_ind,
			    int iconserv, int state_ind, const int *bc,
			    FArrayBox &vol );

  /* 
     Linear advection: Valid for single-phase flow (model = 0, 1)
  */
  virtual void AdvectStateLin(const Box &grd, const Real *dx, Real dt, 
			      FArrayBox &areax, FArrayBox &uedge, 
			      FArrayBox &xflux,
			      FArrayBox &areay, FArrayBox &vedge, 
			      FArrayBox &yflux, 
#if (BL_SPACEDIM == 3 )                               
			      FArrayBox &areaz, FArrayBox &wedge, 
			      FArrayBox &zflux, 
#endif
			      FArrayBox &S, FArrayBox &S_new, FArrayBox &tforces, 
			      int fab_ind,
			      FArrayBox &aofs, int aofs_ind,
			      FArrayBox &rock_phi, int state_ind, const int *bc,
			      FArrayBox &vol, int nscal);

  /* 
     Nonlinear advection: Valid for two-phase flow (model = 2)
  */
  virtual void AdvectStateRmn(const Box &grd, const Real *dx, Real dt, 
			      FArrayBox &areax, FArrayBox &uedge, 
			      FArrayBox &xflux, FArrayBox &kappax, 
			      FArrayBox &areay, FArrayBox &vedge, 
			      FArrayBox &yflux, FArrayBox &kappay,  
#if (BL_SPACEDIM == 3 )                               
			      FArrayBox &areaz, FArrayBox &wedge, 
			      FArrayBox &zflux, FArrayBox &kappaz,
#endif
			      FArrayBox &S, FArrayBox &S_new, FArrayBox &tforces, 
			      FArrayBox &divu, int fab_ind,
			      FArrayBox &aofs, int aofs_ind,
			      FArrayBox &rock_phi, FArrayBox &kappa,
			      FArrayBox &lambda_cc, FArrayBox &dlambda_cc,
			      FArrayBox& kr_coef, const int n_kr_coef,
			      int iconserv, int state_ind, const int *bc,
			      FArrayBox &vol, int nscal);
  /* 
     Nonlinear advection:  Valid for two-phase flow with capillary (model = 2).  
  */
  virtual void AdvectStateCpl(const Box &grd, const Real *dx, Real dt, 
			      FArrayBox &areax, FArrayBox &uedge, 
			      FArrayBox &xflux, FArrayBox &kappax, FArrayBox &lambdax, 
			      FArrayBox &areay, FArrayBox &vedge, 
			      FArrayBox &yflux, FArrayBox &kappay, FArrayBox &lambday, 
#if (BL_SPACEDIM == 3 )                               
			      FArrayBox &areaz, FArrayBox &wedge, 
			      FArrayBox &zflux, FArrayBox &kappaz, FArrayBox &lambdaz, 
#endif
			      FArrayBox &S, FArrayBox &S_new, FArrayBox &tforces, 
			      FArrayBox &divu, int fab_ind,
			      FArrayBox &aofs, int aofs_ind,
			      FArrayBox &rock_phi, FArrayBox &kappa, FArrayBox &pc,
			      FArrayBox &lambda_cc, FArrayBox &dlambda_cc,
			      FArrayBox& kr_coef, const int n_kr_coef,
			      int iconserv, int state_ind, const int *bc,
			      FArrayBox &vol, int nscal);

  virtual void AdvectStatePmr(const Box &grd, const Real *dx, Real dt, 
			      FArrayBox &areax, FArrayBox &uedge, 
			      FArrayBox &xflux, FArrayBox &kappax, 
			      FArrayBox &areay, FArrayBox &vedge, 
			      FArrayBox &yflux, FArrayBox &kappay,  
#if (BL_SPACEDIM == 3 )                               
			      FArrayBox &areaz, FArrayBox &wedge, 
			      FArrayBox &zflux, FArrayBox &kappaz,
#endif
			      FArrayBox &S, FArrayBox &S_new, FArrayBox &tforces, 
			      FArrayBox &divu, int fab_ind,
			      FArrayBox &aofs, int aofs_ind,
			      FArrayBox &rock_phi,
			      FArrayBox &kappa,
			      int iconserv, int state_ind, const int *bc,
			      FArrayBox &vol, int nscal, Real gravity, Real* eigmax);    

  virtual void AdvectTracer (const Box& grd, const Real* dx, Real dt, 
			     FArrayBox& areax, FArrayBox& uedge, FArrayBox& xflux,
			     FArrayBox& areay, FArrayBox& vedge, FArrayBox& yflux, 
#if (BL_SPACEDIM == 3)                               
			     FArrayBox& areaz, FArrayBox& wedge, FArrayBox&  zflux,
#endif
			     FArrayBox& S, FArrayBox& S_new,
			     FArrayBox& St, FArrayBox& St_new,
			     FArrayBox&  tforces, 
			     FArrayBox&  divu, int fab_ind, 
			     FArrayBox&  aofs, int aofs_ind,
			     FArrayBox&  rock_phi, int iconserv,
			     int         state_ind, const int*  bc,
			     FArrayBox&  vol, int nscal);

  virtual void AdvectSyncRmn(const Box &grd, const Real *dx, Real dt, 
			     FArrayBox &areax, FArrayBox &uedge, FArrayBox &ucorr, 
			     FArrayBox &xflux, FArrayBox &kappax, 
			     FArrayBox &areay, FArrayBox &vedge, FArrayBox &vcorr, 
			     FArrayBox &yflux, FArrayBox &kappay,  
#if (BL_SPACEDIM == 3 )                               
			     FArrayBox &areaz, FArrayBox &wedge, FArrayBox &wcorr, 
			     FArrayBox &zflux, FArrayBox &kappaz,
#endif
			     FArrayBox &S, FArrayBox &S_new, FArrayBox &tforces, 
			     FArrayBox &divu, int fab_ind,
			     FArrayBox &aofs, int aofs_ind,
			     FArrayBox &rock_phi, FArrayBox &kappa,
			     FArrayBox& lambda_cc, FArrayBox& dlambda_cc,
			     FArrayBox&  kr_coef, const int   n_kr_coef,
			     int iconserv, int state_ind, const int *bc,
			     FArrayBox &vol, int nscal);

  //Compute the gradient
  void Getdfdn(FArrayBox& S, int state_ind, int ncomps, int dir, int idx);
  void Getsandc(FArrayBox& S, int state_ind, int ncomps, int idx);

  // compute the advective derivative from fluxes
  void ComputeAofs( const Box &grd, 
		    FArrayBox &areax, FArrayBox &uedge, FArrayBox &xflux,  
		    FArrayBox &areay, FArrayBox &vedge, FArrayBox &yflux,  
#if (BL_SPACEDIM == 3 )
		    FArrayBox &areaz, FArrayBox &wedge, FArrayBox &zflux,
#endif
		    FArrayBox &vol,
		    FArrayBox &aofs,  int aofs_ind, int iconserv );

  // compute the advective derivative from fluxes
  void ComputeAofsRmn(const Box &grd, 
		      FArrayBox &areax, FArrayBox &uedge, FArrayBox &xflux,  
		      FArrayBox &areay, FArrayBox &vedge, FArrayBox &yflux,  
#if (BL_SPACEDIM == 3 )
		      FArrayBox &areaz, FArrayBox &wedge, FArrayBox &zflux,
#endif
		      FArrayBox &vol,FArrayBox &aofs, int aofs_ind, int nscal);

  // sync advect a state component
  virtual void SyncAdvect( const Box &grd, const Real *dx, Real dt, int level,
			   const FArrayBox &areax, FArrayBox &uedge,
			   FArrayBox &ucorr, FArrayBox &xflux,
			   const FArrayBox &areay, FArrayBox &vedge,
			   FArrayBox &vcorr, FArrayBox &yflux,
#if (BL_SPACEDIM == 3 )
			   const FArrayBox &areaz, FArrayBox &wedge,
			   FArrayBox &wcorr, FArrayBox &zflux,
#endif
			   FArrayBox &S,
			   FArrayBox &tforces, FArrayBox& divu, 
			   FArrayBox &rock_phi,
			   int ind,
			   FArrayBox &sync, int sync_ind,
			   int iconserv, int state_ind, const int *bc,
			   const FArrayBox &vol );
    
  // compute the advective derivative of corrective fluxes for the mac sync
  void ComputeSyncAofs( const Box &grd,
			const FArrayBox &areax, FArrayBox &ucorr, FArrayBox &xflux,  
			const FArrayBox &areay, FArrayBox &vcorr, FArrayBox &yflux,  
#if (BL_SPACEDIM == 3 )
			const FArrayBox &areaz, FArrayBox &wcorr, FArrayBox &zflux,
#endif                     
			const FArrayBox &vol, FArrayBox &sync,
			int sync_ind, int iconserv );

  void ComputeSyncAofsRmn( const Box &grd,
			   const FArrayBox &areax, FArrayBox &ucorr, FArrayBox &xflux,  
			   const FArrayBox &areay, FArrayBox &vcorr, FArrayBox &yflux,  
#if (BL_SPACEDIM == 3 )
			   const FArrayBox &areaz, FArrayBox &wcorr, FArrayBox &zflux,
#endif                     
			   const FArrayBox &vol, 
			   const FArrayBox& kr_coef, const int n_kr_coef,
			   FArrayBox &sync, int sync_ind, int nscal);
    
  // correct a conservatively-advected scalar for under-over shoots
  void ConservativeScalMinMax( FArrayBox &Sold, FArrayBox &Snew, 
			       int ind_old_s, int ind_old_rho,
			       const int *bc, const Box &grd );

  // correct a convectively-advected scalar for under-over shoots
  void ConvectiveScalMinMax( FArrayBox &Sold, FArrayBox &Snew, 
			     int ind_old, int ind_new,
			     const int *bc, const Box &grd );

  // =============================================================
  // Diagnostic functions follow
  // =============================================================
  void esteig (const Box&  grd,
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
	       Real*       eigmax);

  void esteig_lin (const Box&  grd,
		   FArrayBox&  uedge,
		   FArrayBox&  vedge,
#if (BL_SPACEDIM == 3)               
		   FArrayBox&  wedge,
#endif
		   FArrayBox&  rock_phi,
		   Real*       eigmax);

  void esteig_trc (const Box&  grd,
		   FArrayBox&  uedge,
		   FArrayBox&  vedge,
#if (BL_SPACEDIM == 3)               
		   FArrayBox&  wedge,
#endif
		   FArrayBox&  S,
		   const int   nc,
		   FArrayBox&  rock_phi,
		   Real*       eigmax);

  void esteig_cpl(const Box&  grd,
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
		  Real*       eigmax);

  static bool are_any(const Array<AdvectionForm>& advectionType,
		      const AdvectionForm         testForm,
		      const int                   sComp,
		      const int                   nComp);

  static int how_many(const Array<AdvectionForm>& advectionType,
		      const AdvectionForm         testForm,
		      const int                   sComp,
		      const int                   nComp);


  // =============================================================
  // Source term functions follow
  // =============================================================
    

  // compute the update rule
  // psi^n+1 = psi^n - dt*aofs + dt*tforces
  void Add_aofs_tf( FArrayBox &Sold,
		    FArrayBox &Snew,    int start_ind, int num_comp,
		    FArrayBox &Aofs,    int aofs_ind,
		    FArrayBox &tforces, int tf_ind,
		    FArrayBox &Rockphi,
		    const Box &grd,     Real dt );

  void Add_aofs_tracer( FArrayBox &Sold,
			FArrayBox &Snew,    int start_ind, int num_comp,
			FArrayBox &Aofs,    int aofs_ind,
			FArrayBox &tforces, int tf_ind,
			FArrayBox &Rockphi,
			const Box &grd, const Array<int>& idx_total, Real dt );

  // iconserv==1   => tforces = tforces + visc - divU*S
  // iconserv==0   => tforces = (tforces+visc)/rho
  void Sum_tf_divu_visc( FArrayBox &S, FArrayBox &tforces,
			 int s_ind, int num_comp,
			 FArrayBox &visc,    int v_ind,
			 FArrayBox &divu, int iconserv );
    
  //================ for debugging purposes this is made public
  // internal FABS
  FArrayBox utmp;
  FArrayBox utmpn;
  FArrayBox work; // work space for the Godunov Box
  FArrayBox work_rmn; // work space for the Godunov Box with riemann solver
  FArrayBox xlo, ylo, xhi, yhi;
#if (BL_SPACEDIM == 3)
  FArrayBox zlo,zhi;
#endif

  // internal parameters
  int max_1d;    // largest dimension size allowed
  int scr_size;  // actual length of the scratch space
    
  // 1D arrays used in computing slopes and edges states
  Real *stxlo; 
  Real *stxhi; 
  Real *slxscr;
  Real *stylo;
  Real *styhi;
  Real *slyscr;
#if (BL_SPACEDIM == 3 )    
  Real *stzlo;
  Real *stzhi;
  Real *slzscr;
#endif
    
protected:   // protected work functions
    
  // manipulate 1d scratch space
  void ZeroScratch();
  void SetBogusScratch();
  void SetScratch( int max_size );
  void RemScratch();

protected:   // protected data members

  // sizing boxes
  Box work_bx;
  Box xflux_bx;
  Box yflux_bx;
#if (BL_SPACEDIM == 3 )    
  Box zflux_bx;
#endif

  // flags

  static int verbose;
  static int use_forces_in_trans;

  // algorithm control parameters

  static int slope_order;

};

#endif

