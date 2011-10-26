#ifndef _MacProj_H_
#define _MacProj_H_
//
// $Id: MacProj.H,v 1.33 2011-08-28 05:54:53 gpau Exp $
//
#include <BC_TYPES.H>
#include <BC_TYPES_PM.H>
#include <BCRec.H>
#include <Amr.H>
#include <AmrLevel.H>
#include <FluxRegister.H>
#include <PressBndry.H>
#include "Godunov.H"

class MacProj
{
public:
    
  MacProj (Amr*   parent,
	   int    _finest_level,
	   BCRec* _phys_bc, 
	   int    _do_any_diffuse);
    
  virtual ~MacProj ();

  void install_level (int                   level, 
		      AmrLevel*             level_data,
		      MultiFab&             volume,
		      MultiFab*             area);
  //
  // The level advance mac projection.
  //
  void mac_project (int              level,
		    MultiFab*        u_mac,
		    MultiFab*        mac_coef,
		    MultiFab*        RhoD,
		    MultiFab*        forces,
		    MultiFab*        pressure,
		    const PressBndry&  mac_bndry,
		    const BCRec&     p_bc);

  void print_min_max (MultiFab* u_mac);

  void contribute_to_mac_reg (int       level,
			      MultiFab* u_mac);
  void contribute_to_mac_reg_rhoD (int       level,
				   MultiFab* rhod);
  //
  // The sync solve.
  //
  void mac_sync_solve (int       level,
		       const BCRec&    p_bc,
		       MultiFab* mac_coef,
		       MultiFab* mac_sync_phi,
		       MultiFab* u_mac_sync,
		       IntVect&  fine_ratio);
  //
  // Computing the sync tendency.
  //
  void mac_sync_compute (int                   level,
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
			 Real                  be_cn_theta);

  void set_dirichlet_bcs (int                level,
			  MultiFab*          mac_phi,
			  const MultiFab&    RhoG,
			  const BCRec&       p_bc,
			  const Array<Real>& press_lo,
			  const Array<Real>& press_hi);

  void computeRhoG(FArrayBox* RhoG, FArrayBox* phi,
		   const Geometry& geom, Orientation* outFaces,
		   int numOutFlowFaces, Real gravity,
		   const int* lo_bc, const int* hi_bc,
		   const Real* press_lo, const Real* press_hi);

  static void Initialize ();
  static void Finalize ();

  static void GetOutFlowFaces (bool&        haveOutFlow,
			       Orientation* outFaces,
			       BCRec*       _phys_bc,
			       const BCRec& p_bc,
			       int&         numOutFlowFaces);

  static void GetInFlowFaces  (bool&        haveInFlow,
			       Orientation* InFaces,
			       BCRec*       _phys_bc,
			       const BCRec& p_bc,
			       int&         numInFlowFaces);

  static void GetDirichletFaces  (bool&        haveInFlow,
				  Orientation* InFaces,
				  const BCRec& p_bc,
				  int&         numInFlowFaces);
  //
  // Test the divergence constraint.
  //
  void check_div_cond (int       level,
		       MultiFab& divu,
		       MultiFab  U_edge[],
		       MultiFab  RhoD[]);

  void check_div_cond (int      level,
		       MultiFab U_edge[],
		       MultiFab RhoD[]) const;

protected:
  
  //
  // Boundary conditions.
  //
  void BuildPhiBC (int level);

  //
  // Pointers to amr,amrlevel.
  //
  Amr*             parent;
  PArray<AmrLevel> LevelData;
  //
  // Boundary condition objects.
  //
  BCRec*                phys_bc;
  Array< Array<BCRec> > phi_bcs;
  //
  // MAC sync correction and solution.
  //
  PArray<FluxRegister> mac_reg;
  PArray<FluxRegister> mac_reg_rhoD;
  //
  // Volume and area fractions.
  //
  PArray<MultiFab> volume;
  Array<MultiFab*> area;
  //
  // Parameters.
  //
  int        finest_level;
  int        finest_level_allocated;

  static int  verbose;
  static int  do_any_diffuse; 
  static bool use_cg_solve;
  static Real mac_tol;
  static Real mac_abs_tol;
  static Real mac_sync_tol;
};

#endif
