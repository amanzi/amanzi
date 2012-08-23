#include <winstd.H>

#include <algorithm>
#include <vector>
#include <cmath>

#include <ErrorList.H>
#include <Interpolater.H>
#include <MultiGrid.H>
#include <ArrayLim.H>
#include <Profiler.H>
#include <TagBox.H>
#include <DataServices.H>
#include <AmrData.H>
#include <time.h> 
#include <PMAmr.H>

#include <Godunov.H>
#include <PorousMedia.H>
#include <PROB_PM_F.H>
#include <POROUS_F.H>
#include <POROUSMEDIA_F.H>
#include <VISCOPERATOR_F.H>

#ifdef _OPENMP
#include "omp.h"
#endif

#ifdef AMANZI
#include "exceptions.hh"
#include "errors.hh"
#include "simple_thermo_database.hh"
#include "chemistry_verbosity.hh"
#include "chemistry_exception.hh"
#include "chemistry_output.hh"
extern amanzi::chemistry::ChemistryOutput* amanzi::chemistry::chem_out;
#endif

#define GEOM_GROW   1
#define HYP_GROW    3
#define PRESS_GROW  1

#define SHOWVALARR(val)                        \
{                                              \
    std::cout << #val << " = ";                \
    for (int i=0;i<val.size();++i)             \
    {                                          \
        std::cout << val[i] << " " ;           \
    }                                          \
    std::cout << std::endl;                    \
}                                             
#define SHOWVALARRA(val) { SHOWVALARR(val); BoxLib::Abort();}
#define SHOWVAL(val) { std::cout << #val << " = " << val << std::endl;}
#define SHOWVALA(val) { SHOWVAL(val); BoxLib::Abort();}

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)	\
  const int* fablo = (fab).loVect();		\
  const int* fabhi = (fab).hiVect();		\
  Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)	\
  const int* fablo = (fab).loVect();		\
  const int* fabhi = (fab).hiVect();		\
  const Real* fabdat = (fab).dataPtr();

//
// Static objects.
//
ErrorList PorousMedia::err_list;
BCRec     PorousMedia::phys_bc;
BCRec     PorousMedia::pres_bc;
MacProj*  PorousMedia::mac_projector = 0;
Godunov*  PorousMedia::godunov       = 0;
static double richard_time;
static double richard_time_min = 1.e6;

PM_Error_Value::PM_Error_Value (Real min_time, Real max_time, int max_level, 
                                const PArray<Region>& regions)
    : pmef(0), value(0), min_time(min_time), max_time(max_time), max_level(max_level)
{
    set_regions(regions);
}

PM_Error_Value::PM_Error_Value (PMEF pmef,
                                Real value, Real min_time,
                                Real max_time, int max_level, 
                                const PArray<Region>& regions)
    : pmef(pmef), value(value), min_time(min_time), max_time(max_time), max_level(max_level)
{
    set_regions(regions);
}

void
PM_Error_Value::set_regions(const PArray<Region>& regions_)
{
    regions.clear();
    int nregions=regions_.size();

    // Get a copy of the pointers to regions in a structure that wont 
    //   remove them when it leaves scope
    regions.resize(nregions,PArrayNoManage);
    for (int i=0; i<nregions; ++i)
    {
        Region& r = const_cast<Region&>(regions_[i]);
        regions.set(i,&(r));
    }
}

void
PM_Error_Value::tagCells(int* tag, D_DECL(const int& tlo0,const int& tlo1,const int& tlo2),
                         D_DECL(const int& thi0,const int& thi1,const int& thi2),
                         const int* tagval, const int* clearval,
                         const Real* data, 
                         D_DECL(const int& dlo0,const int& dlo1,const int& dlo2),
                         D_DECL(const int& dhi0,const int& dhi1,const int& dhi2),
                         const Real* mask, 
                         D_DECL(const int& mlo0,const int& mlo1,const int& mlo2),
                         D_DECL(const int& mhi0,const int& mhi1,const int& mhi2),
                         const int* lo, const int* hi, const int* nvar,
                         const int* domain_lo, const int* domain_hi,
                         const Real* dx, const Real* xlo,
                         const Real* prob_lo, const Real* time,
                         const int* level) const
{
    BL_ASSERT(pmef);

    pmef(tag,D_DECL(tlo0,tlo1,tlo2),D_DECL(thi0,thi1,thi2),tagval,clearval,
         data,D_DECL(dlo0,dlo1,dlo2),D_DECL(dhi0,dhi1,dhi2),
         mask,D_DECL(mlo0,mlo1,mlo2),D_DECL(mhi0,mhi1,mhi2),
         lo, hi, nvar,domain_lo,domain_hi,dx,xlo,prob_lo,time,
         level,&value);
}

static Real BL_ONEATM = 101325.0;

namespace
{
  const std::string solid("Solid");
  const std::string absorbed("Absorbed");
  const std::string ctotal("Total");
}

static bool initialized = false;
void
PorousMedia::CleanupStatics ()
{
    regions.clear();
    rocks.clear();
    observations.clear();
    ic_array.clear();
    bc_array.clear();
    tic_array.clear();
    tbc_array.clear();
    initialized = false;
#ifdef AMANZI
    delete amanzi::chemistry::chem_out;
#endif
}

void
PorousMedia::variableCleanUp ()
{
  desc_lst.clear();
  derive_lst.clear();
  err_list.clear();

  delete kappadata;
  kappadata = 0;
  delete phidata;
  phidata = 0;

  delete mac_projector;
  mac_projector = 0;

  delete godunov;
  godunov = 0;

  model_list.clear();
  phase_list.clear();
  comp_list.clear();
  tracer_list.clear();

  source_array.clear();

#ifdef AMANZI
  if (do_chem>0)
    {
      chemSolve.clear();
      components.clear();
      parameters.clear();
    }
#endif
}

PorousMedia::PorousMedia ()
{
    if (!initialized) {
        BoxLib::ExecOnFinalize(PorousMedia::CleanupStatics);
        initialized = true;
    }

  Ssync        = 0;
  advflux_reg  = 0;
  viscflux_reg = 0;
  u_mac_prev   = 0;
  u_macG_prev  = 0;
  u_mac_curr   = 0;
  u_macG_curr  = 0;
  u_macG_trac  = 0;
  u_corr       = 0;
  kappa        = 0;
  kpedge       = 0;
  kr_coef      = 0;
  cpl_coef     = 0;
  lambda       = 0;
  lambda_cc    = 0;
  lambdap1_cc  = 0;
  dlambda_cc   = 0;
  rock_phi     = 0;
  aofs         = 0;
  diffusion    = 0;
  dt_eig       = 0;
  rhs_RhoD     = 0;
	
}

// FIXME: These are assigned in the translator also, should be same list of data
static std::string RpurposeDEF[7] = {"xlobc", "ylobc", "zlobc", "xhibc", "yhibc", "zhibc", "all"};


void
PorousMedia::setup_bound_desc()
{
    bc_descriptor_map.clear();
    const Real* dx   = geom.CellSize();
    const Box& domain = geom.Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Array<Orientation> Faces;
    const BCRec& bc = desc_lst[State_Type].getBC(0);
    getDirichletFaces(Faces,State_Type,bc);

    if (do_tracer_transport) {
        tbc_descriptor_map.resize(ntracers);
    }

    for (int iface = 0; iface < Faces.size(); iface++)
    {
        const Orientation& face = Faces[iface];
        if (PorousMedia::grids_on_side_of_domain(grids,domain,face)) 
        {
	    Box ccBndBox  = BoxLib::adjCell(domain,face,1);

            if (ccBndBox.ok()) {

                // Find BCs for this face
                int idx = face.coordDir() + 3*face.isHigh();
                const std::string& purpose = RpurposeDEF[idx];

                const PArray<RegionData>& bcs = PorousMedia::BCs();
                Array<int> myBCs;
                for (int i=0; i<bcs.size(); ++i) {
                    const PArray<Region>& regions = bcs[i].Regions();
                    int found = 0;
                    for (int j=0; j<regions.size(); ++j) {
                        if (regions[j].purpose == purpose) {
                            found++;
                        }
                    }

                    if (found) {
                        myBCs.push_back(i);
                    }
                }
                if (myBCs.size() > 0) 
                {
                    bc_descriptor_map[face] = BCDesc(ccBndBox,myBCs);
                }
                else {
                    std::cerr << "No BCs responsible for filling face: " << Faces[iface] << std::endl;
                    BoxLib::Abort();
                }

                if (do_tracer_transport) {
                    for (int n=0; n<ntracers; ++n) {
                        const PArray<RegionData>& tbcs = PorousMedia::TBCs(n);
                        Array<int> myTBCs;
                        for (int i=0; i<tbcs.size(); ++i) {
                            const PArray<Region>& tregions = tbcs[i].Regions();
                            int tfound = 0;
                            for (int j=0; j<tregions.size(); ++j) {
                                if (tregions[j].purpose == purpose) {
                                    tfound++;
                                }
                            }
                            
                            if (tfound) {
                                myTBCs.push_back(i);
                            }
                        }
                        if (myTBCs.size() > 0) 
                        {
                            tbc_descriptor_map[n][face] = BCDesc(ccBndBox,myTBCs);
                        }
                        else {
                            std::cerr << "No tracer BCs responsible for filling face: " << Faces[iface] << std::endl;
                            BoxLib::Abort();
                        }
                    }

                }

            }
        }
    }

    
    // setup boundary descriptor for the pressure
    pbc_descriptor_map.clear();
    Faces.clear();
    const BCRec& pbc = desc_lst[Press_Type].getBC(0);
    getDirichletFaces(Faces,Press_Type,pbc);

    for (int iface = 0; iface < Faces.size(); iface++)
    {
        const Orientation& face = Faces[iface];
        if (PorousMedia::grids_on_side_of_domain(grids,domain,face)) 
        {
            Box ccBndBox  = BoxLib::adjCell(domain,face,1);
            if (ccBndBox.ok()) {

                // Find BCs for this face
                int idx = face.coordDir() + 3*face.isHigh();
                const std::string& purpose = RpurposeDEF[idx];

                const PArray<RegionData>& bcs = PorousMedia::BCs();
                Array<int> myBCs;
                for (int i=0; i<bcs.size(); ++i) {
                    const PArray<Region>& regions = bcs[i].Regions();
                    int found = 0;
                    for (int j=0; j<regions.size(); ++j) {
                        if (regions[j].purpose == purpose) {
                            found++;
                        }
                    }

                    if (found) {
                        myBCs.push_back(i);
                    }
                }
                if (myBCs.size() > 0) 
                {
                    pbc_descriptor_map[face] = BCDesc(ccBndBox,myBCs);
                }
                else {
                    std::cerr << "No BCs responsible for filling face: " << Faces[iface] << std::endl;
                    BoxLib::Abort();
                }
            }
        }
    }
}

PorousMedia::PorousMedia (Amr&            papa,
                          int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          Real            time)
  :
  AmrLevel(papa,lev,level_geom,bl,time),
  //
  // Make room for ncomps+ntracers in aux_boundary_data_old.
  // With AMANZI we only use the ntracers parts.  But by using ncomps+ntracers
  // we don't need to worry about the case when ntracers==0.
  //
  aux_boundary_data_old(bl,HYP_GROW,ncomps+ntracers,level_geom),
  FillPatchedOldState_ok(true)
{
    if (!initialized) {
        BoxLib::ExecOnFinalize(CleanupStatics);
        initialized = true;
    }

  //
  // Build metric coefficients for RZ calculations.
  //
  buildMetrics();

  //
  // Set up reflux registers.
  //
  advflux_reg  = 0;
  viscflux_reg = 0;
  if (level > 0 && do_reflux)
    {
      advflux_reg  = new FluxRegister(grids,crse_ratio,level,NUM_SCALARS);
      viscflux_reg = new FluxRegister(grids,crse_ratio,level,NUM_SCALARS);
    }

  //
  // Initialize work multifabs.
  //
  Ssync        = 0;
  u_mac_prev   = 0;
  u_macG_prev  = 0;
  u_mac_curr   = 0;
  u_macG_curr  = 0;
  u_macG_trac  = 0;
  rhs_RhoD     = 0;
  u_corr       = 0;
  kappa        = 0;
  kpedge       = 0;
  kr_coef      = 0;
  cpl_coef     = 0;
  lambda       = 0;
  lambda_cc    = 0;
  lambdap1_cc  = 0;
  dlambda_cc   = 0;
  rock_phi     = 0;
  aofs         = 0;

  //
  // Set up the godunov box.
  //
  SetGodunov();
  //
  // Set up diffusion.
  //
  diffusion = new Diffusion(parent,this,
			    (level > 0) ? getLevel(level-1).diffusion : 0,
			    ndiff,viscflux_reg,volume,area,
			    is_diffusive,visc_coef);
  
  // Allocate space for variable diffusion coefficients
  diffn_cc   = 0;
  diffnp1_cc = 0;
  if (variable_scal_diff) 
    {
      diffn_cc   = new MultiFab(grids, ndiff, 1);
      diffnp1_cc = new MultiFab(grids, ndiff, 1);
    }

  // Allocate space for the capillary pressure diffusive term
  pcn_cc   = 0;
  pcnp1_cc = 0;
  if (have_capillary)
    {
      pcn_cc     = new MultiFab(grids, 1, 2);
      pcnp1_cc   = new MultiFab(grids, 1, 2);
      (*pcn_cc).setVal(0.);
      (*pcnp1_cc).setVal(0.);
    }

  //
  // Set up the mac projector.
  //
  if (mac_projector == 0)
    {
      mac_projector = new MacProj(parent,parent->finestLevel(),
				  &phys_bc,do_any_diffuse);
    }
  mac_projector->install_level(level,this,volume,area);

  //
  // Alloc MultiFab to hold rock quantities
  //
  BL_ASSERT(kappa == 0);
  kappa = new MultiFab(grids,1,3);

  BL_ASSERT(rock_phi == 0);
  rock_phi = new MultiFab(grids,1,3);

  if (model != model_list["single-phase"] &&
      model != model_list["single-phase-solid"] &&
      model != model_list["steady-saturated"])
    {
      BL_ASSERT(kr_coef == 0);
      kr_coef = new MultiFab(grids,5,1);
      (*kr_coef).setVal(0.);

      BL_ASSERT(cpl_coef == 0);
      cpl_coef = new MultiFab(grids,5,1);
      (*cpl_coef).setVal(0.);

      BL_ASSERT(lambda_cc == 0);
      lambda_cc = new MultiFab(grids,ncomps,1);
      (*lambda_cc).setVal(1.);
    
      BL_ASSERT(lambdap1_cc == 0);
      lambdap1_cc = new MultiFab(grids,ncomps,1);
      (*lambdap1_cc).setVal(1.);

      BL_ASSERT(dlambda_cc == 0);
      dlambda_cc = new MultiFab(grids,3,1);
      (*dlambda_cc).setVal(0.);
    }

  BL_ASSERT(lambda == 0);
  lambda = new MultiFab[BL_SPACEDIM];
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      BoxArray edge_grids(grids);
      edge_grids.surroundingNodes(dir);
      lambda[dir].define(edge_grids,1,0,Fab_allocate);
      lambda[dir].setVal(1.e40);
    }
  //
  // Alloc MultiFab to hold u_mac
  //
  BL_ASSERT(u_mac_prev  == 0);
  BL_ASSERT(u_mac_curr  == 0);
  BL_ASSERT(u_macG_trac == 0);
  BL_ASSERT(rhs_RhoD == 0);
  u_mac_prev  = new MultiFab[BL_SPACEDIM];
  u_mac_curr  = new MultiFab[BL_SPACEDIM];
  u_macG_trac = new MultiFab[BL_SPACEDIM];
  rhs_RhoD    = new MultiFab[BL_SPACEDIM];
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      BoxArray edge_grids(grids);
      edge_grids.surroundingNodes(dir);
      u_mac_prev[dir].define(edge_grids,1,0,Fab_allocate);
      u_mac_prev[dir].setVal(1.e40);
      u_mac_curr[dir].define(edge_grids,1,0,Fab_allocate);
      u_mac_curr[dir].setVal(1.e40);
      rhs_RhoD[dir].define(edge_grids,1,0,Fab_allocate);
      rhs_RhoD[dir].setVal(1.e40);
      edge_grids.grow(1);
      u_macG_trac[dir].define(edge_grids,1,0,Fab_allocate);
      u_macG_trac[dir].setVal(1.e40);	
    }
  BL_ASSERT(kpedge == 0);
  kpedge     = new MultiFab[BL_SPACEDIM];
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      BoxArray edge_gridskp(grids);
      edge_gridskp.surroundingNodes(dir).grow(1);
      kpedge[dir].define(edge_gridskp,1,0,Fab_allocate);
      kpedge[dir].setVal(1.e40);
    }

  // Must initialize to zero because we test on zero in estDt.
  dt_eig = 0;

  // Set up boundary condition work
  setup_bound_desc();
}

PorousMedia::~PorousMedia ()
{
  delete Ssync;
  delete advflux_reg;
  delete viscflux_reg;
  delete [] u_mac_prev;
  delete [] u_mac_curr;
  delete [] u_macG_prev;
  delete [] u_macG_curr;
  delete [] u_macG_trac;
  delete [] u_corr;

  u_macG_prev = 0;
  u_macG_curr = 0;
  u_macG_trac = 0;
  u_corr      = 0;
      
  delete [] rhs_RhoD;
  delete [] kpedge;
  delete [] lambda;
  delete kappa;
  delete rock_phi;

  if (kr_coef)
    delete kr_coef;
  if (cpl_coef)
    delete cpl_coef;
  if (lambda_cc)
    delete lambda_cc;
  if (lambdap1_cc)
    delete lambdap1_cc;
  if (dlambda_cc)
    delete dlambda_cc;

 
  // Remove the arrays for variable viscosity and diffusivity
  // and delete the Diffusion object
  if (variable_scal_diff)
    {
      delete diffn_cc;
      delete diffnp1_cc;
    }
  if (have_capillary)
    {
      delete pcn_cc;
      delete pcnp1_cc;
    }
  delete diffusion;
}

void
PorousMedia::allocOldData ()
{
  for (int k = 0; k < num_state_type; k++)
    {
      state[k].allocOldData();
    }
}

void
PorousMedia::removeOldData ()
{
  AmrLevel::removeOldData();
}

void
PorousMedia::SetGodunov()
{
  if (godunov == 0)
    godunov = new Godunov();
}

void
PorousMedia::restart (Amr&          papa,
                      std::istream& is,
                      bool          bReadSpecial)
{
  AmrLevel::restart(papa,is,bReadSpecial);
  is >> dt_eig;

  //int finest_level = parent->finestLevel();
  //for (int k = 0; k <= finest_level; k++)
  //{
  //    Real dt = parent->dtLevel()[k];
  //    Real strt_time =  static_cast<const PMAmr*>(parent)->StartTime();
  //    getLevel(k).setTimeLevel(strt_time,dt,dt);
  //}

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "Estimated time step = " << dt_eig << '\n';
  //
  // Make room for ncomps+ntracers in aux_boundary_data_old.
  // With AMANZI we only use the ntracers parts.  But by using ncomps+ntracers
  // we don't need to worry about the case when ntracers==0.
  //
  aux_boundary_data_old.initialize(grids,HYP_GROW,ncomps+ntracers,Geom());

  FillPatchedOldState_ok = true;

  set_overdetermined_boundary_cells(state[State_Type].curTime());
  //
  // Set the godunov box.
  //
  SetGodunov();
    
  if (mac_projector == 0)
    {
      mac_projector = new MacProj(parent,parent->finestLevel(),
				  &phys_bc,do_any_diffuse);
    }
  mac_projector->install_level(level,this,volume,area );

  //
  // Build metric coefficients for RZ calculations.
  //
  buildMetrics();

  BL_ASSERT(advflux_reg == 0);
  if (level > 0 && do_reflux)
    {
      advflux_reg = new FluxRegister(grids,crse_ratio,level,NUM_SCALARS);
    }
  BL_ASSERT(viscflux_reg == 0);
  if (level > 0 && do_reflux)
    {
      viscflux_reg = new FluxRegister(grids,crse_ratio,level,NUM_SCALARS);
    }

  BL_ASSERT(Ssync == 0);
  if (level < parent->finestLevel())
    Ssync = new MultiFab(grids,NUM_SCALARS,1);

  diffusion = new Diffusion(parent, this,
			    (level > 0) ? getLevel(level-1).diffusion : 0,
			    ndiff, viscflux_reg, volume, area,
			    is_diffusive, visc_coef);
  //
  // Allocate the storage for variable diffusivity
  //
  diffn_cc   = 0;
  diffnp1_cc = 0;    
  if (variable_scal_diff) 
    {
      diffn_cc   = new MultiFab(grids, ndiff, 1);
      diffnp1_cc = new MultiFab(grids, ndiff, 1);
    }
  //
  // Allocate the storage for capillary pressure
  //
  pcn_cc     = 0;
  pcnp1_cc   = 0;    
  if (have_capillary) 
    {
      pcn_cc     = new MultiFab(grids, 1, 2);
      pcnp1_cc   = new MultiFab(grids, 1, 2);
      (*pcn_cc).setVal(0.);
      (*pcnp1_cc).setVal(0.);
    }

  is_first_step_after_regrid = false;
  old_intersect_new          = grids;

  //
  // Alloc MultiFab to hold rock quantities
  //
  BL_ASSERT(kappa == 0);
  kappa = new MultiFab(grids,1,3); 

  BL_ASSERT(rock_phi == 0);
  rock_phi = new MultiFab(grids,1,3);

  if (model != model_list["single-phase"] &&
      model != model_list["single-phase-solid"] &&
      model != model_list["steady-saturated"])
    {
      BL_ASSERT(kr_coef == 0);
      kr_coef = new MultiFab(grids,5,1);
      (*kr_coef).setVal(0.);

      BL_ASSERT(cpl_coef == 0);
      cpl_coef = new MultiFab(grids,5,1);
      (*cpl_coef).setVal(0.);

      BL_ASSERT(lambda_cc == 0);
      lambda_cc = new MultiFab(grids,ncomps,1);
      (*lambda_cc).setVal(1.);
    
      BL_ASSERT(lambdap1_cc == 0);
      lambdap1_cc = new MultiFab(grids,ncomps,1);
      (*lambdap1_cc).setVal(1.);

      BL_ASSERT(dlambda_cc == 0);
      dlambda_cc = new MultiFab(grids,3,1);
      (*dlambda_cc).setVal(0.);
    }

  BL_ASSERT(lambda == 0);
  lambda = new MultiFab[BL_SPACEDIM];
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      BoxArray edge_gridskp(grids);
      edge_gridskp.surroundingNodes(dir);
      lambda[dir].define(edge_gridskp,1,0,Fab_allocate);
      lambda[dir].setVal(1.e40);
    }

  BL_ASSERT(kpedge == 0);
  kpedge     = new MultiFab[BL_SPACEDIM];
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      BoxArray edge_gridskp(grids);
      edge_gridskp.surroundingNodes(dir).grow(1);
      kpedge[dir].define(edge_gridskp,1,0,Fab_allocate);
      kpedge[dir].setVal(1.e40);
    }
  
  init_rock_properties();

  //
  // Alloc MultiFab to hold u_mac
  //
  u_mac_prev  = new MultiFab[BL_SPACEDIM];
  u_mac_curr  = new MultiFab[BL_SPACEDIM];
  u_macG_trac = new MultiFab[BL_SPACEDIM];
  rhs_RhoD    = new MultiFab[BL_SPACEDIM];
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      BoxArray edge_grids(grids);
      edge_grids.surroundingNodes(dir);
      u_mac_prev[dir].define(edge_grids,1,0,Fab_allocate);
    }

  std::string Level = BoxLib::Concatenate("Level_", level, 1);
  std::string FullPath = papa.theRestartFile();
  if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
    FullPath += '/';
  FullPath += Level;
        
  std::string uxfile = "/umac_x";
  std::string uyfile = "/umac_y";
  uxfile = FullPath + uxfile;
  uyfile = FullPath + uyfile;
  VisMF::Read(u_mac_curr[0],uxfile);
  VisMF::Read(u_mac_curr[1],uyfile);

#if (BL_SPACEDIM == 3)
  std::string uzfile = "/umac_z";
  uzfile = FullPath + uzfile;
  VisMF::Read(u_mac_curr[2],uzfile);
#endif

  std::string utxfile = "/umact_x";
  std::string utyfile = "/umact_y";
  utxfile = FullPath + utxfile;
  utyfile = FullPath + utyfile;
  VisMF::Read(u_macG_trac[0],utxfile);
  VisMF::Read(u_macG_trac[1],utyfile);

#if (BL_SPACEDIM == 3)
  std::string utzfile = "/umact_z";
  utzfile = FullPath + utzfile;
  VisMF::Read(u_macG_trac[2],utzfile);
#endif
  
#ifdef MG_USE_FBOXLIB
  if (model != model_list["richard"] &&
      model != model_list["steady-saturated"])
    {
      std::string rxfile = "/rhs_RhoD_x";
      std::string ryfile = "/rhs_RhoD_y";
      rxfile = FullPath + rxfile;
      ryfile = FullPath + ryfile;
      VisMF::Read(rhs_RhoD[0],rxfile);
      VisMF::Read(rhs_RhoD[1],ryfile);
      
#if (BL_SPACEDIM == 3)
      std::string rzfile = "/rhs_RhoD_z";
      rzfile = FullPath + rzfile;
      VisMF::Read(rhs_RhoD[2],rzfile);
#endif
    }
#endif

  is_grid_changed_after_regrid = true;
  if (grids == papa.getLevel(level).boxArray())
    is_grid_changed_after_regrid = false;

  // Set up boundary condition work
  setup_bound_desc();
}

void
PorousMedia::buildMetrics ()
{
  //
  // Build volume and face area arrays.
  //
  geom.GetVolume(volume,grids,GEOM_GROW);
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      geom.GetFaceArea(area[dir],grids,dir,GEOM_GROW);
    }
}

//
// Reset the time levels to time (time) and timestep dt.
// This is done at the start of the timestep in the pressure iteration section.
//

void
PorousMedia::resetState (Real time,
                         Real dt_old,
                         Real dt_new)
{
  for (int k = 0; k < num_state_type; k++)
    {
      state[k].reset();
      state[k].setTimeLevel(time,dt_old,dt_new);
    }
}

//
// Set the time levels to time (time) and timestep dt.
//
void
PorousMedia::setTimeLevel (Real time,
                           Real dt_old,
                           Real dt_new)
{
  for (int k = 0; k < num_state_type; k++)
    state[k].setTimeLevel(time,dt_old,dt_new);
}

void
PorousMedia::set_vel_from_bcs(Real      time,
			      MultiFab* vel)
{
  FArrayBox inflow;
  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation face = oitr();
    if (get_inflow_velocity(face,inflow,time)) {
      // NOTE: Without doing a pressure solve formally, it is not
      //  obvious how to do this in a way that makes sense.  The hack
      //  for the moment is to just do a setval over the entire 
      //  field....
      Real inflow_val = inflow(inflow.box().smallEnd(),0);
      for (int d=0; d<BL_SPACEDIM; ++d) {
	if (d==face.coordDir()) {
	  vel[d].setVal(inflow_val);
	} else {
	  vel[d].setVal(0);
	}
      }
    }
  }
}


//
// This function initializes the all relevant data.  
// It calls subroutines in PROB_$D.F
//
void
PorousMedia::initData ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::initData()");
    
    if (verbose > 1 && ParallelDescriptor::IOProcessor())
        std::cout << "Initializing data ...\n";
    
    // 
    // Initialize rock properties
    //
    init_rock_properties();
    //
    // Initialize the state and the pressure.
    //
    const Real* dx       = geom.CellSize();
    MultiFab&   S_new    = get_new_data(State_Type);
    MultiFab&   P_new    = get_new_data(Press_Type);
    MultiFab&   U_vcr    = get_new_data(  Vcr_Type);
    
    const Real  cur_time = state[State_Type].curTime();
    S_new.setVal(0.);
    P_new.setVal(0.);
    
    //
    // Initialized only based on solutions at the current level
    //
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        BL_ASSERT(grids[mfi.index()] == mfi.validbox());
        
        FArrayBox& sdat = S_new[mfi];
	FArrayBox& pdat = P_new[mfi];
        DEF_LIMITS(sdat,s_ptr,s_lo,s_hi);
        DEF_LIMITS(pdat,p_ptr,p_lo,p_hi);
        
        for (int i=0; i<ic_array.size(); ++i)
        {
            const RegionData& ic = ic_array[i];
            const PArray<Region>& ic_regions = ic.Regions();
            const std::string& type = ic.Type();
            
            if (type == "file") 
	    {
                std::cerr << "Initialization of initial condition based on "
                          << "a file has not been implemented yet.\n";
                BoxLib::Abort("PorousMedia::initData()");
	    }
            else if (type == "scalar") 
	    {
                Array<Real> vals = ic();
                for (int jt=0; jt<ic_regions.size(); ++jt) {
                    regions[jt].setVal(S_new[mfi],vals,
                                       dx,0,0,ncomps);
		}

		// set pressure to single value
		if (model==model_list["richard"]) {
		  const int idx = mfi.index();
		  Array<int> s_bc = getBCArray(State_Type,idx,0,1);
		  calcCapillary((*pcnp1_cc)[mfi],S_new[mfi],(*rock_phi)[mfi],
				(*kappa)[mfi],(*cpl_coef)[mfi],grids[idx],
				s_bc);
		}
	    }
            else if (type == "pressure") 
	    {
	      BL_ASSERT(model == model_list["steady-saturated"]);
	      Array<Real> vals = ic();
	      for (int jt=0; jt<ic_regions.size(); ++jt) {
		regions[jt].setVal(P_new[mfi],vals,
				   dx,0,0,ncomps);
	      }
	      for (int n=0; n<ncomps; ++n) {
		S_new[mfi].setVal(density[n],mfi.validbox(),n,1);
	      }
	    }
            else if (type == "hydrostatic") 
	    {
                Array<Real> vals = ic();
                BL_ASSERT(model >= 2);
                FArrayBox& cdat = (*cpl_coef)[mfi];
                const int n_cpl_coef = cpl_coef->nComp();
                DEF_CLIMITS(cdat,c_ptr,c_lo,c_hi);
                FORT_HYDRO(s_ptr, ARLIM(s_lo),ARLIM(s_hi), 
                           density.dataPtr(),&ncomps, 
                           c_ptr, ARLIM(c_lo),ARLIM(c_hi), &n_cpl_coef,
                           dx, vals.dataPtr(), &gravity);

		// set pressure to phydrostatic
		FORT_HYDRO_PRESSURE(p_ptr, ARLIM(p_lo),ARLIM(p_hi), 
				    density.dataPtr(),&ncomps, 
				    dx, vals.dataPtr(), &gravity);
		

	    }
            else if (type == "zero_total_velocity")
	    {	     
                BL_ASSERT(model != model_list["single-phase"] && 
                          model != model_list["single-phase-solid"]);
                int nc = 1;
		Array<Real> vals = ic();
		FArrayBox& cdat = (*cpl_coef)[mfi];
		const int n_cpl_coef = cpl_coef->nComp();
		DEF_CLIMITS(cdat,c_ptr,c_lo,c_hi);
		FArrayBox& kdat = (*kr_coef)[mfi];
		FArrayBox& kpdat = kpedge[BL_SPACEDIM-1][mfi];
		const int n_kr_coef = kr_coef->nComp();
		DEF_CLIMITS(kdat,k_ptr,k_lo,k_hi);
		DEF_CLIMITS(kpdat,kp_ptr,kp_lo,kp_hi);
		FORT_STEADYSTATE(s_ptr, ARLIM(s_lo),ARLIM(s_hi), 
				 density.dataPtr(),muval.dataPtr(),&ncomps,
				 kp_ptr, ARLIM(kp_lo),ARLIM(kp_hi), 
				 k_ptr, ARLIM(k_lo),ARLIM(k_hi), &n_kr_coef,
				 dx, &vals[0], &nc, &gravity);
		
		// set pressure
		if (model==model_list["richard"]) {
		  const int idx = mfi.index();
		  Array<int> s_bc = getBCArray(State_Type,idx,0,1);
		  calcCapillary((*pcnp1_cc)[mfi],S_new[mfi],(*rock_phi)[mfi],
				(*kappa)[mfi],(*cpl_coef)[mfi],grids[idx],
				s_bc);
		  P_new[mfi].copy((*pcnp1_cc)[mfi]);
		  P_new[mfi].mult(-1.0);
		}
	    }
            else
	    {
                FORT_INITDATA(&level,&cur_time,
                              s_ptr, ARLIM(s_lo),ARLIM(s_hi), 
                              density.dataPtr(), &ncomps, dx);

		// set pressure to single value
		if (model==model_list["richard"]) {
		  const int idx = mfi.index();
		  Array<int> s_bc = getBCArray(State_Type,idx,0,1);
		  calcCapillary((*pcnp1_cc)[mfi],S_new[mfi],(*rock_phi)[mfi],
				(*kappa)[mfi],(*cpl_coef)[mfi],grids[idx],
				s_bc);
		}
	    }
	}
        
        if (ntracers > 0)
        {		
            for (int iTracer=0; iTracer<ntracers; ++iTracer)
            {
                const PArray<RegionData>& rds = tic_array[iTracer];

                for (int i=0; i<rds.size(); ++i)
                {
                    const RegionData& tic = rds[i];
                    const PArray<Region>& tic_regions = tic.Regions();
                    const std::string& tic_type = tic.Type();
                    
                    if (tic_type == "file") 
                    {
                        std::cerr << "Initialization of initial condition based on "
                                  << "a file has not been implemented yet.\n";
                        BoxLib::Abort("PorousMedia::initData()");
                    }
                    else if (tic_type == "concentration") 
                    {
                        Array<Real> val = tic();
                        for (int jt=0; jt<tic_regions.size(); ++jt) {
                            BL_ASSERT(val.size()>=1);
                            BL_ASSERT(sdat.nComp()>ncomps+iTracer);
                            BL_ASSERT(tic_regions.size()>jt);

                            tic_regions[jt].setVal(sdat,val[0],ncomps+iTracer,dx,0);
                        }
                    }
                    else {
                        std::string m = "Unrecognized tracer ic type: " + tic_type;
                        BoxLib::Abort(m.c_str());
                    }
                }
            }
        }
    }

    if (do_chem>0) 
    {
        if (ntracers>0)
        {
            const Real* dx = geom.CellSize();
            MultiFab& AuxChem_new = get_new_data(Aux_Chem_Type);
            for (MFIter mfi(AuxChem_new); mfi.isValid(); ++mfi)
            {
                FArrayBox& fab = AuxChem_new[mfi];
                
                for (ChemICMap::iterator it=sorption_isotherm_ics.begin(); it!=sorption_isotherm_ics.end(); ++it)
                {
                    const std::string& material_name = it->first;
                    const PArray<Region>& rock_regions = find_rock(material_name).regions;

                    ICLabelParmPair& solute_to_pp = it->second; 
                    for (ICLabelParmPair::iterator it1=solute_to_pp.begin(); it1!=solute_to_pp.end(); ++it1) 
                    {
                        const std::string& solute_name = it1->first;
                        ICParmPair& parm_pairs = it1->second;
                        for (ICParmPair::iterator it2=parm_pairs.begin(); it2!=parm_pairs.end(); ++it2) 
                        {
                            const std::string& parameter = it2->first;
                            int comp = sorption_isotherm_label_map[solute_name][parameter];
                            Real value = it2->second;
                            if (false) {
                                std::cout << material_name << "   "
                                          << solute_name << "    "
                                          << parameter << "    "
                                          << value << std::endl;
                            }
                            for (int j=0; j<rock_regions.size(); ++j) {
                                rock_regions[j].setVal(fab,value,comp,dx,0);
                            }
                        }
                    }
                }
            }
        }

        if (nminerals>0) {
            const Real* dx = geom.CellSize();
            MultiFab& AuxChem_new = get_new_data(Aux_Chem_Type);
            for (MFIter mfi(AuxChem_new); mfi.isValid(); ++mfi)
            {
                FArrayBox& fab = AuxChem_new[mfi];
                
                for (ChemICMap::iterator it=mineralogy_ics.begin(); it!=mineralogy_ics.end(); ++it)
                {
                    const std::string& material_name = it->first;
                    const PArray<Region>& rock_regions = find_rock(material_name).regions;
                    
                    ICLabelParmPair& mineral_to_pp = it->second; 
                    for (ICLabelParmPair::iterator it1=mineral_to_pp.begin(); it1!=mineral_to_pp.end(); ++it1) 
                    {
                        const std::string& mineral_name = it1->first;
                        ICParmPair& parm_pairs = it1->second;
                        for (ICParmPair::iterator it2=parm_pairs.begin(); it2!=parm_pairs.end(); ++it2) 
                        {
                            const std::string& parameter = it2->first;
                            int comp = mineralogy_label_map[mineral_name][parameter];
                            Real value = it2->second;
                            for (int j=0; j<rock_regions.size(); ++j) {
                                rock_regions[j].setVal(fab,value,comp,dx,0);
                            }
                        }
                    }
                }
            }
        }
        
       if (nsorption_sites>0) {
            const Real* dx = geom.CellSize();
            MultiFab& AuxChem_new = get_new_data(Aux_Chem_Type);
            for (MFIter mfi(AuxChem_new); mfi.isValid(); ++mfi)
            {
                FArrayBox& fab = AuxChem_new[mfi];
                
                for (ChemICMap::iterator it=surface_complexation_ics.begin(); it!=surface_complexation_ics.end(); ++it)
                {
                    const std::string& material_name = it->first;
                    const Rock& rock = find_rock(material_name);
                    ICLabelParmPair& sorption_site_to_pp = it->second; 
                    for (ICLabelParmPair::iterator it1=sorption_site_to_pp.begin(); it1!=sorption_site_to_pp.end(); ++it1) 
                    {
                        const std::string& sorption_site_name = it1->first;
                        ICParmPair& parm_pairs = it1->second;
                        for (ICParmPair::iterator it2=parm_pairs.begin(); it2!=parm_pairs.end(); ++it2) 
                        {
                            const std::string& parameter = it2->first;
                            int comp = surface_complexation_label_map[sorption_site_name][parameter];
                            Real value = it2->second;
                            
                            const PArray<Region>& rock_regions = rock.regions;
                            for (int j=0; j<rock_regions.size(); ++j) {
                                rock_regions[j].setVal(fab,value,comp,dx,0);
                            }
                        }
                    }
                }
            }
        }
        
        if (ncation_exchange>0) {
            const Real* dx = geom.CellSize();
            MultiFab& AuxChem_new = get_new_data(Aux_Chem_Type);
            for (MFIter mfi(AuxChem_new); mfi.isValid(); ++mfi)
            {
                FArrayBox& fab = AuxChem_new[mfi];
                
                for (ICParmPair::iterator it=cation_exchange_ics.begin(); it!=cation_exchange_ics.end(); ++it)
                {
                    const std::string& material_name = it->first;
                    const Rock& rock = find_rock(material_name);
                    Real value = it->second;
                    int comp = cation_exchange_label_map.begin()->second;
                    
                    const PArray<Region>& rock_regions = rock.regions;
                    for (int j=0; j<rock_regions.size(); ++j) {
                        rock_regions[j].setVal(fab,value,comp,dx,0);
                    }
                }
            }
        }

        // These always get set if do_chem>0
        {
            const Real* dx = geom.CellSize();
            MultiFab& AuxChem_new = get_new_data(Aux_Chem_Type);
            for (MFIter mfi(AuxChem_new); mfi.isValid(); ++mfi)
            {
                FArrayBox& fab = AuxChem_new[mfi];
                
                for (ChemICMap::iterator it=solute_chem_ics.begin(); it!=solute_chem_ics.end(); ++it)
                {
                    const std::string& material_name = it->first;
                    const Rock& rock = find_rock(material_name);
                    ICLabelParmPair& sorption_site_to_pp = it->second; 
                    for (ICLabelParmPair::iterator it1=sorption_site_to_pp.begin(); it1!=sorption_site_to_pp.end(); ++it1) 
                    {
                        const std::string& tracer_name = it1->first;
                        ICParmPair& parm_pairs = it1->second;
                        for (ICParmPair::iterator it2=parm_pairs.begin(); it2!=parm_pairs.end(); ++it2) 
                        {
                            const std::string& parameter = it2->first;
                            int comp = solute_chem_label_map[tracer_name][parameter];
                            Real value = it2->second;
                            
                            const PArray<Region>& rock_regions = rock.regions;
                            for (int j=0; j<rock_regions.size(); ++j) {
                                rock_regions[j].setVal(fab,value,comp,dx,0);
                            }
                        }
                    }
                }
            }
        }

        if (sorption_chem_ics.size()>0)
        {
            const Real* dx = geom.CellSize();
            MultiFab& AuxChem_new = get_new_data(Aux_Chem_Type);
            for (MFIter mfi(AuxChem_new); mfi.isValid(); ++mfi)
            {
                FArrayBox& fab = AuxChem_new[mfi];

                for (ICLabelParmPair::iterator it=sorption_chem_ics.begin(); it!=sorption_chem_ics.end(); ++it) 
                {
                    const std::string& solute_name = it->first;
                    ICParmPair& parm_pair = it->second;                    
                    for (ICParmPair::iterator it1=parm_pair.begin(); it1!=parm_pair.end(); ++it1) 
                    {
                        const std::string& parameter = it1->first;
                        int comp = sorption_chem_label_map[solute_name][parameter];
                        Real value = it1->second;
                        for (int i=0; i<rocks.size(); ++i)
                        {	  
                            const Rock& rock = rocks[i];
                            const PArray<Region>& rock_regions = rock.regions;
                            for (int j=0; j<rock_regions.size(); ++j) {
                                rock_regions[j].setVal(fab,value,comp,dx,0);
                            }
                        }
                    }
                }
            }
        }
    }

        
    if (model == model_list["richard"]) {
      calcInvPressure(S_new,P_new); // Set sat from p

      // FIXME: Why are these needed here?
      FillStateBndry(cur_time,State_Type,0,ncomps);
      FillStateBndry(cur_time,Press_Type,0,1);
      
      if (do_tracer_transport) {
	FillStateBndry(cur_time,State_Type,ncomps,ntracers);
      }
    }

    U_vcr.setVal(0.);
    if (have_capillary) calcCapillary(cur_time);
    //
    // compute lambda
    //
    if (model != model_list["steady-saturated"]) {
      calcLambda(cur_time);
    }

    //
    // Initialize u_mac_curr 
    //
    if (model == model_list["steady-saturated"])
      {
	set_vel_from_bcs(cur_time,u_mac_curr);
      }
    else
      {
	if (model != model_list["richard"])
	  {
	    mac_project(u_mac_curr,rhs_RhoD,cur_time);
	  }

      }
    
    is_grid_changed_after_regrid = false;
        
#ifdef AMANZI
    // Call chemistry to relax initial data to equilibrium
    bool chem_relax_ics = false;
    if (do_chem>0  &&  ic_chem_relax_dt>0)
      {
	MultiFab& Fcnt = get_new_data(FuncCount_Type);
	Fcnt.setVal(1);
	int nGrow = 0;
	strang_chem(cur_time,ic_chem_relax_dt,nGrow);
      }
#endif

    is_first_step_after_regrid = true;
    old_intersect_new          = grids;
}

void
RichardNLSdata::SetNLIterationsTaken(int iters){
  nl_iterations_taken = iters;
}

RichardNLSdata::RichardNLSdata(int slev, int nlevs, PMAmr* _pm_amr)
    : start_level(slev), end_level(slev+nlevs-1),
      num_Jacobian_reuses_remaining(nlevs,0)
{
    pm_amr = _pm_amr;
    BL_ASSERT(pm_amr);

    // Set solver defaults
    max_num_Jacobian_reuse = 0;
    max_nl_iterations = 20;
    max_nl_residual_norm = -1;    
    max_num_consecutive_success = 0;
    max_num_consecutive_failures_1 = 3;
    max_num_consecutive_failures_2 = 4;
    max_num_consecutive_increases = 15;
    consecutive_increase_reduction_factor = 0.4;
    min_nl_iterations_for_dt = 6;
    min_nl_iterations_for_dt_2 = 3;
    max_nl_iterations_for_dt = 10;
    time_step_increase_factor = 1.5;
    time_step_increase_factor_2 = 2.0;
    time_step_reduction_factor = 0.8;
    time_step_retry_factor = 0.5;
    time_step_retry_factor_2 = 0.1;
    time_step_retry_factor_f = 0.01;
    max_time_step_size = 1.e10;

    num_consecutive_success = 0;
    num_consecutive_failures_1 = 0;
    num_consecutive_failures_2 = 0;
    num_consecutive_increases = 0;
    first = true;

    ResetCounters();
    ResetJacobianCounter();

    // Allocate necessary memory and setup pointers
    Build();
}

void RichardNLSdata::SetMaxConsecutiveFails(int max_num) {max_num_consecutive_failures_1=max_num;}
void RichardNLSdata::SetDtRetryFactor(Real factor) {time_step_retry_factor = factor;}
void RichardNLSdata::SetMaxConsecutiveFails2(int max_num) {max_num_consecutive_failures_2=max_num;}
void RichardNLSdata::SetDtRetryFactor2(Real factor) {time_step_retry_factor_2 = factor;}
void RichardNLSdata::SetDtRetryFactorF(Real factor) {time_step_retry_factor_f = factor;}
void RichardNLSdata::SetMaxConsecutiveErrIncrease(int max_incr) {max_num_consecutive_increases=max_incr;}
void RichardNLSdata::SetConsecutiveErrIncreaseDtReduction(Real redux) {consecutive_increase_reduction_factor=redux;}
void RichardNLSdata::SetMaxConsecutiveSuccess(int max_num) {max_num_consecutive_success=max_num;}
void RichardNLSdata::SetMaxNewtonIterations(int max_iter) {max_nl_iterations=max_iter;}
void RichardNLSdata::SetMaxJacobianReuse(int max_num_reuse) {max_num_Jacobian_reuse=max_num_reuse;}
void RichardNLSdata::ResetJacobianCounter(int lev) {num_Jacobian_reuses_remaining[lev]=max_num_Jacobian_reuse;}
void RichardNLSdata::SetMaxNewtonIterationsForDt(int max_iter) {max_nl_iterations_for_dt=max_iter;}
void RichardNLSdata::SetMinNewtonIterationsForDt(int min_iter) {min_nl_iterations_for_dt=min_iter;}
void RichardNLSdata::SetMinNewtonIterationsForDt2(int min_iter) {min_nl_iterations_for_dt_2=min_iter;}
void RichardNLSdata::SetDtIncreaseFactor(Real factor) {time_step_increase_factor=factor;}
void RichardNLSdata::SetDtIncreaseFactor2(Real factor) {time_step_increase_factor_2=factor;}
void RichardNLSdata::SetDtReductionFactor(Real factor) {time_step_reduction_factor=factor;}
void RichardNLSdata::SetMaxDt(Real dt_max) {max_time_step_size=dt_max;}


void
RichardNLSdata::ResetCounters()
{
    nl_iterations_taken = 0;
    nl_residual_norm = -1; 
    //num_consecutive_success = 0;
    //num_consecutive_failures_1 = 0;
    //num_consecutive_failures_2 = 0;
    //num_consecutive_increases = 0;
    last_chance = false;;
    prev_abs_err = -1;
    //first = true;
}

void
RichardNLSdata::ResetJacobianCounter()
{
    int nlevs = end_level - start_level +1;
    for (int lev=0; lev<nlevs; ++lev) {
        ResetJacobianCounter(lev);
    }
}

bool
RichardNLSdata::UpdateJacobian(int lev)
{
    bool do_Jacobian_eval = false;
    num_Jacobian_reuses_remaining[lev]--;
    if (num_Jacobian_reuses_remaining[lev] <= 0) {
        do_Jacobian_eval = true;
        num_Jacobian_reuses_remaining[lev] = max_num_Jacobian_reuse;
    }
    return do_Jacobian_eval;
}
void
RichardNLSdata::Build()
{
  int nLevs = pm_amr->finestLevel() + 1;
  velPhase.resize(nLevs);
  initialState.resize(nLevs,PArrayManage);
  for (int lev = 0; lev <nLevs; lev++) {
    BoxArray grids = pm_amr->getLevel(start_level+lev).boxArray();
    initialState.set(lev,new MultiFab(grids,1,1));
    velPhase[lev] = dynamic_cast<PorousMedia*>(&pm_amr->getLevel(start_level+lev))->UMac_Curr();
  }

  Hcoeffs.resize(BL_SPACEDIM);
  Jacobian.resize(BL_SPACEDIM);
  for (int dir = 0; dir < BL_SPACEDIM; dir++) 
    {
      Hcoeffs[dir].resize(nLevs,PArrayManage);
      Jacobian[dir].resize(nLevs,PArrayManage);
      for (int lev = 0; lev < nLevs; lev++) 
	{
	  BoxArray grids = BoxArray(pm_amr->getLevel(start_level+lev).boxArray()).surroundingNodes(dir);
	  Hcoeffs[dir].set(lev, new MultiFab(grids,1,0));
	  Jacobian[dir].set(lev, new MultiFab(grids,3,0));
	}
    }
  DAlpha.resize(nLevs,PArrayManage);
  for (int lev = 0; lev <nLevs; lev++) {
    BoxArray grids = pm_amr->getLevel(start_level+lev).boxArray();
    DAlpha.set(lev, new MultiFab(grids,1,1)); // Why is this grow cell required/filled?
  }
}

bool
RichardNLSdata::AdjustDt(Real                dt, 
			 RichardNLSdata::Reason nl_solver_status, 
			 Real                abs_err, 
			 Real                rel_err, 
			 Real&               dt_new) // Note: return bool for whether run should stop
{
    dt_new = dt;
    if (first) {
        //num_consecutive_increases = 0;
        //num_consecutive_success = 0;
        first = false;
        prev_abs_err = -1;
    }

    if (nl_solver_status == RICHARD_SUCCESS)
    {
        last_chance = false;

        // "success" is when the error is reduced using small number of iters
        // In this case, increment counter for this event, reset "increase" counter
        // If this keeps happening, increase dt and reset the counter for these events
        //  (when we do, if the problem was  particularly easy, increase dt dramatically)
        if (nl_iterations_taken < min_nl_iterations_for_dt && 
            (prev_abs_err <= 0  ||   prev_abs_err > abs_err) ) {
            
            num_consecutive_success++;
            num_consecutive_increases = 0;
            
            if (num_consecutive_success >= max_num_consecutive_success)
            {
                num_consecutive_success = 0;
                Real fac = time_step_increase_factor;
                if (nl_iterations_taken < min_nl_iterations_for_dt_2) {
                    fac = time_step_increase_factor_2;
                }
                dt_new = dt * fac;
            }
        }

        // "increase" is when the error grows and used large number of iters
        // In this case, increment counter for this event, guarantee recalc of J, and reset "success" counter
        // If this keeps happening, reduce dt and reset the counter for these events
        if (nl_iterations_taken > max_nl_iterations_for_dt && 
            (prev_abs_err <= 0  ||   prev_abs_err < abs_err) )
        {
            ResetJacobianCounter();
            num_consecutive_increases++;
            num_consecutive_success = 0;        
        
            if (nl_iterations_taken > max_nl_iterations_for_dt)
            {
                ResetJacobianCounter();
                dt_new = dt * time_step_reduction_factor;
            }
        }

        num_consecutive_failures_1 = 0;
        num_consecutive_failures_2 = 0;
        prev_abs_err = abs_err;

    }
    else {

        // step was rejected
        num_consecutive_failures_1++;

        if (num_consecutive_failures_1 <= max_num_consecutive_failures_1)
        {
            dt_new = dt * time_step_retry_factor;
#if 0
            // If the last increase was immediately undone, cut back on the dt adjustment knobs...
            //  FIXME: needs more tweaking
            if (num_consecutive_success == 0) {
                time_step_increase_factor = 0.5*(1 + time_step_increase_factor);
                time_step_retry_factor = 0.5*(1 + time_step_retry_factor);
            }
#endif
        }
        else
        {
            num_consecutive_failures_2++;

            if (num_consecutive_failures_2 <= max_num_consecutive_failures_2)
            {
                dt_new = dt * time_step_retry_factor_2;
            }
            else
            {
                if (last_chance)  return false;
                dt_new = dt * time_step_retry_factor_f;
                last_chance = true;
            }
        }

        num_consecutive_success = 0;
        ResetJacobianCounter();
    }
    dt_new = std::min(max_time_step_size,dt_new);
    return true;
}

bool
RichardNLSdata::AdjustDt(Real                dt, 
			 RichardNLSdata::Reason nl_solver_status, 
			 Real&               dt_new) // Note: return bool for whether run should stop
{
    dt_new = dt;
    if (nl_solver_status == RICHARD_SUCCESS)
    {
        last_chance = false;

        // "success" is when the error is reduced using small number of iters
        // In this case, increment counter for this event, reset "increase" counter
        // If this keeps happening, increase dt and reset the counter for these events
        //  (when we do, if the problem was  particularly easy, increase dt dramatically)
        if (nl_iterations_taken < min_nl_iterations_for_dt ) {
            
            num_consecutive_success++;
            num_consecutive_increases = 0;
            
            if (num_consecutive_success >= max_num_consecutive_success)
            {
                num_consecutive_success = 0;
                Real fac = time_step_increase_factor;
                if (nl_iterations_taken < min_nl_iterations_for_dt_2) {
                    fac = time_step_increase_factor_2;
                }
                dt_new = dt * fac;
            }
        }

        // "increase" is when large number of iters
        // In this case, increment counter for this event, guarantee recalc of J, 
	// and reset "success" counter
        // If this keeps happening, reduce dt and reset the counter for these events
        if (nl_iterations_taken > max_nl_iterations_for_dt  )
        {
            ResetJacobianCounter();
            num_consecutive_increases++;
            num_consecutive_success = 0;        
        
            if (nl_iterations_taken > max_nl_iterations_for_dt)
            {
                ResetJacobianCounter();
                dt_new = dt * time_step_reduction_factor;
            }
        }

        num_consecutive_failures_1 = 0;
        num_consecutive_failures_2 = 0;

    }
    else {

        // step was rejected
        num_consecutive_failures_1++;

        if (num_consecutive_failures_1 <= max_num_consecutive_failures_1)
        {
            dt_new = dt * time_step_retry_factor;
#if 0
            // If the last increase was immediately undone, cut back on the dt adjustment knobs...
            //  FIXME: needs more tweaking
            if (num_consecutive_success == 0) {
                time_step_increase_factor = 0.5*(1 + time_step_increase_factor);
                time_step_retry_factor = 0.5*(1 + time_step_retry_factor);
            }
#endif
        }
        else
        {
            num_consecutive_failures_2++;

            if (num_consecutive_failures_2 <= max_num_consecutive_failures_2)
            {
                dt_new = dt * time_step_retry_factor_2;
            }
            else
            {
                if (last_chance)  return false;
                dt_new = dt * time_step_retry_factor_f;
                last_chance = true;
            }
        }

        num_consecutive_success = 0;
        ResetJacobianCounter();
    }

    dt_new = std::min(max_time_step_size,dt_new);
    return true;
}

RichardNLSdata
PorousMedia::BuildInitNLS()
{
    int nlevs = parent->finestLevel() - level + 1;
    PMAmr* pm_amr = dynamic_cast<PMAmr*>(parent);
    if (!pm_amr)
      BoxLib::Abort("Bad cast in PorousMedia::BuildInitNLS");
    RichardNLSdata nld(0,nlevs,pm_amr);
    nld.SetMaxJacobianReuse(0); // Currently switched off because it didnt seem to buy anything
    
    nld.SetMaxConsecutiveFails(steady_max_consecutive_failures_1);
    nld.SetDtRetryFactor(steady_time_step_retry_factor_1);
    
    nld.SetMaxConsecutiveFails2(steady_max_consecutive_failures_2);
    nld.SetDtRetryFactor2(steady_time_step_retry_factor_2);
    nld.SetDtRetryFactorF(steady_time_step_retry_factor_f);
    
    nld.SetMinNewtonIterationsForDt(steady_min_iterations);
    nld.SetDtIncreaseFactor(steady_time_step_increase_factor);
    nld.SetMinNewtonIterationsForDt2(steady_min_iterations_2);
    nld.SetDtIncreaseFactor2(steady_time_step_increase_factor_2);
    
    nld.SetMaxNewtonIterationsForDt(steady_max_iterations);
    nld.SetDtReductionFactor(steady_time_step_reduction_factor);
    
    nld.SetMaxNewtonIterations(steady_limit_iterations);
    
    nld.SetMaxConsecutiveErrIncrease(steady_max_num_consecutive_increases);
    nld.SetConsecutiveErrIncreaseDtReduction(steady_consecutive_increase_reduction_factor);
    
    nld.SetMaxConsecutiveSuccess(steady_max_num_consecutive_success);

    nld.SetMaxDt(steady_max_time_step_size);
    return nld;
}

#include <RichardSolver.H>

std::map<int,std::string> PETSc_Reasons;
static std::string
GetPETScReason(int flag) 
{
  PETSc_Reasons[2] = "SNES_CONVERGED_FNORM_ABS     ";
  PETSc_Reasons[3] = "SNES_CONVERGED_FNORM_RELATIVE"; // ||F|| < atol 
  PETSc_Reasons[4] = "SNES_CONVERGED_SNORM_RELATIVE"; // Newton computed step size small; || delta x || < stol 
  PETSc_Reasons[5] = "SNES_CONVERGED_ITS           "; // maximum iterations reached 
  PETSc_Reasons[7] = "SNES_CONVERGED_TR_DELTA      ";
  PETSc_Reasons[-1] = "SNES_DIVERGED_FUNCTION_DOMAIN"; // the new x location passed the function is not in the domain of F
  PETSc_Reasons[-2] = "SNES_DIVERGED_FUNCTION_COUNT ";
  PETSc_Reasons[-3] = "SNES_DIVERGED_LINEAR_SOLVE   "; // the linear solve failed
  PETSc_Reasons[-4] = "SNES_DIVERGED_FNORM_NAN      ";
  PETSc_Reasons[-5] = "SNES_DIVERGED_MAX_IT         ";
  PETSc_Reasons[-6] = "SNES_DIVERGED_LINE_SEARCH    "; // the line search failed 
  PETSc_Reasons[-7] = "SNES_DIVERGED_INNER          "; // inner solve failed
  PETSc_Reasons[-8] = "SNES_DIVERGED_LOCAL_MIN      "; // || J^T b || is small, implies converged to local minimum of F()
  PETSc_Reasons[0]  = "SNES_CONVERGED_ITERATING     ";
  if (PETSc_Reasons.find(flag)==PETSc_Reasons.end()) {
    BoxLib::Abort("Unknown PETSc return flag");
  }
  return PETSc_Reasons[flag];
}

void
PorousMedia::richard_init_to_steady()
{
#ifdef MG_USE_FBOXLIB
    //
    // Richard initialization
    //
    if (model == model_list["richard"])
    {
        std::string tag = "Richard Init-to-Steady";
        if (richard_init_to_steady_verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << tag << std::endl;
        }        

        bool do_brute_force = false;
        //do_brute_force = true;
        
        if (level == 0) {
            if (do_brute_force)
                richard_eqb_update(u_mac_curr);
            else
            {
                int old_richard_solver_verbose = richard_solver_verbose;
		richard_solver_verbose = richard_init_to_steady_verbose;
		initial_iter = 1;
                
                Real cur_time = state[State_Type].curTime();
                Real prev_time = state[State_Type].prevTime();
                int  finest_level = parent->finestLevel();
                Array<Real> dt_save(finest_level+1);
                Array<int> nc_save(finest_level+1);
                int  n_factor;
                for (int k = 0; k <= finest_level; k++)
                {
                    nc_save[k] = parent->nCycle(k);
                    dt_save[k] = parent->dtLevel()[k];
                    dt_save[k] = parent->getLevel(k).get_state_data(0).curTime()
                        - getLevel(k).get_state_data(0).prevTime();
                }

                Real dt_init = steady_init_time_step;
                Real t_max = (execution_mode=="init_to_steady" ? switch_time : steady_max_psuedo_time);
                int k_max = steady_max_time_steps;
                Real t_eps = 1.e-8*dt_init;

                MultiFab tmp(grids,1,1);
                MultiFab tmpP(grids,1,1);
                int nc = 0; // Component of water in state

                Real prev_abs_err, init_abs_err;
                Real rel_err = -1;
                Real abs_err = -1;

                Real dt = dt_init;
                Real t = 0;
                bool first = true;
                int k = 0;
                bool solved = false;
                bool continue_iterations = (!solved)  &&  (k < k_max)  &&  (t < t_max);

                RichardNLSdata::Reason ret;
                RichardNLSdata nld = BuildInitNLS();

		RichardSolver* rs = 0;
                if (steady_use_PETSc_snes) {
		  RichardSolver::RSParams rsparams;
		  rsparams.max_ls_iterations = richard_max_ls_iterations;
		  rsparams.min_ls_factor = richard_min_ls_factor;
		  rsparams.ls_acceptance_factor = richard_ls_acceptance_factor;
		  rsparams.ls_reduction_factor = richard_ls_reduction_factor;
		  rsparams.monitor_line_search = richard_monitor_line_search;
		  rsparams.errfd = richard_perturbation_scale_for_J;
		  rsparams.maxit = steady_limit_iterations;
		  rsparams.maxf = steady_limit_function_evals;
		  rsparams.atol = steady_abs_tolerance;
		  rsparams.rtol = steady_rel_tolerance;
		  rsparams.stol = steady_abs_update_tolerance;
		  rsparams.use_fd_jac = richard_use_fd_jac;
		  rsparams.use_dense_Jacobian = richard_use_dense_Jacobian;
		  rsparams.upwind_krel = richard_upwind_krel;
		  rsparams.pressure_maxorder = richard_pressure_maxorder;

		  rs = new RichardSolver(*(PMParent()),rsparams);
                }

                int total_num_Newton_iterations = 0;
                int total_rejected_Newton_steps = 0;
                while (continue_iterations)
		{

		  // Advance the state data structures
		  for (int lev=0;lev<= finest_level;lev++)
		    {
		      PorousMedia& pm = getLevel(lev);
		      for (int i = 0; i < num_state_type; i++)
			{
			  pm.state[i].allocOldData();
			  pm.state[i].swapTimeLevels(dt);
			}
		    }
		  // FIXME:
		  // If execution mode is Initialize To Steady, time increments as we evolve here.
		  // If mode is Transient, the loop here does not advance time
		  cur_time = state[Press_Type].curTime();
		  prev_time = state[Press_Type].prevTime();

		  for (int lev=0;lev<= finest_level;lev++)
		    {
		      PorousMedia& pm = getLevel(lev);
		      for (int i = 0; i < num_state_type; i++)
			{
			  MultiFab& od = pm.get_old_data(i);
			  MultiFab& nd = pm.get_new_data(i);
			  MultiFab::Copy(nd,od,0,0,od.nComp(),0);  // Guess for next time step
			}
		    }

		  if (steady_use_PETSc_snes) {
		    rs->ResetRhoSat();
		  }

		  MultiFab& S_new = get_new_data(State_Type);
		  MultiFab::Copy(tmp,get_new_data(Press_Type),nc,0,1,0);
		  tmp.mult(-1.0);

		  if (richard_init_to_steady_verbose && ParallelDescriptor::IOProcessor()) {
                        std::cout << tag << "  t=" << t 
                                  << ", n=" << k << ", dt=" << dt << '\n';
                    }

                    if (do_multilevel_full) {
                        nld.ResetCounters();
                        nld.ResetJacobianCounter();	

			// Save the initial state 
			for (int lev=0;lev<= finest_level;lev++)
			  {
			    PorousMedia&    fine_lev   = getLevel(lev);
			    if (do_richard_sat_solve)
			      {
				MultiFab& S_lev = fine_lev.get_new_data(State_Type);
				MultiFab::Copy(nld.initialState[lev],S_lev,0,0,1,1);
			      }
			    else
			      {
				MultiFab& P_lev = fine_lev.get_new_data(Press_Type);
				MultiFab::Copy(nld.initialState[lev],P_lev,0,0,1,1);
			      }
			  }

			if (steady_use_PETSc_snes) 
			  {
 			      int retCode = rs->Solve(t+dt, dt, k, nld);
                              if (retCode >= 0) {
                                  ret = RichardNLSdata::RICHARD_SUCCESS;
                              } 
                              else {
                                  if (ret == -3) {
                                      ret = RichardNLSdata::RICHARD_LINEAR_FAIL;
                                  }
                                  else {
                                      ret = RichardNLSdata::RICHARD_NONLINEAR_FAIL;
                                      if (richard_solver_verbose>1 && ParallelDescriptor::IOProcessor())
                                          std::cout << "     **************** Newton failed: " << GetPETScReason(retCode) << '\n';
                                  }
                              }
			  }
			else
			  {
			    ret = richard_composite_update(dt,nld);
			  }
			total_num_Newton_iterations += nld.NLIterationsTaken();
                    }
                    else {
                        int curr_nwt_iter = steady_limit_iterations;
                        ret = richard_scalar_update(dt,curr_nwt_iter,u_mac_curr);
                        total_num_Newton_iterations += curr_nwt_iter;
                    }

                    if (ret == RichardNLSdata::RICHARD_SUCCESS)
                    {
                        prev_abs_err = abs_err;
                        MultiFab::Add(tmp,get_new_data(Press_Type),nc,0,1,0);
                        abs_err = tmp.norm2(0);
                        if (first){
                            init_abs_err = abs_err;
                            first = false;
                        }
                        else {
                            rel_err = abs_err / init_abs_err;
                        }
                    }
                    else {

                        if (steady_abort_on_psuedo_timestep_failure)
                        {
                            BoxLib::Abort("Aborting as instructed when timestep fails");
                        }

                        total_rejected_Newton_steps++;
			for (int lev=0;lev<= finest_level;lev++)
			  {
			    PorousMedia&    fine_lev   = getLevel(lev);

			    for (int k = 0; k < num_state_type; k++)
			      {
				fine_lev.state[k].reset();
			      }

			    if (do_richard_sat_solve)
			      {
				MultiFab& S_lev = fine_lev.get_new_data(State_Type);
				MultiFab::Copy(S_lev,nld.initialState[lev],0,0,1,1);
			      }
			    else
			      {
				MultiFab& P_lev = fine_lev.get_new_data(Press_Type);
				MultiFab::Copy(P_lev,nld.initialState[lev],0,0,1,1);
			      }
			  }
                    }

                    Real dt_new;
                    bool cont = nld.AdjustDt(dt,ret,abs_err,rel_err,dt_new); 

                    if (ret == RichardNLSdata::RICHARD_SUCCESS)
                    {
                        k++;
                        t += dt;
			if (execution_mode=="init_to_steady")
			  {
			    solved = false;
			  }
			else
			  {
			    solved = ((abs_err <= steady_abs_update_tolerance) 
				      || ((rel_err>0)  && (rel_err <= steady_rel_update_tolerance)) );
			  }
                        continue_iterations = cont && (!solved)  &&  (k < k_max)  &&  (t < t_max);

                        if (richard_init_to_steady_verbose>1 && ParallelDescriptor::IOProcessor()) {
                            std::cout << tag << "   Step successful, Niters=" << nld.NLIterationsTaken();
                            if ( std::abs(dt - dt_new) > t_eps) {
                                std::cout << "   dt factor = " << dt_new / dt;
                            }
                            std::cout << std::endl;
                        }
                    }
                    else 
                    {
                        if (richard_init_to_steady_verbose>1 && ParallelDescriptor::IOProcessor()) {
                            std::cout << tag << "   Step failed ";
                            if (ret==RichardNLSdata::RICHARD_NONLINEAR_FAIL) {
                                std::cout << "(NL failure)";
                            }
                            else {
                                std::cout << "(L failure) ";
                            }

                            if ( std::abs(dt - dt_new) > t_eps) {
                                std::cout << "   dt factor = " << dt_new / dt;
                            }
                            std::cout << std::endl;
                        }
                    }

                    dt = std::min(dt_new, t_max-t);
                }

                delete rs;
                if (richard_init_to_steady_verbose && ParallelDescriptor::IOProcessor()) {
                    std::cout << tag << " Total psuedo-time advanced: " << t << " in " << k << " steps" << std::endl;
                    std::cout << tag << "      Newton iters: " << total_num_Newton_iterations << std::endl;
                    std::cout << tag << "      Rejected steps: " << total_rejected_Newton_steps << std::endl;
                }

                if (richard_init_to_steady_verbose && ParallelDescriptor::IOProcessor()) {
                    if (solved) {
                        std::cout << tag << " Success!  Steady solution found" << std::endl;
                    }
                    else {
                        std::cout << tag << " Warning: solution is not steady.  Continuing..." << std::endl;
                    }
                }

                richard_solver_verbose = old_richard_solver_verbose;
		initial_iter = 0;

                //
                // Re-instate timestep.
                //	
		PMAmr* p = dynamic_cast<PMAmr*>(parent); BL_ASSERT(p);
		if (execution_mode=="init_to_steady") {
		  p->SetStartTime(t_max);
		  p->SetCumTime(t_max);
		}
		Real time_after_init = p->StartTime();
                parent->setDtLevel(dt_save);
                parent->setNCycle(nc_save);
                for (int k = 0; k <= finest_level; k++)
                {
		  getLevel(k).setTimeLevel(time_after_init,dt_save[k],dt_save[k]);
                }
            }
        }
    }
#endif
}

//
// Fills a new level n with best level n and coarser data available.
//

void
PorousMedia::init (AmrLevel& old)
{
  init_rock_properties();

  PorousMedia*  oldns     = (PorousMedia*) &old;
  const Real    dt_new    = parent->dtLevel(level);
  const Real    cur_time  = oldns->state[State_Type].curTime();
  const Real    prev_time = oldns->state[State_Type].prevTime();
  const Real    dt_old    = cur_time - prev_time;

  MultiFab&     S_new     = get_new_data(State_Type);
  MultiFab&     P_new     = get_new_data(Press_Type);
  MultiFab&     U_cor     = get_new_data(  Vcr_Type);

  U_cor.setVal(0.);

  dt_eig = oldns->dt_eig;
    
  setTimeLevel(cur_time,dt_old,dt_new);
    
  //
  // Get best state data: from old. 
  //
  for (FillPatchIterator fpi(old,S_new,0,cur_time,State_Type,0,NUM_SCALARS);
       fpi.isValid();
       ++fpi)
  {
    S_new[fpi.index()].copy(fpi());
  }
  //
  // Subsequent pressure solve will give the correct pressure.
  // 
  for (FillPatchIterator fpi(old,P_new,0,cur_time,Press_Type,0,1);
       fpi.isValid();
       ++fpi)
  {
    P_new[fpi.index()].copy(fpi());
  }

  //
  // Get best edge-centered velocity data: from old.
  //
  const BoxArray& old_grids = oldns->grids;
  is_grid_changed_after_regrid = true;
  if (old_grids == grids)
    {
      for (int dir=0; dir<BL_SPACEDIM; ++dir)
	{
	  u_mac_curr[dir].copy(oldns->u_mac_curr[dir]);
	  rhs_RhoD[dir].copy(oldns->rhs_RhoD[dir]);
	  u_macG_trac[dir].copy(oldns->u_macG_trac[dir]);
	}
      is_grid_changed_after_regrid = false;
    }


  /*  if (!is_grid_changed_after_regrid && model == model_list["richard"])
    {
      MultiFab P_tmp(grids,1,1);
      MultiFab::Copy(P_tmp,P_new,0,0,1,1);
      P_tmp.mult(-1.0);
      calcInvCapillary(S_new,P_tmp);
      }*/

  //
  // Get best cell-centered velocity data: from old.
  //

#ifdef AMANZI
  if (do_chem>0)
    {
      MultiFab& FC_new  = get_new_data(FuncCount_Type); 
      
      for (FillPatchIterator fpi(old,FC_new,FC_new.nGrow(),cur_time,FuncCount_Type,0,1);
	   fpi.isValid();
           ++fpi)
	{
	  FC_new[fpi.index()].copy(fpi());
	}
    }
#endif
    
  old_intersect_new          = BoxLib::intersect(grids,oldns->boxArray());
  is_first_step_after_regrid = true;
}

void
PorousMedia::init ()
{
  BL_ASSERT(level > 0);
    
  MultiFab& S_new = get_new_data(State_Type);
  MultiFab& P_new = get_new_data(Press_Type);
  MultiFab& U_cor = get_new_data(  Vcr_Type);
   
  const Array<Real>& dt_amr = parent->dtLevel();
  Array<Real>        dt_new(level+1);

  for (int lev = 0; lev < level; lev++)
    dt_new[lev] = dt_amr[lev];
  //
  // Guess new dt from new data (interpolated from coarser level).
  //
  const Real dt = dt_new[level-1]/Real(parent->MaxRefRatio(level-1));
  dt_new[level] = dt;
  parent->setDtLevel(dt_new);

  //
  // Compute dt based on old data.
  //
  PorousMedia& old       = getLevel(level-1);
  const Real   cur_time  = old.state[State_Type].curTime();
  const Real   prev_time = old.state[State_Type].prevTime();
  const Real   dt_old    = (cur_time-prev_time)/Real(parent->MaxRefRatio(level-1));

  setTimeLevel(cur_time,dt_old,dt);
  //
  // Get best coarse state, pressure and velocity data.
  //
  FillCoarsePatch(S_new,0,cur_time,State_Type,0,ncomps);
  if (ntracers>0) {
      if (do_tracer_transport) {
          FillCoarsePatch(S_new,0,cur_time,State_Type,ncomps,ntracers);
      }
      else {
          // FIXME
          // If !do_tracer_transport, we dont have the bc descriptors, so we punt
          S_new.setVal(0.0,ncomps,ntracers);
      }
  }
  FillCoarsePatch(P_new,0,cur_time,Press_Type,0,1);
  U_cor.setVal(0.);

#ifdef AMANZI
  if (do_chem>0)
    {
      FillCoarsePatch(get_new_data(FuncCount_Type),0,cur_time,FuncCount_Type,0,1);
    }
#endif

  init_rock_properties();
  old_intersect_new = grids;
}

//
// ADVANCE FUNCTIONS
//

//
// This function ensures that the multifab registers and boundary
// flux registers needed for syncing the composite grid
//
//     u_mac, umacG, Ssync, fr_adv, fr_visc
//
// are initialized to zero.  These quantities and the  
// advective velocity registers (mac_reg) are compiled by first
// setting them to the coarse value acquired during a coarse timestep
// and then incrementing in the fine values acquired during the
// subcycled fine timesteps.  This compilation procedure occurs in
// different parts for different quantities
//
// * u_mac is set in mac_project.
// * fr_adv, fr_visc are set in scalar_advect
// * Ssync is set in subcycled calls to post_timestep
// * mac_reg is set in mac_project
//
// After these quantities have been compiled during a coarse
// timestep and subcycled fine timesteps.  The post_timestep function
// uses them to sync the fine and coarse levels.  If the coarse level
// is not the base level, post_timestep modifies the next coarsest levels
// registers appropriately.
//
// Note :: There is a little ambiguity as to which level owns the
// boundary flux registers.  The Multifab registers are quantities
// sized by the coarse level BoxArray and belong to the coarse level.
// The fine levels own the boundary registers, since they are sized by
// the boundaries of the fine level BoxArray.
//

void
PorousMedia::advance_setup (Real time,
                            Real dt,
                            int  iteration,
                            int  ncycle)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::advance_setup()");

  const int finest_level = parent->finestLevel();

  if (level < finest_level)
    {
      if (Ssync == 0)
	Ssync = new MultiFab(grids,NUM_SCALARS,1);
      Ssync->setVal(0);
    }

  //
  // Set reflux registers to zero.
  //
  if (do_reflux && level < finest_level)
    {
      getAdvFluxReg(level+1).setVal(0);
      getViscFluxReg(level+1).setVal(0);
    }

  //
  // Alloc space for edge velocities (normal comp only).
  //
  if (u_macG_prev == 0)
    {
      u_macG_prev = new MultiFab[BL_SPACEDIM];

      for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
	  BoxArray edge_grids(grids);
	  edge_grids.surroundingNodes(dir).grow(1);
	  u_macG_prev[dir].define(edge_grids,1,0,Fab_allocate);
	  u_macG_prev[dir].setVal(1.e40);
        }
    }
  if (u_macG_curr == 0)
    {
      u_macG_curr = new MultiFab[BL_SPACEDIM];

      for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
	  BoxArray edge_grids(grids);
	  edge_grids.surroundingNodes(dir).grow(1);
	  u_macG_curr[dir].define(edge_grids,1,0,Fab_allocate);
	  u_macG_curr[dir].setVal(1.e40);
        }
    }
  //
  // Set up state multifabs for the advance.
  //
  for (int k = 0; k < num_state_type; k++)
    {
        if (k!=FuncCount_Type && k!=Aux_Chem_Type) {
            state[k].allocOldData();
            state[k].swapTimeLevels(dt);
        }
    }
  //
  // Alloc MultiFab to hold advective update terms.
  //
  BL_ASSERT(aofs == 0);
  aofs = new MultiFab(grids,NUM_SCALARS,0);
  
  if (model != model_list["steady-saturated"]) {
    //
    // Compute lambda at cell centers
    //

      // FIXME: this test is messed up!
    if (model != model_list["single-phase"] || 
	model != model_list["single-phase-solid"]) 
      {
	calcLambda(time); 
#ifdef MG_USE_FBOXLIB
	if (model != model_list["richard"])
#endif
	  calcDLambda(time);
	MultiFab::Copy(*lambdap1_cc,*lambda_cc,0,0,ncomps,1);
      }
    //
    // Compute diffusion coefficients
    //

    // FIXME: Is this appropriate for richard?
    if (variable_scal_diff)
      {
	calcDiffusivity(time,0,ncomps);
	MultiFab::Copy(*diffnp1_cc,*diffn_cc,0,0,ndiff,1);
      }
    //
    // Compute capillary diffusive coefficients
    //
    if (have_capillary)
      {
	calcCapillary(time);
	MultiFab::Copy(*pcnp1_cc,*pcn_cc,0,0,1,(*pcnp1_cc).nGrow());
      }  
    //
    // If we are not doing a full advection scheme, u_mac_curr 
    // must be recomputed if grid has changed after a timestep.
    //
#ifdef MG_USE_FBOXLIB
    if (model != model_list["richard"])
#endif
      {
	if (do_simple == 0 && (full_cycle == 1 || no_corrector == 1))
	  {
	    if (n_pressure_interval == 0)
	      mac_project(u_mac_curr,rhs_RhoD,time);
	    else
	      {
		if (level == 0)   it_pressure += 1;
		
		if (it_pressure == n_pressure_interval &&
		    parent->levelSteps(level)%parent->nCycle(level)==parent->nCycle(level)-1)
		  {
		    mac_project(u_mac_curr,rhs_RhoD,time);
		    if (level == parent->finestLevel()) it_pressure = 0;
		  }	    
	      }
	  }
	else if (is_grid_changed_after_regrid)
	  {
	    mac_project(u_mac_curr,rhs_RhoD,time);
	  }
      }
    //
    // Alloc MultiFab to hold correction velocity.
    //
    if (u_corr == 0)
      {
	u_corr = new MultiFab[BL_SPACEDIM];
	for (int dir = 0; dir < BL_SPACEDIM; dir++)
	  {
	    BoxArray edge_grids(grids);
	    edge_grids.surroundingNodes(dir).grow(1);
	    u_corr[dir].define(edge_grids,1,0,Fab_allocate);
	    u_corr[dir].setVal(0.);
	  }
      }
  }
    
  //
  // Swap the time levels of u_mac
  //
  MultiFab* dummy = u_mac_curr;
  u_mac_curr = u_mac_prev;
  u_mac_prev = dummy;

#ifdef AMANZI
  if (do_chem>0)
    {
      aux_boundary_data_old.setVal(1.e30);
    }
#endif

  //
  // Copy cell-centered correction velocity computed in 
  // previous timestep to current timestep.
  //
  MultiFab& Uc_old = get_old_data(Vcr_Type);
  MultiFab& Uc_new = get_new_data(Vcr_Type);
  MultiFab::Copy(Uc_new,Uc_old,0,0,BL_SPACEDIM,Uc_new.nGrow());
}

//
// Clean up after the advance function.
//
void
PorousMedia::advance_cleanup (Real dt,
                              int  iteration,
                              int  ncycle)
{
  delete aofs;
  aofs = 0;
}

//
// Compute a timestep at a level. Return largest safe timestep.
//
Real
PorousMedia::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  ncycle)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::advance()");
  if (do_multilevel_full) 
  {
      if (level == 0) {
          if (verbose > 1 && ParallelDescriptor::IOProcessor())
          {
              std::cout << "Advancing all levels:"
                        << " starting time = " << time
                        << " with dt = "               << dt << '\n';
          }
          multilevel_advance(time,dt,iteration,ncycle);
      }
  }
  else
  {
    if (verbose > 1 && ParallelDescriptor::IOProcessor())
      {
	std::cout << "Advancing grids at level " << level
		  << " : starting time = "       << time
		  << " with dt = "               << dt << '\n';
      }

    advance_setup(time,dt,iteration,ncycle);

    FillPatchedOldState_ok = true;
    //
    // Advance the old state for a Strang-split dt/2.  Include grow cells in
    // advance, and squirrel these away for diffusion and Godunov guys to
    // access for overwriting non-advanced fill-patched grow cell data.
    //
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);
    
    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& P_old = get_old_data(Press_Type);
    
    MultiFab::Copy(S_new,S_old,0,0,NUM_SCALARS,S_old.nGrow()); 
    MultiFab::Copy(P_new,P_old,0,0,1,P_old.nGrow());
    
    const Real pcTime = state[State_Type].curTime();
    
    FillStateBndry (pcTime,State_Type,0,ncomps+ntracers);
    FillStateBndry (pcTime,Press_Type,0,1);

    if (do_chem>0)
      {
	if (do_full_strang)
	  {
	    if (verbose>2 && ParallelDescriptor::IOProcessor())
	      std::cout << "... advancing 1/2 strang step for chemistry\n";
	    
	    // Old state is chem-advanced in-place.  Hook set to get grow data
	    //  from squirreled away data rather than via a vanilla fillpatch
	    strang_chem(time,dt/2,HYP_GROW);
	    FillPatchedOldState_ok = false;
	  }
      }

#ifdef MG_USE_FBOXLIB
    if (model == model_list["richard"]) 
      {
	advance_richard(time,dt);
      }
    else
#endif
      {
	// 
	// FIXME: Should we leave this 
	// do_simple: 2 ==> Only solve the tracer equations; assume steady state.
	//            1 ==> Only solve the pressure equation at time 0.
	//            0 ==> Solve the pressure equation at every timestep.
	//
	if (do_simple == 2 && !is_grid_changed_after_regrid)
	  advance_tracer(time,dt);
	else if (do_simple == 1  && !is_grid_changed_after_regrid)
	  advance_simple(time,dt);
	else
	  advance_incompressible(time,dt);
      }
    
    is_grid_changed_after_regrid = false;
    
    // second half of the strang splitting
    if (do_chem>0)
      {      
	if (do_full_strang)
	  {
	    if (verbose>2 && ParallelDescriptor::IOProcessor())
	      std::cout << "Second 1/2 Strang step of chemistry\n";
	    
	    // New state is chem-advanced in-place.  Fillpatch hook unset
	    strang_chem(pcTime,dt/2.0,0);
	    FillPatchedOldState_ok = true;
	  }
	else
	  {
	    if (n_chem_interval == 0)
	      {
		if (verbose>2 && ParallelDescriptor::IOProcessor())
		  std::cout << "... advancing full strang step for chemistry\n";
		strang_chem(pcTime,dt,0);
	      }
	    else
	      {
		if (level == 0)
		  {
		    it_chem += 1;
		    dt_chem += dt;
		  }
		
		if (it_chem == n_chem_interval &&
		    parent->levelSteps(level)%parent->nCycle(level)==parent->nCycle(level)-1 &&
		    level == parent->finestLevel())
		  {
		    if (verbose>2 && ParallelDescriptor::IOProcessor())
		      std::cout << "... advancing full strang step for chemistry with dt ="
				<< dt_chem << "\n";
		    
		    strang_chem(pcTime,dt_chem,0);
		    it_chem = 0;
		    dt_chem = 0;		    
		  }
	      }
	  }
      }
    
    // 
    // Check sum of components
    //
    if (verbose>3) check_sum();
    
    //
    // Clean up after the predicted value at t^n+1.
    //
    advance_cleanup(dt,iteration,ncycle);
  }

  // Dummy value : not used for determining time step.
  Real dt_test = 1.e20; 
  return dt_test; 
}

void
PorousMedia::multilevel_advance (Real time,
				 Real dt,
				 int  iteration,
				 int  ncycle)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::multilevel_advance()");

  BL_ASSERT(do_multilevel_full);

  if (level == 0) 
  {
    for (int lev = 0; lev <= parent->finestLevel(); lev++)
    {
      PorousMedia&    pm_lev   = getLevel(lev);
      
      pm_lev.advance_setup(time,dt,iteration,ncycle);
      pm_lev.FillPatchedOldState_ok = true;

      MultiFab& S_new = pm_lev.get_new_data(State_Type);
      MultiFab& S_old = pm_lev.get_old_data(State_Type);

      MultiFab& P_new = pm_lev.get_new_data(Press_Type);
      MultiFab& P_old = pm_lev.get_old_data(Press_Type);

      S_new.setVal(0.);
      P_new.setVal(0.);
      MultiFab::Copy(S_new,S_old,0,0,ncomps+ntracers,S_old.nGrow()); 
      MultiFab::Copy(P_new,P_old,0,0,1,P_old.nGrow());

      Real pcTime = pm_lev.state[State_Type].curTime();
      int ncomp_fill_bndry = do_tracer_transport ? ncomps+ntracers : ncomps;
      pm_lev.FillStateBndry (pcTime,State_Type,0,ncomp_fill_bndry);
      pm_lev.FillStateBndry (pcTime,Press_Type,0,1);

      if (have_capillary)	pm_lev.calcCapillary(pcTime);

    }
    
    // If do_chem==0, then no reaction.
    // Otherwise, type of reactions depends on magnitude of have_corereact.
    if (do_chem>0)
      {
	if (do_full_strang)
	  {
	    if (verbose>2 && ParallelDescriptor::IOProcessor())
	      std::cout << "... advancing 1/2 strang step for chemistry\n";
	    
	    for (int lev = 0; lev <= parent->finestLevel(); lev++)
	      {
		PorousMedia& pm_lev = getLevel(lev);
		pm_lev.strang_chem(time,dt/2,HYP_GROW);
		// FIXME: Have no code for chem-advancing grow cells
		//FillPatchedOldState_ok = false;
	      }

	    if (verbose>2 && ParallelDescriptor::IOProcessor())
	      std::cout << "PorousMedia::advance(): end of first 1/2 Strang step\n";
	  }
      }
  }
#ifdef MG_USE_FBOXLIB
  if (model == model_list["richard"]) 
  {
      if (verbose>2 && do_chem && ParallelDescriptor::IOProcessor())
          std::cout << "... advancing full dt for flow and transport\n";
      advance_multilevel_richard(time,dt);
    }
#endif
  if (model == model_list["steady-saturated"])
    {
      advance_multilevel_saturated(time,dt);
    }
  //
  // second half of the strang splitting
  //
  if (do_chem>0)
    {
      if (do_full_strang)
	{
	  if (verbose>2 && ParallelDescriptor::IOProcessor())
	    std::cout << "Second 1/2 Strang step of chemistry\n";
	  
	  for (int lev = 0; lev <= parent->finestLevel(); lev++)
	    {
	      PorousMedia&  pm_lev = getLevel(lev);
	      pm_lev.strang_chem(time+dt,dt/2.0,0);
	    }
	}
      else
	{
	  if (n_chem_interval == 0)
	    {
	      if (verbose>2 && ParallelDescriptor::IOProcessor())
		std::cout << "... advancing full strang step for chemistry\n";
	      
	      for (int lev = 0; lev <= parent->finestLevel(); lev++)
		{
		  PorousMedia&  pm_lev = getLevel(lev);
		  pm_lev.strang_chem(time+dt,dt,0);
		}
	    }
	  else
	    {
	      it_chem += 1;
	      dt_chem += dt;
	      
	      if (it_chem == n_chem_interval)
		{
		  if (verbose>2 && ParallelDescriptor::IOProcessor())
		    std::cout << "... advancing full strang step for chemistry with dt ="
			      << dt_chem << "\n";
		  for (int lev = 0; lev <= parent->finestLevel(); lev++)
		    {
		      PorousMedia&  pm_lev = getLevel(lev);
		      pm_lev.strang_chem(time+dt,dt_chem,0);      
                      pm_lev.FillPatchedOldState_ok = true;
		    }
		  
		  it_chem = 0;
		  dt_chem = 0;
		}
	    }	    
	}
    }

  for (int lev = parent->finestLevel(); lev >= 0; lev--)
  {
    PorousMedia&  pm_lev     = getLevel(lev);

    pm_lev.avgDown();

    if (verbose>3) pm_lev.check_sum();      

    pm_lev.advance_cleanup(dt,iteration,ncycle);
  }

  if (verbose>2 && ParallelDescriptor::IOProcessor())
    std::cout << "PorousMedia::advance(): end of multilevel advance\n";

}

void
PorousMedia::advance_incompressible (Real time,
				     Real dt)
{
  // 
  // Time stepping for incompressible flow.  
  // For single-phase constant-density problem, use advance_simple.
  //
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::advance_incompressible()");

  const Real cur_time = state[State_Type].curTime();
  MultiFab& S_new     = get_new_data(State_Type);
  int lscalar         = ncomps - 1; 

  if (n_pressure_interval !=0)
    check_minmax(u_mac_prev);

  MultiFab* rhod_tmp = 0;
  if (do_any_diffuse)
    {
      rhod_tmp = new MultiFab[BL_SPACEDIM];
      for (int dir =0; dir < BL_SPACEDIM; dir++) 
	{
	  BoxArray edge_grids(grids);
	  edge_grids.surroundingNodes(dir);
	  rhod_tmp[dir].define(edge_grids,1,0,Fab_allocate);
	  rhod_tmp[dir].setVal(0.0);
	  rhod_tmp[dir].plus(rhs_RhoD[dir],0,1,0);
	}
    }
     
  if (level == 0) 
    create_umac_grown(u_mac_prev,u_macG_prev);
  else 
    {
      PArray<MultiFab> u_macG_crse(BL_SPACEDIM,PArrayManage);
      GetCrseUmac(u_macG_crse,time);
      create_umac_grown(u_mac_prev,u_macG_crse,u_macG_prev); 
    }
  
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    MultiFab::Copy(u_macG_trac[dir],u_macG_prev[dir],0,0,1,0);
  
  //
  // Predictor: Advance the component conservation equations
  //
  int corrector = 0;
    
  if (no_corrector == 1)
    {
      corrector = 1;

      // copy u_mac_prev to u_mac_curr since we are not solving for u_mac_curr
      for (int dir = 0; dir < BL_SPACEDIM; dir++)
	{
	  MultiFab::Copy(u_mac_curr[dir],u_mac_prev[dir],0,0,1,0);
	  MultiFab::Copy(u_macG_curr[dir],u_macG_prev[dir],0,0,1,0);
	}

      // contribute to velocity register
      mac_projector->contribute_to_mac_reg(level,u_mac_prev);
      if (do_any_diffuse)
	{
	  for (int dir = 0; dir < BL_SPACEDIM; dir++) 
	    rhod_tmp[dir].mult(-1.0);
	  mac_projector->contribute_to_mac_reg_rhoD(level,rhod_tmp);
	}

      // Compute the advective term
      scalar_advection(u_macG_trac,dt,0,lscalar,true);

      // Add the advective and other terms to get scalars at t^{n+1}.
      scalar_update(dt,0,ncomps,corrector,u_macG_trac);

      if (do_tracer_transport && ntracers > 0)
	{
	  int ltracer = ncomps+ntracers-1;
	  tracer_advection(u_macG_trac,dt,ncomps,ltracer,true);
	}

      predictDT(u_macG_prev);
    }

  else
    {
      // Compute the advective term
      scalar_advection(u_macG_trac,dt,0,lscalar,false);

      // Add the advective and other terms to get scalars at t^{n+1}.
      scalar_update(dt,0,ncomps,corrector);

      if (do_chem>0)
	{
	  if (do_full_strang)
	    strang_chem(time,dt/2.0,HYP_GROW);
	}
      
      //
      // Corrector Step
      //    
      if (model > 1)
	calcLambda(cur_time);

      // Do a MAC projection to define edge velocities at time t^(n+1)
      mac_project(u_mac_curr,rhs_RhoD,cur_time);
    
      if (do_any_diffuse)
	{
	  for (int dir = 0; dir < BL_SPACEDIM; dir++) 
	    {
	      rhod_tmp[dir].plus(rhs_RhoD[dir],0,1,0);
	      rhod_tmp[dir].mult(-0.5);
	    }
	  mac_projector->contribute_to_mac_reg_rhoD(level,rhod_tmp);
	}

      if (level == 0) 
	create_umac_grown(u_mac_curr,u_macG_curr);
      else 
	{
	  PArray<MultiFab> u_macG_crse(BL_SPACEDIM,PArrayManage);
	  GetCrseUmac(u_macG_crse,time+dt);
	  create_umac_grown(u_mac_curr,u_macG_crse,u_macG_curr);
	}

      // Create velocity at time t^{n+1/2}.
      MultiFab* u_mac_nph  = new MultiFab[BL_SPACEDIM];
      MultiFab* u_macG_nph = new MultiFab[BL_SPACEDIM];
      for (int dir = 0; dir < BL_SPACEDIM; dir++)
	{
	  BoxArray edge_grids(grids);
	  edge_grids.surroundingNodes(dir);
	  u_mac_nph[dir].define(edge_grids,1,0,Fab_allocate);
	  MultiFab::Copy(u_mac_nph[dir],u_mac_prev[dir],0,0,1,0);
	  u_mac_nph[dir].plus(u_mac_curr[dir],0,1,0);
	  u_mac_nph[dir].mult(0.5);
	  if (do_any_diffuse)
	    u_mac_nph[dir].plus(rhod_tmp[dir],0,1,0);

	  edge_grids.grow(1);
	  u_macG_nph[dir].define(edge_grids,1,0,Fab_allocate);
	  MultiFab::Copy(u_macG_nph[dir],u_macG_prev[dir],0,0,1,0);
	  u_macG_nph[dir].plus(u_macG_curr[dir],0,1,0);
	  u_macG_nph[dir].mult(0.5);

	  MultiFab::Copy(u_macG_trac[dir],u_macG_nph[dir],0,0,1,0);
	}

      mac_projector->contribute_to_mac_reg(level,u_mac_nph);

      // Re-advect component equations 
      corrector = 1;
      if (variable_scal_diff)
	calcDiffusivity (cur_time,0,ncomps);

      scalar_advection(u_macG_trac,dt,0,lscalar,true);
    
      scalar_update(dt,0,ncomps,corrector,u_macG_trac);

      if (do_tracer_transport && ntracers > 0)
	{
	  int ltracer = ncomps+ntracers-1;
	  tracer_advection(u_macG_trac,dt,ncomps,ltracer,true);
	}

      // predict the next time step.  
      predictDT(u_macG_curr);

      delete [] u_mac_nph;
      delete [] u_macG_nph;
    }
  
  if (do_any_diffuse) delete [] rhod_tmp;

  //
  // Check the divergence conditions of v_1 (water)
  //
  MultiFab divutmp(grids,1,0);
  divutmp.setVal(0);
  mac_projector->check_div_cond(level,divutmp,u_macG_trac,rhs_RhoD);
  MultiFab::Copy(S_new,divutmp,0,ncomps+ntracers,1,0);
  if (have_capillary) MultiFab::Copy(S_new,*pcnp1_cc,0,ncomps+ntracers+1,1,1);
  
}

void
PorousMedia::advance_simple (Real time,
			     Real dt)
{
  // 
  // Time stepping for incompressible single-phase single-density flow.
  //
  if (level == 0) 
    create_umac_grown(u_mac_prev,u_macG_prev);
  else 
    {
      PArray<MultiFab> u_macG_crse(BL_SPACEDIM,PArrayManage);
      GetCrseUmac(u_macG_crse,time);
      create_umac_grown(u_mac_prev,u_macG_crse,u_macG_prev); 
    }

  //
  // Single advance the component conservation equations
  //
  int corrector = 0;

  // Compute the coefficients for diffusion operators.
  if (variable_scal_diff)
    {
      calcDiffusivity(time,0,ncomps);
      MultiFab::Copy(*diffnp1_cc,*diffn_cc,0,0,ndiff,diffn_cc->nGrow());
    }

  // Compute the advective term
  scalar_advection(u_macG_prev,dt,0,ncomps,true);
    
  // Add the advective and other terms to get scalars at t^{n+1}.
  scalar_update(dt,0,ncomps,corrector);

  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    u_mac_curr[dir].copy(u_mac_prev[dir]);
}

#ifdef MG_USE_FBOXLIB
void
PorousMedia::advance_richard (Real time,
			      Real dt)
{
  std::string tag = "Richard Time Step: ";
  // 
  // Time stepping for richard's equation
  //
  int curr_nwt_iter = 20;
  RichardNLSdata::Reason ret = richard_scalar_update(dt,curr_nwt_iter,u_mac_curr);
  if (ParallelDescriptor::IOProcessor())
    {
      std::cout << tag;
      if (ret == RichardNLSdata::RICHARD_LINEAR_FAIL) {
	std::cout << " - linear solver failure ";
      }
      else if ( ret == RichardNLSdata::RICHARD_NONLINEAR_FAIL) {
	std::cout << tag << " - nonlinear solver failure ";
      }
      std::cout << std::endl;
    }

  BL_ASSERT(ret == RichardNLSdata::RICHARD_SUCCESS);

  compute_vel_phase(u_mac_curr,0,time+dt);
    
  if (level == 0) {
    create_umac_grown(u_mac_curr,u_macG_trac);
  } else {
      PArray<MultiFab> u_macG_crse(BL_SPACEDIM,PArrayManage);
      GetCrseUmac(u_macG_crse,time);
      create_umac_grown(u_mac_curr,u_macG_crse,u_macG_trac); 
  }

  if (do_tracer_transport && ntracers > 0)
    {
      int ltracer = ncomps+ntracers-1;
      tracer_advection(u_macG_trac,dt,ncomps,ltracer,true);
    }

  // predict the next time step. 
  Real dt_nwt = dt; 
  predictDT(u_macG_trac);
  if (curr_nwt_iter < 10 && dt_grow_max>=1) 
    dt_nwt *= dt_grow_max;
  //if (richard_time < 3.0*richard_time_min)
  //{
  //  dt_nwt = dt_nwt*change_max;
      //if (curr_nwt_iter <= richard_iter && curr_nwt_iter < 4 && dt_nwt < richard_max_dt)
      //	dt_nwt = dt_nwt*1.1;
      //else if (curr_nwt_iter > 5)
      //	dt_nwt = dt_nwt*0.75;
      //else if (curr_nwt_iter < 2) 
      //	dt_nwt = dt_nwt*1.1;
      //richard_iter = curr_nwt_iter;
  //}
  dt_eig = std::min(dt_eig,dt_nwt); 
}

void
PorousMedia::advance_multilevel_richard (Real time,
					 Real dt)
{
  // 
  // Time stepping for richard's equation
  // We will do subcycling if the time step is too big and 
  // and reset the subsequent time step to the smller one.
  //

  if (level == 0) {
    int finest_level = parent->finestLevel();
    int nlevs = finest_level + 1;
    PMAmr* pm_amr = dynamic_cast<PMAmr*>(parent);
    if (!pm_amr) 
      BoxLib::Abort("Bad cast in PorousMedia::advance_multilevel_richard");
    RichardNLSdata nld(0,nlevs,pm_amr);
    RichardSolver* rs = 0;
    if (steady_use_PETSc_snes) {
      RichardSolver::RSParams rsparams;
      rsparams.max_ls_iterations = richard_max_ls_iterations;
      rsparams.min_ls_factor = richard_min_ls_factor;
      rsparams.ls_acceptance_factor = richard_ls_acceptance_factor;
      rsparams.ls_reduction_factor = richard_ls_reduction_factor;
      rsparams.monitor_line_search = richard_monitor_line_search;
      rsparams.errfd = richard_perturbation_scale_for_J;
      rsparams.maxit = steady_limit_iterations;
      rsparams.maxf = steady_limit_function_evals;
      rsparams.atol = steady_abs_tolerance;
      rsparams.rtol = steady_rel_tolerance;
      rsparams.stol = steady_abs_update_tolerance;
      rsparams.use_fd_jac = richard_use_fd_jac;
      rsparams.use_dense_Jacobian = richard_use_dense_Jacobian;
      rsparams.upwind_krel = richard_upwind_krel;
      rsparams.pressure_maxorder = richard_pressure_maxorder;
      
      rs = new RichardSolver(*(PMParent()),rsparams);
    }
    
    Array<Real> dt_save(nlevs);
    Array<int> nc_save(nlevs);
    int n_factor;
  
    for (int k = 0; k <= finest_level; k++)
      {
	nc_save[k] = parent->nCycle(k);
	dt_save[k] = parent->dtLevel()[k];
	dt_save[k] = parent->getLevel(k).get_state_data(0).curTime()
	  - getLevel(k).get_state_data(0).prevTime();
      }
    
    Real dt_iter = dt;
    Real t_max = time+dt;
    int  k_max = steady_max_time_steps;
    Real t_eps = 1.e-8*dt_iter;

    Real t = time;
    int  k = 0;
    bool continue_subtimestep = (k < k_max) && (t < t_max);
    bool continue_iterations  = (dt_iter > t_eps);

    RichardNLSdata::Reason ret;
    Real dt_new = dt_iter;

    // Lazily build structure to save state at "time".  If we must subcycle, the
    // algorithm will overwrite old_time data as it goes.  This saved_state
    // must include all state types involved in this subcycle; we make a set
    // of ids and set them manually to minimize the overhead of this
    std::set<int> types_advanced;
    types_advanced.insert(State_Type);
    types_advanced.insert(Press_Type);
    Array<PArray<MultiFab> > saved_states;

    while(continue_subtimestep) 
      {
	if (t != time  &&  saved_states.size()==0) 
	  {
	    saved_states.resize(types_advanced.size());
	    for (std::set<int>::const_iterator it=types_advanced.begin(), End=types_advanced.end(); it!=End; ++it) 
	      {
		saved_states[*it].resize(nlevs,PArrayManage);
		for (int lev=0; lev<nlevs; ++lev) 
		  {
		    const MultiFab& old = parent->getLevel(lev).get_old_data(*it);
		    saved_states[*it].set(lev, new MultiFab(old.boxArray(), old.nComp(), old.nGrow()));
		    MultiFab::Copy(saved_states[*it][lev],old,0,0,old.nComp(),old.nGrow());
		  }
	      }
	  }

	// Solve the richard equation for the given time step.  
	// Try with progressively smaller timestep if it fails.
	// Save the initial state to enable retry
	for (int lev=0;lev<nlevs;lev++)
	  {
	    PorousMedia&    fine_lev   = getLevel(lev);
	    if (do_richard_sat_solve)
	      {
		MultiFab& S_lev = fine_lev.get_new_data(State_Type);
		MultiFab::Copy(nld.initialState[lev],S_lev,0,0,1,1);
	      }
	    else
	      {
		MultiFab& P_lev = fine_lev.get_new_data(Press_Type);
		MultiFab::Copy(nld.initialState[lev],P_lev,0,0,1,1);
	      }
	  }

	continue_iterations = (dt_iter > t_eps);
	int iter = 0;
	int total_num_Newton_iterations = 0;
	int total_rejected_Newton_steps = 0;
	while (continue_iterations)
	  {
	    nld.ResetCounters();
	    nld.ResetJacobianCounter();

	    if (verbose > 2 && ParallelDescriptor::IOProcessor())
	      std::cout << "Attempting multi-level flow advance: Inner Iter " << iter++ 
			<< " From Time: " << t << " Time Step: " << dt_iter  << std::endl;

	    // Advance the state data structures
	    for (std::set<int>::const_iterator it=types_advanced.begin(), End=types_advanced.end(); it!=End; ++it) 
	      {
		for (int lev=0;lev<= finest_level;lev++)
		  {
		    PorousMedia& pm = getLevel(lev);
		    pm.state[*it].allocOldData();
		    pm.state[*it].swapTimeLevels(dt_iter);
		  }
	      }

	    if (steady_use_PETSc_snes) 
	      {
		rs->ResetRhoSat();
		int retCode = rs->Solve(t+dt_iter, dt_iter, k, nld);
		if (retCode >= 0) {
		  ret = RichardNLSdata::RICHARD_SUCCESS;
		} 
		else {
		  if (ret == -3) {
		    ret = RichardNLSdata::RICHARD_LINEAR_FAIL;
		  }
		  else {
		    ret = RichardNLSdata::RICHARD_NONLINEAR_FAIL;
		    if (richard_solver_verbose>1 && ParallelDescriptor::IOProcessor())
		      std::cout << "     **************** Newton failed: " << GetPETScReason(retCode) << '\n';
		  }
		}
	      }
	    else
	      {
		ret = richard_composite_update(dt_iter,nld);
	      }
	    total_num_Newton_iterations += nld.NLIterationsTaken();

	    bool cont = nld.AdjustDt(dt_iter,ret,dt_new); 

	    if (ret == RichardNLSdata::RICHARD_SUCCESS) 
	      {
		continue_iterations = false;

		for (int lev=0; lev<nlevs; lev++)
		  {
		    PorousMedia&    fine_lev   = getLevel(lev);  
		    
		    // This is now done in the solvers
		    //fine_lev.compute_vel_phase(fine_lev.u_mac_curr,0,t_max);
		    
		    if (lev == 0) 
		      fine_lev.create_umac_grown(fine_lev.u_mac_curr,
						 fine_lev.u_macG_trac);
		    else 
		      {
			PArray<MultiFab> u_macG_crse(BL_SPACEDIM,PArrayManage);
			fine_lev.GetCrseUmac(u_macG_crse,time);
			fine_lev.create_umac_grown(fine_lev.u_mac_curr,u_macG_crse,
						   fine_lev.u_macG_trac); 
		      }
		    
		    if (do_tracer_transport && ntracers > 0)
		      {
			int ltracer = ncomps+ntracers-1;
			
			bool cont_tracer_advection = true;
			Real t_tracer = 0.;
			
			fine_lev.predictDT(fine_lev.u_macG_trac);
			dt_eig = std::min(dt_eig,dt_iter);
			
			
			MultiFab& S_new = fine_lev.get_new_data(State_Type);
			MultiFab& S_old = fine_lev.get_old_data(State_Type);
			MultiFab Stmp(grids,1,1); 
			for (int i=ncomps;i<=ltracer;i++)
			  {
			    MultiFab::Copy(Stmp,S_old,i,0,1,1);
			    MultiFab::Multiply(Stmp,S_old,0,0,1,1);
			    Stmp.divide(S_new,0,1,1);
			    MultiFab::Copy(S_old,Stmp,0,i,1,1);
			  }
			MultiFab::Copy(S_old,S_new,0,0,1,1);
			
			while (cont_tracer_advection)
			  {
			    t_tracer += dt_eig;
			    fine_lev.tracer_advection(fine_lev.u_macG_trac,dt_eig,ncomps,ltracer,true);
			    
			    if (t_tracer < dt_iter)
			      {
				dt_eig = std::min(dt_eig,dt_iter-t_tracer);	
				MultiFab& S_new = fine_lev.get_new_data(State_Type);
				MultiFab& S_old = fine_lev.get_old_data(State_Type);
				MultiFab::Copy(S_old,S_new,ncomps,ncomps,ntracers,1);
			      }
			    else
			      cont_tracer_advection = false;
			  }
		      }
		  }
	      }
	    else
	      {
		continue_iterations = cont  &&  dt_new > t_eps  &&  iter < steady_limit_iterations;

		if (continue_iterations) 
		  {
		    dt_iter = dt_new;
		    for (std::set<int>::const_iterator it=types_advanced.begin(), End=types_advanced.end(); it!=End; ++it) 
		      {
			for (int lev=0; lev<=finest_level; lev++)
			  {
			    getLevel(lev).state[*it].reset();
			  }
		      }
		  }
		else 
		  {
		    if (verbose > 2 && ParallelDescriptor::IOProcessor())
		      std::cout << "MULTILEVEL ADVANCE INNER FLOW STEP failed" << std::endl;
		  }
	      }

	    if (verbose > 2 && ParallelDescriptor::IOProcessor())
	      std::cout << "MULTILEVEL ADVANCE INNER FLOW STEP: Iter " << iter++ 
			<< " Solver Status: " << ret 
			<< " New Time Step: " << dt_iter  << std::endl;
	    
	    if (continue_iterations) 
	      {
		k++;
		dt_iter = std::min(dt_new, t_max - t);

		if (verbose > 2 &&  ParallelDescriptor::IOProcessor())
		  std::cout << "MULTILEVEL ADVANCE SUBCYCLE: Iter: " << k 
			    << " Current time: " << t 
			    << " Current time step: " << dt_iter 
			    << " Next time step " <<  dt_new << std::endl;

		if (k >= k_max  ||  dt_iter <= t_eps) 
		  {
		    delete rs;
		    if (verbose > 2 &&  ParallelDescriptor::IOProcessor())
		      {
			if (k >= k_max) 
			  {
			    std::cout << "MULTILEVEL ADVANCE SUBCYCLE: Too many subcycled steps required."
				      << " Exiting at time = " << t << " (failed to reach " 
				      << time+dt << ")" << std::endl;
			  }
			else
			  {
			    std::cout << "MULTILEVEL ADVANCE SUBCYCLE: Subcycled time step too small."
				      << " Exiting at time = " << t << " (failed to reach " 
				      << time+dt << ")" << std::endl;
			  }
		      }
		    return;
		  }
	      }
	    else
	      {
		// Sub-time step successful
		continue_subtimestep = t + dt_iter < time + dt - t_eps;
		t = std::min(t + dt_iter, time + dt);
		dt_new = std::min(dt_new, time + dt - t);
	      }
	  }
      }
    dt_eig = dt_new;
    delete rs;

    if (saved_states.size()>0) 
      {
	// Restore "old" state
	for (std::set<int>::const_iterator it=types_advanced.begin(), End=types_advanced.end(); it!=End; ++it) 
	  {
	    for (int lev=0; lev<nlevs; ++lev) 
	      {
		PorousMedia& pm = getLevel(lev);
		MultiFab& old = pm.get_old_data(*it);
		MultiFab::Copy(old,saved_states[*it][lev],0,0,old.nComp(),old.nGrow());
		pm.state[*it].setOldTimeLevel(time);
		pm.state[*it].setNewTimeLevel(time + dt);
	      }
	  }
      }
  }
}
#endif

void
PorousMedia::advance_multilevel_saturated (Real time,
					   Real dt)
{
  // 
  // Time stepping for saturated flow with subcycling
  //
  if (level == 0  &&  do_tracer_transport && ntracers > 0)
  {
    int finest_level = parent->finestLevel();
    int nlevs = finest_level + 1;
    
    // Lazily build structure to save state at "time".  If we must subcycle, the
    // algorithm will overwrite old_time data as it goes.  This saved_state
    // must include all state types involved in this subcycle; we make a set
    // of ids and set them manually to minimize the overhead of this
    std::set<int> types_advanced;
    types_advanced.insert(State_Type);
    Array<PArray<MultiFab> > saved_states;
    
    // Subcycle each level independently
    for (int lev=0; lev<nlevs; lev++)
      {
	PorousMedia& fine_lev = getLevel(lev);	

	// Set velocity (u_mac_curr) from bc at t+dt
	fine_lev.set_vel_from_bcs(time+dt,fine_lev.u_mac_curr);

	if (lev == 0) 
	  fine_lev.create_umac_grown(fine_lev.u_mac_curr,
				     fine_lev.u_macG_trac);
	else 
	  {
	    PArray<MultiFab> u_macG_crse(BL_SPACEDIM,PArrayManage);
	    fine_lev.GetCrseUmac(u_macG_crse,time);
	    fine_lev.create_umac_grown(fine_lev.u_mac_curr,u_macG_crse,
				       fine_lev.u_macG_trac); 
	  }
	
        int ltracer = ncomps+ntracers-1;
	
        bool cont_tracer_advection = true;
        Real t_tracer = time;
	
        fine_lev.predictDT(fine_lev.u_macG_trac); // updates dt_eig
        Real dt_sub = std::min(dt_eig,dt);
        Real t_eps = 1.e-8*dt_sub;
        
        MultiFab& S_new = fine_lev.get_new_data(State_Type);
        MultiFab& S_old = fine_lev.get_old_data(State_Type);
        MultiFab Stmp(grids,1,1); 
        for (int i=ncomps;i<=ltracer;i++)
        {
            MultiFab::Copy(Stmp,S_old,i,0,1,1);
            MultiFab::Multiply(Stmp,S_old,0,0,1,1);
            Stmp.divide(S_new,0,1,1);
            MultiFab::Copy(S_old,Stmp,0,i,1,1);
        }
        MultiFab::Copy(S_old,S_new,0,0,1,1);
        
        int num_transport_subcycles = 0;
        std::map<int,MultiFab*> saved_states;
        while (cont_tracer_advection)
        {
            if (num_transport_subcycles > 0) {
                for (std::set<int>::const_iterator it=types_advanced.begin(), End=types_advanced.end(); it!=End; ++it) 
                {
                    const MultiFab& old = fine_lev.get_old_data(*it);                          
                    saved_states[*it] = new MultiFab(old.boxArray(), old.nComp(), old.nGrow());
                    MultiFab::Copy(*saved_states[*it],old,0,0,old.nComp(),old.nGrow());
                }
            }


            PMAmr* p = dynamic_cast<PMAmr*>(parent); BL_ASSERT(p);
            p->SetCumTime(t_tracer);

            fine_lev.setTimeLevel(t_tracer+dt_sub,dt_sub,dt_sub);
            fine_lev.tracer_advection(fine_lev.u_macG_trac,dt_sub,ncomps,ltracer,true);
            t_tracer += dt_sub;
            
            if (std::abs(time + dt - t_tracer) < t_eps) {
                t_tracer = time + dt;
            }
            
            if (t_tracer < time + dt)
            {
                
                // Advance the state data structures
                for (std::set<int>::const_iterator it=types_advanced.begin(), End=types_advanced.end(); it!=End; ++it) 
                {
                    fine_lev.state[*it].allocOldData();
                    fine_lev.state[*it].swapTimeLevels(dt_sub);
                }
                
                dt_sub = std::min(dt_eig,time + dt - t_tracer);
                BL_ASSERT(dt_sub > 0);
                MultiFab& S_new = fine_lev.get_new_data(State_Type);
                MultiFab& S_old = fine_lev.get_old_data(State_Type);
                MultiFab::Copy(S_old,S_new,ncomps,ncomps,ntracers,1);
                
            }
            else {
                cont_tracer_advection = false;
            }
            num_transport_subcycles++;
            
        }
        if (num_transport_subcycles > 1) {
            for (std::set<int>::const_iterator it=types_advanced.begin(), End=types_advanced.end(); it!=End; ++it) 
            {
                MultiFab& old = get_old_data(*it);                          
                MultiFab::Copy(old,*saved_states[*it],0,0,old.nComp(),old.nGrow());
                delete saved_states[*it];
            }
            
            if (verbose > 0 &&  ParallelDescriptor::IOProcessor())
                std::cout << "MULTILEVEL ADVANCE TRANSPORT: # SUBCYCLES: "
                          << num_transport_subcycles << std::endl;
            
        }
        
        // Bring all states up to current time, and reinstate original dt info
        for (std::set<int>::const_iterator it=types_advanced.begin(), End=types_advanced.end(); it!=End; ++it) 
        {
            fine_lev.state[*it].setTimeLevel(time+dt,dt,dt);
        }
      }
  }
  PMAmr* p = dynamic_cast<PMAmr*>(parent); BL_ASSERT(p);
  p->SetCumTime(time); // reset to start of subcycled time period
}

void
PorousMedia::advance_tracer (Real time,
			     Real dt)
{
  // 
  // Time stepping for tracers, assuming steady-state condition. 
  //

  BL_ASSERT(do_tracer_transport);
  BL_ASSERT(ntracers > 0);
    
  int ltracer = ncomps+ntracers-1;
  tracer_advection(u_macG_trac,dt,ncomps,ltracer,true); 
}

void
PorousMedia::create_lambda (Real time) 
{
  // 
  // lambda_T is evaluated at edges.  
  // 

  if (model == model_list["single-phase"] || 
      model == model_list["single-phase-rock"]) 
    {
      for (int dir=0; dir<BL_SPACEDIM; dir++)
	{
	  for (MFIter mfi(lambda[dir]); mfi.isValid(); ++mfi)
	    {
	      const Box& ebox = lambda[dir][mfi].box();
	      lambda[dir][mfi].copy(kpedge[dir][mfi],ebox,0,ebox,0,1);
	    }
	}
    }
  else
    {
      MultiFab& S = get_new_data(State_Type);

      const TimeLevel whichTime = which_time(State_Type,time);
      BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);    
      MultiFab* lcc = (whichTime == AmrOldTime) ? lambda_cc : lambdap1_cc;

      const int*  domlo    = geom.Domain().loVect();
      const int*  domhi    = geom.Domain().hiVect();

      for (FillPatchIterator S_fpi(*this,S,1,time,State_Type,0,ncomps);
	   S_fpi.isValid();
           ++S_fpi)
	{

	  const int i = S_fpi.index();
	  BL_ASSERT(grids[i] == S_fpi.validbox());

	  const int* lo     = S_fpi.validbox().loVect();
	  const int* hi     = S_fpi.validbox().hiVect();

	  const Real* ldat  = (*lcc)[i].dataPtr();
	  const int* l_lo   = (*lcc)[i].loVect();
	  const int* l_hi   = (*lcc)[i].hiVect();

	  const int* lx_lo  = lambda[0][i].loVect();
	  const int* lx_hi  = lambda[0][i].hiVect();
	  const Real* lxdat = lambda[0][i].dataPtr();

	  const int* ly_lo  = lambda[1][i].loVect();
	  const int* ly_hi  = lambda[1][i].hiVect();
	  const Real* lydat = lambda[1][i].dataPtr();

#if(BL_SPACEDIM==3)
	  const int* lz_lo  = lambda[2][i].loVect();
	  const int* lz_hi  = lambda[2][i].hiVect();
	  const Real* lzdat = lambda[2][i].dataPtr();
#endif

	  const int* kx_lo  = kpedge[0][i].loVect();
	  const int* kx_hi  = kpedge[0][i].hiVect();
	  const Real* kxdat = kpedge[0][i].dataPtr();

	  const int* ky_lo  = kpedge[1][i].loVect();
	  const int* ky_hi  = kpedge[1][i].hiVect();
	  const Real* kydat = kpedge[1][i].dataPtr();

#if(BL_SPACEDIM==3)
	  const int* kz_lo  = kpedge[2][i].loVect();
	  const int* kz_hi  = kpedge[2][i].hiVect();
	  const Real* kzdat = kpedge[2][i].dataPtr();
#endif

	  Array<int> bc;
	  bc = getBCArray(State_Type,i,0,1);

	  FORT_MK_MACCOEF (lxdat,ARLIM(lx_lo),ARLIM(lx_hi),
			   lydat,ARLIM(ly_lo),ARLIM(ly_hi),
#if (BL_SPACEDIM==3)
			   lzdat,ARLIM(lz_lo),ARLIM(lz_hi),
#endif
			   kxdat,ARLIM(kx_lo),ARLIM(kx_hi),
			   kydat,ARLIM(ky_lo),ARLIM(ky_hi),
#if (BL_SPACEDIM==3)
			   kzdat,ARLIM(kz_lo),ARLIM(kz_hi),
#endif
			   ldat,ARLIM(l_lo),ARLIM(l_hi),
			   lo,hi,domlo,domhi,bc.dataPtr());
	}
    }
}

void
PorousMedia::mac_project (MultiFab* u_mac, MultiFab* RhoD, Real time)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::mac_project()");

  if (verbose>3 && ParallelDescriptor::IOProcessor())
    std::cout << "... mac_projection at level " << level 
	      << " at time " << time << '\n';
  
  create_lambda(time);

  MultiFab RhoG(grids,1,1); 
  RhoG.setVal(0);
  for (int dir=0; dir < BL_SPACEDIM; dir ++)
    {
      RhoD[dir].setVal(0.0);
      u_mac[dir].setVal(0.0);
    }

  initialize_umac(u_mac,RhoG,RhoD,time);

  const TimeLevel whichTime = which_time(State_Type,time);
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

  MultiFab* phi = 0;
  if (whichTime == AmrOldTime) 
    phi = &get_old_data(Press_Type);
  else if (whichTime == AmrNewTime) 
    phi = &get_new_data(Press_Type);

  // Always start with an initial guess of zero in the interior
  phi->setVal(0.);

  const BCRec& p_bc = desc_lst[Press_Type].getBC(0);

  // Set the boundary conditions *before* we define mac_bndry
  // so the values will end up in mac_bndry
  mac_projector->set_dirichlet_bcs(level, phi, RhoG, p_bc, 
				   press_lo, press_hi);
  phi->FillBoundary();

  PressBndry mac_bndry(grids,1,geom);
  const int src_comp   = 0;
  const int dest_comp  = 0;
  const int num_comp   = 1;
  if (level == 0)
    {
      mac_bndry.setBndryValues(*phi,src_comp,dest_comp,num_comp,p_bc);
    }
  else
    {
      MultiFab CPhi;
      GetCrsePressure(CPhi,time);
      BoxArray crse_boxes = BoxArray(grids).coarsen(crse_ratio);
      const int in_rad     = 0;
      const int out_rad    = 1;
      const int extent_rad = 2;
      BndryRegister crse_br(crse_boxes,in_rad,out_rad,extent_rad,num_comp);
      crse_br.copyFrom(CPhi,extent_rad,src_comp,dest_comp,num_comp);
      mac_bndry.setBndryValues(crse_br,src_comp,*phi,src_comp,
			       dest_comp,num_comp,crse_ratio,p_bc);
    }
  //
  // get source term
  //
  int do_rho_scale = 1;

  MultiFab* forces = 0;

  if (do_source_term)
    {
      forces = new MultiFab(grids,ncomps,0);
      forces->setVal(0.);
      for (MFIter mfi(*forces); mfi.isValid(); ++mfi)
	{
	  int i = mfi.index();
	  getForce((*forces)[mfi],i,0,0,ncomps,time,do_rho_scale);
	}
    }
  const Real strt_time = ParallelDescriptor::second();
  mac_projector->mac_project(level,u_mac,lambda,RhoD,forces,
			     phi,mac_bndry,p_bc);
  
  if (do_source_term)
    delete forces;
    

  if (model != model_list["single-phase"] || 
      model != model_list["single-phase-solid"]) 
    {
      MultiFab* u_phase = new MultiFab[BL_SPACEDIM];
      for (int dir = 0; dir < BL_SPACEDIM; dir++)
	{
	  BoxArray edge_grids(grids);
	  edge_grids.surroundingNodes(dir);
	  u_phase[dir].define(edge_grids,1,0,Fab_allocate);
	  u_phase[dir].setVal(1.e40);
	}

      compute_vel_phase(u_phase,u_mac,time);
      umac_cpy_edge_to_cen(u_phase,Vcr_Type,1);

      delete [] u_phase;

    }

  // compute time spend in mac_project()
  if (verbose > 2)
    {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;
      ParallelDescriptor::ReduceRealMax(run_time,IOProc);
      if (ParallelDescriptor::IOProcessor())
        std::cout << "PorousMedia:mac_project(): lev: " << level
                  << ", time: " << run_time << '\n';
    }
}

void
PorousMedia::initialize_umac (MultiFab* u_mac, MultiFab& RhoG, 
			      MultiFab* RhoD, Real time) 
{

  //
  // u_mac is initilized such that its divergence is 
  //   \nabla \rho g
  // RhoG is initialized such that p + RhoG*\Delta x is 
  //   the hydrostatic pressure 
  // RhoD is initizlized such that its divergence is 
  //   the diffusive term due to variable density formulation.
  //
  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();
  const Real* dx       = geom.CellSize();

  MultiFab& S = get_new_data(State_Type);
  Array<Real> const_diff_coef(ncomps);
  for (int i=0;i<ncomps;i++)
    const_diff_coef[i] = visc_coef[i];
  const TimeLevel whichTime = which_time(State_Type,time);
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
  MultiFab* pc;
  if (have_capillary)
    pc = (whichTime == AmrOldTime) ? pcn_cc : pcnp1_cc;
  else
    {
      pc = new MultiFab(grids,1,1);
      (*pc).setVal(0.);
    }
  MultiFab* lbd = (whichTime == AmrOldTime) ? lambda_cc : lambdap1_cc;
  for (FillPatchIterator S_fpi(*this,S,1,time,State_Type,0,ncomps);
       S_fpi.isValid();
       ++S_fpi)
    {
      const int  i   = S_fpi.index();
      const int* lo  = grids[i].loVect();
      const int* hi  = grids[i].hiVect();
	
      const int* lx_lo  = lambda[0][i].loVect();
      const int* lx_hi  = lambda[0][i].hiVect();
      const Real* lxdat = lambda[0][i].dataPtr();

      const int* ly_lo  = lambda[1][i].loVect();
      const int* ly_hi  = lambda[1][i].hiVect();
      const Real* lydat = lambda[1][i].dataPtr();

      const int* kx_lo  = kpedge[0][i].loVect();
      const int* kx_hi  = kpedge[0][i].hiVect();
      const Real* kxdat = kpedge[0][i].dataPtr();

      const int* ky_lo  = kpedge[1][i].loVect();
      const int* ky_hi  = kpedge[1][i].hiVect();
      const Real* kydat = kpedge[1][i].dataPtr();

      FArrayBox& Sfab   = S_fpi();
      const Real* ndat  = Sfab.dataPtr(); 
      const int*  n_lo  = Sfab.loVect();
      const int*  n_hi  = Sfab.hiVect();
	
      Box bx_mac(u_mac[0][i].box());
      const int* umlo   = bx_mac.loVect();
      const int* umhi   = bx_mac.hiVect();
      const Real* umdat = u_mac[0][i].dataPtr();

      Box by_mac(u_mac[1][i].box());
      const int* vmlo   = by_mac.loVect();
      const int* vmhi   = by_mac.hiVect();
      const Real* vmdat = u_mac[1][i].dataPtr();

      const int* rglo   = RhoG[i].loVect();
      const int* rghi   = RhoG[i].hiVect();
      const Real* rgdat = RhoG[i].dataPtr();

      const int* pclo   = (*pc)[i].loVect();
      const int* pchi   = (*pc)[i].hiVect();
      const Real* pcdat = (*pc)[i].dataPtr();

      const int* lbdlo   = (*lbd)[i].loVect();
      const int* lbdhi   = (*lbd)[i].hiVect();
      const Real* lbddat = (*lbd)[i].dataPtr();

      Box rx_mac(RhoD[0][i].box());
      const int* rxlo   = rx_mac.loVect();
      const int* rxhi   = rx_mac.hiVect();
      const Real* rxdat = RhoD[0][i].dataPtr();

      Box ry_mac(RhoD[1][i].box());
      const int* rylo   = ry_mac.loVect();
      const int* ryhi   = ry_mac.hiVect();
      const Real* rydat = RhoD[1][i].dataPtr();

      const int* p_lo   = (*rock_phi)[i].loVect();
      const int* p_hi   = (*rock_phi)[i].hiVect();
      const Real* pdat  = (*rock_phi)[i].dataPtr();

      Array<int> s_bc;
      s_bc = getBCArray(State_Type,i,0,1);

      Array<int> press_bc;
      press_bc = getBCArray(Press_Type,i,0,1);

#if (BL_SPACEDIM == 2)	
      FORT_INIT_UMAC (umdat,ARLIM(umlo),ARLIM(umhi),
		      vmdat,ARLIM(vmlo),ARLIM(vmhi),
		      pcdat,ARLIM(pclo),ARLIM(pchi),
		      lbddat,ARLIM(lbdlo),ARLIM(lbdhi),
		      lxdat,ARLIM(lx_lo),ARLIM(lx_hi),
		      lydat,ARLIM(ly_lo),ARLIM(ly_hi),
		      kxdat,ARLIM(kx_lo),ARLIM(kx_hi),
		      kydat,ARLIM(ky_lo),ARLIM(ky_hi),
		      rgdat,ARLIM(rglo),ARLIM(rghi),
		      rxdat,ARLIM(rxlo),ARLIM(rxhi),
		      rydat,ARLIM(rylo),ARLIM(ryhi),
		      ndat ,ARLIM(n_lo),ARLIM(n_hi),
		      pdat ,ARLIM(p_lo),ARLIM(p_hi),
		      const_diff_coef.dataPtr(),
		      s_bc.dataPtr(),press_bc.dataPtr(),
		      domain_lo,domain_hi,dx,lo,hi,
		      &wt_lo, &wt_hi,
		      inflow_bc_lo.dataPtr(),inflow_bc_hi.dataPtr());

#elif (BL_SPACEDIM == 3)
      Box bz_mac(u_mac[2][i].box());
      const int* wmlo  = bz_mac.loVect();
      const int* wmhi  = bz_mac.hiVect();
      const Real* wmdat = u_mac[2][i].dataPtr();

      const int* lz_lo  = lambda[2][i].loVect();
      const int* lz_hi  = lambda[2][i].hiVect();
      const Real* lzdat = lambda[2][i].dataPtr();

      const int* kz_lo  = kpedge[2][i].loVect();
      const int* kz_hi  = kpedge[2][i].hiVect();
      const Real* kzdat = kpedge[2][i].dataPtr();

      Box rz_mac(RhoD[2][i].box());
      const int* rzlo  = rz_mac.loVect();
      const int* rzhi  = rz_mac.hiVect();
      const Real* rzdat = RhoD[2][i].dataPtr();

      FORT_INIT_UMAC (umdat,ARLIM(umlo),ARLIM(umhi),
		      vmdat,ARLIM(vmlo),ARLIM(vmhi),
		      wmdat,ARLIM(wmlo),ARLIM(wmhi),
		      pcdat,ARLIM(pclo),ARLIM(pchi),
		      lbddat,ARLIM(lbdlo),ARLIM(lbdhi),
		      lxdat,ARLIM(lx_lo),ARLIM(lx_hi),
		      lydat,ARLIM(ly_lo),ARLIM(ly_hi),
		      lzdat,ARLIM(lz_lo),ARLIM(lz_hi),
		      kxdat,ARLIM(kx_lo),ARLIM(kx_hi),
		      kydat,ARLIM(ky_lo),ARLIM(ky_hi),
		      kzdat,ARLIM(kz_lo),ARLIM(kz_hi),
		      rgdat,ARLIM(rglo),ARLIM(rghi),
		      rxdat,ARLIM(rxlo),ARLIM(rxhi),
		      rydat,ARLIM(rylo),ARLIM(ryhi),
		      rzdat,ARLIM(rzlo),ARLIM(rzhi),
		      ndat,ARLIM(n_lo),ARLIM(n_hi),
		      pdat ,ARLIM(p_lo),ARLIM(p_hi),
		      const_diff_coef.dataPtr(),
		      s_bc.dataPtr(),press_bc.dataPtr(),
		      domain_lo,domain_hi,dx,lo,hi,
		      &wt_lo, &wt_hi,
		      inflow_bc_lo.dataPtr(),inflow_bc_hi.dataPtr());
#endif
    }
    
  RhoG.FillBoundary();

  FArrayBox inflow;
  for (OrientationIter oitr; oitr; ++oitr) {
      Orientation face = oitr();
      if (get_inflow_velocity(face,inflow,time)) {
          int shift = ( face.isHigh() ? -1 : +1 );
          inflow.shiftHalf(face.coordDir(),shift);
          for (MFIter mfi(u_mac[face.coordDir()]); mfi.isValid(); ++mfi) {
              FArrayBox& u = u_mac[face.coordDir()][mfi];
              Box ovlp = inflow.box() & u.box();
	      if (ovlp.ok()) {
                  u.copy(inflow);
              }
          }
      }
  }


  if (!have_capillary)
    delete pc;

}

bool
PorousMedia::get_inflow_velocity(const Orientation& face,
                                 FArrayBox&         ccBndFab,
                                 Real               time)
{
    bool ret = false;
    if (bc_descriptor_map.find(face) != bc_descriptor_map.end()) 
    {
        const Box domain = geom.Domain();
        const int* domhi = domain.hiVect();
        const int* domlo = domain.loVect();
        const Real* dx   = geom.CellSize();

        const BCDesc& bc_desc = bc_descriptor_map[face];
        const Box bndBox = bc_desc.first;
        const Array<int>& face_bc_idxs = bc_desc.second;

        ccBndFab.resize(bndBox,1); ccBndFab.setVal(0);

        for (int i=0; i<face_bc_idxs.size(); ++i) {
            const RegionData& face_bc = bc_array[face_bc_idxs[i]];

            if (face_bc.Type() == "zero_total_velocity") {
                ret = true;
		Array<Real> inflow_tmp = face_bc(time);
		Real inflow_vel = inflow_tmp[0];
                const PArray<Region>& regions = face_bc.Regions();
                for (int j=0; j<regions.size(); ++j)
                {
                    regions[j].setVal(ccBndFab,inflow_vel,0,dx,0);
                }
            }
	}
    }    
    return ret;
}

void
PorousMedia::get_inflow_density(const Orientation& face,
				const RegionData&  face_bc,
				FArrayBox&         ccBndFab,
				Real               time)
{
  if (model == model_list["steady-saturated"]) {
    for (int n=0; n<ncomps; ++n) {
      ccBndFab.setVal(density[n],ccBndFab.box(),n,1);
    }
  }
  else {
    const Real* dx   = geom.CellSize();
    if (face_bc.Type() == "zero_total_velocity") {
        Array<Real> inflow_vel = face_bc(time);
	FArrayBox cdat, ktdat, kdat, vdat;
	const int n_kr_coef = kr_coef->nComp();
	cdat.resize(ccBndFab.box(),n_kr_coef);
	ktdat.resize(ccBndFab.box(),BL_SPACEDIM);
	vdat.resize(ccBndFab.box(),1);
	for (int i=0; i<rocks.size(); ++i)
	  {	  
	    rocks[i].set_constant_krval(cdat,dx);
	    rocks[i].set_constant_kval(ktdat,dx);
	  }
	kdat.resize(ccBndFab.box(),1);
	kdat.copy(ktdat,ccBndFab.box(),face.coordDir(),
		  ccBndFab.box(),0,1);

	const PArray<Region>& regions = face_bc.Regions();
	for (int j=0; j<regions.size(); ++j)
	  regions[j].setVal(vdat,inflow_vel[0],0,dx,0);

	DEF_LIMITS(ccBndFab,s_ptr,s_lo,s_hi);
	DEF_CLIMITS(vdat,v_ptr,v_lo,v_hi);
	DEF_CLIMITS(cdat,c_ptr,c_lo,c_hi);
	DEF_CLIMITS(kdat,k_ptr,k_lo,k_hi);

	int nc = 1;
	FORT_STEADYSTATE_FAB(s_ptr, ARLIM(s_lo),ARLIM(s_hi), 
			     density.dataPtr(),muval.dataPtr(),&ncomps,
			     k_ptr, ARLIM(k_lo),ARLIM(k_hi), 
			     c_ptr, ARLIM(c_lo),ARLIM(c_hi), &n_kr_coef,
			     v_ptr, ARLIM(v_lo),ARLIM(v_hi), 
			     dx, &nc, &gravity);
    }
  }
}


void
PorousMedia::compute_vel_phase (MultiFab* u_phase, MultiFab* u_mac,
				Real time) 
{
  //
  // The phase velocity of component 1 is given by 
  //   v_1 = \lambda_1/\lambda_T ( v_T + \lambda_2 \nabla pc)
  //

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();
  const Real* dx       = geom.CellSize();

  MultiFab& S = get_data(State_Type,time);

  const TimeLevel whichTime = which_time(State_Type,time);
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

  MultiFab* pc;
  if (have_capillary)
    pc = (whichTime == AmrOldTime) ? pcn_cc : pcnp1_cc;
  else
    {
      pc = new MultiFab(grids,1,1);
      (*pc).setVal(0.);
    }

  MultiFab* lbd = (whichTime == AmrOldTime) ? lambda_cc : lambdap1_cc;
    
  for (FillPatchIterator S_fpi(*this,S,1,time,State_Type,0,ncomps);
       S_fpi.isValid();
       ++S_fpi)
    {
      const int  i   = S_fpi.index();
      const int* lo  = grids[i].loVect();
      const int* hi  = grids[i].hiVect();

      const int* kx_lo  = kpedge[0][i].loVect();
      const int* kx_hi  = kpedge[0][i].hiVect();
      const Real* kxdat = kpedge[0][i].dataPtr();

      const int* ky_lo  = kpedge[1][i].loVect();
      const int* ky_hi  = kpedge[1][i].hiVect();
      const Real* kydat = kpedge[1][i].dataPtr();
	
      Box bx_mac(u_mac[0][i].box());
      const int* umlo   = bx_mac.loVect();
      const int* umhi   = bx_mac.hiVect();
      const Real* umdat = u_mac[0][i].dataPtr();

      Box by_mac(u_mac[1][i].box());
      const int* vmlo   = by_mac.loVect();
      const int* vmhi   = by_mac.hiVect();
      const Real* vmdat = u_mac[1][i].dataPtr();

      Box ax_mac(u_phase[0][i].box());
      const int* uplo   = ax_mac.loVect();
      const int* uphi   = ax_mac.hiVect();
      const Real* updat = u_phase[0][i].dataPtr();

      Box ay_mac(u_phase[1][i].box());
      const int* vplo   = ay_mac.loVect();
      const int* vphi   = ay_mac.hiVect();
      const Real* vpdat = u_phase[1][i].dataPtr();

      const int* pclo   = (*pc)[i].loVect();
      const int* pchi   = (*pc)[i].hiVect();
      const Real* pcdat = (*pc)[i].dataPtr();

      const int* lbdlo   = (*lbd)[i].loVect();
      const int* lbdhi   = (*lbd)[i].hiVect();
      const Real* lbddat = (*lbd)[i].dataPtr();

      Array<int> s_bc;
      s_bc = getBCArray(State_Type,i,0,1);

#if (BL_SPACEDIM == 2)	
      FORT_UPHASE (updat,ARLIM(uplo),ARLIM(uphi),
		   vpdat,ARLIM(vplo),ARLIM(vphi),
		   umdat,ARLIM(umlo),ARLIM(umhi),
		   vmdat,ARLIM(vmlo),ARLIM(vmhi),
		   pcdat,ARLIM(pclo),ARLIM(pchi),
		   lbddat,ARLIM(lbdlo),ARLIM(lbdhi),
		   kxdat,ARLIM(kx_lo),ARLIM(kx_hi),
		   kydat,ARLIM(ky_lo),ARLIM(ky_hi),
		   s_bc.dataPtr(),
		   domain_lo,domain_hi,dx,lo,hi);

#elif (BL_SPACEDIM == 3)
      Box bz_mac(u_mac[2][i].box());
      const int* wmlo  = bz_mac.loVect();
      const int* wmhi  = bz_mac.hiVect();
      const Real* wmdat = u_mac[2][i].dataPtr();

      Box az_mac(u_phase[2][i].box());
      const int* wplo  = az_mac.loVect();
      const int* wphi  = az_mac.hiVect();
      const Real* wpdat = u_phase[2][i].dataPtr();

      const int* kz_lo  = kpedge[2][i].loVect();
      const int* kz_hi  = kpedge[2][i].hiVect();
      const Real* kzdat = kpedge[2][i].dataPtr();

      FORT_UPHASE (updat,ARLIM(uplo),ARLIM(uphi),
		   vpdat,ARLIM(vplo),ARLIM(vphi),
		   wpdat,ARLIM(wplo),ARLIM(wphi),
		   umdat,ARLIM(umlo),ARLIM(umhi),
		   vmdat,ARLIM(vmlo),ARLIM(vmhi),
		   wmdat,ARLIM(wmlo),ARLIM(wmhi),
		   pcdat,ARLIM(pclo),ARLIM(pchi),
		   lbddat,ARLIM(lbdlo),ARLIM(lbdhi),
		   kxdat,ARLIM(kx_lo),ARLIM(kx_hi),
		   kydat,ARLIM(ky_lo),ARLIM(ky_hi),
		   kzdat,ARLIM(kz_lo),ARLIM(kz_hi),
		   s_bc.dataPtr(),
		   domain_lo,domain_hi,dx,lo,hi);
#endif
    }

  FArrayBox inflow;
  for (OrientationIter oitr; oitr; ++oitr) {
      Orientation face = oitr();
      if (get_inflow_velocity(face,inflow,time)) {
          int shift = ( face.isHigh() ? -1 : +1 );
          inflow.shiftHalf(face.coordDir(),shift);
          for (MFIter mfi(u_phase[face.coordDir()]); mfi.isValid(); ++mfi) {
              FArrayBox& u = u_phase[face.coordDir()][mfi];
              Box ovlp = inflow.box() & u.box();
	      if (ovlp.ok()) {
                  u.copy(inflow);
              }
          }
      }
  }

  if (!have_capillary)
    delete pc;
}

void
PorousMedia::compute_vel_phase (MultiFab* u_phase, 
				int nc,
				Real time) 
{
  //
  // The phase velocity of component n is given by 
  //   v_n = \lambda_n ( \nabla p_n - \rho \gvec)
  // We are going to assume the p stored in PRESS_TYPE 
  // correspond to p_n.  
  //

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();
  const Real* dx       = geom.CellSize();

  MultiFab& S = get_data(State_Type,time);
  MultiFab& P = get_data(Press_Type,time);

  const TimeLevel whichTime = which_time(State_Type,time);
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

  MultiFab* lbd = (whichTime == AmrOldTime) ? lambda_cc : lambdap1_cc;
    
  for (FillPatchIterator S_fpi(*this,S,1,time,State_Type,0,ncomps);
       S_fpi.isValid();
       ++S_fpi)
    {
      const int  i   = S_fpi.index();
      const int* lo  = grids[i].loVect();
      const int* hi  = grids[i].hiVect();
      
      Box ax_mac(u_phase[0][i].box());
      const int* uplo   = ax_mac.loVect();
      const int* uphi   = ax_mac.hiVect();
      const Real* updat = u_phase[0][i].dataPtr();

      Box ay_mac(u_phase[1][i].box());
      const int* vplo   = ay_mac.loVect();
      const int* vphi   = ay_mac.hiVect();
      const Real* vpdat = u_phase[1][i].dataPtr();

      const int* plo   = P[i].loVect();
      const int* phi   = P[i].hiVect();
      const Real* pdat = P[i].dataPtr();

      const int* lbdlo   = (*lbd)[i].loVect();
      const int* lbdhi   = (*lbd)[i].hiVect();
      const Real* lbddat = (*lbd)[i].dataPtr();

      const int* kx_lo  = kpedge[0][i].loVect();
      const int* kx_hi  = kpedge[0][i].hiVect();
      const Real* kxdat = kpedge[0][i].dataPtr();

      const int* ky_lo  = kpedge[1][i].loVect();
      const int* ky_hi  = kpedge[1][i].hiVect();
      const Real* kydat = kpedge[1][i].dataPtr();

      Array<int> bc = getBCArray(Press_Type,i,0,1);
#if (BL_SPACEDIM == 2)	
      FORT_UPHASE_P (updat,ARLIM(uplo),ARLIM(uphi),
		     vpdat,ARLIM(vplo),ARLIM(vphi),
		     lbddat,ARLIM(lbdlo),ARLIM(lbdhi),
		     pdat,ARLIM(plo),ARLIM(phi),
		     kxdat,ARLIM(kx_lo),ARLIM(kx_hi),
		     kydat,ARLIM(ky_lo),ARLIM(ky_hi),
		     lo,hi,domain_lo,domain_hi,dx,bc.dataPtr());

#elif (BL_SPACEDIM == 3)
      Box az_mac(u_phase[2][i].box());
      const int* wplo  = az_mac.loVect();
      const int* wphi  = az_mac.hiVect();
      const Real* wpdat = u_phase[2][i].dataPtr();

      const int* kz_lo  = kpedge[2][i].loVect();
      const int* kz_hi  = kpedge[2][i].hiVect();
      const Real* kzdat = kpedge[2][i].dataPtr();

      FORT_UPHASE_P (updat,ARLIM(uplo),ARLIM(uphi),
		     vpdat,ARLIM(vplo),ARLIM(vphi),
		     wpdat,ARLIM(wplo),ARLIM(wphi),
		     lbddat,ARLIM(lbdlo),ARLIM(lbdhi),
		     pdat,ARLIM(plo),ARLIM(phi),
		     kxdat,ARLIM(kx_lo),ARLIM(kx_hi),
		     kydat,ARLIM(ky_lo),ARLIM(ky_hi),
		     kzdat,ARLIM(kz_lo),ARLIM(kz_hi),
		     lo,hi,domain_lo,domain_hi,dx,bc.dataPtr());
#endif
    }

  FArrayBox inflow;
  for (OrientationIter oitr; oitr; ++oitr) {
      Orientation face = oitr();
      if (get_inflow_velocity(face,inflow,time)) {
          int shift = ( face.isHigh() ? -1 : +1 );
          inflow.shiftHalf(face.coordDir(),shift);
          for (MFIter mfi(u_phase[face.coordDir()]); mfi.isValid(); ++mfi) {
              FArrayBox& u = u_phase[face.coordDir()][mfi];
              Box ovlp = inflow.box() & u.box();
              if (ovlp.ok()) {
		u.copy(inflow);
              }
          }
      }
  }

}

// =========================================
// Functions related to advection equations.
// =========================================

//
// scalar_advection advects the scalars based on Godunov scheme.
//
void
PorousMedia::scalar_advection (MultiFab* u_macG,
                               Real dt,
                               int  fscalar,
                               int  lscalar,
                               bool reflux_on_this_call)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scalar_advection()");

  if (verbose>3 && ParallelDescriptor::IOProcessor())
  {
    if (reflux_on_this_call) 
      std::cout << "... advect scalars with contribution to refluxing \n";
    else 
      std::cout << "... advect scalars\n";
  }
  
  //
  // Get simulation parameters.
  //
  const Real* dx        = geom.CellSize();
  const Real  prev_time = state[State_Type].prevTime();
  const Real  curr_time = state[State_Type].curTime();
  int nscal             = lscalar - fscalar + 1;

  //
  // Get the viscous terms.
  // For model 0-1 => diffusion term
  //             2 => capillary pressure term
  //
  MultiFab visc_terms(grids,nscal,1);
  visc_terms.setVal(0);
  int do_visc_terms = 1;
  if (be_cn_theta != 1.0 && (do_visc_terms || have_capillary) && !do_cpl_advect)
      getViscTerms(visc_terms,fscalar,nscal,prev_time);

  //
  // Divergence of velocity: set to zero for now.
  //
  MultiFab* divu_fp = new MultiFab(grids,1,1);
  (*divu_fp).setVal(0.);
  //
  // Set up the grid loop.
  //
  FArrayBox flux[BL_SPACEDIM], tforces, pctmp, phitmp, kappatmp;

  Array<int> state_bc;

  // S_new is only used as a container to hold
  // time t^{n+1} inflow boundary conditions
  MultiFab& S_new = get_new_data(State_Type); 
  for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
    setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,fscalar,fscalar,nscal);
  }

  MultiFab fluxes[BL_SPACEDIM];
  
  if (reflux_on_this_call && do_reflux && level < parent->finestLevel())
    {
      for (int i = 0; i < BL_SPACEDIM; i++)
	{
	  BoxArray ba = grids;
	  ba.surroundingNodes(i);
	  fluxes[i].define(ba, nscal, 0, Fab_allocate);
	}
    }
  
  for (FillPatchIterator S_fpi(*this,get_old_data(State_Type),HYP_GROW,
			       prev_time,State_Type,fscalar,nscal);
       S_fpi.isValid();
       ++S_fpi)
    {
      const int i = S_fpi.index();
      
      getForce(tforces,i,1,fscalar,nscal,curr_time);
      godunov->Setup(grids[i], flux[0], flux[1], 
#if (BL_SPACEDIM == 3)  
		     flux[2], 
#endif
		     nscal,model);	   
      
      Real eigmax_m[BL_SPACEDIM] = {D_DECL(-1.e20,-1.e20,-1.e20)};
      
      int state_ind = 0;
      int use_conserv_diff = (advectionType[state_ind] == Conservative);
      
      godunov->Sum_tf_divu_visc(S_fpi(),tforces,state_ind,nscal,
				visc_terms[i],state_ind,
				(*divu_fp)[i],use_conserv_diff);
      
      state_bc = getBCArray(State_Type,i,state_ind,1);
      
      //
      // Polymer model.
      //
      if (model == model_list["polymer"]) 
	{ 
	  godunov->AdvectStatePmr(grids[i], dx, dt, 
				  area[0][i], u_macG[0][i], flux[0], kpedge[0][i],
				  area[1][i], u_macG[1][i], flux[1], kpedge[1][i],
#if (BL_SPACEDIM == 3)                        
				  area[2][i], u_macG[2][i], flux[2], kpedge[2][i],
#endif
				  S_fpi(),S_new[i],tforces,
				  (*divu_fp)[i] , state_ind,
				  (*aofs)[i]    , state_ind,
				  (*rock_phi)[i], (*kappa)[i],
				  use_conserv_diff,
				  state_ind,state_bc.dataPtr(),volume[i],
				  nscal,gravity,eigmax_m);
	}

      //
      // Single phase  model.
      //
      else if (model == model_list["single-phase"] || 
	       model == model_list["single-phase-solid"])
	{
	  godunov->AdvectStateLin(grids[i], dx, dt, 
				  area[0][i], u_macG[0][i], flux[0],
				  area[1][i], u_macG[1][i], flux[1], 
#if (BL_SPACEDIM == 3)                        
				  area[2][i], u_macG[2][i], flux[2],
#endif
				  S_fpi(),S_new[i],tforces, state_ind,
				  (*aofs)[i]    , state_ind,
				  (*rock_phi)[i], state_ind,
				  state_bc.dataPtr(),volume[i],nscal);	
	}
      //
      // Two-phase two-component model.
      //
      else if (model == model_list["two-phase"])
	{
	  const int n_kr_coef = kr_coef->nComp();
	  if (do_cpl_advect) 
	    {
	      Box box = (*pcn_cc)[i].box();
	      pctmp.resize(box,1);
	      pctmp.copy((*pcn_cc)[i],box,0,box,0,1);
	      pctmp.plus((*pcnp1_cc)[i],box,0,0,1);
	      pctmp.mult(0.5);
	      godunov->AdvectStateCpl(grids[i], dx, dt, 
				      area[0][i], u_macG[0][i], flux[0], kpedge[0][i], lambda[0][i],
				      area[1][i], u_macG[1][i], flux[1], kpedge[1][i], lambda[1][i],
#if (BL_SPACEDIM == 3)                        
				      area[2][i], u_macG[2][i], flux[2], kpedge[2][i], lambda[2][i],
#endif
				      S_fpi(), S_new[i], tforces,
				      (*divu_fp)[i] , state_ind,
				      (*aofs)[i]    , state_ind,
				      (*rock_phi)[i], (*kappa)[i],  pctmp,
				      (*lambda_cc)[i],(*dlambda_cc)[i], 
				      (*kr_coef)[i],n_kr_coef,
				      use_conserv_diff,
				      state_ind,state_bc.dataPtr(),volume[i],nscal);
	    }
	  else
	    godunov->AdvectStateRmn(grids[i], dx, dt, 
				    area[0][i], u_macG[0][i], flux[0], kpedge[0][i],
				    area[1][i], u_macG[1][i], flux[1], kpedge[1][i],
#if (BL_SPACEDIM == 3)                        
				    area[2][i], u_macG[2][i], flux[2], kpedge[2][i],
#endif
				    S_fpi(),S_new[i],tforces,
				    (*divu_fp)[i] , state_ind,
				    (*aofs)[i]    , state_ind,
				    (*rock_phi)[i], (*kappa)[i], 
				    (*lambda_cc)[i],(*dlambda_cc)[i], 
				    (*kr_coef)[i],n_kr_coef,
				    use_conserv_diff,
				    state_ind,state_bc.dataPtr(),volume[i],nscal);
	
	}

      //
      // Set aofs of components in solid phase to zero.
      //
      if ((model == model_list["single-phase-solid"]) && (nphases > 1)) 
	{
	  for (int ii=0; ii<ncomps; ii++) 
	    {
	      if (solid.compare(pNames[pType[ii]]) == 0)
		(*aofs)[i].setVal(0.,ii);	      
	    }
	}
      
      if (reflux_on_this_call)
        {
	  if (do_reflux)
	    {
	      if (level < parent->finestLevel())
		{
		  for (int d = 0; d < BL_SPACEDIM; d++)
		    fluxes[d][i].copy(flux[d]);
		}
	      
	      if (level > 0)
		{
		  for (int d = 0; d < BL_SPACEDIM; d++)
		    advflux_reg->FineAdd(flux[d],d,i,0,state_ind,nscal,dt);
		}
	    }
        }
      
      //
      // Allocate the eigenvalues into scalar array.
      //
      
      if (model == model_list["two-phase"])
	{
	  S_new[i].setVal(0.0,ncomps+ntracers);
	  S_new[i].setVal(0.0,ncomps+ntracers+1);
	  godunov->Getdfdn(S_new[i],ncomps+ntracers,ncomps,0,4);
	  godunov->Getdfdn(S_new[i],ncomps+ntracers+1,ncomps,1,4);
	}
    } 
  for (int d = 0; d < BL_SPACEDIM; d++)
    lambda[d].FillBoundary();

  //MultiFab::Copy(S_new,*aofs,0,ncomps+ntracers,1,0);
  //MultiFab::Copy(S_new,visc_terms,0,ncomps+ntracers+1,1,0);

  delete divu_fp;
  
  if (do_reflux && level < parent->finestLevel() && reflux_on_this_call)
    {
      for (int d = 0; d < BL_SPACEDIM; d++)
	getAdvFluxReg(level+1).CrseInit(fluxes[d],d,0,0,nscal,-dt);
    }
}

void
PorousMedia::scalar_update (Real      dt,
                            int       first_scalar,
                            int       num_comp,
			    int       corrector,
			    MultiFab* u_mac)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scalar_update()");
  if (verbose > 3 && ParallelDescriptor::IOProcessor())
    std::cout << "... update scalars\n";

  int last_scalar = num_comp-1;

  scalar_advection_update(dt, first_scalar, last_scalar, corrector);
  if (do_any_diffuse) 
    scalar_diffusion_update(dt, first_scalar, last_scalar, corrector);
  if (have_capillary)
    {
      if (!do_cpl_advect)
	scalar_capillary_update(dt, corrector, u_mac);
      else
	{
	  if (do_cpl_advect == 2)
	    {
	      Real pcTime = state[State_Type].curTime();
	      calcCapillary(pcTime);
	    }
	  else
	    diff_capillary_update(dt, corrector, u_mac);
	}
    }

  if (idx_dominant > -1)
    scalar_adjust_constraint(first_scalar,last_scalar);
}

void
PorousMedia::scalar_advection_update (Real dt,
                                      int  first_scalar,
                                      int  last_scalar,
				      int  corrector)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scalar_advection_update()");

  MultiFab&  S_old    = get_old_data(State_Type);
  MultiFab&  S_new    = get_new_data(State_Type);
  MultiFab&  Aofs     = *aofs;
  MultiFab&  Rockphi  = *rock_phi;
  FArrayBox  tforces;
    
  // Compute inviscid estimate of scalars.
  // component first_scalar -> last_scalar: N
  Real pcTime = state[State_Type].curTime();
  int nscal = last_scalar - first_scalar + 1;
  for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
      const int i = mfi.index();
      getForce(tforces,i,0,first_scalar,nscal,pcTime);
      godunov->Add_aofs_tf(S_old[i],S_new[i],first_scalar,nscal,
			   Aofs[i],first_scalar,tforces,0,Rockphi[i],grids[i],dt);
    }


  FillStateBndry(pcTime,State_Type,first_scalar,nscal);
  S_new.FillBoundary();
  if (idx_dominant > -1 && last_scalar < ncomps)
    scalar_adjust_constraint(first_scalar,last_scalar);

  //
  // Write out the min and max of each component of the new state.
  //
  if (corrector || verbose > 3) check_minmax();
}

void
PorousMedia::tracer_advection_update (Real dt,
                                      int  first_scalar,
                                      int  last_scalar,
				      int  corrector)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::tracer_advection_update()");

  BL_ASSERT(do_tracer_transport);

  MultiFab&  S_old    = get_old_data(State_Type);
  MultiFab&  S_new    = get_new_data(State_Type);
  MultiFab&  Aofs     = *aofs;
  MultiFab&  Rockphi  = *rock_phi;
  FArrayBox  tforces;
    
  int nscal = ncomps + ntracers;
    
  //
  // Advect only the Total
  //
  const Array<int>& idx_total = group_map["Total"];
  if (ntracers > 0) 
  {
      if (idx_total.size()) 
      {
          Real pcTime = state[State_Type].curTime();
          for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
          {
              const int i = mfi.index();
              getForce_Tracer(tforces,i,0,0,ntracers,pcTime);
              
              godunov->Add_aofs_tracer(S_old[i],S_new[i],0,nscal,
                                       Aofs[i],0,tforces,0,Rockphi[i],grids[i],
                                       idx_total,dt);
              
          }
      }
      else {
          MultiFab::Copy(S_new,S_old,ncomps,ncomps,ntracers,0);
      }
  }
  S_new.FillBoundary();
  //
  // Write out the min and max of each component of the new state.
  //
  if (corrector || verbose > 3) check_minmax(first_scalar,last_scalar);
}

void
PorousMedia::scalar_diffusion_update (Real dt,
				      int first_scalar,
				      int last_scalar,
				      int corrector)

{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scalar_diffusion_update()");

  if (verbose > 2 && ParallelDescriptor::IOProcessor())
    std::cout << "... diffuse scalars\n";

  BL_ASSERT(model == model_list["single-phase"] || model == model_list["single-phase-solid"]);

  const Real strt_time = ParallelDescriptor::second();

  // Build single component edge-centered array of MultiFabs for fluxes
  MultiFab** fluxSCn;
  MultiFab** fluxSCnp1;

  diffusion->allocFluxBoxesLevel(fluxSCn,  0,1);
  diffusion->allocFluxBoxesLevel(fluxSCnp1,0,1);


  MultiFab* rho;
  MultiFab& S_new = get_new_data(State_Type);
  rho = new MultiFab(grids,1,1);
  MultiFab::Copy(*rho,S_new,0,0,1,1);

  for (int kk = 1; kk<ncomps; kk++)
    {
      if (solid.compare(pNames[pType[kk]]) != 0) 
	MultiFab::Add(*rho,S_new,kk,0,1,1);
    }

  diffusion->set_rho(rho);

  for (int kk = first_scalar; kk <= last_scalar; kk++)
    {
      if (is_diffusive[kk])
        {
	  MultiFab*  delta_rhs   = 0;
	  MultiFab*  alpha       = 0;
	  MultiFab** cmp_diffn   = 0;
	  MultiFab** cmp_diffnp1 = 0;

	  alpha     = new MultiFab(grids, 1, 1);
	  MultiFab::Copy(*alpha,*rock_phi,0,0,1,alpha->nGrow());

	  if (variable_scal_diff)
            {
	      Real diffTime = state[State_Type].prevTime();
	      diffusion->allocFluxBoxesLevel(cmp_diffn, 0, 1);
	      getDiffusivity(cmp_diffn, diffTime, kk, 0, 1);

	      diffTime = state[State_Type].curTime();
	      diffusion->allocFluxBoxesLevel(cmp_diffnp1, 0, 1);
	      getDiffusivity(cmp_diffnp1, diffTime, kk, 0, 1);
            }
	    
	  diffusion->diffuse_scalar(dt,kk,be_cn_theta,
				    fluxSCn,fluxSCnp1,0,delta_rhs,
				    alpha,cmp_diffn,cmp_diffnp1);

	  if (variable_scal_diff)
            {
	      diffusion->removeFluxBoxesLevel(cmp_diffn);
	      diffusion->removeFluxBoxesLevel(cmp_diffnp1);
            }
	    
	  delete delta_rhs;
	  delete alpha;

	  //
	  // Increment the viscous flux registers
	  //
	  if (do_reflux && corrector)
            {
	      FArrayBox fluxtot;

	      for (int d = 0; d < BL_SPACEDIM; d++)
                {
		  MultiFab fluxes;

		  if (level < parent->finestLevel())
		    fluxes.define((*fluxSCn[d]).boxArray(), 1, 0, Fab_allocate);

		  for (MFIter fmfi(*fluxSCn[d]); fmfi.isValid(); ++fmfi)
                    {
		      const Box& ebox = (*fluxSCn[d])[fmfi].box();

		      fluxtot.resize(ebox,1);
		      fluxtot.copy((*fluxSCn[d])[fmfi],ebox,0,ebox,0,1);
		      fluxtot.plus((*fluxSCnp1[d])[fmfi],ebox,0,0,1);

		      if (level < parent->finestLevel())
			fluxes[fmfi].copy(fluxtot);

		      if (level > 0)
			getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),0,kk,1,dt);
                    }

		  if (level < parent->finestLevel())
		    getLevel(level+1).getViscFluxReg().CrseInit(fluxes,d,0,kk,1,-dt);
                }
            }
        }
    }

  delete rho;

  diffusion->removeFluxBoxesLevel(fluxSCn);
  diffusion->removeFluxBoxesLevel(fluxSCnp1);
    
  // Make sure values on bc is correct
  Real pcTime = state[State_Type].curTime();
  FillStateBndry(pcTime,State_Type,0,ncomps);

  if (verbose > 2)
  {
    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "PorousMedia::scalar_diffusion_update(): time: " << run_time << '\n';
  }
    
  // Write out the min and max of each component of the new state
  if (corrector && verbose>3) check_minmax();
    
}

void
PorousMedia::diffuse_adjust_dominant(MultiFab&              Phi_new,
				     int                    sComp,
				     Real                   dt,
				     MultiFab**             fluxn,
				     MultiFab**             fluxnp1)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::diffuse_adjust_dominant()");

  FArrayBox update;
  FArrayBox tmpfab;
  int nscal = 1;
  for (MFIter mfi(Phi_new); mfi.isValid(); ++mfi)
    {
      int iGrid = mfi.index();

      const Box& box = mfi.validbox();
      const int* lo    = mfi.validbox().loVect();
      const int* hi    = mfi.validbox().hiVect();

      update.resize(box,1);
      tmpfab.resize(box,1);
      tmpfab.setVal(0.);

      const int* p_lo  = (*rock_phi)[mfi].loVect();
      const int* p_hi  = (*rock_phi)[mfi].hiVect();
      const Real* pdat = (*rock_phi)[mfi].dataPtr();


      FORT_RECOMP_UPDATE(lo, hi,
			 update.dataPtr(),
			 ARLIM(update.loVect()),ARLIM(update.hiVect()),
			 pdat, ARLIM(p_lo), ARLIM(p_hi),
			 (*fluxn[0])[iGrid].dataPtr(),
			 ARLIM((*fluxn[0])[iGrid].loVect()),
			 ARLIM((*fluxn[0])[iGrid].hiVect()),
			 (*fluxn[1])[iGrid].dataPtr(),
			 ARLIM((*fluxn[1])[iGrid].loVect()),
			 ARLIM((*fluxn[1])[iGrid].hiVect()),
#if BL_SPACEDIM == 3
			 (*fluxn[2])[iGrid].dataPtr(),
			 ARLIM((*fluxn[2])[iGrid].loVect()),
			 ARLIM((*fluxn[2])[iGrid].hiVect()),
#endif
			 volume[iGrid].dataPtr(),
			 ARLIM(volume[iGrid].loVect()),ARLIM(volume[iGrid].hiVect()),
			 &nscal);

      update.mult(dt,box,0,1);
      tmpfab.plus(update,box,0,0,1);

      if (fluxnp1 != 0) 
	{
	  FORT_RECOMP_UPDATE(lo,hi,
			     update.dataPtr(),
			     ARLIM(update.loVect()),ARLIM(update.hiVect()),
			     pdat, ARLIM(p_lo), ARLIM(p_hi),
			     (*fluxnp1[0])[iGrid].dataPtr(),
			     ARLIM((*fluxnp1[0])[iGrid].loVect()),
			     ARLIM((*fluxnp1[0])[iGrid].hiVect()),
			     (*fluxnp1[1])[iGrid].dataPtr(),
			     ARLIM((*fluxnp1[1])[iGrid].loVect()),
			     ARLIM((*fluxnp1[1])[iGrid].hiVect()),
#if BL_SPACEDIM == 3
			     (*fluxnp1[2])[iGrid].dataPtr(),
			     ARLIM((*fluxnp1[2])[iGrid].loVect()),
			     ARLIM((*fluxnp1[2])[iGrid].hiVect()),
#endif
			     volume[iGrid].dataPtr(),
			     ARLIM(volume[iGrid].loVect()),ARLIM(volume[iGrid].hiVect()),
			     &nscal);

	  update.mult(dt,box,0,1);
	  tmpfab.plus(update,box,0,0,1);
	}

      tmpfab.plus(Phi_new[iGrid],box,sComp,0,1);
      Phi_new[mfi].copy(tmpfab,box,0,box,sComp,1);
    }

}

//
// This routine advects the scalars
//
void
PorousMedia::tracer_advection (MultiFab* u_macG,
                               Real dt,
                               int  fscalar,
                               int  lscalar,
                               bool reflux_on_this_call)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::tracer_advection()");

  BL_ASSERT(do_tracer_transport);

  if (verbose > 2 && ParallelDescriptor::IOProcessor())
  {
    std::cout << "... advect tracers\n";
  }

  //
  // Get simulation parameters.
  //
  const Real* dx        = geom.CellSize();
  const Real  prev_time = state[State_Type].prevTime();
  const Real  cur_time  = state[State_Type].curTime();
  int nscal             = ntracers;
    
  //
  // Get the viscous terms.
  //
  MultiFab visc_terms(grids,nscal,1);
  visc_terms.setVal(0);

  //
  // Set up the grid loop.
  //
  FArrayBox flux[BL_SPACEDIM], tforces;

  FArrayBox sat, satn;

  Array<int> state_bc;

  MultiFab* divu_fp = new MultiFab(grids,1,1);
  (*divu_fp).setVal(0.);

  MultiFab fluxes[BL_SPACEDIM];
  if (reflux_on_this_call && do_reflux && level < parent->finestLevel())
    {
      for (int i = 0; i < BL_SPACEDIM; i++)
        {
	  BoxArray ba = grids;
	  ba.surroundingNodes(i);
	  fluxes[i].define(ba, nscal, 0, Fab_allocate);
        }
    }
 
  for (FillPatchIterator S_fpi(*this,get_old_data(State_Type),HYP_GROW,
			       prev_time,State_Type,fscalar,nscal),
	 Sn_fpi(*this,get_new_data(State_Type),HYP_GROW,
		cur_time,State_Type,fscalar,nscal),
	 St_fpi(*this,get_old_data(State_Type),HYP_GROW,  
		prev_time,State_Type,0,ncomps),
	 Stn_fpi(*this,get_new_data(State_Type),HYP_GROW,  
		 cur_time,State_Type,0,ncomps);
       S_fpi.isValid() && Sn_fpi.isValid() && St_fpi.isValid() && Stn_fpi.isValid(); 
       ++S_fpi,++Sn_fpi,++St_fpi,++Stn_fpi)
    {
      const int i = S_fpi.index();
      getForce_Tracer(tforces,i,1,fscalar,nscal,cur_time);

      godunov->Setup_tracer(grids[i], flux[0], flux[1],
#if (BL_SPACEDIM == 3)  
			    flux[2], 	    
#endif		 
			    nscal);

      int aofs_ind  = ncomps;
      int state_ind = 0;
      int use_conserv_diff = (advectionType[state_ind] == Conservative);

      godunov->Sum_tf_divu_visc(S_fpi(),tforces,state_ind,nscal,
				visc_terms[i],state_ind,
				(*divu_fp)[i],use_conserv_diff);
      
      state_bc = getBCArray(State_Type,i,state_ind,1);

      sat.resize(BoxLib::grow(grids[i],HYP_GROW),1);
      satn.resize(BoxLib::grow(grids[i],HYP_GROW),1);
      sat.copy(St_fpi(),0,0,1);
      satn.copy(Stn_fpi(),0,0,1);
      sat.mult(1.0/density[0]);
      satn.mult(1.0/density[0]);
      godunov->AdvectTracer(grids[i], dx, dt, 
			    area[0][i], u_macG[0][i], flux[0], 
			    area[1][i], u_macG[1][i], flux[1], 
#if (BL_SPACEDIM == 3)                        
			    area[2][i], u_macG[2][i], flux[2], 
#endif
			    S_fpi(), Sn_fpi(), sat, satn, tforces,
			    (*divu_fp)[i] , state_ind,
			    (*aofs)[i]    , aofs_ind,
			    (*rock_phi)[i], 
			    use_conserv_diff,
			    state_ind,state_bc.dataPtr(),volume[i],
			    nscal);

      if (reflux_on_this_call)
	{
	  if (do_reflux)
	    {
	      if (level < parent->finestLevel())
		{
		  for (int d = 0; d < BL_SPACEDIM; d++)
		    fluxes[d][i].copy(flux[d]);
		}

	      if (level > 0)
		{
		  for (int d = 0; d < BL_SPACEDIM; d++)
		    advflux_reg->FineAdd(flux[d],d,i,0,fscalar,nscal,dt);
		}
	    }
	}
    }

  delete divu_fp;

  if (do_reflux && level < parent->finestLevel() && reflux_on_this_call)
    {
      for (int d = 0; d < BL_SPACEDIM; d++)
	getAdvFluxReg(level+1).CrseInit(fluxes[d],d,0,fscalar,nscal,-dt);
    }

  int corrector = 1;
  tracer_advection_update (dt, fscalar, lscalar, corrector);
}

DistributionMapping
PorousMedia::getFuncCountDM (const BoxArray& bxba, int ngrow)
{
  //
  // Sometimes "mf" is the valid region of the State.
  // Sometimes it's the region covered by AuxBoundaryData.
  // When ngrow>0 were doing AuxBoundaryData with nGrow()==ngrow.
  // Taken from LMC/HeatTransfer.cpp
  //
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getFuncCountDM()");

  DistributionMapping rr;
  rr.RoundRobinProcessorMap(bxba.size(),ParallelDescriptor::NProcs());

  MultiFab fctmpnew;
  fctmpnew.define(bxba, 1, 0, rr, Fab_allocate);
  fctmpnew.setVal(1);

  if (ngrow == 0)
    {
      //
      // Working on valid region of state.
      //
      fctmpnew.copy(get_new_data(FuncCount_Type));  // Parallel copy.
    }
  else
    {
      //
      // Can't directly use a parallel copy from FuncCount_Type to fctmpnew.
      //
      MultiFab& FC = get_new_data(FuncCount_Type);

      BoxArray ba = FC.boxArray();
      ba.grow(ngrow);
      MultiFab grownFC(ba, 1, 0);
      grownFC.setVal(1);
                
      for (MFIter mfi(FC); mfi.isValid(); ++mfi)
	grownFC[mfi].copy(FC[mfi]);

      fctmpnew.copy(grownFC);  // Parallel copy.
    }

  int count = 0;
  Array<long> vwrk(bxba.size());
  for (MFIter mfi(fctmpnew); mfi.isValid(); ++mfi)
    vwrk[count++] = static_cast<long>(fctmpnew[mfi].sum(0));

  fctmpnew.clear();

#if BL_USE_MPI
  const int IOProc = ParallelDescriptor::IOProcessorNumber();

  Array<int> nmtags(ParallelDescriptor::NProcs(),0);
  Array<int> offset(ParallelDescriptor::NProcs(),0);

  for (int i = 0; i < vwrk.size(); i++)
    nmtags[rr.ProcessorMap()[i]]++;

  BL_ASSERT(nmtags[ParallelDescriptor::MyProc()] == count);

  for (int i = 1; i < offset.size(); i++)
    offset[i] = offset[i-1] + nmtags[i-1];

  Array<long> vwrktmp = vwrk;

  MPI_Gatherv(vwrk.dataPtr(),
	      count,
	      ParallelDescriptor::Mpi_typemap<long>::type(),
	      vwrktmp.dataPtr(),
	      nmtags.dataPtr(),
	      offset.dataPtr(),
	      ParallelDescriptor::Mpi_typemap<long>::type(),
	      IOProc,
	      ParallelDescriptor::Communicator());

  if (ParallelDescriptor::IOProcessor())
    {
      //
      // We must now assemble vwrk in the proper order.
      //
      std::vector< std::vector<int> > table(ParallelDescriptor::NProcs());

      for (int i = 0; i < vwrk.size(); i++)
	table[rr.ProcessorMap()[i]].push_back(i);

      int idx = 0;
      for (int i = 0; i < table.size(); i++)
	for (int j = 0; j < table[i].size(); j++)
	  vwrk[table[i][j]] = vwrktmp[idx++]; 
    }
  //
  // Send the properly-ordered vwrk to all processors.
  //
  ParallelDescriptor::Bcast(vwrk.dataPtr(), vwrk.size(), IOProc);
#endif

  DistributionMapping res;
  //
  // This call doesn't invoke the MinimizeCommCosts() stuff.
  //
  res.KnapSackProcessorMap(vwrk,ParallelDescriptor::NProcs());

  return res;
}

#if defined(AMANZI)
static
void
TagUnusedGrowCells(MultiFab&    state, 
		   int          state_idx,
		   const BCRec& bc,
		   PorousMedia& pm, 
		   int          ngrow, 
		   Real         tagVal,
		   int          comp,
		   int          nComp)
{
  const BoxArray& ba = state.boxArray();
  const Geometry& geom = pm.Geom();
  const Box& domain = geom.Domain();
  Array<Orientation> Faces;
  pm.getDirichletFaces(Faces,state_idx,bc);

  BoxList bl_unused;
  for (int iface = 0; iface < Faces.size(); iface++)
    {
      const Orientation& face = Faces[iface];
      Box bnd = BoxLib::adjCell(domain,face,ngrow);
      int dir = face.coordDir();
      for (int d=0; d<BL_SPACEDIM; ++d) 
	{
	  bnd.grow(d,ngrow);
	}
      bl_unused.join(BoxLib::boxDiff(bnd,BoxLib::adjCell(domain,face,1)));
    }
  BoxLib::removeOverlap(bl_unused);
  BoxArray ba_unused(bl_unused);

  for (MFIter mfi(state); mfi.isValid(); ++mfi) 
    {
      FArrayBox& fab = state[mfi];
      const Box& box = fab.box();
      std::vector< std::pair<int,Box> > isects = ba_unused.intersections(box);
      for (int ii = 0, N = isects.size(); ii < N; ii++)
	{
	  fab.setVal(tagVal,isects[ii].second,comp,nComp);
	}
    }
}

static
BoxArray
ChemistryGrids (const MultiFab& state,
                const Amr*      parent,
                int             level,
                int             ngrow)
{
    //
    // Let's chop the grids up a bit.
    //
    // We want to try and level out the chemistry work.
    //
    const int NProcs = ParallelDescriptor::NProcs();

    BoxArray ba = state.boxArray();

    if (ngrow>0) {
        BoxList bl = BoxList(ba).accrete(ngrow);
        ba = BoxArray(BoxLib::removeOverlap(bl));
    }

    bool done = false;

    for (int cnt = 1; !done; cnt *= 2)
    {
        const int ChunkSize = parent->maxGridSize(level)/cnt;

        if (ChunkSize < 16)
            //
            // Don't let grids get too small. 
            //
            break;

        IntVect chunk(D_DECL(ChunkSize,ChunkSize,ChunkSize));

        for (int j = 0; j < BL_SPACEDIM && ba.size() < 3*NProcs; j++)
        {
            chunk[j] /= 2;

            ba.maxSize(chunk);

            if (ba.size() >= 3*NProcs) done = true;
        }
    }

    return ba;
}
#endif

//
// ODE-solve for chemistry: cell-by-cell
//
void
PorousMedia::strang_chem (Real time,
			  Real dt,
                          int  ngrow)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::strang_chem()");
  const Real strt_time = ParallelDescriptor::second();

  const TimeLevel whichTime = which_time(State_Type,time);
  MultiFab& S    = (whichTime == AmrOldTime ? get_old_data(State_Type)     : get_new_data(State_Type));
  MultiFab& Fcnt = (whichTime == AmrOldTime ? get_old_data(FuncCount_Type) : get_new_data(FuncCount_Type));
  //MultiFab& Aux  = (whichTime == AmrOldTime ? get_old_data(Aux_Chem_Type)  : get_new_data(Aux_Chem_Type));
  MultiFab& Aux  = get_new_data(Aux_Chem_Type);

  // FIXME: Curently have no code to set Aux_Chem_Type data at the physical boundaries, so we just punt
  ngrow = 0;

#if defined(AMANZI)
  //
  // ngrow == 0 -> we're working on the valid region of state.
  //
  // ngrow > 0  -> we're working on aux_boundary_data_old with that many grow cells.
  //
  int tnum = 1;
#ifdef _OPENMP
  tnum = omp_get_max_threads();
#endif

  BL_ASSERT(S.nComp() >= ncomps+ntracers);
  for (int ithread = 0; ithread < tnum; ithread++)
    {
      BL_ASSERT(components[ithread].mineral_volume_fraction.size() == nminerals);
      BL_ASSERT(components[ithread].mineral_specific_surface_area.size() == nminerals);
      BL_ASSERT(components[ithread].total.size() == ntracers);
      BL_ASSERT(components[ithread].free_ion.size() == ntracers);
      if (using_sorption) {
	BL_ASSERT(components[ithread].total_sorbed.size() == ntracers);
      }
      if (nsorption_isotherms > 0) {
	BL_ASSERT(components[ithread].isotherm_kd.size() == ntracers);
	BL_ASSERT(components[ithread].isotherm_langmuir_b.size() == ntracers);
	BL_ASSERT(components[ithread].isotherm_freundlich_n.size() == ntracers);
      }
      BL_ASSERT(components[ithread].ion_exchange_sites.size() == 0);
    }

  //
  // Assume we are always doing funccount.
  //
  BoxArray            ba = ChemistryGrids(S, parent, level, ngrow);
  DistributionMapping dm = getFuncCountDM(ba,ngrow);

  if (verbose > 2 && ParallelDescriptor::IOProcessor())
    {
      if (ngrow == 0)
	std::cout << "*** strang_chem: FABs in tmp MF covering valid region: " << ba.size() << std::endl;
      else
	std::cout << "*** strang_chem: FABs in tmp MF covering aux_boundary_data_old: " << ba.size() << std::endl;
    }
  
  MultiFab stateTemp, phiTemp, volTemp, fcnCntTemp, auxTemp;

  stateTemp.define(ba, S.nComp(), 0, dm, Fab_allocate);
  auxTemp.define(ba, Aux.nComp(), 0, dm, Fab_allocate);

  stateTemp.copy(S,0,0,ncomps+ntracers);  // Parallel copy.
  auxTemp.copy(Aux,0,0,Aux.nComp());  // Parallel copy.

  Real tagVal = -1;
  if (ngrow>0) {
    for (int n=0; n<ncomps+ntracers; ++n) 
      {      
	const BCRec& theBC = AmrLevel::desc_lst[State_Type].getBC(n);
	TagUnusedGrowCells(S,State_Type,theBC,*this,ngrow,tagVal,n,1);
      }
  }
  
  phiTemp.define(ba, 1, 0, dm, Fab_allocate);

  if (ngrow == 0)
    {
      phiTemp.copy(*rock_phi,0,0,1);
    }
  else
    {
      BL_ASSERT(rock_phi->nGrow() >= ngrow);
      MultiFab phiGrow(BoxArray(rock_phi->boxArray()).grow(ngrow), 1, 0);
      for (MFIter mfi(*rock_phi); mfi.isValid(); ++mfi)
	phiGrow[mfi].copy((*rock_phi)[mfi],0,0,1);
      phiTemp.copy(phiGrow,0,0,1);  // Parallel copy.
    }
  //
  // This gets set by the chemistry solver.
  //
  fcnCntTemp.define(ba, 1, 0, dm, Fab_allocate);
  //
  // It's cheaper to just build a new volume than doing a parallel copy
  // from the existing one.  Additionally this also works when ngrow > 0.
  //
  volTemp.define(ba, 1, 0, dm, Fab_allocate);
  for (MFIter mfi(volTemp); mfi.isValid(); ++mfi)
    geom.GetVolume(volTemp[mfi], volTemp.boxArray(), mfi.index(), 0);
  
  // Make some handy indices, do some more checking to be sure all is set up correctly
  Array<int> guess_comp(ntracers,-1);
  Array<int> activity_comp(ntracers,-1);
  Array<int> sorbed_comp(ntracers,-1);
  Array<int> kd_comp(ntracers,-1);
  std::map<int,int> freundlich_n_comp;
  std::map<int,int> langmuir_b_comp;
  Array<int> volume_fraction_comp(nminerals,-1);
  Array<int> specific_surface_area_comp(nminerals,-1);
  int cation_exchange_capacity_comp = -1;

  int nAux = Aux.nComp();
  for (int i = 0; i < ntracers; ++i) {
    const std::string& name = tNames[i];

    BL_ASSERT(sorption_isotherm_label_map.find(name)!=sorption_isotherm_label_map.end());
    guess_comp[i] = solute_chem_label_map[name]["Free_Ion_Guess"];
    BL_ASSERT(guess_comp[i] >= 0  &&  guess_comp[i] < nAux);
    //std::cout << "Free Ion Guess for " << name << " in comp: " << guess_comp[i] << std::endl;

    activity_comp[i] = solute_chem_label_map[name]["Activity_Coefficient"];
    BL_ASSERT(activity_comp[i] >= 0  &&  activity_comp[i] < nAux);
    //std::cout << "Activity Coefficient for " << name << " in comp: " << activity_comp[i] << std::endl;

    if (using_sorption) {
      sorbed_comp[i] = sorption_chem_label_map[name]["Total_Sorbed"];
      BL_ASSERT(sorbed_comp[i] >= 0  &&  sorbed_comp[i] < nAux);
      //std::cout << "Total Sorbed for " << name << " in comp: " << sorbed_comp[i] << std::endl;
    }

    if (nsorption_isotherms > 0) {
      kd_comp[i] = sorption_isotherm_label_map[name]["Kd"];
      //std::cout << "Kd for " << name << " in comp: " << kd_comp[i] << std::endl;
      BL_ASSERT(kd_comp[i] >= 0  &&  kd_comp[i] < nAux);
      if (sorption_isotherm_label_map[name].count("Freundlich_n")) {
          freundlich_n_comp[i] = sorption_isotherm_label_map[name]["Freundlich_n"];
          BL_ASSERT(freundlich_n_comp[i] >= 0  &&  freundlich_n_comp[i] < nAux);
          //std::cout << "Freundlich n for " << name << " in comp: " << freundlich_n_comp[i] << std::endl;
      }
      else if (sorption_isotherm_label_map[name].count("Langmuir_b")) {
          langmuir_b_comp[i] = sorption_isotherm_label_map[name]["Langmuir_b"];
          BL_ASSERT(langmuir_b_comp[i] >= 0  &&  langmuir_b_comp[i] < nAux);
          //std::cout << "Langmuir b for " << name << " in comp: " << langmuir_b_comp[i] << std::endl;
      }
    }
  }

  for (int i = 0; i < nminerals; ++i) {
    const std::string& name = minerals[i];

    BL_ASSERT(mineralogy_label_map.find(name)!=mineralogy_label_map.end());
    volume_fraction_comp[i] = mineralogy_label_map[name]["Volume_Fraction"];
    BL_ASSERT(volume_fraction_comp[i] >= 0  &&  volume_fraction_comp[i] < nAux);
    //std::cout << "Volume Fraction for " << name << " in comp: " << volume_fraction_comp[i] << std::endl;

    specific_surface_area_comp[i] = mineralogy_label_map[name]["Specific_Surface_Area"];
    BL_ASSERT(specific_surface_area_comp[i] >= 0  &&  specific_surface_area_comp[i] < nAux);
    //std::cout << "Specific Surface Area for " << name << " in comp: " << specific_surface_area_comp[i] << std::endl;
  }

  if (ncation_exchange > 0) {
    cation_exchange_capacity_comp = cation_exchange_label_map["Cation_Exchange_Capacity"];
    BL_ASSERT(cation_exchange_capacity_comp >= 0  &&  cation_exchange_capacity_comp < nAux);      
    //std::cout << "Cation Exchange Capacity in comp: " << cation_exchange_capacity_comp << std::endl;
  }


//  HACK...should be unnecessary
  for (MFIter mfi(stateTemp); mfi.isValid(); ++mfi) {
    setPhysBoundaryValues(stateTemp[mfi],State_Type,time,0,0,ncomps+ntracers);
  }


  // Grab the auxiliary chemistry data
  for (MFIter mfi(stateTemp); mfi.isValid(); ++mfi)
    {
      FArrayBox& fab     = stateTemp[mfi];
      FArrayBox& phi_fab = phiTemp[mfi];
      FArrayBox& vol_fab = volTemp[mfi];
      FArrayBox& fct_fab = fcnCntTemp[mfi];
      FArrayBox& aux_fab = auxTemp[mfi];
      const int* lo      = fab.loVect();
      const int* hi      = fab.hiVect();

#if (BL_SPACEDIM == 2)
      for (int iy=lo[1]; iy<=hi[1]; iy++)
        {
	  int threadid = 0;
	  amanzi::chemistry::SimpleThermoDatabase&     TheChemSolve = chemSolve[threadid];
	  amanzi::chemistry::Beaker::BeakerComponents& TheComponent = components[threadid];
	  amanzi::chemistry::Beaker::BeakerParameters& TheParameter = parameters[threadid];

	  for (int ix=lo[0]; ix<=hi[0]; ix++)
            {
#elif (BL_SPACEDIM == 3)
      
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) 
#endif
      for (int iz = lo[2]; iz<=hi[2]; iz++)
	{		
	  int threadid = 0;
#ifdef _OPENMP
	  threadid = omp_get_thread_num();
#endif
	  amanzi::chemistry::SimpleThermoDatabase&     TheChemSolve = chemSolve[threadid];
	  amanzi::chemistry::Beaker::BeakerComponents& TheComponent = components[threadid];
	  amanzi::chemistry::Beaker::BeakerParameters& TheParameter = parameters[threadid];
	  
	  for (int iy = lo[1]; iy<=hi[1]; iy++)
	    {
	      for (int ix = lo[0]; ix<=hi[0]; ix++)
		{
#else
#error "We only support 2 or 3-D"
#endif
		  int idx_minerals = 0, idx_sorbed = 0, idx_total = 0;
		  
		  IntVect iv(D_DECL(ix,iy,iz));
		  
		  bool allzero = true;
		  bool skip_cell = fab(iv,0) < tagVal+0.1;
		  for (int i = 0; i < ntracers && !skip_cell; ++i)
		    {
		      allzero = allzero && (fab(iv, ncomps+i) == 0);
		    }
		  skip_cell |= allzero;
		  
		  if (!skip_cell) 
		    {
		      for (int i = 0; i < ntracers; ++i)
			{
			  TheComponent.total[i] = fab(iv,ncomps+i);
			  TheComponent.free_ion[i] = aux_fab(iv,guess_comp[i]);

			  if (using_sorption) {
			    TheComponent.total_sorbed[i] = aux_fab(iv, sorbed_comp[i]);
			  }

			  if (nsorption_isotherms > 0) {
			    TheComponent.isotherm_kd[i] = aux_fab(iv,kd_comp[i]);
                            if (sorption_isotherm_label_map[tNames[i]].count("Freundlich_n")) {
                                TheComponent.isotherm_freundlich_n[i] = aux_fab(iv,freundlich_n_comp[i]);
                            }
                            else {
                                TheComponent.isotherm_langmuir_b[i] = aux_fab(iv,langmuir_b_comp[i]);
                            }
			  }
			}
		      
		      // FIXME:
		      //if (ncation_exchange_capacities > 0) {
		      //TheComponent.cation_exchange_capacity = aux_fab(iv,cation_exchange_capacity_comp);
		      //}

		      // TODO: loop over minerals
		      // TODO: loop over surface complexation sites
		      
		      TheParameter.porosity   = phi_fab(iv,0);
		      TheParameter.saturation = std::min(1., std::max(0., fab(iv,0) / density[0]));
		      TheParameter.volume     = vol_fab(iv,0);
		      TheParameter.water_density = density[0];
		      
		      amanzi::chemistry::Beaker::SolverStatus stat;
		      
		      try
			{
			  //TheComponent.Display("-- before rxn step: \n");
			  TheChemSolve.ReactionStep(&TheComponent,TheParameter,dt);
			  
			  stat = TheChemSolve.status();
			  
			  fct_fab(iv,0) = use_funccount ? stat.num_rhs_evaluations : 1;
			}
		      catch (const amanzi::chemistry::ChemistryException& geochem_error)
			{
			  std::cout << iv << " : ";
			  for (int icmp = 0; icmp < ntracers; icmp++)
			    std::cout << fab(iv,icmp+ncomps) << ' ';
			  std::cout << std::endl;
			  TheComponent.Display("components: ");
			  BoxLib::Abort(geochem_error.what());	    
			}
		      //
		      // After calculating the change in the tracer species,
		      // update the state variables.
		      //
		      
		      // loop over the solutes and store the calculated
		      // values in the appropriate place in the fabs
		      for (int i = 0; i < ntracers; ++i)
			{
			  fab(iv, ncomps+i) = TheComponent.total[i];
			  aux_fab(iv,guess_comp[i]) = TheComponent.free_ion[i];
			  if (using_sorption) {
			    aux_fab(iv,sorbed_comp[i]) = TheComponent.total_sorbed[i];
			  }
			  if (nsorption_isotherms) {
			    aux_fab(iv,kd_comp[i]) = TheComponent.isotherm_kd[i];
			    aux_fab(iv,freundlich_n_comp[i]) = TheComponent.isotherm_freundlich_n[i];
			    aux_fab(iv,langmuir_b_comp[i]) = TheComponent.isotherm_langmuir_b[i];
			  }
			}

		      // FIXME:
		      //if (ncation_exchange > 0) {
		      //  aux_fab(iv,cation_exchange_capacity_comp) = TheComponent.cation_exchange_capacity;
		      //}

		      // TODO: loop over minerals, porosity, etc
		    }
		}
	    }
#if (BL_SPACEDIM == 3)
	}
#endif
    }
    
  phiTemp.clear();
  volTemp.clear();
    
  S.copy(stateTemp,ncomps,ncomps,ntracers); // Parallel copy, tracers only
  stateTemp.clear();
    
  Aux.copy(auxTemp,0,0,Aux.nComp()); // Parallel copy, everything.
  auxTemp.clear();
    
  if (ngrow == 0)
    {
      Fcnt.copy(fcnCntTemp,0,0,1); // Parallel copy.
      fcnCntTemp.clear();
	
      S.FillBoundary();
      Aux.FillBoundary();
      Fcnt.FillBoundary();
	
      geom.FillPeriodicBoundary(S,true);
      geom.FillPeriodicBoundary(Aux,true);
      geom.FillPeriodicBoundary(Fcnt,true);
    }
  else
    {
      //
      // Can't directly use a parallel copy to update FuncCount_Type.
      //
      MultiFab grownFcnt(BoxArray(Fcnt.boxArray()).grow(ngrow), 1, 0);
      grownFcnt.setVal(1);
      for (MFIter mfi(Fcnt); mfi.isValid(); ++mfi)
	grownFcnt[mfi].copy(Fcnt[mfi]);
	
      grownFcnt.copy(fcnCntTemp); // Parallel copy.
      fcnCntTemp.clear();
	
      for (MFIter mfi(grownFcnt); mfi.isValid(); ++mfi)
	Fcnt[mfi].copy(grownFcnt[mfi]);
    }
    
#else /* Not AMANZI */
    
  if (do_chem==0)
    {
      MultiFab tmp;
      tmp.define(S.boxArray(),ncomps,0,Fab_allocate);
      tmp.copy(S,0,0,ncomps);
	
      for (MFIter mfi(S);mfi.isValid();++mfi)
	{
	  const int* s_lo  = tmp[mfi].loVect();
	  const int* s_hi  = tmp[mfi].hiVect();
	  const Real* sdat = tmp[mfi].dataPtr();
	    
	  if (ncomps == 4) 
	    FORT_CHEM_DUMMY(sdat,ARLIM(s_lo),ARLIM(s_hi),&dt,&ncomps);
	}     
      S.copy(tmp,0,0,ncomps);
    }
#endif
    
  if (verbose > 2 && ParallelDescriptor::IOProcessor())
    {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;
      ParallelDescriptor::ReduceRealMax(run_time,IOProc);
	
      std::cout << "PorousMedia::strang_chem time: " << run_time << '\n';
    }
}
    
void
PorousMedia::set_preferred_boundary_values (MultiFab& S,
					    int       state_index,
					    int       src_comp,
					    int       dst_comp,
					    int       num_comp,
					    Real      time) const
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::set_preferred_boundary_values()");

  if (state_index == State_Type)
  {
      const TimeLevel whichTime = which_time(State_Type,time);
      //
      // To get chem-advanced data instead of FP'd data at old time.
      //
      // For AMANZI the chem-advanced data are the tracers.
      //
      if (!FillPatchedOldState_ok && whichTime == AmrOldTime)
      {
          if (src_comp == ncomps && num_comp == ntracers)
          {
              aux_boundary_data_old.copyTo(S, src_comp, dst_comp, num_comp);
          }
      }
  }
}

//
// Compute capillary update.  This assumes there are only 2 phases and
// incompressible.  We only solve for component 1, and solution to 
// component 2 are deduced from component 1.
//
void
PorousMedia::scalar_capillary_update (Real      dt,
				      int       corrector,
				      MultiFab* u_mac)

{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scalar_capillary_update()");

  BL_ASSERT(nphases == 2);
  BL_ASSERT(have_capillary == 1);

  const Real strt_time = ParallelDescriptor::second();

  // Build single component edge-centered array of MultiFabs for fluxes
  MultiFab** fluxSCn;
  MultiFab** fluxSCnp1;
  const int nGrow = 0;
  const int nComp = 1;
  diffusion->allocFluxBoxesLevel(fluxSCn,  nGrow,nComp);
  diffusion->allocFluxBoxesLevel(fluxSCnp1,nGrow,nComp);

  int nc = 0; 
  int nd = 1;
  MultiFab*  delta_rhs = 0;
  MultiFab*  alpha     = 0;
  MultiFab** cmp_pcn   = 0;
  MultiFab** cmp_pcnp1 = 0;
  MultiFab** cmp_pcnp1_dp = 0;
  MultiFab&  S_new = get_new_data(State_Type);

  MultiFab* sat_res_mf = new MultiFab(grids,1,1);
  sat_res_mf->setVal(1.);
  for (MFIter mfi(*sat_res_mf); mfi.isValid();++mfi)
    {
      const Box& box = (*sat_res_mf)[mfi].box();
      (*sat_res_mf)[mfi].minus((*cpl_coef)[mfi],box,3,0,1);
    }
  sat_res_mf->mult(density[nc]);
  diffusion->set_rho(sat_res_mf); 

  MultiFab* S_nwt = new MultiFab(grids,1,1);
  MultiFab::Copy(*S_nwt,S_new,nc,0,nComp,1);

  alpha = new MultiFab(grids, 1, 1);
  MultiFab::Copy(*alpha,*rock_phi,0,0,1,alpha->nGrow());
  
  // Newton method.
  // initialization
  Real pcTime = state[State_Type].prevTime();
  diffusion->allocFluxBoxesLevel(cmp_pcn,0,1);
  calcCapillary(pcTime);
  calcDiffusivity_CPL(cmp_pcn,lambda_cc); 
  diffusion->diffuse_init_CPL(dt,nc,be_cn_theta,
			      fluxSCn,0,delta_rhs,
			      alpha,cmp_pcn,pcn_cc,S_nwt);
  pcTime = state[State_Type].curTime();
  FillStateBndry (pcTime,State_Type,0,ncomps);
  diffusion->allocFluxBoxesLevel(cmp_pcnp1,0,1);
  diffusion->allocFluxBoxesLevel(cmp_pcnp1_dp,0,1);
  calcCapillary(pcTime);
  calcLambda(pcTime);
  calcDiffusivity_CPL(cmp_pcnp1,lambdap1_cc);
  calcDiffusivity_CPL_dp(cmp_pcnp1_dp,lambdap1_cc,pcTime,1);

  int  max_itr_nwt = 20;
#if (BL_SPACEDIM == 3)
  Real max_err_nwt = 1e-8;
#else
  Real max_err_nwt = 1e-8;
#endif
  int  itr_nwt = 0;
  Real err_nwt = 1e10;
  Real be_theta = be_cn_theta;
  while ((itr_nwt < max_itr_nwt) && (err_nwt > max_err_nwt)) 
    {
      diffusion->diffuse_iter_CPL(dt,nc,ncomps,be_theta,
				  0,alpha,cmp_pcnp1,cmp_pcnp1_dp,
				  pcnp1_cc,S_nwt,&err_nwt);

      if (verbose > 3 && ParallelDescriptor::IOProcessor())
	std::cout << "Newton iteration " << itr_nwt 
	          << " : Error = "       << err_nwt << "\n"; 

      //scalar_adjust_constraint(0,ncomps-1);
      //FillStateBndry(pcTime,State_Type,0,ncomps);
      //calcCapillary(pcTime);
      //calcLambda(pcTime);
      calcDiffusivity_CPL(cmp_pcnp1,lambdap1_cc);
      calcDiffusivity_CPL_dp(cmp_pcnp1_dp,lambdap1_cc,pcTime,1);
      itr_nwt += 1;

      if (verbose > 3)
	check_minmax();
    }
    
  diffusion->compute_flux(nc,dt,be_cn_theta,fluxSCnp1,pcnp1_cc,cmp_pcnp1);

  if (verbose > 3 && ParallelDescriptor::IOProcessor())
    {
      if (itr_nwt < max_itr_nwt)
	std::cout << "Newton converged at iteration " << itr_nwt
		  << " with error " << err_nwt << '\n';
      else
	std::cout << "Newton failed to converged: termination error is "
		  <<  err_nwt << '\n'; 
    }

  //
  // add to phase velocity
  //
  if (u_mac != 0) {
      
    FArrayBox fluxtot;

    for (int d = 0; d < BL_SPACEDIM; d++) 
      {
	for (MFIter fmfi(*fluxSCn[d]); fmfi.isValid(); ++fmfi) 
	  {
	    const Box& ebox = (*fluxSCn[d])[fmfi].box();
	    fluxtot.resize(ebox,nComp);
	    fluxtot.copy((*fluxSCn[d])[fmfi],ebox,0,ebox,0,nComp);
	    if (no_corrector == 1)
	      fluxtot.mult(2.0);
	    else
	      fluxtot.plus((*fluxSCnp1[d])[fmfi],ebox,0,0,nComp);

	    fluxtot.mult(-1.0/density[nc]);
	    fluxtot.divide(area[d][fmfi],0,0,1);
	    u_mac[d][fmfi].plus(fluxtot,ebox,0,0,nComp);
	  }
	u_mac[d].FillBoundary();
      }
  }

  //
  // Increment the viscous flux registers
  // The fluxes are - beta \nabla p_c. We accummulate flux assuming 
  // it is on the LHS.  Thus, we need to multiply by -dt due to the sign change.
  // 

  if (do_reflux && corrector) {

    FArrayBox fluxtot;
	  
    for (int d = 0; d < BL_SPACEDIM; d++) 
      {
	MultiFab fluxes;

	if (level < parent->finestLevel())
	  fluxes.define((*fluxSCn[d]).boxArray(), ncomps, 0, Fab_allocate);

	for (MFIter fmfi(*fluxSCn[d]); fmfi.isValid(); ++fmfi)
	  {
	    // for component nc
	    const Box& ebox = (*fluxSCn[d])[fmfi].box();

	    fluxtot.resize(ebox,ncomps);
	    fluxtot.copy((*fluxSCn[d])[fmfi],ebox,0,ebox,nc,1);
	    fluxtot.plus((*fluxSCnp1[d])[fmfi],ebox,0,nc,1);

	    (*fluxSCn[d])[fmfi].mult(-density[nd]/density[nc]);
	    (*fluxSCnp1[d])[fmfi].mult(-density[nd]/density[nc]);
	    fluxtot.copy((*fluxSCn[d])[fmfi],ebox,0,ebox,nd,1);
	    fluxtot.plus((*fluxSCnp1[d])[fmfi],ebox,0,nd,1);

	    if (level < parent->finestLevel())
	      fluxes[fmfi].copy(fluxtot);

	    if (level > 0)
	      getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),0,0,ncomps,-dt);
	  }

	if (level < parent->finestLevel())
	  getLevel(level+1).getViscFluxReg().CrseInit(fluxes,d,0,0,ncomps,dt);
		  
      }
  }
    
  //     nc = 0; 
  //     MultiFab::Copy(*S_nwt,S_new,nc,0,nComp,1);
  //     diffusion->check_consistency(dt,nc,mone,be_cn_theta,
  // 				 rho_flag,0,alpha,
  // 				 cmp_pcn,cmp_pcnp1,
  // 				 pcn_cc, pcnp1_cc,
  // 				 S_nwt,&err_nwt);

  delete delta_rhs;
  delete S_nwt;
  delete sat_res_mf;
  delete alpha;

  diffusion->removeFluxBoxesLevel(cmp_pcn);
  diffusion->removeFluxBoxesLevel(cmp_pcnp1);
  diffusion->removeFluxBoxesLevel(cmp_pcnp1_dp);
	   
  diffusion->removeFluxBoxesLevel(fluxSCn);
  diffusion->removeFluxBoxesLevel(fluxSCnp1);
    
  if (verbose > 3)
    {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      if (ParallelDescriptor::IOProcessor())
        std::cout << "PorousMedia::scalar_CPL_update(): time: " << run_time << '\n';
    }
  
  //
  // Write out the min and max of each component of the new state.
  //
  if (verbose > 3) check_minmax();

}

//
// Compute capillary update.  This assumes there are only 2 phases and
// incompressible.  We only solve for component 1, and solution to 
// component 2 are deduced from component 1.
//
void
PorousMedia::diff_capillary_update (Real      dt,
				    int       corrector,
				    MultiFab* u_mac)
  
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::diff_capillary_update()");

  BL_ASSERT(nphases == 2);
  BL_ASSERT(have_capillary == 1);

  const Real strt_time = ParallelDescriptor::second();

  // Build single component edge-centered array of MultiFabs for fluxes
  MultiFab** fluxSCn;
  MultiFab** fluxSCnp1;
  const int nGrow = 0;
  const int nComp = 1;
  diffusion->allocFluxBoxesLevel(fluxSCn,  nGrow,nComp);
  diffusion->allocFluxBoxesLevel(fluxSCnp1,nGrow,nComp);

  int nc = 0; 
  int nd = 1;
  MultiFab*  delta_rhs    = 0;
  MultiFab*  alpha        = 0;
  MultiFab** tmp          = 0;
  MultiFab** cmp_pcnp1_dp = 0;
  MultiFab*  S_nwt = 0;
  MultiFab&  S_new = get_new_data(State_Type);
  diffusion->allocFluxBoxesLevel(cmp_pcnp1_dp,0,1);

  MultiFab* sat_res_mf = new MultiFab(grids,1,1);
  sat_res_mf->setVal(1.);
  for (MFIter mfi(*sat_res_mf); mfi.isValid();++mfi)
    {
      const Box& box = (*sat_res_mf)[mfi].box();
      (*sat_res_mf)[mfi].minus((*cpl_coef)[mfi],box,3,0,1);
    }
  sat_res_mf->mult(density[nc]);
  diffusion->set_rho(sat_res_mf); 

  S_nwt = new MultiFab(grids,1,1);
  MultiFab::Copy(*S_nwt,S_new,nc,0,nComp,1);

  alpha = new MultiFab(grids, 1, 1);
  MultiFab::Copy(*alpha,*rock_phi,0,0,1,alpha->nGrow());

  tmp = new MultiFab* [BL_SPACEDIM];
  for (int d=0; d<BL_SPACEDIM; d++)
    tmp[d] = &lambda[d];
  
  MultiFab* Stmp = new MultiFab(grids,1,1);
  MultiFab::Copy(*Stmp,*pcn_cc,0,0,1,1);
  MultiFab::Add(*Stmp,*pcnp1_cc,0,0,1,1);
  (*Stmp).mult(0.5);
  
  // Newton method.
  // initialization
  diffusion->diffuse_init_CPL(dt,nc,-be_cn_theta,
			      fluxSCn,0,delta_rhs,
			      alpha,tmp,Stmp,S_nwt);

  Real pcTime = state[State_Type].prevTime();
  
  Stmp->setVal(0);
  calcCapillary(pcTime);
  calcLambda(pcTime);
  calcDiffusivity_CPL(tmp,lambda_cc);

  diffusion->diffuse_init_CPL(dt,nc,be_cn_theta,
			      fluxSCnp1,0,delta_rhs,
			      alpha,tmp,pcn_cc,Stmp);

  MultiFab::Add(*S_nwt,*Stmp,0,0,1,0);
  delete Stmp;

  for (int d=0; d<BL_SPACEDIM; d++)
    MultiFab::Add(*fluxSCn[d],*fluxSCnp1[d],0,0,1,0);

  pcTime = state[State_Type].curTime();
  calcCapillary(pcTime);
  calcLambda(pcTime);
  calcDiffusivity_CPL(tmp,lambdap1_cc);
  calcDiffusivity_CPL_dp(cmp_pcnp1_dp,lambdap1_cc,pcTime,1);

  int  max_itr_nwt = 20;
#if (BL_SPACEDIM == 3)
  Real max_err_nwt = 1e-8;
#else
  Real max_err_nwt = 1e-8;
#endif
  int  itr_nwt = 0;
  Real err_nwt = 1e10;
  Real be_theta = be_cn_theta;

  while ((itr_nwt < max_itr_nwt) && (err_nwt > max_err_nwt)) 
    {
      diffusion->diffuse_iter_CPL(dt,nc,ncomps,be_theta,
				  0,alpha,tmp,cmp_pcnp1_dp,
				  pcnp1_cc,S_nwt,&err_nwt);

      if (verbose > 3 && ParallelDescriptor::IOProcessor())
	std::cout << "Newton iteration " << itr_nwt 
	          << " : Error = "       << err_nwt << "\n"; 

      //scalar_adjust_constraint(0,ncomps-1);
      //FillStateBndry(pcTime,State_Type,0,ncomps);
      //calcCapillary(pcTime);
      //calcLambda(pcTime);
      calcDiffusivity_CPL(tmp,lambdap1_cc);
      calcDiffusivity_CPL_dp(cmp_pcnp1_dp,lambdap1_cc,pcTime,1);
      itr_nwt += 1;

      if (verbose > 3)
	check_minmax();
    }

  diffusion->compute_flux(nc,dt,be_cn_theta,fluxSCnp1,pcnp1_cc,tmp);
    
  if (verbose > 3 && ParallelDescriptor::IOProcessor())
    {
      if (itr_nwt < max_itr_nwt)
	std::cout << "Newton converged at iteration " << itr_nwt
		  << " with error " << err_nwt << '\n';
      else
	std::cout << "Newton failed to converged: termination error is "
		  <<  err_nwt << '\n'; 
    }

  //
  // add to phase velocity
  //
  if (u_mac != 0) {
      
    FArrayBox fluxtot;

    for (int d = 0; d < BL_SPACEDIM; d++) 
      {
	for (MFIter fmfi(*fluxSCn[d]); fmfi.isValid(); ++fmfi) 
	  {
	    const Box& ebox = (*fluxSCn[d])[fmfi].box();
	    fluxtot.resize(ebox,nComp);
	    fluxtot.copy((*fluxSCn[d])[fmfi],ebox,0,ebox,0,nComp);
	    if (no_corrector == 1)
	      fluxtot.mult(2.0);
	    else
	      fluxtot.plus((*fluxSCnp1[d])[fmfi],ebox,0,0,nComp);

	    fluxtot.mult(-1.0/density[nc]);
	    fluxtot.divide(area[d][fmfi],0,0,1);
	    u_mac[d][fmfi].plus(fluxtot,ebox,0,0,nComp);
	  }
	u_mac[d].FillBoundary();
      }
  }

  //
  // Increment the viscous flux registers
  // The fluxes are - beta \nabla p_c. We accummulate flux assuming 
  // it is on the LHS.  Thus, we need to multiply by -dt due to the sign change.
  // 

  if (do_reflux && corrector) {

    FArrayBox fluxtot;
	  
    for (int d = 0; d < BL_SPACEDIM; d++) 
      {
	MultiFab fluxes;
	
	if (level < parent->finestLevel())
	  fluxes.define((*fluxSCn[d]).boxArray(), ncomps, 0, Fab_allocate);
	
	for (MFIter fmfi(*fluxSCn[d]); fmfi.isValid(); ++fmfi)
	  {
	    // for component nc
	    const Box& ebox = (*fluxSCn[d])[fmfi].box();
	    
	    fluxtot.resize(ebox,ncomps);
	    fluxtot.copy((*fluxSCn[d])[fmfi],ebox,0,ebox,nc,1);
	    fluxtot.plus((*fluxSCnp1[d])[fmfi],ebox,0,nc,1);

	    (*fluxSCn[d])[fmfi].mult(-density[nd]/density[nc]);
	    (*fluxSCnp1[d])[fmfi].mult(-density[nd]/density[nc]);
	    fluxtot.copy((*fluxSCn[d])[fmfi],ebox,0,ebox,nd,1);
	    fluxtot.plus((*fluxSCnp1[d])[fmfi],ebox,0,nd,1);

	    if (level < parent->finestLevel())
	      fluxes[fmfi].copy(fluxtot);

	    if (level > 0)
	      getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),0,0,ncomps,-dt);
	  }

	if (level < parent->finestLevel())
	  getLevel(level+1).getViscFluxReg().CrseInit(fluxes,d,0,0,ncomps,dt);
		  
      }
  }

  delete delta_rhs;
  delete S_nwt;
  delete alpha;
  delete sat_res_mf;
  delete [] tmp;

  diffusion->removeFluxBoxesLevel(cmp_pcnp1_dp);
  diffusion->removeFluxBoxesLevel(fluxSCn);
  diffusion->removeFluxBoxesLevel(fluxSCnp1);
    
  if (verbose > 3)
    {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      if (ParallelDescriptor::IOProcessor())
        std::cout << "PorousMedia::diff_CPL_update(): time: " << run_time << '\n';
    }

  //
  // Write out the min and max of each component of the new state.
  //
  if (verbose > 3) check_minmax();
}

#ifdef MG_USE_FBOXLIB
//
// Richard equation: Equilibrium solver
//
void
PorousMedia::richard_eqb_update (MultiFab* u_mac)

{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::richards_update()");
  BL_ASSERT(nphases == 1);
  BL_ASSERT(have_capillary == 1);

  const Real strt_time = ParallelDescriptor::second();

  // Build single component edge-centered array of MultiFabs for fluxes
  MultiFab** fluxSC;
  const int nGrow = 0;
  const int nComp = 1;
  diffusion->allocFluxBoxesLevel(fluxSC,nGrow,nComp);

  int nc = 0; 
  MultiFab** cmp_pcp1    = 0;
  MultiFab** cmp_pcp1_dp = 0;
  MultiFab sat_res_mf(grids,1,1);
  sat_res_mf.setVal(1.);
  for (MFIter mfi(sat_res_mf); mfi.isValid();++mfi)
    {
      const Box& box = sat_res_mf[mfi].box();
      sat_res_mf[mfi].minus((*cpl_coef)[mfi],box,3,0,1);
    }
  //sat_res_mf.mult(density[nc]);
  diffusion->set_rho(&sat_res_mf); 

  // Compute first res_fix = \nabla v_inflow.  
  // Its value does not change.
  MultiFab res_fix(grids,1,0);
  res_fix.setVal(0.);
  calc_richard_velbc(res_fix,u_mac);

  // Newton method.
  // initialization
  int do_upwind = 1;
  int  max_itr_nwt = 10;
  Real max_err_nwt = 1e-8;
  int  itr_nwt = 0;
  Real err_nwt = 1e10;
  Real pcTime = state[State_Type].curTime();
  MultiFab& P_new = get_new_data(Press_Type);
  FillStateBndry(pcTime,State_Type,0,ncomps);
  diffusion->allocFluxBoxesLevel(cmp_pcp1,0,1);
  diffusion->allocFluxBoxesLevel(cmp_pcp1_dp,0,3);
  calcCapillary(pcTime);
  calcLambda(pcTime);
  MultiFab::Copy(P_new,*pcnp1_cc,0,0,1,1);
  P_new.mult(-1.0,1);
  calc_richard_coef(cmp_pcp1,lambdap1_cc,u_mac,0,do_upwind,pcTime);
  MultiFab* dalpha=0;
  bool do_n = true;
  Real dt=-1;
  calc_richard_jac (cmp_pcp1_dp,dalpha,lambdap1_cc,u_mac,pcTime,dt,0,do_upwind,do_n);
  while ((itr_nwt < max_itr_nwt) && (err_nwt > max_err_nwt)) 
    {
      diffusion->richard_iter_eqb(nc,gravity,density,res_fix,
				  cmp_pcp1,cmp_pcp1_dp,
				  u_mac,do_upwind,&err_nwt);      
      if (verbose > 3 && ParallelDescriptor::IOProcessor())
	std::cout << "Newton iteration " << itr_nwt 
	          << " : Error = "       << err_nwt << "\n"; 
      scalar_adjust_constraint(0,ncomps-1);
      FillStateBndry(pcTime,State_Type,0,ncomps);
      calcCapillary(pcTime);
      calcLambda(pcTime);
      MultiFab::Copy(P_new,*pcnp1_cc,0,0,1,1);
      P_new.mult(-1.0,1);
      compute_vel_phase(u_mac,0,pcTime);
      calc_richard_coef(cmp_pcp1,lambdap1_cc,u_mac,0,do_upwind,pcTime);
      calc_richard_jac (cmp_pcp1_dp,dalpha,lambdap1_cc,u_mac,pcTime,dt,0,do_upwind,do_n);
      itr_nwt += 1;

      if (verbose > 3)	check_minmax();
    }
    
  diffusion->compute_flux(nc,1.0,1.0,fluxSC,pcnp1_cc,cmp_pcp1);

  if (verbose > 3 && ParallelDescriptor::IOProcessor())
    {
      if (itr_nwt < max_itr_nwt)
	std::cout << "Newton converged at iteration " << itr_nwt
		  << " with error " << err_nwt << '\n';
      else
	std::cout << "Newton failed to converged: termination error is "
		  <<  err_nwt << '\n'; 
    }

  /*
  //
  // Increment the viscous flux registers
  // The fluxes are - beta \nabla p_c. We accummulate flux assuming 
  // it is on the LHS.  Thus, we need to multiply by -dt due to the sign change.
  // 

  if (do_reflux && corrector) {

      FArrayBox fluxtot;
      for (int d = 0; d < BL_SPACEDIM; d++) 
	{
	  MultiFab fluxes;
	  
	  if (level < parent->finestLevel())
	    fluxes.define((*fluxSC[d]).boxArray(), nComp, 0, Fab_allocate);
	  
	  for (MFIter fmfi(*fluxSC[d]); fmfi.isValid(); ++fmfi)
	    {
	      // for component nc
	      const Box& ebox = (*fluxSC[d])[fmfi].box();
	      
	      fluxtot.resize(ebox,nComp);
	      fluxtot.copy((*fluxSC[d])[fmfi],ebox,0,ebox,0,1);

	      if (level < parent->finestLevel())
		fluxes[fmfi].copy(fluxtot);
	      
	      if (level > 0)
		getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),0,0,nComp,-dt);
	    }
	  
	  if (level < parent->finestLevel())
	    getLevel(level+1).getViscFluxReg().CrseInit(fluxes,d,0,0,nComp,dt);
      }
  }
  */

  diffusion->removeFluxBoxesLevel(cmp_pcp1);
  diffusion->removeFluxBoxesLevel(cmp_pcp1_dp);
  diffusion->removeFluxBoxesLevel(fluxSC);
    
  if (verbose > 3)
    {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      

      if (ParallelDescriptor::IOProcessor())
        std::cout << "PorousMedia::richard_update(): time: " 
		  << run_time << '\n';
    }
  //
  // Write out the min and max of each component of the new state.
  //
  if (verbose > 3) check_minmax();

}

//
// Richard equation: Time-dependent solver.  Only doing a first-order implicit scheme
//
RichardNLSdata::Reason
PorousMedia::richard_scalar_update (Real dt, int& total_nwt_iter, MultiFab* u_mac)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::richards_update()");
  BL_ASSERT(have_capillary == 1);

  const Real strt_time = ParallelDescriptor::second();

  const int nGrow = 0;
  int nc = 0; 
  MultiFab** cmp_pcp1    = 0;
  MultiFab** cmp_pcp1_dp = 0;
  MultiFab sat_res_mf(grids,1,1);
  sat_res_mf.setVal(1.);
  for (MFIter mfi(sat_res_mf); mfi.isValid();++mfi)
    {
      const Box& box = sat_res_mf[mfi].box();
      sat_res_mf[mfi].minus((*cpl_coef)[mfi],box,3,0,1);
    }
  diffusion->set_rho(&sat_res_mf);

  MultiFab& S_new = get_new_data(State_Type);
  MultiFab& P_new = get_new_data(Press_Type);
  MultiFab* alpha = new MultiFab(grids,1,1);
  MultiFab::Copy(*alpha,*rock_phi,0,0,1,alpha->nGrow());

  // Compute first res_fix = -\phi * n^k + dt*\nabla v_inflow.  
  // Its value does not change.
  Real pcTime = state[State_Type].curTime();
  MultiFab res_fix(grids,1,0);
  MultiFab::Copy(res_fix,S_new,nc,0,1,0);
  for (MFIter mfi(res_fix); mfi.isValid(); ++mfi)
    res_fix[mfi].mult((*alpha)[mfi],mfi.validbox(),0,0,1);
  res_fix.mult(-1.0);
  compute_vel_phase(u_mac,0,pcTime);
  calc_richard_velbc(res_fix,u_mac,dt*density[0]);

  // Newton method.
  // initialization
  int do_upwind = 1;
  int  max_itr_nwt = total_nwt_iter;
  Real max_err_nwt = 1e-6;
  int  itr_nwt = 0;
  Real err_nwt = 1e10;
  FillStateBndry(pcTime,State_Type,0,ncomps);
  FillStateBndry(pcTime,Press_Type,0,1);
  diffusion->allocFluxBoxesLevel(cmp_pcp1,0,1);
  diffusion->allocFluxBoxesLevel(cmp_pcp1_dp,0,3);

  calcLambda(pcTime);
  calc_richard_coef(cmp_pcp1,lambdap1_cc,u_mac,0,do_upwind,pcTime);

  MultiFab* dalpha = 0;
  if (!do_richard_sat_solve) {
      dalpha = new MultiFab(grids,1,1,Fab_allocate); // Note: requires a delete
  }
  calc_richard_jac (cmp_pcp1_dp,dalpha,lambdap1_cc,u_mac,pcTime,dt,0,do_upwind,do_richard_sat_solve);

  Diffusion::NewtonStepInfo linear_status;
  linear_status.success = true;
  linear_status.max_ls_iterations = richard_max_ls_iterations;
  linear_status.min_ls_factor = richard_min_ls_factor;
  linear_status.ls_acceptance_factor = richard_ls_acceptance_factor;
  linear_status.ls_reduction_factor = richard_ls_reduction_factor;
  linear_status.monitor_linear_solve = richard_monitor_linear_solve;
  linear_status.monitor_line_search = richard_monitor_line_search;

  if (do_richard_sat_solve)
    {
      calcCapillary(pcTime);
      MultiFab::Copy(P_new,*pcnp1_cc,0,0,1,1);
      P_new.mult(-1.0,1);
      while ((itr_nwt < max_itr_nwt) && (err_nwt > max_err_nwt) && (linear_status.success)) 
	{
            itr_nwt++;
            diffusion->richard_iter(dt,nc,gravity,density,res_fix,
                                    alpha,cmp_pcp1,cmp_pcp1_dp,
                                    u_mac,do_upwind,linear_status);
            
            err_nwt = linear_status.residual_norm_post_ls;

            if (linear_status.success) {
                if (richard_solver_verbose>0 && ParallelDescriptor::IOProcessor()) {
                    std::cout << "     Iteration (n) " << itr_nwt 
                              << " : Error = "       << err_nwt << " (tol: " << max_err_nwt << ")\n";
                }
                if (model != model_list["richard"])
                    scalar_adjust_constraint(0,ncomps-1);
                FillStateBndry(pcTime,State_Type,0,ncomps);
                calcCapillary(pcTime);
                calcLambda(pcTime);
                MultiFab::Copy(P_new,*pcnp1_cc,0,0,1,1);
                P_new.mult(-1.0,1);
                compute_vel_phase(u_mac,0,pcTime);
                calc_richard_coef(cmp_pcp1,lambdap1_cc,u_mac,0,do_upwind,pcTime);
                calc_richard_jac (cmp_pcp1_dp,dalpha,lambdap1_cc,u_mac,pcTime,dt,0,do_upwind,do_richard_sat_solve);
                if (verbose > 3) check_minmax();
            }
	}
    }
  else
    {
        //MultiFab dalpha(grids,1,1);
        //calc_richard_alpha(&dalpha,pcTime);

      while ((itr_nwt < max_itr_nwt) && (err_nwt > max_err_nwt)  && (linear_status.success)) 
	{
            itr_nwt++;
            diffusion->richard_iter_p(dt,nc,gravity,density,res_fix,
                                      alpha,dalpha,cmp_pcp1,cmp_pcp1_dp,
                                      u_mac,do_upwind,linear_status);

	    err_nwt = linear_status.residual_norm_post_ls;

            if (linear_status.success) {
                if (richard_solver_verbose>0 && ParallelDescriptor::IOProcessor()) {
                    std::cout << "     Iteration (p) " << itr_nwt 
                              << " : Error = "       << err_nwt << " (tol: " << max_err_nwt << ")\n";
                }
                calcInvPressure (S_new,P_new); 
                calcLambda(pcTime);
                compute_vel_phase(u_mac,0,pcTime);
                calc_richard_coef(cmp_pcp1,lambdap1_cc,u_mac,0,do_upwind,pcTime);
                calc_richard_jac (cmp_pcp1_dp,dalpha,lambdap1_cc,u_mac,pcTime,dt,0,do_upwind,do_richard_sat_solve);
                //calc_richard_alpha(&dalpha,pcTime);
                if (verbose > 3)  check_minmax();
            }
	}
    }

  RichardNLSdata::Reason retVal = RichardNLSdata::RICHARD_SUCCESS;
  if (!linear_status.success) {
      retVal = RichardNLSdata::RICHARD_LINEAR_FAIL;
  }
  if (itr_nwt >= max_itr_nwt) {
      retVal = RichardNLSdata::RICHARD_NONLINEAR_FAIL;
      if (richard_solver_verbose>1 && ParallelDescriptor::IOProcessor())
          std::cout << "     **************** Newton failed in richard_scalar_update: too many iterations (max = " << max_itr_nwt << '\n'; 
  }

  MultiFab** fluxSC;
  const int nComp = 1;
  if (retVal == RichardNLSdata::RICHARD_SUCCESS) {
      diffusion->allocFluxBoxesLevel(fluxSC,nGrow,nComp);
      diffusion->richard_flux(nc,-1.0,gravity,density,fluxSC,pcnp1_cc,cmp_pcp1);
  }

  delete alpha;
  if (!do_richard_sat_solve) {
      delete dalpha;
  }
  diffusion->removeFluxBoxesLevel(cmp_pcp1);
  diffusion->removeFluxBoxesLevel(cmp_pcp1_dp);

  if (retVal != RichardNLSdata::RICHARD_SUCCESS) {
      return retVal;
  }

  if (richard_solver_verbose>1 && ParallelDescriptor::IOProcessor() && retVal == RichardNLSdata::RICHARD_SUCCESS) {
      std::cout << "     Newton converged in " << itr_nwt << " iterations (max = "
                << total_nwt_iter << ") with err: " 
                << err_nwt << " (tol = " << max_err_nwt << ")\n";
  }
  total_nwt_iter = itr_nwt;
  
  //
  // Increment the viscous flux registers
  // The fluxes are - beta \nabla p_c. We accummulate flux assuming 
  // it is on the LHS.  Thus, we need to multiply by -dt.
  // 
  if (do_reflux) 
    {
      FArrayBox fluxtot;
      for (int d = 0; d < BL_SPACEDIM; d++) 
	{
	  MultiFab fluxes;
	  
	  if (level < parent->finestLevel())
	    fluxes.define((*fluxSC[d]).boxArray(), nComp, 0, Fab_allocate);
	  
	  for (MFIter fmfi(*fluxSC[d]); fmfi.isValid(); ++fmfi)
	    {
	      // for component nc
	      const Box& ebox = (*fluxSC[d])[fmfi].box();
	      
	      fluxtot.resize(ebox,nComp);
	      fluxtot.copy((*fluxSC[d])[fmfi],ebox,0,ebox,0,1);
	      
	      if (level < parent->finestLevel())
		fluxes[fmfi].copy(fluxtot);
	      
	      if (level > 0)
		getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),0,0,nComp,-dt);
	    }
	  
	  if (level < parent->finestLevel())
	    getLevel(level+1).getViscFluxReg().CrseInit(fluxes,d,0,0,nComp,dt);
      }
  }
      
  diffusion->removeFluxBoxesLevel(fluxSC);    

 
  Real run_time = ParallelDescriptor::second() - strt_time;
  richard_time = run_time;
  ParallelDescriptor::ReduceRealMax(richard_time);
  richard_time_min = std::min(richard_time_min,richard_time);
  
  if (verbose > 3)
    { 
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      ParallelDescriptor::ReduceRealMax(run_time,IOProc);


      if (ParallelDescriptor::IOProcessor())
        std::cout << "PorousMedia::richard_update(): time: " 
		  << run_time << '\n';
    }
  //
  // Write out the min and max of each component of the new state.
  //
  if (verbose > 3) check_minmax();

  return retVal;
}



//
// Richard equation: Time-dependent solver.  Only doing a first-order implicit scheme
//
RichardNLSdata::Reason 
PorousMedia::richard_composite_update (Real dt, RichardNLSdata& nl_data)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::richards_composite_update()");
  BL_ASSERT(have_capillary == 1);

  const Real strt_time = ParallelDescriptor::second();

  int nlevs = parent->finestLevel() - level + 1;
  int nc = 0;

  Array<MultiFab*>& u_mac_local = nl_data.velPhase;
  Array < PArray <MultiFab> >& cmp_pcp1 = nl_data.Hcoeffs;
  Array < PArray <MultiFab> >& cmp_pcp1_dp = nl_data.Jacobian;
  PArray <MultiFab>& dalpha = nl_data.DAlpha;

  // Create a nlevs-level array for the coefficients
  PArray <MultiFab> alpha(nlevs,PArrayManage);
  PArray <MultiFab> res_fix(nlevs,PArrayManage);
    
  int do_upwind = 1;
  int  max_itr_nwt = nl_data.MaxNLIterations();
  Real max_err_nwt = 1e-6;
  Real err_nwt = 1e10;
  Real pcTime = state[State_Type].curTime();

  for (int lev=0; lev<nlevs; lev++)
  {
      PorousMedia&    fine_lev   = getLevel(lev);
      const BoxArray& fine_grids = fine_lev.boxArray();      
      MultiFab& S_lev = fine_lev.get_new_data(State_Type);
      MultiFab& P_lev = fine_lev.get_new_data(Press_Type);

      alpha.set(lev,new MultiFab(fine_grids,1,1));
      MultiFab::Copy(alpha[lev],*(fine_lev.rock_phi),0,0,1,1);

      res_fix.set(lev,new MultiFab(fine_grids,1,1));
      if (do_richard_sat_solve) {
	MultiFab::Copy(res_fix[lev],nl_data.initialState[lev],nc,0,1,0);
      }
      else {
	fine_lev.FillStateBndry(pcTime,Press_Type,0,1);
	fine_lev.calcInvPressure(res_fix[lev],nl_data.initialState[lev]);
      }

      for (MFIter mfi(res_fix[lev]); mfi.isValid(); ++mfi)
	res_fix[lev][mfi].mult(alpha[lev][mfi],mfi.validbox(),0,0,1);
      res_fix[lev].mult(-1.0);

      fine_lev.compute_vel_phase(u_mac_local[lev],0,pcTime);
      fine_lev.calc_richard_velbc(res_fix[lev],u_mac_local[lev],dt*density[0]);

      if (do_richard_sat_solve) 
	{
            // FIXME: pulled from above after calcLambda
	  fine_lev.calcCapillary(pcTime);
            // FIXME: in the scalar version, this is done inside the do_richard_sat_solve loop below
            //fine_lev.calcLambda(pcTime);
	  MultiFab::Copy(P_lev,*(fine_lev.pcnp1_cc),0,0,1,1);
	  P_lev.mult(-1.0,1);
	}
  }

  Diffusion::NewtonStepInfo linear_status;
  linear_status.success = true;
  linear_status.max_ls_iterations = richard_max_ls_iterations;
  linear_status.min_ls_factor = richard_min_ls_factor;
  linear_status.ls_acceptance_factor = richard_ls_acceptance_factor;
  linear_status.ls_reduction_factor = richard_ls_reduction_factor;
  linear_status.monitor_linear_solve = richard_monitor_linear_solve;
  linear_status.monitor_line_search = richard_monitor_line_search;

  Real solve_time = 0;

  if (do_richard_sat_solve)
    {
        while ((nl_data.NLIterationsTaken() < nl_data.MaxNLIterations()) && (err_nwt > max_err_nwt) && (linear_status.success)) 
	{
          nl_data++;
	  diffusion->richard_composite_iter(dt,nlevs,nc,gravity,density,res_fix,
					    alpha,cmp_pcp1,cmp_pcp1_dp,u_mac_local,
					    do_upwind,linear_status); 

          err_nwt = linear_status.residual_norm_post_ls;

	  for (int lev=0; lev<nlevs; lev++)
          {
	      PorousMedia&  fine_lev = getLevel(lev);   
	      MultiFab& P_lev        = fine_lev.get_new_data(Press_Type);

	      fine_lev.FillStateBndry(pcTime,State_Type,0,ncomps);
	      fine_lev.calcCapillary(pcTime);
	      fine_lev.calcLambda(pcTime);
	      MultiFab::Copy(P_lev,*(fine_lev.pcnp1_cc),0,0,1,1);
	      P_lev.mult(-1.0,1);

	      MultiFab* tmp_cmp_pcp1[BL_SPACEDIM];
	      MultiFab* tmp_cmp_pcp1_dp[BL_SPACEDIM];
	      for (int dir=0;dir<BL_SPACEDIM;dir++)
		{
		  tmp_cmp_pcp1[dir] = &cmp_pcp1[dir][lev];
		  tmp_cmp_pcp1_dp[dir] = &cmp_pcp1_dp[dir][lev];
		}

	      fine_lev.compute_vel_phase(u_mac_local[lev],0,pcTime);
	      fine_lev.calc_richard_coef(tmp_cmp_pcp1,fine_lev.lambdap1_cc,
					 u_mac_local[lev],0,do_upwind,pcTime);
              if (nl_data.UpdateJacobian(lev)) {
                  fine_lev.calc_richard_jac(tmp_cmp_pcp1_dp,0,fine_lev.lambdap1_cc,
                                            u_mac_local[lev],pcTime,dt,0,do_upwind,do_richard_sat_solve);
              }
 	    }
        }
    }

  else
    {
      while ((nl_data.NLIterationsTaken() < nl_data.MaxNLIterations()) && (err_nwt > max_err_nwt) && (linear_status.success)) 
	{
          nl_data++;

	  bool update_jac = false;
	  for (int lev=0;lev<nlevs;lev++) update_jac = nl_data.UpdateJacobian(lev) || update_jac;
	  const Real tmp_time = ParallelDescriptor::second();
	  diffusion->richard_composite_iter_p(dt,nlevs,nc,gravity,density,res_fix,
                                              alpha,dalpha,cmp_pcp1,cmp_pcp1_dp,
                                              u_mac_local,update_jac,do_upwind,linear_status); 
	  solve_time+= ParallelDescriptor::second() - tmp_time;
          err_nwt = linear_status.residual_norm_post_ls;

          if (richard_solver_verbose>1 && ParallelDescriptor::IOProcessor()) {
              std::cout << "  " << nl_data.NLIterationsTaken() << "  Amanzi-S Newton Function norm " << err_nwt << "\n"; 
          }
	}
    }

  RichardNLSdata::Reason retVal = RichardNLSdata::RICHARD_SUCCESS;
  if (!linear_status.success) {
      retVal = RichardNLSdata::RICHARD_LINEAR_FAIL;
  }
  if (nl_data.NLIterationsTaken() >= nl_data.MaxNLIterations()) {
      retVal = RichardNLSdata::RICHARD_NONLINEAR_FAIL;
      if (richard_solver_verbose>1 && ParallelDescriptor::IOProcessor())
          std::cout << "     **************** Newton failed in richard_composite_update: too many iterations (max = "
                    << nl_data.MaxNLIterations() << ")\n"; 
  }

  if (retVal != RichardNLSdata::RICHARD_SUCCESS) {
      return retVal;
  }

  if (richard_solver_verbose>1 && ParallelDescriptor::IOProcessor() && retVal == RichardNLSdata::RICHARD_SUCCESS) {
      std::cout << "     Newton converged in " << nl_data.NLIterationsTaken() << " iterations (max = "
                << nl_data.MaxNLIterations() << ") with err: " 
                << err_nwt << " (tol = " << max_err_nwt << ")\n";
  }


  Real run_time = ParallelDescriptor::second() - strt_time;
  richard_time = run_time;
  ParallelDescriptor::ReduceRealMax(richard_time);
  ParallelDescriptor::ReduceRealMax(solve_time);
  richard_time_min = std::min(richard_time_min,richard_time);
  
  if (verbose > 3)
    { 
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      ParallelDescriptor::ReduceRealMax(run_time,IOProc);


      if (ParallelDescriptor::IOProcessor())
        std::cout << "PorousMedia::richard_update(): time: " 
		  << run_time << ' ' << solve_time  << '\n';
    }
  //
  // Write out the min and max of each component of the new state.
  //
  if (verbose > 3) check_minmax();

  return retVal;
}

#endif

//
// Enforce the constraint sum_i s_i = 1.  This is achieved by adjusting 
// the saturation of the dominant component specified in the input.
//
void
PorousMedia::scalar_adjust_constraint (int  first_scalar,
				       int  last_scalar)
{
  if (idx_dominant > -1) {
      MultiFab&  S_new = get_new_data(State_Type);
      
      MultiFab S_adj(grids,1,S_new.nGrow());
      MultiFab S_div(grids,1,S_new.nGrow());
      S_adj.setVal(1.0);
      
      for (int kk=first_scalar; kk<=last_scalar; kk++)
      {
          if (solid.compare(pNames[pType[kk]]) != 0 && 
              kk != idx_dominant) 
          {
              MultiFab::Copy(S_div,S_new,kk,0,1,S_div.nGrow());
              S_div.mult(1.0/density[kk],S_div.nGrow());
              S_adj.minus(S_div,0,1,S_adj.nGrow());
          }
      }
      S_adj.mult(density[idx_dominant],S_div.nGrow());
  MultiFab::Copy(S_new,S_adj,0,idx_dominant,1,S_new.nGrow());
  S_new.FillBoundary();
  geom.FillPeriodicBoundary(S_new,true);
  }
}

void
coarsenMask(FArrayBox& crse, const FArrayBox& fine, const IntVect& ratio)
{
    const Box& fbox = fine.box();
    const Box cbox = BoxLib::coarsen(fbox,ratio);
    crse.resize(cbox,1); crse.setVal(0);

    Box b1(BoxLib::refine(cbox,ratio));

    const int* flo      = fbox.loVect();
    const int* fhi      = fbox.hiVect();
    IntVect    d_length = fbox.size();
    const int* flen     = d_length.getVect();
    const int* clo      = cbox.loVect();
    IntVect    cbox_len = cbox.size();
    const int* clen     = cbox_len.getVect();
    const int* lo       = b1.loVect();
    int        longlen  = b1.longside();

    const Real* fdat = fine.dataPtr();
    Real* cdat = crse.dataPtr();

    Array<Real> t(longlen,0);

    int klo = 0, khi = 0, jlo = 0, jhi = 0, ilo, ihi;
    D_TERM(ilo=flo[0]; ihi=fhi[0]; ,
           jlo=flo[1]; jhi=fhi[1]; ,
           klo=flo[2]; khi=fhi[2];)

#define IXPROJ(i,r) (((i)+(r)*std::abs(i))/(r) - std::abs(i))
#define IOFF(j,k,lo,len) D_TERM(0, +(j-lo[1])*len[0], +(k-lo[2])*len[0]*len[1])
   
   int ratiox = 1, ratioy = 1, ratioz = 1;
   D_TERM(ratiox = ratio[0];,
          ratioy = ratio[1];,
          ratioz = ratio[2];)

   for (int k = klo; k <= khi; k++)
   {
       const int kc = IXPROJ(k,ratioz);
       for (int j = jlo; j <= jhi; j++)
       {
           const int   jc = IXPROJ(j,ratioy);
           Real*       c = cdat + IOFF(jc,kc,clo,clen);
           const Real* f = fdat + IOFF(j,k,flo,flen);
           //
           // Copy fine grid row of values into tmp array.
           //
           for (int i = ilo; i <= ihi; i++)
               t[i-lo[0]] = f[i-ilo];

           for (int off = 0; off < ratiox; off++)
           {
               for (int ic = 0; ic < clen[0]; ic++)
               {
                   const int i = ic*ratiox + off;
                   c[ic] = std::max(c[ic],t[i]);
               }
           }
       }
   }

#undef IXPROJ
#undef IOFF
}





//
// Tag cells for refinement
//
void
PorousMedia::errorEst (TagBoxArray& tags,
		       int         clearval,
		       int         tagval,
		       Real        time,
		       int         n_error_buf, 
		       int         ngrow)
{
  const int*  domain_lo = geom.Domain().loVect();
  const int*  domain_hi = geom.Domain().hiVect();
  const Real* dx        = geom.CellSize();
  const Real* prob_lo   = geom.ProbLo();
  const Real* prob_hi   = geom.ProbHi();
  Array<int>  itags;

  //
  // Tag cells for refinement
  //
  for (int j = 0; j < err_list.size(); j++)
  {
      const ErrorRec::ErrorFunc& efunc = err_list[j].errFunc();
      const PM_Error_Value* pmfunc = dynamic_cast<const PM_Error_Value*>(&efunc);
      if (pmfunc==0) 
      {
          MultiFab* mf = derive(err_list[j].name(), time, err_list[j].nGrow());
          
          for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
          {
              RealBox     gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
              itags               = tags[mfi.index()].tags();
              int*        tptr    = itags.dataPtr();
              const int*  tlo     = tags[mfi.index()].box().loVect();
              const int*  thi     = tags[mfi.index()].box().hiVect();
              const int*  lo      = mfi.validbox().loVect();
              const int*  hi      = mfi.validbox().hiVect();
              const Real* xlo     = gridloc.lo();
              Real*       dat     = (*mf)[mfi].dataPtr();
              const int*  dlo     = (*mf)[mfi].box().loVect();
              const int*  dhi     = (*mf)[mfi].box().hiVect();
              const int   ncomp   = (*mf)[mfi].nComp();
              
              err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                    &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                    lo,hi, &ncomp, domain_lo, domain_hi,
                                    dx, xlo, prob_lo, &time, &level);
                      
              //
              // Don't forget to set the tags in the TagBox.
              //
              tags[mfi.index()].tags(itags);
          }
          delete mf;
      }
      else {

          Real min_time = pmfunc->MinTime();
          Real max_time = pmfunc->MaxTime();
          int max_level = pmfunc->MaxLevel();

          if ( (max_level<0) || (max_level>parent->maxLevel()) ) {
              max_level = parent->maxLevel();
          }

          if ( ( (min_time>=max_time) || (min_time<=time) && (max_time>=time) )
               && (level<max_level) )
          {
              IntVect cumRatio = IntVect(D_DECL(1,1,1));
              for (int i=level; i<max_level; ++i) {
                  cumRatio *= parent->refRatio()[i];
              }
              
              const Geometry& fgeom = parent->Geom(max_level);
              const Real* dx_fine = fgeom.CellSize();
              const Real* plo = fgeom.ProbLo();

              const PArray<Region>& my_regions = pmfunc->Regions();

              MultiFab* mf = 0;
              const std::string& name = err_list[j].name();

              if (!pmfunc->regionOnly())
                  mf = derive(err_list[j].name(), time, err_list[j].nGrow());

              FArrayBox mask, cmask;
              for (MFIter mfi(tags); mfi.isValid(); ++mfi)
              {
                  TagBox& tagbox = tags[mfi];
                  const Box fine_box = Box(tagbox.box()).refine(cumRatio);
                  
                  mask.resize(fine_box,1); mask.setVal(0);                  
                  for (int j=0; j<my_regions.size(); ++j) {
                      my_regions[j].setVal(mask,1,0,dx_fine,0);
                  }                  
                  coarsenMask(cmask,mask,cumRatio);

                  if (cmask.max()>0)
                  {
                      itags               = tags[mfi.index()].tags();
                      int*        tptr    = itags.dataPtr();
                      const int*  tlo     = tags[mfi.index()].box().loVect();
                      const int*  thi     = tags[mfi.index()].box().hiVect();
                      const Real* mdat    = cmask.dataPtr();
                      const int*  mlo     = cmask.box().loVect();
                      const int*  mhi     = cmask.box().hiVect();

                      if (pmfunc->regionOnly())
                      {
                          const Box& crse_box = cmask.box();
                          BL_ASSERT(crse_box == tagbox.box());
                          int numPts = crse_box.numPts();
                          for (int i=0; i<numPts; ++i) {
                              if (mdat[i]==1) {
                                  tptr[i] = tagval;
                              }
                          }
                      }
                      else {

                          RealBox     gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
                          const int*  lo      = mfi.validbox().loVect();
                          const int*  hi      = mfi.validbox().hiVect();
                          const Real* xlo     = gridloc.lo();
                          Real*       dat     = (*mf)[mfi].dataPtr();
                          const int*  dlo     = (*mf)[mfi].box().loVect();
                          const int*  dhi     = (*mf)[mfi].box().hiVect();
                          const int   ncomp   = (*mf)[mfi].nComp();
                          
                          Real value = pmfunc->Value();
                          pmfunc->tagCells(tptr,ARLIM(tlo),ARLIM(thi),
                                           &tagval, &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                           mdat, ARLIM(mlo), ARLIM(mhi),
                                           lo,hi, &ncomp, domain_lo, domain_hi,
                                           dx, xlo, prob_lo, &time, &level);
                      }
                          
                      //
                      // Don't forget to set the tags in the TagBox.
                      //
                      tags[mfi.index()].tags(itags);
                  }
              }

              delete mf;
          }
      }
  }
#if 0
  //
  // Tag cells for refinement based on permeability values
  //
  if (do_kappa_refine == 1)
    { 
      Real kpset = 1.e-6;
      Array<int> itags;
      
      for (MFIter mfi(*kappa); mfi.isValid(); ++mfi)
	{
	  const int* k_lo  = (*kappa)[mfi].loVect();
	  const int* k_hi  = (*kappa)[mfi].hiVect();
	  const Real* kdat = (*kappa)[mfi].dataPtr();

	  itags            = tags[mfi.index()].tags();
	  const int* t_lo  = tags[mfi.index()].box().loVect();
	  const int* t_hi  = tags[mfi.index()].box().hiVect();
	  const int* tdat  = itags.dataPtr();

	  const int*  lo   = mfi.validbox().loVect();
	  const int*  hi   = mfi.validbox().hiVect();
	
	  FORT_KPERROR(tdat,ARLIM(t_lo),ARLIM(t_hi),
		       kdat,ARLIM(k_lo),ARLIM(k_hi),
		       &tagval, &kpset, dx, prob_lo, prob_hi,
		       lo, hi, domain_lo, domain_hi, &level);
	
	  tags[mfi.index()].tags(itags);
	}
    }
#endif
}

Real
PorousMedia::sumDerive (const std::string& name, Real time)
{
    Real      sum = 0.0;
    MultiFab* mf  = derive(name,time,0);

    BL_ASSERT(!(mf == 0));

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = mf->get(mfi);

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0, N = isects.size(); ii < N; ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }

        sum += fab.sum(0);
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
PorousMedia::volWgtSum (const std::string& name,
			Real           time)
{
  Real        sum     = 0;
  const Real* dx      = geom.CellSize();
  MultiFab*   mf      = derive(name,time,0);

  BoxArray baf;

  if (level < parent->finestLevel())
  {
      baf = parent->boxArray(level+1);
      baf.coarsen(fine_ratio);
  }

  for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
      FArrayBox& fab = (*mf)[mfi];

      if (level < parent->finestLevel())
        {
            if (level < parent->finestLevel())
            {
                std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

                for (int ii = 0, N = isects.size(); ii < N; ii++)
                {
                    fab.setVal(0,isects[ii].second,0,fab.nComp());
                }
            }
        }
      Real        s;
      const Real* dat = fab.dataPtr();
      const int*  dlo = fab.loVect();
      const int*  dhi = fab.hiVect();
      const int*  lo  = grids[mfi.index()].loVect();
      const int*  hi  = grids[mfi.index()].hiVect();

      FORT_SUMMASS(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),dx,&s);

      sum += s;
    }

  delete mf;

  ParallelDescriptor::ReduceRealSum(sum);

  return sum;
}

void
PorousMedia::sum_integrated_quantities ()
{
  const int finest_level = parent->finestLevel();

  Real time = state[State_Type].curTime();
  Real mass = 0.0;

  for (int lev = 0; lev <= finest_level; lev++)
    {
      PorousMedia& ns_level = getLevel(lev);
      mass += ns_level.volWgtSum("Water",time);
    }

  if (verbose > 3 && ParallelDescriptor::IOProcessor())
    {
      const int old_prec = std::cout.precision(12);
      std::cout << "TIME= " << time << " MASS= " << mass << '\n';
      std::cout.precision(old_prec);
    }
}

void
PorousMedia::setPlotVariables()
{
    ParmParse pp("amr");

    // By default, do not add state variables
    if (pp.contains("plot_vars"))
    {
        std::string nm;
      
        int nPltVars = pp.countval("plot_vars");
      
        for (int i = 0; i < nPltVars; i++)
        {
            pp.get("plot_vars", nm, i);

            if (nm == "ALL") 
                parent->fillStatePlotVarList();
            else if (nm == "NONE")
                parent->clearStatePlotVarList();
            else
                parent->addStatePlotVar(nm);
        }
    }

    // Search for "ALL" in list
    bool has_all = false;
    Array<std::string> names_to_derive;

    if (pp.contains("derive_plot_vars"))
    {
        std::string nm;
      
        int nDrvPltVars = pp.countval("derive_plot_vars");
        names_to_derive.resize(nDrvPltVars);
        pp.getarr("derive_plot_vars",names_to_derive,0,nDrvPltVars);
      
        for (int i = 0; i < nDrvPltVars; i++)
        {
            if (names_to_derive[i] == "ALL") 
                has_all = true;
        }
    }

    if (has_all || names_to_derive.size()==0) {
        names_to_derive = UserDerives();
    }

    for (int i=0; i<names_to_derive.size(); ++i) {

        const std::string name = names_to_derive[i];

        if (derive_lst.canDerive(name)) {

            if (derive_lst.get(name)->deriveType() == IndexType::TheCellType())
            {
                parent->addDerivePlotVar(name);
            }
        }
    }
}

std::string
PorousMedia::thePlotFileType () const
{
  //
  // Increment this whenever the writePlotFile() format changes.
  //
  static const std::string the_plot_file_type("PorousMedia-V1.1");

  return the_plot_file_type;
}

void
PorousMedia::writePlotFile (const std::string& dir,
			    std::ostream&  os,
			    VisMF::How     how)
{
  if ( ! Amr::Plot_Files_Output() ) return;

  int i, n;
  //
  // The list of indices of State to write to plotfile.
  // first component of pair is state_type,
  // second component of pair is component # within the state_type
  //
  std::vector<std::pair<int,int> > plot_var_map;

  int noutput = desc_lst.size();
  for (int typ = 0; typ < noutput; typ++)
    for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
      if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
	  desc_lst[typ].getType() == IndexType::TheCellType())
	plot_var_map.push_back(std::pair<int,int>(typ,comp));

  int num_derive = 0;
  std::list<std::string> derive_names;
  const std::list<DeriveRec>& dlist = derive_lst.dlist();

  for (std::list<DeriveRec>::const_iterator it = dlist.begin();
       it != dlist.end();
       ++it)
    {
      if (parent->isDerivePlotVar(it->name()))
	{
	  derive_names.push_back(it->name());
	  num_derive += it->numDerive();
	}
    }

  int n_data_items = plot_var_map.size() + num_derive;
  Real cur_time = state[State_Type].curTime();

  if (level == 0 && ParallelDescriptor::IOProcessor())
    {
      //
      // The first thing we write out is the plotfile type.
      //
      os << thePlotFileType() << '\n';

      if (n_data_items == 0)
	BoxLib::Error("Must specify at least one valid data item to plot");

      os << n_data_items << '\n';

      //
      // Names of variables -- first state, then derived
      //
      for (i =0; i < plot_var_map.size(); i++)
        {
	  int typ  = plot_var_map[i].first;
	  int comp = plot_var_map[i].second;
	  os << desc_lst[typ].name(comp) << '\n';
        }

      for (std::list<std::string>::const_iterator it = derive_names.begin();
	   it != derive_names.end();
	   ++it)
        {
	  const DeriveRec* rec = derive_lst.get(*it);
	  for (i = 0; i < rec->numDerive(); i++)
	    os << rec->variableName(i) << '\n';
        }
      os << BL_SPACEDIM << '\n';
      os << parent->cumTime() << '\n';
      int f_lev = parent->finestLevel();
      os << f_lev << '\n';
      for (i = 0; i < BL_SPACEDIM; i++)
	os << Geometry::ProbLo(i) << ' ';
      os << '\n';
      for (i = 0; i < BL_SPACEDIM; i++)
	os << Geometry::ProbHi(i) << ' ';
      os << '\n';
      for (i = 0; i < f_lev; i++)
	os << parent->refRatio(i)[0] << ' ';
      os << '\n';
      for (i = 0; i <= f_lev; i++)
	os << parent->Geom(i).Domain() << ' ';
      os << '\n';
      for (i = 0; i <= f_lev; i++)
	os << parent->levelSteps(i) << ' ';
      os << '\n';
      for (i = 0; i <= f_lev; i++)
        {
	  for (int k = 0; k < BL_SPACEDIM; k++)
	    os << parent->Geom(i).CellSize()[k] << ' ';
	  os << '\n';
        }
      os << (int) Geometry::Coord() << '\n';
      os << "0\n"; // Write bndry data.
    }
  // Build the directory to hold the MultiFab at this level.
  // The name is relative to the directory containing the Header file.
  //
  static const std::string BaseName = "/Cell";

  std::string Level = BoxLib::Concatenate("Level_", level, 1);
  //
  // Now for the full pathname of that directory.
  //
  std::string FullPath = dir;
  if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
    FullPath += '/';
  FullPath += Level;
  //
  // Only the I/O processor makes the directory if it doesn't already exist.
  //
  if (ParallelDescriptor::IOProcessor())
    if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
      BoxLib::CreateDirectoryFailed(FullPath);
  //
  // Force other processors to wait till directory is built.
  //
  ParallelDescriptor::Barrier();

  if (ParallelDescriptor::IOProcessor())
    {
      os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
      os << parent->levelSteps(level) << '\n';

      for (i = 0; i < grids.size(); ++i)
        {
	  RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
	  for (n = 0; n < BL_SPACEDIM; n++)
	    os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
      //
      // The full relative pathname of the MultiFabs at this level.
      // The name is relative to the Header file containing this name.
      // It's the name that gets written into the Header.
      //
      if (n_data_items > 0)
        {
	  std::string PathNameInHeader = Level;
	  PathNameInHeader += BaseName;
	  os << PathNameInHeader << '\n';
        }
    }

  //
  // We combine all of the multifabs -- state, derived, etc -- into one
  // multifab -- plotMF.
  // NOTE: we are assuming that each state variable has one component,
  // but a derived variable is allowed to have multiple components.
  int       cnt   = 0;
  int       ncomp = 1;
  const int nGrow = 0;
  MultiFab  plotMF(grids,n_data_items,nGrow);
  MultiFab* this_dat = 0;
  //
  // Cull data from state variables -- use no ghost cells.
  //
  for (i = 0; i < plot_var_map.size(); i++)
    {
      int typ  = plot_var_map[i].first;
      int comp = plot_var_map[i].second;
      this_dat = &state[typ].newData();
      MultiFab::Copy(plotMF,*this_dat,comp,cnt,ncomp,nGrow);
      cnt+= ncomp;
    }
  //
  // Cull data from derived variables.
  // 
  Real plot_time;

  if (derive_names.size() > 0)
    {
      for (std::list<std::string>::const_iterator it = derive_names.begin();
	   it != derive_names.end();
	   ++it) 
	{
	  plot_time = cur_time;
	  const DeriveRec* rec = derive_lst.get(*it);
	  ncomp = rec->numDerive();
	  MultiFab* derive_dat = derive(*it,plot_time,nGrow);
	  MultiFab::Copy(plotMF,*derive_dat,0,cnt,ncomp,nGrow);
	  delete derive_dat;
	  cnt += ncomp;
	}
    }
  //
  // Use the Full pathname when naming the MultiFab.
  //
  std::string TheFullPath = FullPath;
  TheFullPath += BaseName;
  VisMF::Write(plotMF,TheFullPath,how,true);
}

Real
PorousMedia::estTimeStep (MultiFab* u_mac)
{
  if (fixed_dt > 0.0)
    {
      Real factor = 1.0;

      if (!(level == 0))
        {
	  int ratio = 1;
	  for (int lev = 1; lev <= level; lev++)
            {
	      ratio *= parent->nCycle(lev);
            }
	  factor = 1.0/double(ratio);
        }

      return factor*fixed_dt;
    }

  Real estdt        = 1.0e+20; // FIXME: need more robust
  const Real cur_time = state[State_Type].curTime();

  if (dt_eig != 0.0)
    {
        if (cfl>0) {
            estdt = cfl * dt_eig;
        }
	else
	  estdt = dt_eig;
    } 
  else 
    {
      int making_new_umac = 0;
      
      // Need to define the MAC velocities in order to define the initial dt 
      if (u_mac == 0) 
	{
	  making_new_umac = 1;

	  u_mac = new MultiFab[BL_SPACEDIM];
	  for (int dir = 0; dir < BL_SPACEDIM; dir++)
	    {
	      BoxArray edge_grids(grids);
	      edge_grids.surroundingNodes(dir);
	      u_mac[dir].define(edge_grids,1,0,Fab_allocate);
	      u_mac[dir].setVal(0.);
	    }
#ifdef MG_USE_FBOXLIB
	  if (model == model_list["richard"]) {
              if (!steady_use_PETSc_snes) {
                  compute_vel_phase(u_mac,0,cur_time);
              }
          }
	  else
#endif
	    {
	      MultiFab* RhoD;
	      RhoD  = new MultiFab[BL_SPACEDIM];
	      for (int dir = 0; dir < BL_SPACEDIM; dir++)
		{
		  BoxArray edge_grids(grids);
		  edge_grids.surroundingNodes(dir);
		  RhoD[dir].define(edge_grids,1,0,Fab_allocate);
		  RhoD[dir].setVal(0.);
		}

	      initial_mac_project(u_mac,RhoD,cur_time);
	      delete [] RhoD;
	    }
	}


      predictDT(u_mac);

      if (cfl>0) {
          estdt = cfl*dt_eig;
      }
      if (making_new_umac)
	delete [] u_mac;
    }

  // 
  // Limit by max_dt
  //
#ifdef MG_USE_FBOXLIB
  if (model == model_list["richard"]) {
      Real richard_max_dt = (initial_iter  ?  steady_richard_max_dt  :  transient_richard_max_dt);
      if (richard_max_dt>0) {
          estdt = std::min(richard_max_dt,estdt);
      }
  }  
#endif

  return estdt;
}

Real
PorousMedia::initialTimeStep (MultiFab* u_mac)
{
    Real dt_0;

    if (dt_init>0) {
        dt_0 = dt_init;
    }
    else {
        dt_0 = estTimeStep(u_mac);
    }

    std::cout << "initialTimeStep:  dt_0: " << dt_0 << ", " << dt_init << std::endl;

    const Real cur_time = state[State_Type].curTime();
    if (stop_time > cur_time) {
        dt_0 = std::min(dt_0, stop_time - cur_time);
    }
    
    if (init_shrink>0) 
        dt_0 *= init_shrink;
    return dt_0;
}

void 
PorousMedia::predictDT (MultiFab* u_macG)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::predictDT()");

  const Real* dx       = geom.CellSize();
  const Real  cur_time = state[State_Type].curTime();

  dt_eig = 1.e20; // FIXME: Need more robust

  Real eigmax[BL_SPACEDIM] = { D_DECL(0,0,0) };
  for (FillPatchIterator S_fpi(*this,get_new_data(State_Type),GEOM_GROW,
			       cur_time,State_Type,0,ncomps);
       S_fpi.isValid();
       ++S_fpi)
    {
      const int i = S_fpi.index();

      Array<int> state_bc;
      state_bc = getBCArray(State_Type,i,0,1);

      Real eigmax_m[BL_SPACEDIM] = {D_DECL(0,0,0)};
      
      if (model == model_list["single-phase"])
	{
	  godunov->esteig_lin (grids[i], u_macG[0][i], u_macG[1][i],
#if (BL_SPACEDIM == 3)    
			       u_macG[2][i],
#endif
			       (*rock_phi)[i], eigmax_m);
	}
      else if (model == model_list["two-phase"])
	{
	  const int n_kr_coef = kr_coef->nComp();
	  if (do_cpl_advect)
	    {
	      godunov->esteig_cpl (grids[i], dx,
				   u_macG[0][i],kpedge[0][i],
				   u_macG[1][i],kpedge[1][i],
#if (BL_SPACEDIM == 3)    
				   u_macG[2][i],kpedge[2][i],
#endif
				   S_fpi(), (*pcnp1_cc)[i],
				   (*rock_phi)[i], 
				   (*kr_coef)[i], n_kr_coef,
				   state_bc.dataPtr(),eigmax_m);
	    }
	  else
	    godunov->esteig (grids[i], dx,
			     u_macG[0][i],kpedge[0][i],
			     u_macG[1][i],kpedge[1][i],
#if (BL_SPACEDIM == 3)    
			     u_macG[2][i],kpedge[2][i],
#endif
			     S_fpi(),(*rock_phi)[i], 
			     (*kr_coef)[i], n_kr_coef,
			     state_bc.dataPtr(),eigmax_m);
	}
    
      if (do_tracer_transport && ntracers > 0)
	{
	  godunov->esteig_trc (grids[i], u_macG[0][i], u_macG[1][i],
#if (BL_SPACEDIM == 3)    
			       u_macG[2][i],
#endif
			       S_fpi(),1,(*rock_phi)[i] ,eigmax_m);
	}

      for (int dir = 0; dir < BL_SPACEDIM; dir++)
	{
	  eigmax[dir] = std::max(eigmax[dir],eigmax_m[dir]);
	  if (eigmax_m[dir] > 1.e-15) 
	    dt_eig = std::min(dt_eig,dx[dir]/eigmax_m[dir]);

	}
    }
  
  ParallelDescriptor::ReduceRealMin(dt_eig);

  if (verbose > 3)
    {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      ParallelDescriptor::ReduceRealMax(&eigmax[0], BL_SPACEDIM, IOProc);

      if (ParallelDescriptor::IOProcessor())
	{
	  for (int dir = 0; dir < BL_SPACEDIM; dir++)
	    std::cout << "Max Eig in dir " << dir << " = " << eigmax[dir] << '\n';
	  std::cout << "Max timestep = " << dt_eig << '\n';
	}
    }
}

void
PorousMedia::computeNewDt (int                   finest_level,
                           int                   sub_cycle,
                           Array<int>&           n_cycle,
                           const Array<IntVect>& ref_ratio,
                           Array<Real>&          dt_min,
                           Array<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
  //
  // We are at the end of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  //
  if (level > 0) return;

  const int max_level = parent->maxLevel();

  n_cycle[0] = 1;
  for (int i = 1; i <= max_level; i++)
    {
      n_cycle[i] = sub_cycle ? parent->MaxRefRatio(i-1) : 1;
    }

  Real dt_0     = 1.0e20;
  int  n_factor = 1;
    
  for (int i = 0; i <= finest_level; i++)
    {
      PorousMedia* pm = dynamic_cast<PorousMedia*>(&parent->getLevel(i));
      dt_min[i] = std::min(dt_min[i],getLevel(i).estTimeStep(pm->u_mac_curr));
    }

  if (fixed_dt <= 0.0) 
    {
      if (post_regrid_flag == 1)
	{
          //
          // Limit dt's by pre-regrid dt
          //
          for (int i = 0; i <= finest_level; i++)
	    {
              dt_min[i] = std::min(dt_min[i],dt_level[i]);
	    }
	}
      else
	{
          //
          // Limit dt's by change_max * old dt if dt_(grow,shrink)_max>0
          //
            for (int i = 0; i <= finest_level; i++)
            {
                if (dt_grow_max >= 1) {
                    dt_min[i] = std::min(dt_min[i],dt_grow_max* dt_level[i]);
                }
                if (dt_shrink_max <= 1) {
                    dt_min[i] = std::max(dt_min[i],dt_shrink_max * dt_level[i]);
                }
	    }
	}
    }

  //
  // Find the minimum over all levels
  //
  for (int i = 0; i <= finest_level; i++)
    {
      n_factor *= n_cycle[i];
      dt_0      = std::min(dt_0,n_factor*dt_min[i]);
    }
  // 
  // Limit by max_dt: redundant?
  //
#ifdef MG_USE_FBOXLIB
  if (model == model_list["richard"]) {
      Real richard_max_dt = (initial_iter  ?  steady_richard_max_dt  :  transient_richard_max_dt);
      if (richard_max_dt>0) {
          dt_0 = std::min(richard_max_dt,dt_0);
      }
  }  
#endif


  n_factor = 1;
  for (int i = 0; i <= max_level; i++)
    {
      n_factor   *= n_cycle[i];
      dt_level[i] = dt_0/( (Real)n_factor );
    }
}

void
PorousMedia::computeInitialDt (int                   finest_level,
                               int                   sub_cycle,
                               Array<int>&           n_cycle,
                               const Array<IntVect>& ref_ratio,
                               Array<Real>&          dt_level, 
                               Real                  stop_time)
{
  //
  // Grids have been constructed, compute dt for all levels.
  //
  if (level > 0)
    return;

  if (verbose>2 && ParallelDescriptor::IOProcessor())
    std::cout << "... computing dt at level 0 only in computeInitialDt\n";

  const int max_level = parent->maxLevel();

  n_cycle[0] = 1;
  for (int i = 1; i <= max_level; i++)
    {
      n_cycle[i] = sub_cycle ? parent->MaxRefRatio(i-1) : 1;
    }

  Real dt_0;
  if (dt_init < 0) {
      dt_0 = 1.0e100;
      int n_factor = 1;
      for (int i = 0; i <= finest_level; i++)
      {
          
          const PorousMedia* pm = dynamic_cast<const PorousMedia*>(&parent->getLevel(i));
          dt_level[i] = getLevel(i).initialTimeStep(pm->u_mac_curr);
          n_factor   *= n_cycle[i];
          dt_0        = std::min(dt_0,n_factor*dt_level[i]);
      }
  }
  else {
      dt_0 = dt_init;
  }

  // FIXME: should be stop_time >= start_time
  if (stop_time >= 0.0)
    {
      const Real eps      = 0.0001*dt_0;
      const Real cur_time = state[State_Type].curTime();
      if ((cur_time + dt_0) > (stop_time - eps))
	dt_0 = stop_time - cur_time;
    }

  if (init_shrink>0) 
      dt_0 *= init_shrink;

  int n_factor = 1;
  for (int i = 0; i <= max_level; i++)
    {
      n_factor   *= n_cycle[i];
      dt_level[i] = dt_0/( (Real)n_factor );
    }
}

//
// This function estimates the initial timesteping used by the model.
//

void
PorousMedia::post_init_estDT (Real&        dt_init_local,
                              Array<int>&  nc_save,
                              Array<Real>& dt_save,
                              Real         stop_time)
{
  const Real strt_time    = state[State_Type].curTime();
  const int  finest_level = parent->finestLevel();

  if (verbose>2 && ParallelDescriptor::IOProcessor())
    std::cout << "... computing dt at all levels in post_init_estDT\n";

  if (dt_init > 0) {

      dt_init_local = dt_init;

  }
  else {

      // FIXME: should be std::max(0,stop_time - start_time)
      dt_init_local = std::abs(stop_time);

      // Create a temporary data structure for this solve -- this u_mac just
      //   used to compute dt.
      
      int n_factor;
      MultiFab* umac = 0;
      for (int k = 0; k <= finest_level; k++)
      {
          dt_save[k] = getLevel(k).initialTimeStep(umac);
          
          n_factor   = 1;
          for (int m = finest_level; m > k; m--) {
              n_factor *= parent->nCycle(m);
          }
          dt_init_local = std::min( dt_init_local, dt_save[k]/((Real) n_factor) );
      }
      
      dt_init_local *= n_factor;

  }

  // Make something workable if stop>=start
  if (stop_time <= 0.0)
  {
      dt_init_local = std::abs(dt_init);
  }

  if (init_shrink>0) 
      dt_init_local *= init_shrink;

  BL_ASSERT(dt_init_local != 0);

  int n_factor = 1;
  dt_save[0] = dt_init_local;
  for (int k = 0; k <= finest_level; k++)
  {
      nc_save[k] = parent->nCycle(k);
      n_factor  *= nc_save[k];
      dt_save[k] = dt_init_local/( (Real) n_factor);
  }

  parent->setDtLevel(dt_save);
  parent->setNCycle(nc_save);
  for (int k = 0; k <= finest_level; k++)
    {
      getLevel(k).setTimeLevel(strt_time,dt_init_local,dt_init_local);
    }
}

//
// Fills in amrLevel okToContinue.
//

int
PorousMedia::okToContinue ()
{
    bool ret = true;
    std::string reason_for_stopping;
    bool successfully_completed = false;

    if (level == 0) {
        if (parent->dtLevel(0) <= dt_cutoff) {
            ret = false; reason_for_stopping = "Dt at level 0 too small";
        }

        if (parent->levelSteps(0) >= max_step) {
            ret = false; reason_for_stopping = "Hit maximum allowed time steps";
            successfully_completed = true;
        }

        if (parent->cumTime() >= stop_time) {
            ret = false; reason_for_stopping = "Hit maximum allowed time";
            successfully_completed = true;
        }

        if (!ret) {
            //
            // Print final solutions
            //
            if (verbose > 3)
            {      
                for (int lev = 0; lev <= parent->finestLevel(); lev++)
                {
                    if (verbose>2 && ParallelDescriptor::IOProcessor())
                        std::cout << "Final solutions at level = " 
                                  << lev << '\n';
                    
                    getLevel(lev).check_minmax(); 
                    
                }
            }
            
            //
            // Compute observations
            //
            Observation::setAmrPtr(parent);
            if (successfully_completed  &&  ParallelDescriptor::IOProcessor()) 
            {
                if (verbose>1)
                {
                    std::cout << "Computed observations:\n";
                    for (int i=0; i<observations.size(); ++i)
                    {
                        const std::map<int,Real> vals = observations[i].vals;
                        for (std::map<int,Real>::const_iterator it=vals.begin();it!=vals.end(); ++it) 
                        {
                            int j = it->first;
                            std::cout << i << " " << observations[i].name << " " 
                                      << j << " " << observations[i].times[j] << " "
                                      << it->second << std::endl;
                        }
                    }
                    std::cout << "\n";

                }
      
                std::ofstream out;
                out.open(obs_outputfile.c_str(),std::ios::out);
                out.precision(16);
                out.setf(std::ios::scientific);
                for (int i=0; i<observations.size(); ++i)
                {
                    const std::map<int,Real>& vals = observations[i].vals;
                    for (std::map<int,Real>::const_iterator it=vals.begin();it!=vals.end(); ++it) 
                    {
                        int j = it->first;
                        out << i << " " << observations[i].name << " " 
                            << j << " "  << observations[i].times[j] << " "
                            << it->second << std::endl;
                    }
                }
                out.close();
            }
        }
    }
    
    if (!ret && ParallelDescriptor::IOProcessor()) {
        std::cout << "Stopping simulation: " << reason_for_stopping << std::endl;
    }
    return ret;
}

//
// THE MAIN HOOKS INTO AMR AND AMRLEVEL
//

//
// Integration cycle on fine level grids is complete .
// post_timestep() is responsible for syncing levels together.
//
// The registers used for level syncing are initialized in the
// coarse level advance and incremented in the fine level advance.
// These quantities are described in comments above advance_setup.
//

void
PorousMedia::post_timestep (int crse_iteration)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_timestep()");

  const int finest_level = parent->finestLevel();

  if (!(do_multilevel_full)  &&  level < finest_level)
    {
      //avgDown();
      
      if (do_reflux) 
	{
	  reflux();
#ifdef MG_USE_FBOXLIB
	  if (model == model_list["richard"])
	    richard_sync();
	  else
#endif
	    mac_sync();
	}
    }

  //
  // Test for conservation.
  //
  if (level==0 && sum_interval>0 && 
      parent->levelSteps(0)%sum_interval == 0)
    sum_integrated_quantities();
    
  old_intersect_new          = grids;
  is_first_step_after_regrid = false;

  if (level == 0)
    {
      Observation::setAmrPtr(parent);
      Real prev_time = state[State_Type].prevTime();
      Real curr_time = state[State_Type].curTime();
      for (int i=0; i<observations.size(); ++i)
          observations[i].process(prev_time, curr_time, parent->levelSteps(0));
    }


}

PMAmr*
PorousMedia::PMParent()
{
    PMAmr* pm_parent = dynamic_cast<PMAmr*>(parent);
    if (!pm_parent) {
        BoxLib::Abort("Bad cast");
    }
    return pm_parent;
}

//
// Build any additional data structures after restart.
//
void PorousMedia::post_restart()
{
  if (level==0)
  {
      PMAmr::GetLayout().Rebuild();
  }

  if (level == 0)
    {
      Observation::setAmrPtr(parent);
      Real prev_time = state[State_Type].prevTime();
      Real curr_time = state[State_Type].curTime();
      for (int i=0; i<observations.size(); ++i)
          observations[i].process(prev_time, curr_time, parent->levelSteps(0));
    }

}

//
// Build any additional data structures after regrid.
//
void
PorousMedia::post_regrid (int lbase,
                          int new_finest)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_regrid()");

  if (lbase==0)
  {
      PMAmr::GetLayout().Rebuild();
  }

  //if (level > lbase)
  {
    //
    // Alloc MultiFab to hold rock quantities
    //
    if (kpedge   == 0) {
      kpedge = new MultiFab[BL_SPACEDIM];
      for (int dir = 0; dir < BL_SPACEDIM; dir++)
	{
	  BoxArray edge_grids(grids);
	  edge_grids.surroundingNodes(dir).grow(1);
	  kpedge[dir].define(edge_grids,1,0,Fab_allocate);
	}
    }	      
  }
}

void 
PorousMedia::init_rock_properties ()
{

  //
  // Determine rock properties.
  //
  const Real* dx = geom.CellSize();
  const int* domain_hi = geom.Domain().hiVect();

  const int max_level = parent->maxLevel();
  const Geometry& fgeom  = parent->Geom(max_level);

  int fratio = fine_ratio[0];
  int twoexp = 1;
  int ng_twoexp = 1;
  for (int ii = 0; ii<max_level; ii++) 
    {
      if (ii >= level) twoexp *= parent->refRatio(ii)[0];
      ng_twoexp *= parent->refRatio(ii)[0];
    }	
  ng_twoexp = ng_twoexp*3;

  int curr_grid_size = parent->maxGridSize(level);
  int new_grid_size  = 4;
  if (twoexp < curr_grid_size)
    new_grid_size  = curr_grid_size/twoexp;


  // permeability
  
  if (permeability_from_fine)
  {
      BoxArray tba(grids);
      tba.maxSize(new_grid_size);
      MultiFab tkappa(tba,1,3);
      tkappa.setVal(1.e40);
      
      MultiFab* tkpedge;
      tkpedge = new MultiFab[BL_SPACEDIM];
      for (int dir = 0; dir < BL_SPACEDIM; dir++)
      {
	  BoxArray tbe(tba);
	  tbe.surroundingNodes(dir).grow(1);
	  tkpedge[dir].define(tbe,1,0,Fab_allocate);
	  tkpedge[dir].setVal(1.e40);
      }
      
      BoxArray ba(tkappa.size());
      BoxArray ba2(tkappa.size());
      for (int i = 0; i < ba.size(); i++)
      {
	  Box bx = tkappa.box(i);
	  bx.refine(twoexp);
	  ba.set(i,bx);
          
	  bx.grow(ng_twoexp);
	  ba2.set(i,bx);
      }
      
      MultiFab mftmp(ba2,BL_SPACEDIM,0);
      mftmp.copy(*kappadata); 
      
      // mfbig has same CPU distribution as kappa
      MultiFab mfbig_kappa(ba,BL_SPACEDIM,ng_twoexp); 
      for (MFIter mfi(mftmp); mfi.isValid(); ++mfi)
          mfbig_kappa[mfi].copy(mftmp[mfi]);
      mftmp.clear();
      mfbig_kappa.FillBoundary();
      fgeom.FillPeriodicBoundary(mfbig_kappa,true);
      
      for (MFIter mfi(tkappa); mfi.isValid(); ++mfi)
      {
	  const int* lo    = mfi.validbox().loVect();
	  const int* hi    = mfi.validbox().hiVect();
          
	  const int* k_lo  = tkappa[mfi].loVect();
	  const int* k_hi  = tkappa[mfi].hiVect();
	  const Real* kdat = tkappa[mfi].dataPtr();

	  const int* kx_lo  = tkpedge[0][mfi].loVect();
	  const int* kx_hi  = tkpedge[0][mfi].hiVect();
	  const Real* kxdat = tkpedge[0][mfi].dataPtr();

	  const int* ky_lo  = tkpedge[1][mfi].loVect();
	  const int* ky_hi  = tkpedge[1][mfi].hiVect();
	  const Real* kydat = tkpedge[1][mfi].dataPtr();

#if(BL_SPACEDIM==3)
	  const int* kz_lo  = tkpedge[2][mfi].loVect();
	  const int* kz_hi  = tkpedge[2][mfi].hiVect();
	  const Real* kzdat = tkpedge[2][mfi].dataPtr();
#endif

	  const int* mf_lo  = mfbig_kappa[mfi].loVect();
	  const int* mf_hi  = mfbig_kappa[mfi].hiVect();
	  const Real* mfdat = mfbig_kappa[mfi].dataPtr();

	  FORT_INITKAPPA2(mfdat,ARLIM(mf_lo),ARLIM(mf_hi),
			  kdat,ARLIM(k_lo),ARLIM(k_hi),
			  kxdat,ARLIM(kx_lo),ARLIM(kx_hi),
			  kydat,ARLIM(ky_lo),ARLIM(ky_hi),
#if(BL_SPACEDIM==3)
			  kzdat,ARLIM(kz_lo),ARLIM(kz_hi),
#endif		      
			  lo,hi,&level,&max_level, &fratio);
      }
      
      mfbig_kappa.clear();
      
      for (int dir = 0; dir < BL_SPACEDIM; dir++)
          kpedge[dir].copy(tkpedge[dir]);
      delete [] tkpedge;
      
      BoxArray tba2(tkappa.boxArray());
      tba2.grow(3);
      MultiFab tmpgrow(tba2,1,0);
    
      for (MFIter mfi(tkappa); mfi.isValid(); ++mfi)
          tmpgrow[mfi].copy(tkappa[mfi]);

      tkappa.clear();

      tba2 = kappa->boxArray();
      tba2.grow(3);
      MultiFab tmpgrow2(tba2,1,0);

      tmpgrow2.copy(tmpgrow);
      tmpgrow.clear();

      for (MFIter mfi(tmpgrow2); mfi.isValid(); ++mfi)
          (*kappa)[mfi].copy(tmpgrow2[mfi]);
  }

  else 
  {
      BoxLib::Abort("rock layers not yet implemented");
#if 0
      int nlayer = rock_array.size();
      Array<Real> kappaval_x(nlayer), kappaval_y(nlayer), kappaval_z(nlayer);
      int mediumtype = 0;
      for (int i=0;i<nlayer;i++)
	{
	  kappaval_x[i] = rock_array[i].permeability[0];
	  kappaval_y[i] = rock_array[i].permeability[1];
#if(BL_SPACEDIM==3)     
	  kappaval_z[i] = rock_array[i].permeability[2];
#endif
	}

      for (MFIter mfi(*kappa); mfi.isValid(); ++mfi)
	{
	  const int* lo    = mfi.validbox().loVect();
	  const int* hi    = mfi.validbox().hiVect();
	  
	  const int* k_lo  = (*kappa)[mfi].loVect();
	  const int* k_hi  = (*kappa)[mfi].hiVect();
	  const Real* kdat = (*kappa)[mfi].dataPtr();
	  
	  const int* kx_lo  = kpedge[0][mfi].loVect();
	  const int* kx_hi  = kpedge[0][mfi].hiVect();
	  const Real* kxdat = kpedge[0][mfi].dataPtr();

	  const int* ky_lo  = kpedge[1][mfi].loVect();
	  const int* ky_hi  = kpedge[1][mfi].hiVect();
	  const Real* kydat = kpedge[1][mfi].dataPtr();

#if(BL_SPACEDIM==3)
	  const int* kz_lo  = kpedge[2][mfi].loVect();
	  const int* kz_hi  = kpedge[2][mfi].hiVect();
	  const Real* kzdat = kpedge[2][mfi].dataPtr();
#endif
	  FORT_INITKAPPA(kdat,ARLIM(k_lo),ARLIM(k_hi),
			 kxdat,ARLIM(kx_lo),ARLIM(kx_hi),
			 kydat,ARLIM(ky_lo),ARLIM(ky_hi),
#if(BL_SPACEDIM==3)
			 kzdat,ARLIM(kz_lo),ARLIM(kz_hi),
#endif		      
			 lo,hi,dx,geom.ProbHi(),
			 &level,&max_level,&mediumtype,
			 kappaval_x.dataPtr(), kappaval_y.dataPtr(),
#if (BL_SPACEDIM==3)
			 kappaval_z.dataPtr(),
#endif
			 &nlayer, &fratio);
        }
#endif
  }

  kappa->FillBoundary();
  (*kpedge).FillBoundary();
   
  // porosity

  if (porosity_from_fine) 
    {      
      BoxArray tba(grids);
      tba.maxSize(new_grid_size);
      MultiFab trock_phi(tba,1,3);
      trock_phi.setVal(1.e40);

      BoxArray ba(trock_phi.size());
      BoxArray ba2(trock_phi.size());
      for (int i = 0; i < ba.size(); i++)
	{
	  Box bx = trock_phi.box(i);
	  bx.refine(twoexp);
	  ba.set(i,bx);
	  bx.grow(ng_twoexp);
	  ba2.set(i,bx);
	}

      MultiFab mftmp(ba2,1,0);      
      mftmp.copy(*phidata);     
      
      // mfbig has same CPU distribution as phi
      MultiFab mfbig_phi(ba,1,ng_twoexp);
      for (MFIter mfi(mftmp); mfi.isValid(); ++mfi)
	mfbig_phi[mfi].copy(mftmp[mfi]);
      mftmp.clear();
      mfbig_phi.FillBoundary();
      fgeom.FillPeriodicBoundary(mfbig_phi,true);

      for (MFIter mfi(trock_phi); mfi.isValid(); ++mfi)
	{
	  const int* lo    = mfi.validbox().loVect();
	  const int* hi    = mfi.validbox().hiVect();

	  const int* p_lo  = trock_phi[mfi].loVect();
	  const int* p_hi  = trock_phi[mfi].hiVect();
	  const Real* pdat = trock_phi[mfi].dataPtr();
	  
	  const int*  mfp_lo = mfbig_phi[mfi].loVect();
	  const int*  mfp_hi = mfbig_phi[mfi].hiVect();
	  const Real* mfpdat = mfbig_phi[mfi].dataPtr();
	  
	  FORT_INITPHI2 (mfpdat, ARLIM(mfp_lo), ARLIM(mfp_hi),
			 pdat,ARLIM(p_lo),ARLIM(p_hi),
			 lo,hi,&level,&max_level, &fratio);
	}
      mfbig_phi.clear();

      BoxArray tba2(trock_phi.boxArray());
      tba2.grow(3);
      MultiFab tmpgrow(tba2,1,0);
      
      for (MFIter mfi(trock_phi); mfi.isValid(); ++mfi)
	tmpgrow[mfi].copy(trock_phi[mfi]);
      
      trock_phi.clear();

      tba2 = rock_phi->boxArray();
      tba2.grow(3);
      MultiFab tmpgrow2(tba2,1,0);

      tmpgrow2.copy(tmpgrow);
      tmpgrow.clear();

      for (MFIter mfi(tmpgrow2); mfi.isValid(); ++mfi)
	(*rock_phi)[mfi].copy(tmpgrow2[mfi]);
    }
  else
    { 
        BoxLib::Abort("!(porosity_from_fine) no yet supported");
#if 0
      int porosity_type = 0;
      (*rock_phi).setVal(rock_array[0].porosity);

      if (porosity_type != 0)
	{
	  int porosity_nlayer = rock_array.size();
	  Array<Real> porosity_val(porosity_nlayer);
	  for (int i=0;i<porosity_nlayer;i++)
	    porosity_val[i]=rock_array[i].porosity;

	  for (MFIter mfi(*rock_phi); mfi.isValid(); ++mfi)
	    {
	      const int* p_lo  = (*rock_phi)[mfi].loVect();
	      const int* p_hi  = (*rock_phi)[mfi].hiVect();
	      const Real* pdat = (*rock_phi)[mfi].dataPtr();
	      
	      FORT_INITPHI (pdat,ARLIM(p_lo),ARLIM(p_hi),
			    domain_hi, dx, &porosity_type,
			    porosity_val.dataPtr(),&porosity_nlayer);
	    }
	}
#endif
    }
  rock_phi->FillBoundary();

  if (model != model_list["single-phase"] && 
      model != model_list["single-phase-solid"] &&
      model != model_list["steady-saturated"])
    {
      bool do_fine_average = true;
      // relative permeability
      FArrayBox tmpfab;
      Real dxf[BL_SPACEDIM];
      for (int i = 0; i<BL_SPACEDIM; i++)
	dxf[i] = dx[i]/twoexp;
      int n_kr_coef = kr_coef->nComp();
      for (MFIter mfi(*kr_coef); mfi.isValid(); ++mfi)
	{
	  if (do_fine_average)
	    {
	      // build data on finest grid
	      Box bx = (*kr_coef)[mfi].box();
	      bx.refine(twoexp);
	      tmpfab.resize(bx,n_kr_coef);
	      tmpfab.setVal(0.);

	      for (int i=0; i<rocks.size(); i++) {
                  rocks[i].set_constant_krval(tmpfab,dxf);
              }

	      // average onto coarse grid
	      const int* p_lo  = (*kr_coef)[mfi].loVect();
	      const int* p_hi  = (*kr_coef)[mfi].hiVect();
	      const Real* pdat = (*kr_coef)[mfi].dataPtr();
	  
	      const int*  mfp_lo = tmpfab.loVect();
	      const int*  mfp_hi = tmpfab.hiVect();
	      const Real* mfpdat = tmpfab.dataPtr();
	      
	      FORT_INITKR (mfpdat, ARLIM(mfp_lo), ARLIM(mfp_hi),
			   pdat,ARLIM(p_lo),ARLIM(p_hi),&n_kr_coef,
			   &level, &max_level, &fratio);
	    }
	  else
	    {
                for (int i=0; i<rocks.size(); i++) {
                    rocks[i].set_constant_krval((*kr_coef)[mfi],dx);
                }
	    }
	}
      // capillary pressure
      int n_cpl_coef = cpl_coef->nComp();
      for (MFIter mfi(*cpl_coef); mfi.isValid(); ++mfi)
	{
	  if (do_fine_average)
	    {
	      // build data on finest grid
	      Box bx = (*cpl_coef)[mfi].box();
	      bx.refine(twoexp);
	      tmpfab.resize(bx,n_cpl_coef);
	      tmpfab.setVal(0.);

	      for (int i=0; i<rocks.size(); i++) {
                  rocks[i].set_constant_cplval(tmpfab,dxf);
              }
	      
	      // average onto coarse grid
	      const int* p_lo  = (*cpl_coef)[mfi].loVect();
	      const int* p_hi  = (*cpl_coef)[mfi].hiVect();
	      const Real* pdat = (*cpl_coef)[mfi].dataPtr();
	      
	      const int*  mfp_lo = tmpfab.loVect();
	      const int*  mfp_hi = tmpfab.hiVect();
	      const Real* mfpdat = tmpfab.dataPtr();
	  
	      FORT_INITKR (mfpdat, ARLIM(mfp_lo), ARLIM(mfp_hi),
			   pdat,ARLIM(p_lo),ARLIM(p_hi),&n_cpl_coef,
			   &level,&max_level, &fratio);
	    }
	  else
	    {
                for (int i=0; i<rocks.size(); i++) {
                  rocks[i].set_constant_cplval((*cpl_coef)[mfi],dx);
                }
	    }
	}
    }
}

//
// Ensure state, and pressure are consistent.
//

void
PorousMedia::post_init (Real stop_time)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_init()");

  if (level > 0)
    //
    // Nothing to sync up at level > 0.
    //
    return;

  const int   finest_level = parent->finestLevel();
  Real        dt_init_local = 0.;
  Array<Real> dt_save(finest_level+1);
  Array<int>  nc_save(finest_level+1);

  //
  // Ensure state is consistent, i.e. velocity field is non-divergent,
  // Coarse levels are fine level averages, pressure is zero.
  // Call initial_mac_project in order to get a good initial dt.
  //
  post_init_state();
  //
  // Estimate the initial timestepping.
  //
  post_init_estDT(dt_init_local, nc_save, dt_save, stop_time);

  const Real strt_time       = state[State_Type].curTime();
  for (int k = 0; k <= finest_level; k++)
    getLevel(k).setTimeLevel(strt_time,dt_save[k],dt_save[k]);

  parent->setDtLevel(dt_save);
  parent->setNCycle(nc_save);

  //
  // Compute the initial estimate of conservation.
  //
  if (sum_interval > 0)
    sum_integrated_quantities();

  if (level == 0)
    {
      Observation::setAmrPtr(parent);
      Real prev_time = state[State_Type].prevTime();
      Real curr_time = state[State_Type].curTime();
      for (int i=0; i<observations.size(); ++i)
          observations[i].process(prev_time, curr_time, parent->levelSteps(0));
    }
}

//
// MULTILEVEL SYNC FUNCTIONS
//


//
// This function ensures that the state is initially consistent
// with respect to the divergence condition and fields are initially consistent
//

void
PorousMedia::post_init_state ()
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_init_state()");

  //
  // Richard initialization
  //
  int  finest_level = parent->finestLevel();
  for (int lev=0;lev<= finest_level;lev++)
    {
      PorousMedia& pm = getLevel(lev);
      for (int i = 0; i < num_state_type; i++)
	{
	  pm.state[i].allocOldData();
	  MultiFab& od = pm.get_old_data(i);
	  MultiFab& nd = pm.get_new_data(i);
	  MultiFab::Copy(od,nd,0,0,nd.nComp(),0);
	}
    }

  if (model == model_list["richard"] && do_richard_init_to_steady) {
      richard_init_to_steady();        
  }

  PorousMedia::initial_step = true;

  //
  // Average scalar and pressure data down from finer levels
  // so that conserved data is consistant between levels.
  //
  for (int k = finest_level-1; k>= 0; k--)
    {
      getLevel(k).avgDown();
    }
}

//
// Compute an initial MAC velocity in order to get a good first dt
//
void
PorousMedia::initial_mac_project (MultiFab* u_mac, MultiFab* RhoD, Real time)
{
  mac_project(u_mac,RhoD,time);
}

//
// Helper function for PorousMedia::SyncInterp().
//

static
void
set_bc_new (int*            bc_new,
            int             n,
            int             src_comp,
            const int*      clo,
            const int*      chi,
            const int*      cdomlo,
            const int*      cdomhi,
            const BoxArray& cgrids,
            int**           bc_orig_qty)
            
{
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
      bc_new[bc_index]             = INT_DIR;
      bc_new[bc_index+BL_SPACEDIM] = INT_DIR;
 
      if (clo[dir] < cdomlo[dir] || chi[dir] > cdomhi[dir])
        {
	  for (int crse = 0; crse < cgrids.size(); crse++)
            {
	      const int* c_lo = cgrids[crse].loVect();
	      const int* c_hi = cgrids[crse].hiVect();

	      if (clo[dir] < cdomlo[dir] && c_lo[dir] == cdomlo[dir])
		bc_new[bc_index] = bc_orig_qty[crse][bc_index];
	      if (chi[dir] > cdomhi[dir] && c_hi[dir] == cdomhi[dir])
		bc_new[bc_index+BL_SPACEDIM] = bc_orig_qty[crse][bc_index+BL_SPACEDIM]; 
            }
        }
    }
}

//
// Interpolate A cell centered Sync correction from a
// coarse level (c_lev) to a fine level (f_lev).
//
// This routine interpolates the num_comp components of CrseSync
// (starting at src_comp) and either increments or puts the result into
// the num_comp components of FineSync (starting at dest_comp)
// The components of bc_orig_qty corespond to the quantities of CrseSync.
//

void
PorousMedia::SyncInterp (MultiFab&      CrseSync,
			 int            c_lev,
			 MultiFab&      FineSync,
			 int            f_lev,
			 IntVect&       ratio,
			 int            src_comp,
			 int            dest_comp,
			 int            num_comp,
			 int            increment,
			 Real           dt_clev, 
			 int**          bc_orig_qty,
			 SyncInterpType which_interp,
			 int            state_comp)
{
  BL_ASSERT(which_interp >= 0 && which_interp <= 5);

  Interpolater* interpolater = 0;

  switch (which_interp)
    {
    case PC_T:           interpolater = &pc_interp;           break;
    case CellCons_T:     interpolater = &cell_cons_interp;    break;
    case CellConsLin_T:  interpolater = &lincc_interp;        break;
    case CellConsProt_T: interpolater = &protected_interp;    break;
    default:
      BoxLib::Abort("PorousMedia::SyncInterp(): how did this happen");
    }

  PorousMedia&   fine_level  = getLevel(f_lev);
  const BoxArray& fgrids     = fine_level.boxArray();
  const Geometry& fgeom      = parent->Geom(f_lev);
  const BoxArray& cgrids     = getLevel(c_lev).boxArray();
  const Geometry& cgeom      = parent->Geom(c_lev);
  const Real*     dx_crse    = cgeom.CellSize();
  Box             cdomain    = BoxLib::coarsen(fgeom.Domain(),ratio);
  const int*      cdomlo     = cdomain.loVect();
  const int*      cdomhi     = cdomain.hiVect();
  int*            bc_new     = new int[2*BL_SPACEDIM*(src_comp+num_comp)];

  BoxArray cdataBA(fgrids.size());

  for (int i = 0; i < fgrids.size(); i++)
    cdataBA.set(i,interpolater->CoarseBox(fgrids[i],ratio));
  //
  // Note: The boxes in cdataBA may NOT be disjoint !!!
  //
  MultiFab cdataMF(cdataBA,num_comp,0);

  cdataMF.setVal(0);

  cdataMF.copy(CrseSync, src_comp, 0, num_comp);
  //
  // Set physical boundary conditions in cdataMF.
  //
  for (MFIter mfi(cdataMF); mfi.isValid(); ++mfi)
    {
      int         i       = mfi.index();
      RealBox     gridloc = RealBox(fine_level.boxArray()[i],
				    fine_level.Geom().CellSize(),
				    fine_level.Geom().ProbLo());
      FArrayBox&  cdata   = cdataMF[mfi];
      const int*  clo     = cdata.loVect();
      const int*  chi     = cdata.hiVect();
      const Real* xlo     = gridloc.lo();

      for (int n = 0; n < num_comp; n++)
        {
	  set_bc_new(bc_new,n,src_comp,clo,chi,cdomlo,cdomhi,cgrids,bc_orig_qty);

	  FORT_FILCC(cdata.dataPtr(n), ARLIM(clo), ARLIM(chi),
		     cdomlo, cdomhi, dx_crse, xlo,
		     &(bc_new[2*BL_SPACEDIM*(n+src_comp)]));
        }
    }
  cgeom.FillPeriodicBoundary(cdataMF, 0, num_comp);
  //
  // Interpolate from cdataMF to fdata and update FineSync.
  // Note that FineSync and cdataMF will have the same distribution
  // since the length of their BoxArrays are equal.
  //
  FArrayBox    fdata;
  Array<BCRec> bc_interp(num_comp);

  MultiFab* fine_stateMF = 0;
  if (interpolater == &protected_interp)
    fine_stateMF = &(getLevel(f_lev).get_new_data(State_Type));

  for (MFIter mfi(cdataMF); mfi.isValid(); ++mfi)
    {
      int        i     = mfi.index();
      FArrayBox& cdata = cdataMF[mfi];
      const int* clo   = cdata.loVect();
      const int* chi   = cdata.hiVect();

      fdata.resize(fgrids[i], num_comp);
      //
      // Set the boundary condition array for interpolation.
      //
      for (int n = 0; n < num_comp; n++)
        {
	  set_bc_new(bc_new,n,src_comp,clo,chi,cdomlo,cdomhi,cgrids,bc_orig_qty);
        }

      for (int n = 0; n < num_comp; n++)
        {
	  for (int dir = 0; dir < BL_SPACEDIM; dir++)
            {
	      int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
	      bc_interp[n].setLo(dir,bc_new[bc_index]);
	      bc_interp[n].setHi(dir,bc_new[bc_index+BL_SPACEDIM]);
            }
        }

      interpolater->interp(cdata,0,fdata,0,num_comp,fgrids[i],ratio,
			   cgeom,fgeom,bc_interp,src_comp,State_Type);

      if (increment)
        {
	  fdata.mult(dt_clev);

	  if (interpolater == &protected_interp) {

	    cdata.mult(dt_clev);
	    FArrayBox& fine_state = (*fine_stateMF)[i];
	    interpolater->protect(cdata,0,fdata,0,fine_state,state_comp,
				  num_comp,fgrids[i],ratio,
				  cgeom,fgeom,bc_interp);
	    Real dt_clev_inv = 1./dt_clev;
	    cdata.mult(dt_clev_inv);

	  }
            
	  FineSync[i].plus(fdata,0,dest_comp,num_comp);
        }
      else
        {
	  FineSync[i].copy(fdata,0,dest_comp,num_comp);
        }
    }

  delete [] bc_new;
}

//
// Interpolate sync pressure correction to a finer level.
//

void
PorousMedia::SyncProjInterp (MultiFab& phi,
			     int       c_lev,
			     MultiFab& P_new,
			     MultiFab& P_old,
			     int       f_lev,
			     IntVect&  ratio,
			     bool      first_crse_step_after_initial_iters,
			     Real      cur_crse_pres_time,
			     Real      prev_crse_pres_time)
{
  const Geometry& fgeom   = parent->Geom(f_lev);
  const BoxArray& P_grids = P_new.boxArray();
  const Geometry& cgeom   = parent->Geom(c_lev);

  BoxArray crse_ba(P_grids.size());

  for (int i = 0; i < P_grids.size(); i++)
    {
      crse_ba.set(i,node_bilinear_interp.CoarseBox(P_grids[i],ratio));
    }

  Array<BCRec> bc(BL_SPACEDIM);
  MultiFab     crse_phi(crse_ba,1,0);

  crse_phi.setVal(1.e200);
  crse_phi.copy(phi,0,0,1);

  FArrayBox     fine_phi;
  PorousMedia& fine_lev            = getLevel(f_lev);
  const Real    cur_fine_pres_time  = fine_lev.state[Press_Type].curTime();
  const Real    prev_fine_pres_time = fine_lev.state[Press_Type].prevTime();

  if (state[Press_Type].descriptor()->timeType() == 
      StateDescriptor::Point && first_crse_step_after_initial_iters)
    {
      const Real time_since_zero  = cur_crse_pres_time - prev_crse_pres_time;
      const Real dt_to_prev_time  = prev_fine_pres_time - prev_crse_pres_time;
      const Real dt_to_cur_time   = cur_fine_pres_time - prev_crse_pres_time;
      const Real cur_mult_factor  = dt_to_cur_time / time_since_zero;
      const Real prev_mult_factor = dt_to_prev_time / dt_to_cur_time;

      for (MFIter mfi(crse_phi); mfi.isValid(); ++mfi)
        {
	  fine_phi.resize(P_grids[mfi.index()],1);
	  fine_phi.setVal(1.e200);
	  node_bilinear_interp.interp(crse_phi[mfi],0,fine_phi,0,1,
				      fine_phi.box(),ratio,cgeom,fgeom,bc,
				      0,Press_Type);
	  fine_phi.mult(cur_mult_factor);
	  P_new[mfi.index()].plus(fine_phi);
	  fine_phi.mult(prev_mult_factor);
	  P_old[mfi.index()].plus(fine_phi);
        }
    }
  else 
    {
      for (MFIter mfi(crse_phi); mfi.isValid(); ++mfi)
        {
	  fine_phi.resize(P_grids[mfi.index()],1);
	  fine_phi.setVal(1.e200);
	  node_bilinear_interp.interp(crse_phi[mfi],0,fine_phi,0,1,
				      fine_phi.box(),ratio,cgeom,fgeom,bc,
				      0,Press_Type);
	  P_new[mfi.index()].plus(fine_phi);
	  P_old[mfi.index()].plus(fine_phi);

        }
    }
}

//
// Averages a multifab of fine data down onto a multifab of coarse data.
//
// This should be an Amrlevel or Multifab function
//
void
PorousMedia::avgDown (MultiFab* s_crse,
		      int c_lev,
		      MultiFab* s_fine, 
		      int f_lev) 
{
    PorousMedia&   fine_lev = getLevel(f_lev);
    PorousMedia&   crse_lev = getLevel(c_lev);
    const BoxArray& fgrids  = fine_lev.grids;
    MultiFab&       fvolume = fine_lev.volume;
    const BoxArray& cgrids  = crse_lev.grids;
    MultiFab&       cvolume = crse_lev.volume;
    IntVect         ratio   = parent->refRatio(c_lev);

    int nc = (*s_crse).nComp();
    avgDown(cgrids,fgrids,*s_crse,*s_fine,cvolume,fvolume,c_lev,f_lev,0,nc,ratio);
}

void
PorousMedia::avgDown (const BoxArray& cgrids,
		      const BoxArray& fgrids,
		      MultiFab&       S_crse,
		      MultiFab&       S_fine,
		      MultiFab&       cvolume,
		      MultiFab&       fvolume,
		      int             c_level,
		      int             f_level,
		      int             scomp,
		      int             ncomp,
		      const IntVect&  fratio)
{
  BL_ASSERT(cgrids == S_crse.boxArray());
  BL_ASSERT(fgrids == S_fine.boxArray());
  BL_ASSERT(cvolume.boxArray() == cgrids);
  BL_ASSERT(fvolume.boxArray() == fgrids);
  BL_ASSERT(S_crse.nComp() == S_fine.nComp());
  BL_ASSERT(fvolume.nComp() == 1 && cvolume.nComp() == 1);

  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::avgDown()");
  //
  // Coarsen() the fine stuff on processors owning the fine data.
  //
  BoxArray crse_S_fine_BA(fgrids.size());

  for (int i = 0; i < fgrids.size(); ++i)
    {
      crse_S_fine_BA.set(i,BoxLib::coarsen(fgrids[i],fratio));
    }

  MultiFab crse_S_fine(crse_S_fine_BA,ncomp,0);
  MultiFab crse_fvolume(crse_S_fine_BA,1,0);

  crse_fvolume.copy(cvolume);

  for (MFIter mfi(S_fine); mfi.isValid(); ++mfi)
    {
      const int i = mfi.index();

      avgDown(S_fine[i],crse_S_fine[i],fvolume[i],crse_fvolume[i],
	      f_level,c_level,crse_S_fine_BA[i],scomp,ncomp,fratio);
    }

  S_crse.copy(crse_S_fine,0,scomp,ncomp);
}

//
// Average fine down to coarse in the ovlp intersection.
//

void
PorousMedia::avgDown (const FArrayBox& fine_fab,
		      const FArrayBox& crse_fab, 
		      const FArrayBox& fine_vol,
		      const FArrayBox& crse_vol,
		      int              f_level,
		      int              c_level,
		      const Box&       ovlp,
		      int              scomp,
		      int              ncomp,
		      const IntVect&   fratio)
{
  avgDown_doit(fine_fab,crse_fab,fine_vol,crse_vol,
	       f_level,c_level,ovlp,scomp,ncomp,fratio);
}



//
// Actually average the data down (this is static)
//

void
PorousMedia::avgDown_doit (const FArrayBox& fine_fab,
			   const FArrayBox& crse_fab, 
			   const FArrayBox& fine_vol,
			   const FArrayBox& crse_vol,
			   int              f_level,
			   int              c_level,
			   const Box&       ovlp,
			   int              scomp,
			   int              ncomp,
			   const IntVect&   fratio)
{
  //
  //  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
  //        because the crse fab is a temporary which was made starting at comp 0, it is
  //        not the actual state data.
  //
  const int*  ovlo   = ovlp.loVect();
  const int*  ovhi   = ovlp.hiVect();
  const int*  flo    = fine_fab.loVect();
  const int*  fhi    = fine_fab.hiVect();
  const Real* f_dat  = fine_fab.dataPtr(scomp);
  const int*  fvlo   = fine_vol.loVect();
  const int*  fvhi   = fine_vol.hiVect();
  const Real* fv_dat = fine_vol.dataPtr();
  const int*  clo    = crse_fab.loVect();
  const int*  chi    = crse_fab.hiVect();
  const Real* c_dat  = crse_fab.dataPtr();
  const int*  cvlo   = crse_vol.loVect();
  const int*  cvhi   = crse_vol.hiVect();
  const Real* cv_dat = crse_vol.dataPtr();

  FORT_AVGDOWN(c_dat,ARLIM(clo),ARLIM(chi),&ncomp,
	       f_dat,ARLIM(flo),ARLIM(fhi),
	       cv_dat,ARLIM(cvlo),ARLIM(cvhi),
	       fv_dat,ARLIM(fvlo),ARLIM(fvhi),
	       ovlo,ovhi,fratio.getVect());
}

static
void
SyncMacAcrossPeriodicEdges (MultiFab&       u_mac_crse_in_dir,
                            const MultiFab& crse_src,
                            const Geometry& cgeom,
                            int             dir,
                            int             nc)
{
  if (cgeom.isPeriodic(dir))
    {
      const Box cdmn = BoxLib::surroundingNodes(cgeom.Domain(),dir);

      const int N = 2;
      const int L = cdmn.length(dir) - 1;

      Box sides[N] = {cdmn,cdmn};

      sides[0].shift(dir, +L); // The hi end.
      sides[1].shift(dir, -L); // The lo end.

      const IntVect ZeroVector(D_DECL(0,0,0));

      IntVect shifts[N] = {ZeroVector,ZeroVector};

      shifts[0][dir] = -L; // How to shift hi -> lo
      shifts[1][dir] = +L; // How to shift lo -> hi

      for (int which = 0; which < N; ++which)
        {
	  Array<int> pmap;

	  BoxList bl(cdmn.ixType());

	  std::vector< std::pair<int,Box> > isects;

	  isects = crse_src.boxArray().intersections(sides[which]);

	  for (int i = 0; i < isects.size(); i++)
            {
	      const Box bx = crse_src.boxArray()[isects[i].first] & cdmn;

	      if (bx.ok())
                {
		  bl.push_back(bx);
		  pmap.push_back(crse_src.DistributionMap()[isects[i].first]);
                }
            }

	  if (!bl.isEmpty())
            {
	      pmap.push_back(ParallelDescriptor::MyProc()); // The sentinel.
	      MultiFab mf;
	      mf.define(BoxArray(bl), nc, 0, DistributionMapping(pmap), Fab_allocate);
	      mf.copy(crse_src);
	      mf.shift(shifts[which]);
	      u_mac_crse_in_dir.copy(mf);
            }
        }
    }
}

//
// Average edged values down a level
//
void
PorousMedia::SyncEAvgDown (PArray<MultiFab> u_mac_crse,
			   PArray<MultiFab> u_mac_fine,
			   int c_lev)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::SyncEAvgDown()");

  const Geometry& cgeom = parent->Geom(c_lev);
  IntVect         ratio = parent->refRatio(c_lev);
  int             nc    = u_mac_fine[0].nComp();
  
  for (int n = 0; n < u_mac_fine.size(); ++n)
    {
      //
      // crse_src & fine_src must have same parallel distribution.
      // We'll use the KnapSack distribution for the fine_src_ba.
      // Since fine_src_ba should contain more points, this'll lead
      // to a better distribution.
      //
      BoxArray fine_src_ba = u_mac_fine[n].boxArray();
      BoxArray crse_src_ba = BoxArray(fine_src_ba.size());

      for (int i=0; i<fine_src_ba.size();++i)
	{
	  crse_src_ba.set(i,Box(fine_src_ba[i]).coarsen(ratio));
	  fine_src_ba.set(i,Box(crse_src_ba[i]).refine(ratio));
	}

      std::vector<long> wgts(fine_src_ba.size());
    
      for (unsigned int i = 0; i < wgts.size(); i++)
	wgts[i] = fine_src_ba[i].numPts();
	
      DistributionMapping dm;
      dm.KnapSackProcessorMap(wgts,ParallelDescriptor::NProcs());

      MultiFab crse_src,  fine_src; 
    
      crse_src.define(crse_src_ba, nc, 0, dm, Fab_allocate);
      fine_src.define(fine_src_ba, nc, 0, dm, Fab_allocate);
    
      crse_src.setVal(1.e200);
      fine_src.setVal(1.e200);
	
      fine_src.copy(u_mac_fine[n]);
    
      for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
	{
	  const int  nComp = nc;
	  const Box& box   = crse_src[mfi].box();
	  const int* rat   = ratio.getVect();
	  FORT_EDGE_AVGDOWN(box.loVect(), box.hiVect(), &nComp, rat, &n,
			    crse_src[mfi].dataPtr(),
			    ARLIM(crse_src[mfi].loVect()),
			    ARLIM(crse_src[mfi].hiVect()),
			    fine_src[mfi].dataPtr(),
			    ARLIM(fine_src[mfi].loVect()),
			    ARLIM(fine_src[mfi].hiVect()));
	}

      fine_src.clear();
    
      u_mac_crse[n].copy(crse_src);

      SyncMacAcrossPeriodicEdges(u_mac_crse[n], crse_src, cgeom, n, nc);

    }
}

void
PorousMedia::SyncEAvgDown (MultiFab* u_mac_crse,
			   int c_lev,
			   MultiFab* u_mac_fine, 
			   int f_lev) 
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::SyncEAvgDown()");

  BL_ASSERT(f_lev > 0);

  const Geometry& cgeom      = parent->Geom(c_lev);
  const BoxArray& fgrids     = getLevel(f_lev).grids;
  IntVect    ratio           = IntVect::TheUnitVector();
  ratio                     *= parent->refRatio(c_lev);
  BoxArray f_bnd_ba = fgrids;
  BoxArray c_bnd_ba = BoxArray(f_bnd_ba.size());

  int nc = u_mac_fine[0].nComp();

  for (int i = 0; i < f_bnd_ba.size(); ++i)
    {
      c_bnd_ba.set(i,Box(f_bnd_ba[i]).coarsen(ratio));
      f_bnd_ba.set(i,Box(c_bnd_ba[i]).refine(ratio));
    }

  for (int n = 0; n < BL_SPACEDIM; ++n)
    {
      //
      // crse_src & fine_src must have same parallel distribution.
      // We'll use the KnapSack distribution for the fine_src_ba.
      // Since fine_src_ba should contain more points, this'll lead
      // to a better distribution.
      //
      BoxArray crse_src_ba(c_bnd_ba);
      BoxArray fine_src_ba(f_bnd_ba);

      crse_src_ba.surroundingNodes(n);
      fine_src_ba.surroundingNodes(n);

      std::vector<long> wgts(fine_src_ba.size());

      for (unsigned int i = 0; i < wgts.size(); i++)
	{
	  wgts[i] = fine_src_ba[i].numPts();
	}
      DistributionMapping dm;
      //
      // This call doesn't invoke the MinimizeCommCosts() stuff.
      // There's very little to gain with these types of coverings
      // of trying to use SFC or anything else.
      // This also guarantees that these DMs won't be put into the
      // cache, as it's not representative of that used for more
      // usual MultiFabs.
      //
      dm.KnapSackProcessorMap(wgts,ParallelDescriptor::NProcs());

      MultiFab crse_src,  fine_src; 

      crse_src.define(crse_src_ba, nc, 0, dm, Fab_allocate);
      fine_src.define(fine_src_ba, nc, 0, dm, Fab_allocate);
	    
      crse_src.setVal(1.e200);
      fine_src.setVal(1.e200);
	
      fine_src.copy(u_mac_fine[n]);
        
      for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
	{
	  const int  nComp = nc;
	  const Box& box   = crse_src[mfi].box();
	  const int* rat   = ratio.getVect();
	  FORT_EDGE_AVGDOWN(box.loVect(), box.hiVect(), &nComp, rat, &n,
			    crse_src[mfi].dataPtr(),
			    ARLIM(crse_src[mfi].loVect()),
			    ARLIM(crse_src[mfi].hiVect()),
			    fine_src[mfi].dataPtr(),
			    ARLIM(fine_src[mfi].loVect()),
			    ARLIM(fine_src[mfi].hiVect()));
	}
      fine_src.clear();

      u_mac_crse[n].copy(crse_src);

      SyncMacAcrossPeriodicEdges(u_mac_crse[n], crse_src, cgeom, n, nc);

    }
}

void
PorousMedia::SyncEAvgDown (MultiFab* u_mac_crse[],
			   int c_lev,
			   MultiFab* u_mac_fine[], 
			   int f_lev) 
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::SyncEAvgDown()");

  BL_ASSERT(f_lev > 0);

  const Geometry& cgeom      = parent->Geom(c_lev);
  const BoxArray& fgrids     = getLevel(f_lev).grids;
  IntVect    ratio           = IntVect::TheUnitVector();
  ratio                     *= parent->refRatio(c_lev);
  BoxArray f_bnd_ba = fgrids;
  BoxArray c_bnd_ba = BoxArray(f_bnd_ba.size());

  int nc = (*u_mac_fine[0]).nComp();

  for (int i = 0; i < f_bnd_ba.size(); ++i)
    {
      c_bnd_ba.set(i,Box(f_bnd_ba[i]).coarsen(ratio));
      f_bnd_ba.set(i,Box(c_bnd_ba[i]).refine(ratio));
    }

  for (int n = 0; n < BL_SPACEDIM; ++n)
    {
      //
      // crse_src & fine_src must have same parallel distribution.
      // We'll use the KnapSack distribution for the fine_src_ba.
      // Since fine_src_ba should contain more points, this'll lead
      // to a better distribution.
      //
      BoxArray crse_src_ba(c_bnd_ba);
      BoxArray fine_src_ba(f_bnd_ba);

      crse_src_ba.surroundingNodes(n);
      fine_src_ba.surroundingNodes(n);

      std::vector<long> wgts(fine_src_ba.size());

      for (unsigned int i = 0; i < wgts.size(); i++)
	{
	  wgts[i] = fine_src_ba[i].numPts();
	}
      DistributionMapping dm;
      //
      // This call doesn't invoke the MinimizeCommCosts() stuff.
      // There's very little to gain with these types of coverings
      // of trying to use SFC or anything else.
      // This also guarantees that these DMs won't be put into the
      // cache, as it's not representative of that used for more
      // usual MultiFabs.
      //
      dm.KnapSackProcessorMap(wgts,ParallelDescriptor::NProcs());

      MultiFab crse_src,  fine_src; 

      crse_src.define(crse_src_ba, nc, 0, dm, Fab_allocate);
      fine_src.define(fine_src_ba, nc, 0, dm, Fab_allocate);
	    
      crse_src.setVal(1.e200);
      fine_src.setVal(1.e200);
	
      fine_src.copy(*u_mac_fine[n]);
        
      for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
	{
	  const int  nComp = nc;
	  const Box& box   = crse_src[mfi].box();
	  const int* rat   = ratio.getVect();
	  FORT_EDGE_AVGDOWN(box.loVect(), box.hiVect(), &nComp, rat, &n,
			    crse_src[mfi].dataPtr(),
			    ARLIM(crse_src[mfi].loVect()),
			    ARLIM(crse_src[mfi].hiVect()),
			    fine_src[mfi].dataPtr(),
			    ARLIM(fine_src[mfi].loVect()),
			    ARLIM(fine_src[mfi].hiVect()));
	}

      fine_src.clear();

      u_mac_crse[n]->copy(crse_src);

      SyncMacAcrossPeriodicEdges(*u_mac_crse[n], crse_src, cgeom, n, nc);

    }
}

//
// The Mac Sync correction function
//
void
PorousMedia::mac_sync ()
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::mac_sync()");

  const int  numscal   = ncomps; 
  const Real prev_time = state[State_Type].prevTime();
  const Real curr_time = state[State_Type].curTime();
  const Real dt        = parent->dtLevel(level);
  const BCRec& p_bc    = desc_lst[Press_Type].getBC(0);
        
  //
  // Compute the u_mac for the correction.
  //
  MultiFab* p_corr = new MultiFab(grids,1,1);
  for (int i=0; i < BL_SPACEDIM; i++)
    u_corr[i].setVal(0.);
  create_lambda(curr_time); 
  mac_projector->mac_sync_solve(level,p_bc,lambda,p_corr,u_corr,fine_ratio);
  
  //
  // Assign rock_phi to alpha
  //
  MultiFab* alpha = new MultiFab(grids, 1, 1);
  MultiFab::Copy(*alpha,*rock_phi,0,0,1,alpha->nGrow());
  //
  // Update coarse grid state by adding correction from mac_sync solve
  // the correction is the advective tendency of the new velocities.
  //
  mac_projector->mac_sync_compute(level,u_macG_curr,u_corr,
				  Ssync,lambda,rock_phi,kappa,
				  lambda_cc,dlambda_cc,kr_coef,
				  kpedge,p_corr,
				  level > 0 ? &getAdvFluxReg(level) : 0,
				  advectionType, prev_time, dt,
				  ncomps,be_cn_theta);
  //
  // The following used to be done in mac_sync_compute.  Ssync is
  //   the source for a rate of change to rock_phi*S over the time step, so
  //   Ssync*dt is the source to the actual sync amount.
  //
  MultiFab& S_new = get_new_data(State_Type);

  if (verbose > 3)
    {
      Real tmp = (*Ssync).norm2(0);
      if (ParallelDescriptor::IOProcessor())
	std::cout << "SSYNC NORM  AFTER = " << tmp << '\n';
      Ssync->mult(-dt,Ssync->nGrow());
    
      MultiFab::Copy(S_new,*Ssync,0,ncomps+ntracers+1,1,1);
    }
    
  //
  // Diffusion solve for Ssync
  //    
  bool any_diffusive = false;
  for (int kk  = 0; kk < ncomps; kk++)
    if (is_diffusive[kk])
      any_diffusive = true;
    
  if (any_diffusive)
    {
      MultiFab tmp(grids,1,1);
      MultiFab** fluxSC  = 0;
      diffusion->allocFluxBoxesLevel(fluxSC,0,1);
      
      tmp.setVal(0.);
      for (int i=0; i < BL_SPACEDIM; ++i){
	(*fluxSC[i]).setVal(0.);
      }

      //
      // Set up rho function for diffusive solve
      //
      MultiFab* rho = new MultiFab(grids,1,1);
      MultiFab::Copy(*rho,S_new,0,0,1,1);
      for (int kk = 1; kk<ncomps; kk++)
	{
	  if (solid.compare(pNames[pType[kk]]) != 0) 
	    MultiFab::Add(*rho,S_new,kk,0,1,1);
	}
      diffusion->set_rho(rho);
      delete rho;
      
      for (int kk = 0; kk<ncomps; kk++)
	{
	  if (is_diffusive[kk])
	    {
	      MultiFab** cmp_diffn=0;
	  
	      if (variable_scal_diff)
		{
		  Real diffTime = state[State_Type].curTime();
		  diffusion->allocFluxBoxesLevel(cmp_diffn, 0, 1);
		  getDiffusivity(cmp_diffn, diffTime,kk,0,1);
		}
	      diffusion->diffuse_Ssync(Ssync,kk,dt,be_cn_theta,
				       fluxSC,0,cmp_diffn,alpha);
	      if (variable_scal_diff)
		diffusion->removeFluxBoxesLevel(cmp_diffn);

	      if (level > 0)
		{
		  for (int d = 0; d < BL_SPACEDIM; d++)
		    {
		      Real mult = dt;
		      MultiFab& fluxSCd = *fluxSC[d];
		      for (MFIter fmfi(fluxSCd); fmfi.isValid(); ++fmfi)
			getViscFluxReg().FineAdd(fluxSCd[fmfi],d,
						 fmfi.index(),
						 0,kk,1,mult);
		    }
		}
	    }
	}
      diffusion->removeFluxBoxesLevel(fluxSC);
    }
    
  // 
  // Capillary-solve.  Since capillary function is nonlinear, we cannot
  // do a simple capillary-diffuse solve for Ssync.  A full nonlinear
  // parabolic solve is needed to determine the new solution at the end of 
  // coarse timestep.  
  //
  if  (have_capillary)
    {
      const int nGrow = 0;
      const int nComp = 1;
      MultiFab** fluxSC    = 0;
      MultiFab** fluxSCp1  = 0;
      diffusion->allocFluxBoxesLevel(fluxSC,  nGrow,nComp);
      diffusion->allocFluxBoxesLevel(fluxSCp1,nGrow,nComp);
      
      int nc = 0; 
      int nd = 1;
      MultiFab*  delta_rhs = 0;
      MultiFab** cmp_pcn   = 0;
      MultiFab** cmp_pcnp1 = 0;
      MultiFab** cmp_pcnp1_dp = 0;
      MultiFab*  S_nwt = 0;
      MultiFab&  S_new = get_new_data(State_Type);

      MultiFab* sat_res_mf = new MultiFab(grids,1,1);
      sat_res_mf->setVal(1.);
      for (MFIter mfi(*sat_res_mf); mfi.isValid();++mfi)
	{
	  const Box& box = (*sat_res_mf)[mfi].box();
	  (*sat_res_mf)[mfi].minus((*cpl_coef)[mfi],box,3,0,1);
	}
      sat_res_mf->mult(density[nc]);
      diffusion->set_rho(sat_res_mf); 

      MultiFab S_tmp(grids,ncomps,1);
      MultiFab::Copy(S_tmp,S_new,0,0,ncomps,1);

      S_nwt = new MultiFab(grids,1,1);
      MultiFab::Copy(*S_nwt,S_new,nc,0,nComp,1);
      
      delta_rhs = new MultiFab(grids,1,1);
      MultiFab::Copy(*delta_rhs,*Ssync,nc,0,nComp,1);

      //
      // Newton iteration
      //

      // initialization
      Real pcTime = state[State_Type].prevTime();
      diffusion->allocFluxBoxesLevel(cmp_pcn,0,1);
      calcCapillary(pcTime);
      calcDiffusivity_CPL(cmp_pcn,lambda_cc); 
 
      pcTime = state[State_Type].curTime();
      FillStateBndry (pcTime,State_Type,0,ncomps);
      diffusion->allocFluxBoxesLevel(cmp_pcnp1,0,1);
      diffusion->allocFluxBoxesLevel(cmp_pcnp1_dp,0,1);
      calcCapillary(pcTime);
      calcLambda(pcTime);
      calcDiffusivity_CPL(cmp_pcnp1,lambdap1_cc);
      calcDiffusivity_CPL_dp(cmp_pcnp1_dp,lambdap1_cc,pcTime,1);
      
      int  max_itr_nwt = 20;
#if (BL_SPACEDIM == 3)
      Real max_err_nwt = 1e-8;
#else
      Real max_err_nwt = 1e-8;
#endif
 
      int  itr_nwt = 0;
      Real err_nwt = 1e10;
      diffusion->diffuse_init_CPL(dt,nc,be_cn_theta,
				  fluxSC,0,delta_rhs,
				  alpha,cmp_pcn,pcn_cc,S_nwt);
      while ((itr_nwt < max_itr_nwt) && (err_nwt > max_err_nwt)) 
	{
	  diffusion->diffuse_iter_CPL(dt,nc,ncomps,be_cn_theta,
				      0,alpha,cmp_pcnp1,cmp_pcnp1_dp,
				      pcnp1_cc,S_nwt,&err_nwt);

	  if (verbose > 3 && ParallelDescriptor::IOProcessor())
	    std::cout << "Newton iteration " << itr_nwt 
		      << " : Error = "       << err_nwt << "\n"; 

	  scalar_adjust_constraint(0,ncomps-1);
	  FillStateBndry (pcTime,State_Type,0,ncomps);
	  calcCapillary(pcTime);
	  calcLambda(pcTime);
	  calcDiffusivity_CPL(cmp_pcnp1,lambdap1_cc);
	  calcDiffusivity_CPL_dp(cmp_pcnp1_dp,lambdap1_cc,pcTime,1);
	  itr_nwt += 1;	  

	  if (verbose > 3) 
	    check_minmax();
	}

      diffusion->compute_flux(nc,dt,be_cn_theta,fluxSCp1,pcnp1_cc,cmp_pcnp1);

      if (verbose > 3 && ParallelDescriptor::IOProcessor())
	{
	  if (itr_nwt < max_itr_nwt)
	    std::cout << "Newton converged at iteration " << itr_nwt
		      << " with error " << err_nwt << '\n';
	  else
	    std::cout << "Newton failed to converged: termination error is "
		      <<  err_nwt << '\n'; 
	}
      
      if (level > 0)
	{
	  for (int d = 0; d < BL_SPACEDIM; d++)
	    {
	      Real mult = -dt;
	      MultiFab& fluxSCd = *fluxSCp1[d];
	      for (MFIter fmfi(fluxSCd); fmfi.isValid(); ++fmfi)
		getViscFluxReg().FineAdd(fluxSCd[fmfi],d,
					 fmfi.index(),
					 0,nc,1,mult);
	  
	      fluxSCd.mult(-density[nd]/density[nc]);
	      for (MFIter fmfi(fluxSCd); fmfi.isValid(); ++fmfi)
		getViscFluxReg().FineAdd(fluxSCd[fmfi],d,
					 fmfi.index(),
					 0,nd,1,mult);
	    }
	}
      
      // Determine the corrector after capillary-solve
      for (MFIter mfi(*S_nwt); mfi.isValid();++mfi)
	{
	  const Box& box = mfi.validbox();
	  (*Ssync)[mfi].copy(S_new[mfi],box,0,box,0,ncomps);
	  (*Ssync)[mfi].minus(S_tmp[mfi],box,0,0,ncomps);
	}
	
      
  
      delete delta_rhs;
      delete S_nwt;
      delete sat_res_mf;
      
      diffusion->removeFluxBoxesLevel(fluxSC);
      diffusion->removeFluxBoxesLevel(fluxSCp1);
      diffusion->removeFluxBoxesLevel(cmp_pcn);
      diffusion->removeFluxBoxesLevel(cmp_pcnp1);
      diffusion->removeFluxBoxesLevel(cmp_pcnp1_dp);
    }
    
  delete p_corr;
  delete alpha;

  //
  // Add the sync correction to the state.
  //
  if (have_capillary == 0 && !any_diffusive)
    {
      for (MFIter mfi(*Ssync); mfi.isValid(); ++mfi) 
	{
	  for (int nc = 0; nc < ncomps; nc++)
	    (*Ssync)[mfi].divide((*rock_phi)[mfi],0,nc,1);
	}
    }
    
  if (have_capillary == 0)
    {
      for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
	S_new[mfi].plus((*Ssync)[mfi],mfi.validbox(),
			0,0,numscal);
    }
    
  if (idx_dominant > -1)
    scalar_adjust_constraint(0,ncomps-1);
      
  //
  // Get boundary conditions.
  //
  Array<int*>         sync_bc(grids.size());
  Array< Array<int> > sync_bc_array(grids.size());
      
  for (int i = 0; i < grids.size(); i++)
    {
      sync_bc_array[i] = getBCArray(State_Type,i,0,numscal);
      sync_bc[i]       = sync_bc_array[i].dataPtr();
    }

  //
  // Interpolate the sync correction to the finer levels.
  //
  IntVect    ratio = IntVect::TheUnitVector();
  const Real mult  = 1.0;
  for (int lev = level+1; lev <= parent->finestLevel(); lev++)
    {
      ratio                     *= parent->refRatio(lev-1);
      PorousMedia&     fine_lev  = getLevel(lev);
      const BoxArray& fine_grids = fine_lev.boxArray();
      MultiFab sync_incr(fine_grids,numscal,0);
      sync_incr.setVal(0.0);
      
      SyncInterp(*Ssync,level,sync_incr,lev,ratio,0,0,
		 numscal,1,mult,sync_bc.dataPtr());
      
      MultiFab& S_new = fine_lev.get_new_data(State_Type);
      for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
	S_new[mfi].plus(sync_incr[mfi],fine_grids[mfi.index()],
			0,0,numscal);
    }
}

//
// The Mac Sync correction function
//
#ifdef MG_USE_FBOXLIB
void
PorousMedia::richard_sync ()
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::richard_sync()");

  BL_ASSERT(have_capillary);

  const Real dt = parent->dtLevel(level);
       
  //
  //   Ssync is the source for a rate of change to rock_phi*S over the time step, so
  //   Ssync*dt*density[0] is the source to the actual sync amount.
  //

  if (verbose > 3)
    {
      Real tmp = (*Ssync).norm2(0);
      if (ParallelDescriptor::IOProcessor())
	std::cout << "SSYNC NORM  AFTER = " << tmp << '\n';
    }

  // 
  // Capillary-solve.  Since capillary function is nonlinear, we cannot
  // do a simple capillary-diffuse solve for Ssync.  A full nonlinear
  // parabolic solve is needed to determine the new solution at the end of 
  // coarse timestep.  
  //

  // Build single component edge-centered array of MultiFabs for fluxes
  MultiFab** fluxSC;
  const int nGrow = 0;
  const int nComp = 1;
  diffusion->allocFluxBoxesLevel(fluxSC,nGrow,nComp);
  
  int nc = 0; 
  MultiFab** cmp_pcp1    = 0;
  MultiFab** cmp_pcp1_dp = 0;
  MultiFab sat_res_mf(grids,1,1);
  sat_res_mf.setVal(1.);
  for (MFIter mfi(sat_res_mf); mfi.isValid();++mfi)
    {
      const Box& box = sat_res_mf[mfi].box();
      sat_res_mf[mfi].minus((*cpl_coef)[mfi],box,3,0,1);
    }
  diffusion->set_rho(&sat_res_mf); 

  bool sync_n = true;

  MultiFab& S_new  = get_new_data(State_Type);
  MultiFab& S_old  = get_old_data(State_Type);
  MultiFab& P_new  = get_new_data(Press_Type);
  MultiFab* alpha  = new MultiFab(grids,1,1);
  MultiFab* dalpha = 0;
  MultiFab Tmp(grids,1,1);
  
  if (sync_n)
    MultiFab::Copy(Tmp,S_new,0,0,1,1);
  else
    MultiFab::Copy(Tmp,P_new,0,0,1,1);
  MultiFab::Copy(*alpha,*rock_phi,0,0,1,alpha->nGrow());
  
  if (!do_richard_sat_solve) dalpha = new MultiFab(grids,1,1);
      
  // Compute first res_fix = -\phi * n^k + dt*\nabla v_inflow.  
  // Its value does not change.
  MultiFab res_fix(grids,1,0);
  MultiFab::Copy(res_fix,S_old,nc,0,1,0);
  for (MFIter mfi(res_fix); mfi.isValid(); ++mfi)
    res_fix[mfi].mult((*alpha)[mfi],mfi.validbox(),0,0,1);
  res_fix.mult(-1.0);
  Ssync->mult(-dt*density[0]);
  MultiFab::Add(res_fix,*Ssync,nc,0,1,0);
  calc_richard_velbc(res_fix,u_mac_curr,dt*density[0]);
  // Newton method.
  // initialization
  int do_upwind = 1;
  int  max_itr_nwt = 20;
  Real max_err_nwt = 1e-8;
  int  itr_nwt = 0;
  Real err_nwt = 1e10;
  Real pcTime = state[State_Type].curTime();
  FillStateBndry(pcTime,State_Type,0,ncomps);
  diffusion->allocFluxBoxesLevel(cmp_pcp1,0,1);
  diffusion->allocFluxBoxesLevel(cmp_pcp1_dp,0,3);
  
  calcLambda(pcTime);
  calc_richard_coef(cmp_pcp1,lambdap1_cc,u_mac_curr,0,do_upwind,pcTime);
  calc_richard_jac (cmp_pcp1_dp,dalpha,lambdap1_cc,u_mac_curr,pcTime,dt,0,do_upwind,do_richard_sat_solve);

  
  //if (!do_richard_sat_solve) calc_richard_alpha(dalpha,pcTime);

  Diffusion::NewtonStepInfo linear_status;
  linear_status.status = "";
  linear_status.success = true;
  linear_status.reason = "";
  linear_status.ls_iterations = -1;
  linear_status.max_ls_iterations = 10;
  linear_status.min_ls_factor = 1.e-8;
  linear_status.ls_factor = -1;
  linear_status.ls_acceptance_factor = 1.;
  linear_status.ls_reduction_factor = 0.1;
  linear_status.residual_norm_pre_ls = -1; 
  linear_status.residual_norm_post_ls = -1;
  linear_status.initial_residual_norm = -1;
  linear_status.initial_solution_norm = -1;
  linear_status.monitor_linear_solve = false;
  linear_status.monitor_line_search = false;

  if (do_richard_sat_solve)
    {
      calcCapillary(pcTime);
      MultiFab::Copy(P_new,*pcnp1_cc,0,0,1,1);
      P_new.mult(-1.0,1);
      while ((itr_nwt < max_itr_nwt) && (err_nwt > max_err_nwt) && (linear_status.status!="Finished")) 
	{
            diffusion->richard_iter(dt,nc,gravity,density,res_fix,
                                    alpha,cmp_pcp1,cmp_pcp1_dp,
                                    u_mac_curr,do_upwind,linear_status);    
            
            err_nwt = linear_status.residual_norm_post_ls;

            if (linear_status.success) {
                if (verbose > 3 && ParallelDescriptor::IOProcessor())
                    std::cout << "Newton iteration " << itr_nwt 
                              << " : Error = "       << err_nwt << "\n"; 
                if (model != model_list["richard"])
                    scalar_adjust_constraint(0,ncomps-1);
                FillStateBndry(pcTime,State_Type,0,ncomps);
                calcCapillary(pcTime);
                calcLambda(pcTime);
                MultiFab::Copy(P_new,*pcnp1_cc,0,0,1,1);
                P_new.mult(-1.0,1);
                compute_vel_phase(u_mac_curr,0,pcTime);
                calc_richard_coef(cmp_pcp1,lambdap1_cc,u_mac_curr,0,do_upwind,pcTime);
                calc_richard_jac (cmp_pcp1_dp,0,lambdap1_cc,u_mac_curr,pcTime,dt,0,do_upwind,do_richard_sat_solve);
                itr_nwt += 1;
                if (verbose > 3)  check_minmax();
            }
	}
    }
  else
    {
        //MultiFab dalpha(grids,1,1);
        //calc_richard_alpha(&dalpha,pcTime);

      while ((itr_nwt < max_itr_nwt) && (err_nwt > max_err_nwt) && (linear_status.status!="Finished")) 
	{
            diffusion->richard_iter_p(dt,nc,gravity,density,res_fix,
                                      alpha,dalpha,cmp_pcp1,cmp_pcp1_dp,
                                      u_mac_curr,do_upwind,linear_status);

	    err_nwt = linear_status.residual_norm_post_ls;

            if (linear_status.success) {
                linear_status.success = true;
                if (verbose > 3 && ParallelDescriptor::IOProcessor())
                    std::cout << "Newton iteration " << itr_nwt 
                              << " : Error = "       << err_nwt << "\n"; 
                calcInvPressure (S_new,P_new); 
                if (model != model_list["richard"])
                    scalar_adjust_constraint(0,ncomps-1);
                calcLambda(pcTime);
                compute_vel_phase(u_mac_curr,0,pcTime);
                calc_richard_coef(cmp_pcp1,lambdap1_cc,u_mac_curr,0,do_upwind,pcTime);
                calc_richard_jac (cmp_pcp1_dp,dalpha,lambdap1_cc,u_mac_curr,pcTime,dt,0,do_upwind,do_richard_sat_solve);
                //calc_richard_alpha(&dalpha,pcTime);
                itr_nwt += 1;
                if (verbose > 3)  check_minmax();
            }
        }
    }

  RichardNLSdata::Reason retVal = RichardNLSdata::RICHARD_SUCCESS;
  if (!linear_status.success) {
      retVal = RichardNLSdata::RICHARD_LINEAR_FAIL;
  }

  if (err_nwt > max_err_nwt) {
      retVal = RichardNLSdata::RICHARD_NONLINEAR_FAIL;
      std::cout << "     **************** Newton failed in richard_sync: too many iterations\n"; 
  }

  if (retVal == RichardNLSdata::RICHARD_SUCCESS)
  {  
      diffusion->richard_flux(nc,-1.0,gravity,density,fluxSC,pcnp1_cc,cmp_pcp1);

      if (verbose > 3 && ParallelDescriptor::IOProcessor())
      {
          std::cout << "Newton converged at iteration " << itr_nwt
                    << " with error " << err_nwt << '\n';
      }
  
      if (level > 0) 
      {
          for (int d = 0; d < BL_SPACEDIM; d++) 
          {
              for (MFIter fmfi(*fluxSC[d]); fmfi.isValid(); ++fmfi)
                  getViscFluxReg().FineAdd((*fluxSC[d])[fmfi],d,fmfi.index(),0,nc,nComp,-dt);
          }
      }
  
      // Determine the corrector after capillary-solve
      for (MFIter mfi(*Ssync); mfi.isValid();++mfi)
      {
          const Box& box = mfi.validbox();
          if (sync_n)
          {
              (*Ssync)[mfi].copy(S_new[mfi],box,0,box,0,ncomps);
              (*Ssync)[mfi].minus(Tmp[mfi],box,0,0,ncomps);
          }
          else
          {
              (*Ssync)[mfi].copy(P_new[mfi],box,0,box,0,1);
              (*Ssync)[mfi].minus(Tmp[mfi],box,0,0,1);
          }
      }

      MultiFab::Copy(S_new,*Ssync,0,ncomps+ntracers,1,0);
  
      //
      // Get boundary conditions.
      //
      Array<int*>         sync_bc(grids.size());
      Array< Array<int> > sync_bc_array(grids.size());
      
      for (int i = 0; i < grids.size(); i++)
      {
          sync_bc_array[i] = getBCArray(Press_Type,i,0,1);
          sync_bc[i]       = sync_bc_array[i].dataPtr();
      }

      //
      // Interpolate the sync correction to the finer levels.
      //
      IntVect    ratio = IntVect::TheUnitVector();
      const Real mult  = 1.0;
      for (int lev = level+1; lev <= parent->finestLevel(); lev++)
      {
          ratio                     *= parent->refRatio(lev-1);
          PorousMedia&     fine_lev  = getLevel(lev);
          const BoxArray& fine_grids = fine_lev.boxArray();
          MultiFab sync_incr(fine_grids,1,0);
          sync_incr.setVal(0.0);
      
          SyncInterp(*Ssync,level,sync_incr,lev,ratio,0,0,
                     1,1,mult,sync_bc.dataPtr());

          MultiFab& S_new = fine_lev.get_new_data(Press_Type);
          if (sync_n)
          {
              for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
                  S_new[mfi].plus(sync_incr[mfi],fine_grids[mfi.index()],
                                  0,0,1);
          }
          else
          {	    
              MultiFab& P_new = fine_lev.get_new_data(Press_Type);
              for (MFIter mfi(P_new); mfi.isValid(); ++mfi)
                  P_new[mfi].plus(sync_incr[mfi],fine_grids[mfi.index()],
                                  0,0,1);
              MultiFab P_tmp(fine_grids,1,0);
              MultiFab::Copy(P_tmp,P_new,0,0,1,0);
              P_tmp.mult(-1.0);
              fine_lev.calcInvCapillary(sync_incr,P_tmp);
              MultiFab::Copy(S_new,sync_incr,0,0,1,0);
          }
      }
  }

  delete alpha;
  if (dalpha) delete dalpha;
  diffusion->removeFluxBoxesLevel(cmp_pcp1);
  diffusion->removeFluxBoxesLevel(cmp_pcp1_dp);
  diffusion->removeFluxBoxesLevel(fluxSC);
}
#endif

//
// The reflux function
//
void
PorousMedia::reflux ()
{
  if (level == parent->finestLevel())
    return;

  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::reflux()");

  BL_ASSERT(do_reflux);
  //
  // First do refluxing step.
  //
  FluxRegister& fr_adv  = getAdvFluxReg(level+1);
  FluxRegister& fr_visc = getViscFluxReg(level+1);
  Real          dt_crse = parent->dtLevel(level);
  Real          scale   = 1.0/dt_crse;

  fr_visc.Reflux(*Ssync,volume,scale,0,0,NUM_SCALARS,geom);
  fr_adv.Reflux (*Ssync,volume,scale,0,0,NUM_SCALARS,geom);
  //
  // This is necessary in order to zero out the contribution to any
  // coarse grid cells which underlie fine grid cells.
  //
  BoxArray baf = getLevel(level+1).boxArray();

  baf.coarsen(fine_ratio);

  for (MFIter mfi(*Ssync); mfi.isValid(); ++mfi)
  {
      BL_ASSERT(grids[mfi.index()] == mfi.validbox());

      std::vector< std::pair<int,Box> > isects = baf.intersections(mfi.validbox());

      for (int i = 0, N = isects.size(); i < N; i++)
      {
          (*Ssync)[mfi.index()].setVal(0,isects[i].second,0,NUM_SCALARS);
      }
  }
}

//
// Average fine information from the complete set of state types to coarse.
//

void
PorousMedia::avgDown ()
{
  if (level == parent->finestLevel())
    return;

  PorousMedia&   fine_lev = getLevel(level+1);
  const BoxArray& fgrids  = fine_lev.grids;
  MultiFab&       fvolume = fine_lev.volume;
  //
  // Average down the state at the new time.
  //
  MultiFab& S_crse = get_new_data(State_Type);
  MultiFab& S_fine = fine_lev.get_new_data(State_Type);
  avgDown(grids,fgrids,S_crse,S_fine,volume,fvolume,level,level+1,0,S_crse.nComp(),fine_ratio);

  //
  // Average down the pressure at the new time.
  //
  MultiFab& P_crse = get_new_data(Press_Type);
  MultiFab& P_fine = fine_lev.get_new_data(Press_Type);
  avgDown(grids,fgrids,P_crse,P_fine,volume,fvolume,level,level+1,0,1,fine_ratio);

  if (do_reflux && u_macG_curr != 0)
    SyncEAvgDown(u_macG_curr,level,fine_lev.u_macG_curr,level+1);

  //
  // Average down the cell-centered velocity at the new time.
  //
#ifdef AMANZI
  if (do_chem>0)
    {
      MultiFab& FC_crse = get_new_data(FuncCount_Type);
      MultiFab& FC_fine = fine_lev.get_new_data(FuncCount_Type);
      avgDown(grids,fgrids,FC_crse,FC_fine,volume,fvolume,
	      level,level+1,0,1,fine_ratio);
    }
#endif
}

//
// ACCESS FUNCTIONS FOLLOW
//

//
// Virtual access function for getting the advective flux out of the
// advection routines for diagnostics and refluxing.
//

void
PorousMedia::pullFluxes (int        i,
                         int        start_ind,
                         int        ncomp,
                         FArrayBox& xflux,
                         FArrayBox& yflux,
#if (BL_SPACEDIM == 3)
                         FArrayBox& zflux,
#endif
                         Real       dt)
{
  //
  // Add fluxes into the refluxing counters.
  //
  if (do_reflux)
    {
      if (level < parent->finestLevel())
        {
	  FluxRegister& fr = getAdvFluxReg(level+1);
	  fr.CrseInit(xflux,xflux.box(),0,0,start_ind,ncomp,-dt);
	  fr.CrseInit(yflux,yflux.box(),1,0,start_ind,ncomp,-dt);
#if (BL_SPACEDIM == 3)                              
	  fr.CrseInit(zflux,zflux.box(),2,0,start_ind,ncomp,-dt);
#endif
        }
      if (level > 0)
        {
	  advflux_reg->FineAdd(xflux,0,i,0,start_ind,ncomp,dt);
	  advflux_reg->FineAdd(yflux,1,i,0,start_ind,ncomp,dt);
#if (BL_SPACEDIM == 3)                                
	  advflux_reg->FineAdd(zflux,2,i,0,start_ind,ncomp,dt);
#endif
        }
    }
}

//
// Virtual access function for getting the forcing terms for the 
// pressure and scalars.  
//
void
PorousMedia::getForce (FArrayBox& force,
		       int        gridno,
		       int        ngrow,
		       int        scomp,
		       int        ncomp,
		       const Real time,
		       int        do_rho_scale)
{      

  force.resize(BoxLib::grow(grids[gridno],ngrow),ncomp);

  force.setVal(0);
  if (do_source_term)
    { 
      const Real* dx       = geom.CellSize();

      BoxLib::Abort("FIXME");
#if 0
      for (int i = 0; i< source_array.size(); i++)
	if (!source_array[i].var_type.compare("comp"))
	    source_array[i].setVal(force, region_array, dx); 
      
      if (do_rho_scale)
	{
	  for (int i = 0; i< ncomps; i++)
	    force.mult(1.0/density[i],i);
	}
#endif
    }
}

//
// Virtual access function for getting the forcing terms for the 
// tracers.  
//
void
PorousMedia::getForce_Tracer (FArrayBox& force,
			      int        gridno,
			      int        ngrow,
			      int        scomp,
			      int        ncomp,
			      const Real time)
{      
  force.resize(BoxLib::grow(grids[gridno],ngrow),ncomp);

  force.setVal(0.);
  if (do_source_term)
  {   
      BoxLib::Abort("Sources no longer supported");
#if 0
      const Real* dx = geom.CellSize();
      for (int i = 0; i< source_array.size(); i++)
	if (!source_array[i].var_type.compare("tracer"))
	  source_array[i].setVal(force, region_array, dx); 
#endif
  }
}

//
// Fills ghost cells of states.
//
void
PorousMedia::FillStateBndry (Real time,
                             int  state_idx,
                             int  src_comp, 
                             int  ncomp) 
{
  MultiFab& S = get_data(state_idx,time);

  if (S.nGrow() == 0)
    return;

  for (FillPatchIterator fpi(*this,S,S.nGrow(),time,state_idx,src_comp,ncomp);
       fpi.isValid();
       ++fpi)
    {
      //
      // Fill all ghost cells interior & exterior to valid region.
      //
      BoxList boxes = BoxLib::boxDiff(fpi().box(),grids[fpi.index()]);
      for (BoxList::iterator bli = boxes.begin(); bli != boxes.end(); ++bli)
        {
	  S[fpi.index()].copy(fpi(),*bli,0,*bli,src_comp,ncomp);
        }
    }
}


void 
PorousMedia::getViscTerms (MultiFab& visc_terms,
                           int       src_comp, 
                           int       ncomp,
                           Real      time)
{
  // 
  // Initialize all viscous terms to zero
  //
  const int nGrow = visc_terms.nGrow();
  visc_terms.setVal(0,0,ncomp,nGrow);

  //
  // Get Scalar Diffusive Terms
  //
  const int first_scal = src_comp;
  const int num_scal   = ncomp;

  if (num_scal > 0)

    {
      for (int icomp = first_scal; icomp < first_scal+num_scal; icomp++)
	{
	  if (is_diffusive[icomp])
	    {
	      MultiFab** cmp_diffn = 0;
	  
	      if (variable_scal_diff)
		{
		  diffusion->allocFluxBoxesLevel(cmp_diffn, 0, 1);
		  getDiffusivity(cmp_diffn, time, icomp, 0, 1);
		}
	      diffusion->getViscTerms(visc_terms,src_comp,icomp,
				      time,0,cmp_diffn);
	      if (variable_scal_diff)
		{
		  diffusion->removeFluxBoxesLevel(cmp_diffn);
		}	  
	    }
	}

      //
      // Get Capillary Diffusive Terms at time n
      //
      if (have_capillary)
	{
	  int nc = 0;
	  MultiFab** cmp_pcn = 0;
	  diffusion->allocFluxBoxesLevel(cmp_pcn,0,1);

	  calcCapillary(time);
	  calcDiffusivity_CPL(cmp_pcn,lambda_cc);

	  // multiply by kedge
	  for (int dir = 0; dir < BL_SPACEDIM; dir++)
	    {
	      for (MFIter mfi(*cmp_pcn[dir]); mfi.isValid(); ++mfi)
		(*cmp_pcn[dir])[mfi].mult(kpedge[dir][mfi],0,0,1);
	      (*cmp_pcn[dir]).FillBoundary();
	    }

	  diffusion->getCplViscTerms(visc_terms,nc,time,density.dataPtr(),0,
				     cmp_pcn,pcn_cc);
	  diffusion->removeFluxBoxesLevel(cmp_pcn);
	}	
    }

  //
  // Ensure consistent grow cells
  //    
  if (nGrow > 0)
    {
      for (MFIter mfi(visc_terms); mfi.isValid(); ++mfi)
        {
	  FArrayBox& vt  = visc_terms[mfi];
	  const Box& box = mfi.validbox();
	  FORT_VISCEXTRAP(vt.dataPtr(),ARLIM(vt.loVect()),ARLIM(vt.hiVect()),
			  box.loVect(),box.hiVect(),&ncomp);
        }
      visc_terms.FillBoundary(0,ncomp);
      //
      // Note: this is a special periodic fill in that we want to
      // preserve the extrapolated grow values when periodic --
      // usually we preserve only valid data.  The scheme relies on
      // the fact that there is good data in the "non-periodic" grow cells.
      // ("good" data produced via VISCEXTRAP above)
      //
      geom.FillPeriodicBoundary(visc_terms,0,ncomp,true);
    }
}

//
// Functions for calculating the variable viscosity and diffusivity.
// These default to setting the variable viscosity and diffusivity arrays
// to the values in visc_coef and diff_coef.  These functions would
// need to be replaced in any class derived from PorousMedia that
// wants variable coefficients.
//

void 
PorousMedia::calcDiffusivity (const Real time, 
			      const int  src_comp, 
			      const int  ncomp)
{
  //
  // NOTE:  The component numbers passed into PorousMedia::calcDiffusivity
  //        correspond to the components in the state.  
  //

  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calcDiffusivity()");

  MultiFab& S = get_data(State_Type,time);
  //
  // Select time level to work with (N or N+1)
  //
  const TimeLevel whichTime = which_time(State_Type,time);
    
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
    
  // diffn_cc and diffnp1_cc are in PorousMedia class.
  MultiFab* diff_cc         = (whichTime == AmrOldTime) ? diffn_cc : diffnp1_cc;
  const int nGrow           = 1;

  Array<Real> const_diff_coef(ncomp);
  for (int i=0;i<ncomp;i++)
    const_diff_coef[i] = visc_coef[i];
    
  //
  // Calculate diffusivity
  //
  for (FillPatchIterator fpi(*this,S,nGrow,time,State_Type,0,ncomp);
       fpi.isValid();
       ++fpi)
    {
      const int idx   = fpi.index();
      const Box box   = BoxLib::grow(grids[idx],nGrow);
      const int vflag = -1;

      FArrayBox&  Sfab  = fpi();
      const Real* ndat  = Sfab.dataPtr(); 
      const int*  n_lo  = Sfab.loVect();
      const int*  n_hi  = Sfab.hiVect();

      const Real* ddat  = (*diff_cc)[fpi].dataPtr(); 
      const int*  d_lo  = (*diff_cc)[fpi].loVect();
      const int*  d_hi  = (*diff_cc)[fpi].hiVect();

      const Real* pdat  = (*rock_phi)[fpi].dataPtr();
      const int*  p_lo  = (*rock_phi)[fpi].loVect();
      const int*  p_hi  = (*rock_phi)[fpi].hiVect();

      BL_ASSERT(box == fpi().box());
      FORT_SPECTEMPVISC(box.loVect(),box.hiVect(),
			ndat, ARLIM(n_lo), ARLIM(n_hi),
			ddat, ARLIM(d_lo), ARLIM(d_hi),
			pdat, ARLIM(p_lo),ARLIM(p_hi),
			const_diff_coef.dataPtr(), &ncomp, &vflag);
    }
}

void 
PorousMedia::getDiffusivity (MultiFab*  diffusivity[BL_SPACEDIM],
			     const Real time,
			     const int  state_comp,
			     const int  dst_comp,
			     const int  ncomp)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getDiffusivity()");

  //
  // Pick correct diffusivity component
  //
  int diff_comp = state_comp;

  //
  // Select time level to work with (N or N+1)
  //   
  const TimeLevel whichTime = which_time(State_Type,time);
    
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

  MultiFab* diff_cc  = (whichTime == AmrOldTime) ? diffn_cc : diffnp1_cc;

  //
  // Fill edge-centered diffusivities based on diffn_cc or diffnp1_cc
  //
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      for (MFIter ecMfi(*diffusivity[dir]); ecMfi.isValid(); ++ecMfi)
        {
	  center_to_edge_plain((*diff_cc)[ecMfi],(*diffusivity[dir])[ecMfi],
			       diff_comp,dst_comp,ncomp);
        }
    }
}

void 
PorousMedia::calcDiffusivity_CPL (MultiFab*  diffusivity[BL_SPACEDIM],
				  const Real time)
{

  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calcDiffusivity_CPL()");
  const TimeLevel whichTime = which_time(State_Type,time);
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);    
  MultiFab* lcc = (whichTime == AmrOldTime) ? lambda_cc : lambdap1_cc;
  calcDiffusivity_CPL(diffusivity,lcc);
}

void 
PorousMedia::calcDiffusivity_CPL (MultiFab*        diffusivity[BL_SPACEDIM],
				  const MultiFab*  lbd_cc)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calcDiffusivity_CPL()");   

  const int*  domlo    = geom.Domain().loVect();
  const int*  domhi    = geom.Domain().hiVect();
  const int   ncomp    = (*diffusivity[0]).nComp();
  for (MFIter mfi(*lbd_cc); mfi.isValid(); ++mfi)
    {
      const int idx   = mfi.index();
      const int* lo   = mfi.validbox().loVect();
      const int* hi   = mfi.validbox().hiVect();

      const int* lbd_lo  = (*lbd_cc)[idx].loVect();
      const int* lbd_hi  = (*lbd_cc)[idx].hiVect();
      const Real* lbddat = (*lbd_cc)[idx].dataPtr();

      const int* dfx_lo  = (*diffusivity[0])[idx].loVect();
      const int* dfx_hi  = (*diffusivity[0])[idx].hiVect();
      const Real* dfxdat = (*diffusivity[0])[idx].dataPtr();

      const int* dfy_lo  = (*diffusivity[1])[idx].loVect();
      const int* dfy_hi  = (*diffusivity[1])[idx].hiVect();
      const Real* dfydat = (*diffusivity[1])[idx].dataPtr();

#if(BL_SPACEDIM==3)
      const int* dfz_lo  = (*diffusivity[2])[idx].loVect();
      const int* dfz_hi  = (*diffusivity[2])[idx].hiVect();
      const Real* dfzdat = (*diffusivity[2])[idx].dataPtr();
#endif
      Array<int> bc;
      bc = getBCArray(State_Type,idx,0,1);
      FORT_GETDIFFUSE_CPL(lbddat, ARLIM(lbd_lo), ARLIM(lbd_hi),
			  dfxdat, ARLIM(dfx_lo), ARLIM(dfx_hi),
			  dfydat, ARLIM(dfy_lo), ARLIM(dfy_hi),
#if(BL_SPACEDIM==3)
			  dfzdat, ARLIM(dfz_lo), ARLIM(dfz_hi),
#endif
			  lo,hi,domlo,domhi,bc.dataPtr(),&ncomp);
    }
  // multiply by kedge
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      for (MFIter ecMfi(*diffusivity[dir]); ecMfi.isValid(); ++ecMfi)
        {
	  (*diffusivity[dir])[ecMfi].mult(kpedge[dir][ecMfi],0,0,1);
        }
      (*diffusivity[dir]).FillBoundary();
    }  
}

void 
PorousMedia::calcDiffusivity_CPL_dp (MultiFab* diffusivity[BL_SPACEDIM],
				     const MultiFab* lbd_cc,
				     const Real time,
				     const int ncomp)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calcDiffusivity_CPL_dp()");

  MultiFab& S = get_data(State_Type,time);
  const int nGrow = 1;    

  const int*  domlo    = geom.Domain().loVect();
  const int*  domhi    = geom.Domain().hiVect();
  const int n_cpl_coef = cpl_coef->nComp(); 

  // Calculate diffusivity with the dp/ds term.
  for (FillPatchIterator fpi(*this,S,nGrow,time,State_Type,0,ncomps);
       fpi.isValid();
       ++fpi)
    {
      const int idx   = fpi.index();
      const Box box   = BoxLib::grow(grids[idx],nGrow);

      BL_ASSERT(box == fpi().box());

      FArrayBox Htmp(box,1);
      Htmp.setVal(0.);
      const Real* hdat = Htmp.dataPtr();

      const Real* ndat = fpi().dataPtr(); 
      const int*  n_lo = fpi().loVect();
      const int*  n_hi = fpi().hiVect();

      const Real* lbddat = (*lbd_cc)[fpi].dataPtr();
      const int* lbd_lo  = (*lbd_cc)[fpi].loVect();
      const int* lbd_hi  = (*lbd_cc)[fpi].hiVect();	

      const Real* pdat   = (*rock_phi)[fpi].dataPtr();
      const int* p_lo    = (*rock_phi)[fpi].loVect();
      const int* p_hi    = (*rock_phi)[fpi].hiVect();

      const Real* kdat   = (*kappa)[fpi].dataPtr();
      const int* k_lo    = (*kappa)[fpi].loVect();
      const int* k_hi    = (*kappa)[fpi].hiVect();

      const int* lo      = fpi.validbox().loVect();
      const int* hi      = fpi.validbox().hiVect();

      const int* dfx_lo  = (*diffusivity[0])[idx].loVect();
      const int* dfx_hi  = (*diffusivity[0])[idx].hiVect();
      const Real* dfxdat = (*diffusivity[0])[idx].dataPtr();

      const int* dfy_lo  = (*diffusivity[1])[idx].loVect();
      const int* dfy_hi  = (*diffusivity[1])[idx].hiVect();
      const Real* dfydat = (*diffusivity[1])[idx].dataPtr();

#if(BL_SPACEDIM==3)
      const int* dfz_lo  = (*diffusivity[2])[idx].loVect();
      const int* dfz_hi  = (*diffusivity[2])[idx].hiVect();
      const Real* dfzdat = (*diffusivity[2])[idx].dataPtr();
#endif

      const Real* cpdat  = (*cpl_coef)[fpi].dataPtr(); 
      const int*  cp_lo  = (*cpl_coef)[fpi].loVect();
      const int*  cp_hi  = (*cpl_coef)[fpi].hiVect();

      Array<int> bc;
      bc = getBCArray(State_Type,idx,0,1);

      FORT_GETDIFFUSE_CPL_dp(ndat, hdat, ARLIM(n_lo), ARLIM(n_hi),
			     lbddat, ARLIM(lbd_lo), ARLIM(lbd_hi),
			     dfxdat, ARLIM(dfx_lo), ARLIM(dfx_hi),
			     dfydat, ARLIM(dfy_lo), ARLIM(dfy_hi),
#if(BL_SPACEDIM==3)
			     dfzdat, ARLIM(dfz_lo), ARLIM(dfz_hi),
#endif
			     pdat, ARLIM(p_lo), ARLIM(p_hi),
			     kdat, ARLIM(k_lo), ARLIM(k_hi),
			     cpdat, ARLIM(cp_lo), ARLIM(cp_hi),
			     &n_cpl_coef,
			     lo, hi, domlo, domhi,
			     bc.dataPtr(), &ncomp);
    }
    
  // multiply by kedge
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      for (MFIter ecMfi(*diffusivity[dir]); ecMfi.isValid(); ++ecMfi)
        {
	  (*diffusivity[dir])[ecMfi].mult(kpedge[dir][ecMfi],0,0,1);
        }
      (*diffusivity[dir]).FillBoundary();
    }
}

#ifdef MG_USE_FBOXLIB
void 
PorousMedia::calc_richard_coef (MultiFab*        diffusivity[BL_SPACEDIM],
				const MultiFab*  lbd_cc,
				const MultiFab*  umac,
				const int        nc,
				const int        do_upwind,
				double           time)
{

  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calc_richard_coef()");

  const int*  domlo    = geom.Domain().loVect();
  const int*  domhi    = geom.Domain().hiVect();
  const int ncp1 = nc + 1;

  // Calculate diffusivity for the richard's equation
  for (MFIter mfi(*lbd_cc); mfi.isValid(); ++mfi)
    {

      const int idx      = mfi.index();
      const int* lo      = mfi.validbox().loVect();
      const int* hi      = mfi.validbox().hiVect();
      
      const int* lbd_lo  = (*lbd_cc)[idx].loVect();
      const int* lbd_hi  = (*lbd_cc)[idx].hiVect();
      const Real* lbddat = (*lbd_cc)[idx].dataPtr();

      const int* ux_lo   = umac[0][idx].loVect();
      const int* ux_hi   = umac[0][idx].hiVect();
      const Real* uxdat  = umac[0][idx].dataPtr();

      const int* uy_lo   = umac[1][idx].loVect();
      const int* uy_hi   = umac[1][idx].hiVect();
      const Real* uydat  = umac[1][idx].dataPtr();

      const int* dfx_lo  = (*diffusivity[0])[idx].loVect();
      const int* dfx_hi  = (*diffusivity[0])[idx].hiVect();
      const Real* dfxdat = (*diffusivity[0])[idx].dataPtr();

      const int* dfy_lo  = (*diffusivity[1])[idx].loVect();
      const int* dfy_hi  = (*diffusivity[1])[idx].hiVect();
      const Real* dfydat = (*diffusivity[1])[idx].dataPtr();

#if(BL_SPACEDIM==3)      
      const int* uz_lo   = umac[2][idx].loVect();
      const int* uz_hi   = umac[2][idx].hiVect();
      const Real* uzdat  = umac[2][idx].dataPtr();

      const int* dfz_lo  = (*diffusivity[2])[idx].loVect();
      const int* dfz_hi  = (*diffusivity[2])[idx].hiVect();
      const Real* dfzdat = (*diffusivity[2])[idx].dataPtr();
#endif

      Array<int> bc;
      bc = getBCArray(State_Type,idx,0,1);

      FORT_RICHARD_COEF(lbddat, ARLIM(lbd_lo), ARLIM(lbd_hi),
			dfxdat, ARLIM(dfx_lo), ARLIM(dfx_hi),
			dfydat, ARLIM(dfy_lo), ARLIM(dfy_hi),
#if(BL_SPACEDIM==3)
			dfzdat, ARLIM(dfz_lo), ARLIM(dfz_hi),
#endif
			uxdat, ARLIM(ux_lo), ARLIM(ux_hi),
			uydat, ARLIM(uy_lo), ARLIM(uy_hi),
#if(BL_SPACEDIM==3)
			uzdat, ARLIM(uz_lo), ARLIM(uz_hi),
#endif
			lo,hi,domlo,domhi,bc.dataPtr(),
			&ncp1,&do_upwind);

    }

  // multiply by kedge
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      for (MFIter ecMfi(*diffusivity[dir]); ecMfi.isValid(); ++ecMfi)
	  (*diffusivity[dir])[ecMfi].mult(kpedge[dir][ecMfi],0,0,1);
      (*diffusivity[dir]).FillBoundary();
    }

  FArrayBox inflow;
  for (OrientationIter oitr; oitr; ++oitr) {
      Orientation face = oitr();
      if (get_inflow_velocity(face,inflow,time)) {
	int shift = ( face.isHigh() ? -1 : +1 );
	inflow.setVal(0.);
          inflow.shiftHalf(face.coordDir(),shift);
          for (MFIter mfi(*diffusivity[face.coordDir()]); mfi.isValid(); ++mfi) {
	    FArrayBox& u = (*diffusivity[face.coordDir()])[mfi];
              Box ovlp = inflow.box() & u.box();
              if (ovlp.ok()) {
  		u.copy(inflow);
              }
          }
      }
  }
}

void 
PorousMedia::calc_richard_jac (MultiFab*       diffusivity[BL_SPACEDIM],
                               MultiFab*       dalpha,
			       const MultiFab* lbd_cc,                                
			       const MultiFab* umac,
			       Real            time,
			       Real            dt,
			       int             nc,
			       int             do_upwind,
			       bool            do_richard_sat_solve)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calc_richard_jac()");

  MultiFab& S = get_data(State_Type,time);
  MultiFab& P = get_data(Press_Type,time);
  const int nGrow = 1;    
  //
  // Select time level to work with (N or N+1)
  //
  const TimeLevel whichTime = which_time(State_Type,time);
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
  MultiFab* pc_cc = (whichTime == AmrOldTime) ? pcn_cc : pcnp1_cc;

  const Real* dx       = geom.CellSize();
  const int*  domlo    = geom.Domain().loVect();
  const int*  domhi    = geom.Domain().hiVect();
  const int n_cpl_coef = cpl_coef->nComp(); 
  const int n_kr_coef  = kr_coef->nComp(); 

  //FIXME: analytical jacobian is not working
  bool do_analytic_jac = false;//true;

#ifdef BL_USE_PETSC
  PetscErrorCode ierr;
  Mat& J = PMAmr::GetLayout().Jacobian();
  Vec& JRowScale = PMAmr::GetLayout().JRowScale();
  BaseFab<int> nodeNums;
  const BCRec& theBC = AmrLevel::desc_lst[Press_Type].getBC(0);
  const Layout& layout = PMAmr::GetLayout();

#endif

  for (FillPatchIterator fpi(*this,S,nGrow,time,State_Type,0,ncomps);
       fpi.isValid();
       ++fpi)
    {
      const int idx   = fpi.index();
      const Box box   = BoxLib::grow(grids[idx],nGrow);

      BL_ASSERT(box == fpi().box());

      const Real* ndat = fpi().dataPtr(); 
      const int*  n_lo = fpi().loVect();
      const int*  n_hi = fpi().hiVect();

      const Real* prdat  = P[fpi].dataPtr(); 
      const int*  pr_lo  = P[fpi].loVect();
      const int*  pr_hi  = P[fpi].hiVect();

      const Real* lbddat = (*lbd_cc)[fpi].dataPtr();
      const int* lbd_lo  = (*lbd_cc)[fpi].loVect();
      const int* lbd_hi  = (*lbd_cc)[fpi].hiVect();

      const Real* pcdat  = (*pc_cc)[fpi].dataPtr(); 
      const int*  pc_lo  = (*pc_cc)[fpi].loVect();
      const int*  pc_hi  = (*pc_cc)[fpi].hiVect();

      const Real* pdat   = (*rock_phi)[fpi].dataPtr();
      const int* p_lo    = (*rock_phi)[fpi].loVect();
      const int* p_hi    = (*rock_phi)[fpi].hiVect();

      const Real* kdat   = (*kappa)[fpi].dataPtr();
      const int* k_lo    = (*kappa)[fpi].loVect();
      const int* k_hi    = (*kappa)[fpi].hiVect();

      const int* lo      = fpi.validbox().loVect();
      const int* hi      = fpi.validbox().hiVect();

      const int* ux_lo   = umac[0][idx].loVect();
      const int* ux_hi   = umac[0][idx].hiVect();
      const Real* uxdat  = umac[0][idx].dataPtr();

      const int* uy_lo   = umac[1][idx].loVect();
      const int* uy_hi   = umac[1][idx].hiVect();
      const Real* uydat  = umac[1][idx].dataPtr();

      const int* dfx_lo  = (*diffusivity[0])[idx].loVect();
      const int* dfx_hi  = (*diffusivity[0])[idx].hiVect();
      const Real* dfxdat = (*diffusivity[0])[idx].dataPtr();

      const int* dfy_lo  = (*diffusivity[1])[idx].loVect();
      const int* dfy_hi  = (*diffusivity[1])[idx].hiVect();
      const Real* dfydat = (*diffusivity[1])[idx].dataPtr();

      const int* kpx_lo  = kpedge[0][idx].loVect();
      const int* kpx_hi  = kpedge[0][idx].hiVect();
      const Real* kpxdat = kpedge[0][idx].dataPtr();

      const int* kpy_lo  = kpedge[1][idx].loVect();
      const int* kpy_hi  = kpedge[1][idx].hiVect();
      const Real* kpydat = kpedge[1][idx].dataPtr();

#if(BL_SPACEDIM==3)    
      const int* uz_lo   = umac[2][idx].loVect();
      const int* uz_hi   = umac[2][idx].hiVect();
      const Real* uzdat  = umac[2][idx].dataPtr();

      const int* dfz_lo  = (*diffusivity[2])[idx].loVect();
      const int* dfz_hi  = (*diffusivity[2])[idx].hiVect();
      const Real* dfzdat = (*diffusivity[2])[idx].dataPtr();

      const int* kpz_lo  = kpedge[2][idx].loVect();
      const int* kpz_hi  = kpedge[2][idx].hiVect();
      const Real* kpzdat = kpedge[2][idx].dataPtr();
#endif 
      const Real* krdat  = (*kr_coef)[fpi].dataPtr(); 
      const int*  kr_lo  = (*kr_coef)[fpi].loVect();
      const int*  kr_hi  = (*kr_coef)[fpi].hiVect();
      const Real* cpdat  = (*cpl_coef)[fpi].dataPtr(); 
      const int*  cp_lo  = (*cpl_coef)[fpi].loVect();
      const int*  cp_hi  = (*cpl_coef)[fpi].hiVect();

      if (dalpha) {
          const Real* adat = (*dalpha)[fpi].dataPtr();
          const int*  a_lo = (*dalpha)[fpi].loVect();
          const int*  a_hi = (*dalpha)[fpi].hiVect();

          FORT_RICHARD_ALPHA(adat, ARLIM(a_lo), ARLIM(a_hi),
                             ndat, ARLIM(n_lo), ARLIM(n_hi),
                             pdat, ARLIM(p_lo), ARLIM(p_hi),
                             kdat, ARLIM(k_lo), ARLIM(k_hi),
                             cpdat, ARLIM(cp_lo), ARLIM(cp_hi),
                             &n_cpl_coef, lo, hi);
      }

      Array<int> bc;
      bc = getBCArray(Press_Type,idx,0,1);

      if (do_analytic_jac) 
	FORT_RICHARD_AJAC(ndat, ARLIM(n_lo), ARLIM(n_hi),
			  dfxdat, ARLIM(dfx_lo), ARLIM(dfx_hi),
			  dfydat, ARLIM(dfy_lo), ARLIM(dfy_hi),
#if(BL_SPACEDIM==3)
			  dfzdat, ARLIM(dfz_lo), ARLIM(dfz_hi),
#endif	
			  uxdat, ARLIM(ux_lo), ARLIM(ux_hi),
			  uydat, ARLIM(uy_lo), ARLIM(uy_hi),
#if(BL_SPACEDIM==3)
			  uzdat, ARLIM(uz_lo), ARLIM(uz_hi),
#endif
			  kpxdat, ARLIM(kpx_lo), ARLIM(kpx_hi),
			  kpydat, ARLIM(kpy_lo), ARLIM(kpy_hi),
#if(BL_SPACEDIM==3)
			  kpzdat, ARLIM(kpz_lo), ARLIM(kpz_hi),
#endif
			  lbddat, ARLIM(lbd_lo), ARLIM(lbd_hi),
			  pcdat, ARLIM(pc_lo), ARLIM(pc_hi),
			  pdat, ARLIM(p_lo), ARLIM(p_hi),
			  kdat, ARLIM(k_lo), ARLIM(k_hi),
			  krdat, ARLIM(kr_lo), ARLIM(kr_hi), &n_kr_coef,
			  cpdat, ARLIM(cp_lo), ARLIM(cp_hi), &n_cpl_coef,
			  lo, hi, domlo, domhi, dx, bc.dataPtr(), 
			  rinflow_bc_lo.dataPtr(),rinflow_bc_hi.dataPtr(), 
			  &do_upwind);
      else
	{
	  if (do_richard_sat_solve)
	    FORT_RICHARD_NJAC(ndat,   ARLIM(n_lo), ARLIM(n_hi),
			      dfxdat, ARLIM(dfx_lo), ARLIM(dfx_hi),
			      dfydat, ARLIM(dfy_lo), ARLIM(dfy_hi),
#if(BL_SPACEDIM==3)
			      dfzdat, ARLIM(dfz_lo), ARLIM(dfz_hi),
#endif	
			      uxdat, ARLIM(ux_lo), ARLIM(ux_hi),
			      uydat, ARLIM(uy_lo), ARLIM(uy_hi),
#if(BL_SPACEDIM==3)
			      uzdat, ARLIM(uz_lo), ARLIM(uz_hi),
#endif
			      kpxdat, ARLIM(kpx_lo), ARLIM(kpx_hi),
			      kpydat, ARLIM(kpy_lo), ARLIM(kpy_hi),
#if(BL_SPACEDIM==3)
			      kpzdat, ARLIM(kpz_lo), ARLIM(kpz_hi),
#endif
			      lbddat, ARLIM(lbd_lo), ARLIM(lbd_hi),
			      pcdat, ARLIM(pc_lo), ARLIM(pc_hi),
			      pdat, ARLIM(p_lo), ARLIM(p_hi),
			      kdat, ARLIM(k_lo), ARLIM(k_hi),
			      krdat, ARLIM(kr_lo), ARLIM(kr_hi), &n_kr_coef,
			      cpdat, ARLIM(cp_lo), ARLIM(cp_hi), &n_cpl_coef,
			      lo, hi, domlo, domhi, dx, bc.dataPtr(), 
			      rinflow_bc_lo.dataPtr(),rinflow_bc_hi.dataPtr(), 
			      &richard_perturbation_scale_for_J, &do_upwind);
	  else
	    {
	      FORT_RICHARD_NJAC2(dfxdat, ARLIM(dfx_lo), ARLIM(dfx_hi),
				 dfydat, ARLIM(dfy_lo), ARLIM(dfy_hi),
#if(BL_SPACEDIM==3)
				 dfzdat, ARLIM(dfz_lo), ARLIM(dfz_hi),
#endif	
				 uxdat, ARLIM(ux_lo), ARLIM(ux_hi),
				 uydat, ARLIM(uy_lo), ARLIM(uy_hi),
#if(BL_SPACEDIM==3)
				 uzdat, ARLIM(uz_lo), ARLIM(uz_hi),
#endif
				 kpxdat, ARLIM(kpx_lo), ARLIM(kpx_hi),
				 kpydat, ARLIM(kpy_lo), ARLIM(kpy_hi),
#if(BL_SPACEDIM==3)
				 kpzdat, ARLIM(kpz_lo), ARLIM(kpz_hi),
#endif
				 lbddat, ARLIM(lbd_lo), ARLIM(lbd_hi),
				 prdat, ARLIM(pr_lo), ARLIM(pr_hi),
				 pdat, ARLIM(p_lo), ARLIM(p_hi),
				 kdat, ARLIM(k_lo), ARLIM(k_hi),
				 krdat, ARLIM(kr_lo), ARLIM(kr_hi), &n_kr_coef,
				 cpdat, ARLIM(cp_lo), ARLIM(cp_hi), &n_cpl_coef,
				 lo, hi, domlo, domhi, dx, bc.dataPtr(), 
				 rinflow_bc_lo.dataPtr(),rinflow_bc_hi.dataPtr(), 
				 &richard_perturbation_scale_for_J, &do_upwind);
	    }
	}

#ifdef BL_USE_PETSC
      const Box& vbox = grids[idx];
      Box gbox = Box(vbox).grow(1);
      nodeNums.resize(gbox,1);
      layout.SetNodeIds(nodeNums,level,idx);

      Array<int> cols(1+2*BL_SPACEDIM);
      Array<int> rows(1);
      Array<Real> vals(cols.size());
      FArrayBox* wrk[BL_SPACEDIM];
      for (int d=0; d<BL_SPACEDIM; ++d) {
          wrk[d] = &((*diffusivity[d])[idx]);
      }

      for (IntVect iv(vbox.smallEnd()), iEnd=vbox.bigEnd(); iv<=iEnd; vbox.next(iv))
      {
          cols[0] = nodeNums(iv);
          if (cols[0]>=0) {
              rows[0] = cols[0];
              vals[0] = (dalpha ? (*dalpha)[fpi](iv,0)  :  0);
              Real rdt = (dt>0  ?  density[nc]*dt : 1); // The "b" factor
              int cnt = 1;
              for (int d=0; d<BL_SPACEDIM; ++d) {
                  vals[0] -= rdt * (*wrk[d])(iv,2);
                  IntVect ivp = iv + BoxLib::BASISV(d);
                  int np = nodeNums(ivp,0);
                  if (np>=0) {
                      cols[cnt]  = np; 
                      vals[cnt]  = -rdt * (*wrk[d])(iv,0);
                      cnt++;
                  }
                  else {
                      if (theBC.hi()[d]==FOEXTRAP) {
                          vals[0] -= rdt * (*wrk[d])(iv,0);
                      }
                  }
                  
                  IntVect ivn = iv - BoxLib::BASISV(d);
                  int nn = nodeNums(ivn,0);
                  if (nn>=0) {
                      cols[cnt]  = nn; 
                      vals[cnt]  = -rdt * (*wrk[d])(iv,1);
                      cnt++;
                  }
                  else {
                      if (theBC.lo()[d]==FOEXTRAP) {
                          vals[0] -= rdt * (*wrk[d])(iv,1);
                      }
                  }
              }
              
              // Normalize matrix entries
#if 1
              Real max_abs = 1;
#else
              Real max_abs = 0;
              for (int n=0; n<cnt; ++n) {
                  max_abs = std::max(max_abs,std::abs(vals[n]));
              }
              max_abs = 1/max_abs;
              for (int n=0; n<cnt; ++n) {
                  vals[n] *= max_abs;
              }
#endif
              ierr = MatSetValues(J,rows.size(),rows.dataPtr(),cnt,cols.dataPtr(),vals.dataPtr(),INSERT_VALUES); CHKPETSC(ierr);
              ierr = VecSetValues(JRowScale,1,&(rows[0]),&max_abs,INSERT_VALUES);
          }
      }
#endif

    }

#ifdef BL_USE_PETSC
  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKPETSC(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKPETSC(ierr);
  ierr = VecAssemblyBegin(JRowScale); CHKPETSC(ierr);  
  ierr = VecAssemblyEnd(JRowScale); CHKPETSC(ierr);  
#endif
}

void 
PorousMedia::calc_richard_alpha (MultiFab*     alpha,
				 const Real    time)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calc_richard_alpha()");

  const int nGrow = 1;    
  MultiFab S(grids,ncomps,nGrow);
  const int n_cpl_coef = cpl_coef->nComp(); 
  for (FillPatchIterator fpi(*this,S,nGrow,time,State_Type,0,ncomps);
       fpi.isValid();
       ++fpi)
    {
      S[fpi].copy(fpi());
    }
  calc_richard_alpha(alpha,S);
}

void 
PorousMedia::calc_richard_alpha (MultiFab*       alpha,
				 const MultiFab& S) const
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calc_richard_alpha()");
  int nGrow = S.nGrow();
  const int n_cpl_coef = cpl_coef->nComp(); 
  for (MFIter mfi(S); mfi.isValid(); ++mfi) 
    {
      const int idx   = mfi.index();
      const Box box   = BoxLib::grow(grids[idx],nGrow);

      const int* lo      = mfi.validbox().loVect();
      const int* hi      = mfi.validbox().hiVect();

      const Real* ndat = S[mfi].dataPtr(); 
      const int*  n_lo = S[mfi].loVect();
      const int*  n_hi = S[mfi].hiVect();

      const Real* adat = (*alpha)[mfi].dataPtr();
      const int*  a_lo = (*alpha)[mfi].loVect();
      const int*  a_hi = (*alpha)[mfi].hiVect();

      const Real* pdat   = (*rock_phi)[mfi].dataPtr();
      const int* p_lo    = (*rock_phi)[mfi].loVect();
      const int* p_hi    = (*rock_phi)[mfi].hiVect();

      const Real* kdat   = (*kappa)[mfi].dataPtr();
      const int* k_lo    = (*kappa)[mfi].loVect();
      const int* k_hi    = (*kappa)[mfi].hiVect();

      const Real* cpdat  = (*cpl_coef)[mfi].dataPtr(); 
      const int*  cp_lo  = (*cpl_coef)[mfi].loVect();
      const int*  cp_hi  = (*cpl_coef)[mfi].hiVect();

			  
      FORT_RICHARD_ALPHA(adat, ARLIM(a_lo), ARLIM(a_hi),
			 ndat, ARLIM(n_lo), ARLIM(n_hi),
			 pdat, ARLIM(p_lo), ARLIM(p_hi),
			 kdat, ARLIM(k_lo), ARLIM(k_hi),
			 cpdat, ARLIM(cp_lo), ARLIM(cp_hi),
			 &n_cpl_coef, lo, hi);
    }
}

void 
PorousMedia::calc_richard_velbc (MultiFab& res, 
				 MultiFab* u_phase,
				 const Real dt)  
{ 
  //
  // Add boundary condition to residual
  //
  const int* domlo = geom.Domain().loVect(); 
  const int* domhi = geom.Domain().hiVect();
  const Real* dx   = geom.CellSize();

  for (MFIter mfi(res); mfi.isValid(); ++mfi)
    {
      const int* lo = mfi.validbox().loVect();
      const int* hi = mfi.validbox().hiVect();
	
      FArrayBox& rg       = res[mfi];  
      FArrayBox& ux       = u_phase[0][mfi];
      FArrayBox& uy       = u_phase[1][mfi];
      DEF_LIMITS (rg,rg_dat,rglo,rghi);
      DEF_LIMITS (ux,ux_dat,uxlo,uxhi);
      DEF_LIMITS (uy,uy_dat,uylo,uyhi);

#if (BL_SPACEDIM == 3)
      FArrayBox& uz       = u_phase[2][mfi];
      DEF_LIMITS (uz,uz_dat,uzlo,uzhi);
#endif
      FORT_RICHARD_VELBC (rg_dat, ARLIM(rglo), ARLIM(rghi),
			  ux_dat, ARLIM(uxlo), ARLIM(uxhi),
			  uy_dat, ARLIM(uylo), ARLIM(uyhi),
#if (BL_SPACEDIM == 3)
			  uz_dat, ARLIM(uzlo), ARLIM(uzhi),
#endif
			  lo,hi,domlo,domhi,dx,
			  rinflow_bc_lo.dataPtr(),
			  rinflow_bc_hi.dataPtr(), 
			  &dt);
    }
}
#endif

void 
PorousMedia::calcCapillary (const Real time)
{
  //
  // Calculate the capillary pressure.  
  //
  MultiFab& S = get_data(State_Type,time);
  FillStateBndry(time,State_Type,0,ncomps);
  //
  // Select time level to work with (N or N+1)
  //
  const TimeLevel whichTime = which_time(State_Type,time);
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
  MultiFab* pc_cc = (whichTime == AmrOldTime) ? pcn_cc : pcnp1_cc;
  const int nGrow = 1;

  pc_cc->setVal(-20);

  for (FillPatchIterator fpi(*this,S,nGrow,time,State_Type,0,ncomps);
       fpi.isValid();
       ++fpi)
    {      
      const int idx  = fpi.index();
      const Box box   = BoxLib::grow(grids[idx],nGrow);
      BL_ASSERT(box == fpi().box());

      Array<int> s_bc;
      s_bc = getBCArray(State_Type,idx,0,1);

      calcCapillary ((*pc_cc)[fpi],fpi(),(*rock_phi)[fpi],(*kappa)[fpi],
		     (*cpl_coef)[fpi],grids[idx],s_bc);
    }
 pc_cc->FillBoundary();
}

void 
PorousMedia::calcCapillary (MultiFab* pc,
			    const MultiFab& S)
{
  //
  // Calculate the capillary pressure for a given state.  
  //
  BL_ASSERT(S.nGrow() >=1); // Assumes that boundary cells have been properly filled
  BL_ASSERT(pc->nGrow() >= 0); // Fill boundary cells (in F)
  for (MFIter mfi(S); mfi.isValid(); ++mfi) 
    {
      
      const int idx  = mfi.index();
      Array<int> s_bc;
      s_bc = getBCArray(State_Type,idx,0,1);

      calcCapillary ((*pc)[mfi],S[mfi],(*rock_phi)[mfi],(*kappa)[mfi],
		     (*cpl_coef)[mfi],grids[idx],s_bc);
    }
  pc->FillBoundary();
}

void 
PorousMedia::calcCapillary (FArrayBox&       pc,
			    const FArrayBox& S,
			    const FArrayBox& rockphi,
			    const FArrayBox& rockkappa,
			    const FArrayBox& cplcoef,
			    const Box&       grids,
			    const Array<int>& s_bc)
{
  //
  // Calculate the capillary pressure for a fab
  //
  const int n_cpl_coef = cplcoef.nComp();
  const int* lo  = grids.loVect();
  const int* hi  = grids.hiVect();
  
  const Real* ndat  = S.dataPtr(); 
  const int*  n_lo  = S.loVect();
  const int*  n_hi  = S.hiVect();

  const Real* ddat  = pc.dataPtr(); 
  const int*  d_lo  = pc.loVect();
  const int*  d_hi  = pc.hiVect();

  const Real* pdat = rockphi.dataPtr();
  const int* p_lo  = rockphi.loVect();
  const int* p_hi  = rockphi.hiVect();

  const Real* kdat = rockkappa.dataPtr();
  const int* k_lo  = rockkappa.loVect();
  const int* k_hi  = rockkappa.hiVect();

  const Real* cpdat  = cplcoef.dataPtr(); 
  const int*  cp_lo  = cplcoef.loVect();
  const int*  cp_hi  = cplcoef.hiVect();

  FORT_MK_CPL( ddat, ARLIM(d_lo), ARLIM(d_hi),
	       ndat, ARLIM(n_lo), ARLIM(n_hi),
	       pdat, ARLIM(p_lo), ARLIM(p_hi),
	       kdat, ARLIM(k_lo), ARLIM(k_hi),
	       cpdat, ARLIM(cp_lo), ARLIM(cp_hi),
	       &n_cpl_coef, lo, hi, s_bc.dataPtr());
}

void 
PorousMedia::calcInvCapillary (MultiFab& S,
			       const MultiFab& pc)
{
  //
  // Calculate inverse capillary pressure
  //    
  const int n_cpl_coef = cpl_coef->nComp();
  for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {

      FArrayBox& Sfab   = S[mfi];
      const Real* ndat  = Sfab.dataPtr(); 
      const int*  n_lo  = Sfab.loVect();
      const int*  n_hi  = Sfab.hiVect();

      const Real* ddat  = pc[mfi].dataPtr(); 
      const int*  d_lo  = pc[mfi].loVect();
      const int*  d_hi  = pc[mfi].hiVect();

      const Real* pdat = (*rock_phi)[mfi].dataPtr();
      const int* p_lo  = (*rock_phi)[mfi].loVect();
      const int* p_hi  = (*rock_phi)[mfi].hiVect();

      const Real* kdat = (*kappa)[mfi].dataPtr();
      const int* k_lo  = (*kappa)[mfi].loVect();
      const int* k_hi  = (*kappa)[mfi].hiVect();

      const Real* cpdat  = (*cpl_coef)[mfi].dataPtr(); 
      const int*  cp_lo  = (*cpl_coef)[mfi].loVect();
      const int*  cp_hi  = (*cpl_coef)[mfi].hiVect();

      FORT_MK_INV_CPL( ddat, ARLIM(d_lo), ARLIM(d_hi),
		       ndat, ARLIM(n_lo), ARLIM(n_hi),
		       pdat, ARLIM(p_lo), ARLIM(p_hi),
		       kdat, ARLIM(k_lo), ARLIM(k_hi),
		       cpdat, ARLIM(cp_lo), ARLIM(cp_hi),
		       &n_cpl_coef); 
    }
}

void 
PorousMedia::calcInvPressure (MultiFab& S,
			      const MultiFab& p)
{
  //
  // Calculate inverse pressure
  //    
  const int n_cpl_coef = cpl_coef->nComp();
  for (MFIter mfi(S); mfi.isValid(); ++mfi)
    { 
      FArrayBox pc;
      const Box& fbox = S[mfi].box(); 
      pc.resize(fbox,1);
      pc.copy(p[mfi],fbox,0,fbox,0,1);
      pc.mult(-1.0);

      FArrayBox& Sfab   = S[mfi];
      const Real* ndat  = Sfab.dataPtr(); 
      const int*  n_lo  = Sfab.loVect();
      const int*  n_hi  = Sfab.hiVect();

      const Real* ddat  = pc.dataPtr(); 
      const int*  d_lo  = pc.loVect();
      const int*  d_hi  = pc.hiVect();

      const Real* pdat = (*rock_phi)[mfi].dataPtr();
      const int* p_lo  = (*rock_phi)[mfi].loVect();
      const int* p_hi  = (*rock_phi)[mfi].hiVect();

      const Real* kdat = (*kappa)[mfi].dataPtr();
      const int* k_lo  = (*kappa)[mfi].loVect();
      const int* k_hi  = (*kappa)[mfi].hiVect();

      const Real* cpdat  = (*cpl_coef)[mfi].dataPtr(); 
      const int*  cp_lo  = (*cpl_coef)[mfi].loVect();
      const int*  cp_hi  = (*cpl_coef)[mfi].hiVect();

      FORT_MK_INV_CPL( ddat, ARLIM(d_lo), ARLIM(d_hi),
		       ndat, ARLIM(n_lo), ARLIM(n_hi),
		       pdat, ARLIM(p_lo), ARLIM(p_hi),
		       kdat, ARLIM(k_lo), ARLIM(k_hi),
		       cpdat, ARLIM(cp_lo), ARLIM(cp_hi),
		       &n_cpl_coef); 
    }
}

void 
PorousMedia::smooth_pc (MultiFab* pc)
{
  //
  // Calculate the capillary pressure for a given state.  
  //
  const int n_cpl_coef = cpl_coef->nComp();
  for (MFIter mfi(*pc); mfi.isValid(); ++mfi) 
    {
      const int idx  = mfi.index();
      const int* lo  = grids[idx].loVect();
      const int* hi  = grids[idx].hiVect();

      const Real* ddat  = (*pc)[mfi].dataPtr(); 
      const int*  d_lo  = (*pc)[mfi].loVect();
      const int*  d_hi  = (*pc)[mfi].hiVect();

      const Real* cpdat  = (*cpl_coef)[mfi].dataPtr(); 
      const int*  cp_lo  = (*cpl_coef)[mfi].loVect();
      const int*  cp_hi  = (*cpl_coef)[mfi].hiVect();

      FORT_SMOOTH_CPL( ddat, ARLIM(d_lo), ARLIM(d_hi),
		       cpdat, ARLIM(cp_lo), ARLIM(cp_hi),
		       &n_cpl_coef, lo, hi);
    }
  pc->FillBoundary();
}


void 
PorousMedia::calcLambda (const Real time, MultiFab* lbd_cc)
{
  //
  // Calculate the lambda values at cell-center. 
  //
  MultiFab& S = get_data(State_Type,time);
  FillStateBndry(time,State_Type,0,ncomps);
  MultiFab* lcc;
  if (lbd_cc == 0)
    {
      const TimeLevel whichTime = which_time(State_Type,time);
      BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
      lcc = (whichTime == AmrOldTime) ? lambda_cc : lambdap1_cc;
    }
  else
    lcc = lbd_cc;  

  const int nGrow = 1;
  const int n_kr_coef = kr_coef->nComp();  
  for (FillPatchIterator fpi(*this,S,nGrow,time,State_Type,0,ncomps);
       fpi.isValid();
       ++fpi)
    {
      const int idx   = fpi.index();
      const Box box   = BoxLib::grow(grids[idx],nGrow);
      BL_ASSERT(box == fpi().box());

      const FArrayBox& Sfab = fpi();
      const Real* ndat  = Sfab.dataPtr(); 
      const int*  n_lo  = Sfab.loVect();
      const int*  n_hi  = Sfab.hiVect();

      const Real* ddat  = (*lcc)[fpi].dataPtr(); 
      const int*  d_lo  = (*lcc)[fpi].loVect();
      const int*  d_hi  = (*lcc)[fpi].hiVect();

      const Real* krdat  = (*kr_coef)[fpi].dataPtr(); 
      const int*  kr_lo  = (*kr_coef)[fpi].loVect();
      const int*  kr_hi  = (*kr_coef)[fpi].hiVect();
	
      FORT_MK_LAMBDA( ddat, ARLIM(d_lo), ARLIM(d_hi),
		      ndat, ARLIM(n_lo), ARLIM(n_hi),
		      krdat, ARLIM(kr_lo),ARLIM(kr_hi),
		      &n_kr_coef);
    }
  lcc->FillBoundary();
}

void 
PorousMedia::calcLambda (MultiFab* lbd, const MultiFab& S)
{
  //
  // Calculate the lambda values at cell-center. 
  //   
  const int n_kr_coef = kr_coef->nComp();
  for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
      const FArrayBox& Sfab = S[mfi];
      const Real* ndat  = Sfab.dataPtr(); 
      const int*  n_lo  = Sfab.loVect();
      const int*  n_hi  = Sfab.hiVect();

      const Real* ddat  = (*lbd)[mfi].dataPtr(); 
      const int*  d_lo  = (*lbd)[mfi].loVect();
      const int*  d_hi  = (*lbd)[mfi].hiVect();

      const Real* krdat  = (*kr_coef)[mfi].dataPtr(); 
      const int*  kr_lo  = (*kr_coef)[mfi].loVect();
      const int*  kr_hi  = (*kr_coef)[mfi].hiVect();
	
      FORT_MK_LAMBDA( ddat, ARLIM(d_lo), ARLIM(d_hi),
		      ndat, ARLIM(n_lo), ARLIM(n_hi),
		      krdat, ARLIM(kr_lo),ARLIM(kr_hi),
		      &n_kr_coef);
    }
  lbd->FillBoundary();
}

void 
PorousMedia::calcDLambda (const Real time, MultiFab* dlbd_cc)
{
  //
  // Calculate the lambda values at cell-center. 
  //

  MultiFab& S = get_data(State_Type,time);

  MultiFab* dlcc;
  if (dlbd_cc == 0)
    dlcc = dlambda_cc;
  else
    dlcc = dlbd_cc;

  const int nGrow = 1;    
  const int n_kr_coef = kr_coef->nComp();
  for (FillPatchIterator fpi(*this,S,nGrow,time,State_Type,0,ncomps);
       fpi.isValid();
       ++fpi)
    {
      const int idx   = fpi.index();
      const Box box   = BoxLib::grow(grids[idx],nGrow);

      BL_ASSERT(box == fpi().box());

      FArrayBox& Sfab   = fpi();
      const Real* ndat  = Sfab.dataPtr(); 
      const int*  n_lo  = Sfab.loVect();
      const int*  n_hi  = Sfab.hiVect();

      const Real* ddat  = (*dlcc)[fpi].dataPtr(); 
      const int*  d_lo  = (*dlcc)[fpi].loVect();
      const int*  d_hi  = (*dlcc)[fpi].hiVect();

      const Real* krdat  = (*kr_coef)[fpi].dataPtr(); 
      const int*  kr_lo  = (*kr_coef)[fpi].loVect();
      const int*  kr_hi  = (*kr_coef)[fpi].hiVect();

      FORT_MK_DLAMBDA( ddat, ARLIM(d_lo), ARLIM(d_hi),
		       ndat, ARLIM(n_lo), ARLIM(n_hi), 
		       krdat, ARLIM(kr_lo),ARLIM(kr_hi),
		       &n_kr_coef);
    }

  (*dlcc).FillBoundary();
    
}

void
PorousMedia::set_overdetermined_boundary_cells (Real time)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::set_overdetermined_boundary_cells()");

}


void
PorousMedia::center_to_edge_plain (const FArrayBox& ccfab,
				   FArrayBox&       ecfab,
				   int              sComp,
				   int              dComp,
				   int              nComp)
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::center_to_edge_plain()");

  //
  // This routine fills an edge-centered FAB from a cell-centered FAB.
  // It assumes that the data in all cells of the cell-centered FAB is
  // valid and totally ignores any concept of boundary conditions.  
  // It is assummed that the cell-centered FAB fully contains the 
  // edge-centered FAB.  If anything special needs to be done at boundaries, 
  // a varient of this routine needs to be written.  See 
  // HeatTransfer::center_to_edge_fancy().
  //
  const Box&      ccbox = ccfab.box();
  const Box&      ecbox = ecfab.box();
  const IndexType ixt   = ecbox.ixType();
  //
  // Get direction for interpolation to edges
  //
  int dir = -1;
  for (int d = 0; d < BL_SPACEDIM; d++)
    if (ixt.test(d))
      dir = d;
  //
  // Miscellanious checks
  //
  BL_ASSERT(!(ixt.cellCentered()) && !(ixt.nodeCentered()));
  BL_ASSERT(BoxLib::grow(ccbox,-BoxLib::BASISV(dir)).contains(BoxLib::enclosedCells(ecbox)));
  BL_ASSERT(sComp+nComp <= ccfab.nComp() && dComp+nComp <= ecfab.nComp());
  //
  // Shift cell-centered data to edges
  //
  Box fillBox = ccbox; 
  for (int d = 0; d < BL_SPACEDIM; d++)
    if (d != dir)
      fillBox.setRange(d, ecbox.smallEnd(d), ecbox.length(d));
    
  const int isharm = def_harm_avg_cen2edge;
  FORT_CEN2EDG(fillBox.loVect(), fillBox.hiVect(),
	       ARLIM(ccfab.loVect()), ARLIM(ccfab.hiVect()),
	       ccfab.dataPtr(sComp),
	       ARLIM(ecfab.loVect()), ARLIM(ecfab.hiVect()),
	       ecfab.dataPtr(dComp),
	       &nComp, &dir, &isharm);
}

// ===================
// Boundary Conditions
// ===================

void
PorousMedia::setPhysBoundaryValues (FArrayBox& dest,
                                    int        state_indx,
                                    Real       time,
                                    int        dest_comp,
                                    int        src_comp,
                                    int        num_comp)
{
    // The default behavior of an AmrLevel, for reference:
    state[state_indx].FillBoundary(dest,time,geom.CellSize(),
                                   geom.ProbDomain(),dest_comp,src_comp,num_comp);

    if (state_indx==State_Type) {
        int n_c = std::min(ncomps-1,src_comp+num_comp-1) - src_comp + 1;
        int s_t = src_comp + n_c;
        int n_t = num_comp - n_c;

        if (n_c > 0) {
            dirichletStateBC(dest,time,src_comp,dest_comp,n_c);
        }
        if (n_t > 0) {
            dirichletTracerBC(dest,time,s_t,dest_comp+n_c,n_t);
        }
    }
    else if (state_indx==Press_Type) {
      dirichletPressBC(dest,time);
    }
}

void
PorousMedia::getDirichletFaces (Array<Orientation>& Faces,
				const int           comp_Type,
				const BCRec&        _bc)
{
  Faces.resize(0);
  for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
      if ((comp_Type == Press_Type && _bc.lo(idir) == EXT_DIR) ||
	  (comp_Type == State_Type && _bc.lo(idir) == EXT_DIR))
        {
	  const int len = Faces.size();
	  Faces.resize(len+1);
	  Faces.set(len,Orientation(idir,Orientation::low));
        }
      if ((comp_Type == Press_Type && _bc.hi(idir) == EXT_DIR) ||
	  (comp_Type == State_Type && _bc.hi(idir) == EXT_DIR))
        {
	  const int len = Faces.size();
	  Faces.resize(len+1);
	  Faces.set(len,Orientation(idir,Orientation::high));
        }
    }
}

bool
PorousMedia::grids_on_side_of_domain (const BoxArray&    _grids,
				      const Box&         _domain,
				      const Orientation& _Face) 
{
  // FIOXME: this should use the intersections code
    const int idir = _Face.coordDir();

    if (_Face.isLow())
      {
        for (int igrid = 0; igrid < _grids.size(); igrid++)
	  { 
            if (_grids[igrid].smallEnd(idir) == _domain.smallEnd(idir))
	      return true;
	  }
      }
  
    if (_Face.isHigh())
      {
        for (int igrid = 0; igrid < _grids.size(); igrid++)
	  {
            if (_grids[igrid].bigEnd(idir) == _domain.bigEnd(idir))
	      return true;
	  }
      }

    return false;
}

void
PorousMedia::dirichletStateBC (FArrayBox& fab, Real time,int sComp, int dComp, int nComp)
{
    if (bc_descriptor_map.size()) 
    {
        const Box domain = geom.Domain();
        const int* domhi = domain.hiVect();
        const int* domlo = domain.loVect();
        const Real* dx   = geom.CellSize();
        FArrayBox sdat;

        for (std::map<Orientation,BCDesc>::const_iterator
                 it=bc_descriptor_map.begin(); it!=bc_descriptor_map.end(); ++it) 
        {
            const Box bndBox = it->second.first;
            const Array<int>& face_bc_idxs = it->second.second;
            sdat.resize(bndBox,nComp); sdat.setVal(0);

            for (int i=0; i<face_bc_idxs.size(); ++i) {
                const RegionData& face_bc = bc_array[face_bc_idxs[i]]; 

		if (face_bc.Type() == "zero_total_velocity") {
		  get_inflow_density(it->first,face_bc,sdat,time);
		}
		else {
		  face_bc.apply(sdat,dx,0,nComp,time);
		}
	    }
            Box ovlp = bndBox & fab.box();
            if (ovlp.ok()) {
                fab.copy(sdat,ovlp,0,ovlp,dComp,nComp);
            }
        }
    }    
}  

void
PorousMedia::dirichletTracerBC (FArrayBox& fab, Real time, int sComp, int dComp, int nComp)
{
    BL_ASSERT(do_tracer_transport);
    
    for (int n=0; n<nComp; ++n) 
    {
        int tracer_idx = sComp+n-ncomps;
        if (tbc_descriptor_map[tracer_idx].size()) 
        {
            const Box domain = geom.Domain();
            const int* domhi = domain.hiVect();
            const int* domlo = domain.loVect();
            const Real* dx   = geom.CellSize();
            FArrayBox sdat;
            
            for (std::map<Orientation,BCDesc>::const_iterator
                     it=tbc_descriptor_map[tracer_idx].begin(); it!=tbc_descriptor_map[tracer_idx].end(); ++it) 
            {
                const Box bndBox = it->second.first;
                const Array<int>& face_bc_idxs = it->second.second;
                sdat.resize(bndBox,1); sdat.setVal(0);
                
                for (int i=0; i<face_bc_idxs.size(); ++i) {
                    const RegionData& face_tbc = tbc_array[tracer_idx][face_bc_idxs[i]];
                    face_tbc.apply(sdat,dx,0,1,time);
                }

                Box ovlp = bndBox & fab.box();
                if (ovlp.ok()) {
                    fab.copy(sdat,ovlp,0,ovlp,dComp+n,1);
                }
            }
        }    
    }
}

void
PorousMedia::dirichletPressBC (FArrayBox& fab, Real time)
{
  Array<int> bc(BL_SPACEDIM*2,0);

    if (pbc_descriptor_map.size()) 
    {
        const Box domain = geom.Domain();
        const int* domhi = domain.hiVect();
        const int* domlo = domain.loVect();
        const Real* dx   = geom.CellSize();
        FArrayBox sdat, prdat;
        for (std::map<Orientation,BCDesc>::const_iterator
                 it=pbc_descriptor_map.begin(); it!=pbc_descriptor_map.end(); ++it) 
        {
            const Box bndBox = it->second.first;
            const Array<int>& face_bc_idxs = it->second.second;
            sdat.resize(bndBox,ncomps); sdat.setVal(0);
	    prdat.resize(bndBox,1); prdat.setVal(0);

            for (int i=0; i<face_bc_idxs.size(); ++i) {
                const RegionData& face_bc = bc_array[face_bc_idxs[i]]; 

		if (model==model_list["richard"]) {
		  face_bc.apply(sdat,dx,0,ncomps,time);
		  FArrayBox cpldat, phidat, kpdat, ktdat;
		  const int n_cpl_coef = cpl_coef->nComp();
		  cpldat.resize(bndBox,n_cpl_coef);
		  phidat.resize(bndBox,1);
		  ktdat.resize(bndBox,BL_SPACEDIM);
		  kpdat.resize(bndBox,1);
		  for (int i=0; i<rocks.size(); ++i)
		    {	  
		      rocks[i].set_constant_cplval(cpldat,dx);
		      rocks[i].set_constant_kval(ktdat,dx);
		      rocks[i].set_constant_pval(phidat,dx);
		    }
		  kpdat.copy(ktdat,bndBox,it->first.coordDir(),bndBox,0,1);
		  calcCapillary(prdat,sdat,phidat,kpdat,cpldat,bndBox,bc);
		}
	    }

            Box ovlp = bndBox & fab.box();
            if (ovlp.ok()) {
                fab.copy(prdat,ovlp,0,ovlp,0,1);
            }
        }
    }    
}  

MultiFab*
PorousMedia::derive (const std::string& name,
                     Real               time,
                     int                ngrow)
{
    BL_ASSERT(ngrow >= 0);
    
    if (const DeriveRec* rec = derive_lst.get(name))
    {
        const DeriveRec* rec = derive_lst.get(name);
        BoxArray dstBA(grids);
        MultiFab* mf = new MultiFab(dstBA, rec->numDerive(), ngrow);
        int dcomp = 0;
        derive(name,time,*mf,dcomp);
        return mf;
    }
    else
    {
        //
        // If we got here, cannot derive given name.
        //
        std::string msg("PorousMedia::derive(): unknown variable: ");
        msg += name;
        BoxLib::Error(msg.c_str());
    }
}

void
PorousMedia::derive (const std::string& name,
                     Real               time,
                     MultiFab&          mf,
                     int                dcomp)
{
    const DeriveRec* rec = derive_lst.get(name);

    bool not_found_yet = false;

    if (name=="Material_ID") {
        
        BL_ASSERT(dcomp < mf.nComp());

        const int ngrow = mf.nGrow();
        
        BoxArray dstBA(mf.boxArray());
        BL_ASSERT(rec->deriveType() == dstBA[0].ixType());

        const Real* dx = geom.CellSize();

        mf.setVal(-1,dcomp,1,ngrow);

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = mf[mfi];
            for (int i=0; i<rocks.size(); ++i)
            {
                Real val = (Real)i;
                const PArray<Region>& rock_regions = rocks[i].regions;
                for (int j=0; j<rock_regions.size(); ++j) {
                    rock_regions[j].setVal(fab,val,dcomp,dx,0);
                }
            }
        }

    }
    else if (name=="Grid_ID") {
        
        BL_ASSERT(dcomp < mf.nComp());

        const int ngrow = mf.nGrow();
        
        BoxArray dstBA(mf.boxArray());
        BL_ASSERT(rec->deriveType() == dstBA[0].ixType());

        mf.setVal(-1,dcomp,1,ngrow);

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            mf[mfi].setVal(mfi.index());
        }

    }
    else if (name=="Core_ID") {
        
        BL_ASSERT(dcomp < mf.nComp());

        const int ngrow = mf.nGrow();
        
        BoxArray dstBA(mf.boxArray());
        BL_ASSERT(rec->deriveType() == dstBA[0].ixType());

        mf.setVal(-1,dcomp,1,ngrow);

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            mf[mfi].setVal(ParallelDescriptor::MyProc());
        }

    }
    else if (name=="Cell_ID") {
        
        BL_ASSERT(dcomp < mf.nComp());

        const int ngrow = mf.nGrow();
        
        BoxArray dstBA(mf.boxArray());
        BL_ASSERT(rec->deriveType() == dstBA[0].ixType());

        mf.setVal(-1,dcomp,1,ngrow);
        Layout& layout = PMAmr::GetLayout();

        Layout::IntFab ifab;
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            Box gbox = mf[mfi].box();
            ifab.resize(gbox,1);
            layout.SetNodeIds(ifab,level,mfi.index());
            const int* idat = ifab.dataPtr();
            Real* rdat = mf[mfi].dataPtr();
            int numpts = gbox.numPts();
            for (int i=0; i<numpts; ++i) {
                rdat[i] = Real(idat[i]);
            }
        }
    }
    else if (name=="Capillary_Pressure") {
        
        int ncomp = rec->numDerive();
        if (have_capillary)
        {
            const BoxArray& BA = mf.boxArray();
            BL_ASSERT(rec->deriveType() == BA[0].ixType());

            int ngrow = 1;
            MultiFab S(BA,ncomps,ngrow);
            FillPatchIterator fpi(*this,S,ngrow,time,State_Type,0,ncomps);
            for ( ; fpi.isValid(); ++fpi)
            {
                S[fpi].copy(fpi(),0,0,ncomps);
            }
            
            MultiFab tmpmf(BA,ncomp,1);
            calcCapillary(&tmpmf,S);
            MultiFab::Copy(mf,tmpmf,0,dcomp,ncomp,0);
            mf.mult(BL_ONEATM,dcomp,ncomp,0);
        }
        else if (model == model_list["steady-saturated"]) {
	  mf.setVal(0,dcomp,ncomp);
	}
	else {

	  BoxLib::Abort("PorousMedia::derive: cannot derive Capillary Pressure");
        }
    }
    else if (name=="Volumetric_Water_Content") {
        
        // Note, assumes one comp per phase
        int scomp = -1;
        for (int i=0; i<cNames.size(); ++i) {
            if (cNames[i] == "Water") {
                if (pNames[i] != "Aqueous") {
                    BoxLib::Abort("No Water in the Aqueous phase");
                }
                scomp = i;
            }
        }

        if (scomp>=0)
        {
            const BoxArray& BA = mf.boxArray();
            BL_ASSERT(rec->deriveType() == BA[0].ixType());
            int ngrow = mf.nGrow();
            BL_ASSERT(mf.nGrow()<=3); // rock_phi only has this many

            int ncomp = 1; // Just water
            BL_ASSERT(rec->numDerive()==ncomp);

            FillPatchIterator fpi(*this,mf,ngrow,time,State_Type,scomp,ncomp);
            for ( ; fpi.isValid(); ++fpi)
            {
                mf[fpi].copy(fpi(),0,dcomp,ncomp);
                mf[fpi].mult((*rock_phi)[fpi],0,dcomp,ncomp);
            }
            mf.mult(1/density[scomp],dcomp,ncomp);
        }            
        else {
            BoxLib::Abort("PorousMedia::derive: cannot derive Volumetric_Water_Content");
        }
    }
    else if (name=="Aqueous_Saturation") {

        // Sum all components in the Aqueous phase
        // FIXME: Assumes one comp per phase
        int scomp = -1;
        int naq = 0;
        for (int ip=0; ip<pNames.size(); ++ip) {
            if (pNames[ip] == "Aqueous") {
                scomp = ip;
                naq++;
            }
        }

        if (naq==1)
        {
            const BoxArray& BA = mf.boxArray();
            BL_ASSERT(rec->deriveType() == BA[0].ixType());
            int ngrow = mf.nGrow();
            BL_ASSERT(mf.nGrow()<=1); // state only has this many

            int ncomp = 1; // Just aqueous
            BL_ASSERT(rec->numDerive()==ncomp);
            FillPatchIterator fpi(*this,mf,ngrow,time,State_Type,scomp,ncomp);
            for ( ; fpi.isValid(); ++fpi)
            {
                mf[fpi].copy(fpi(),0,dcomp,ncomp);
            }
            BL_ASSERT(scomp>=0 && scomp<ncomps);
            mf.mult(1./density[scomp],dcomp,ncomp);
        }            
        else {
            BoxLib::Abort("PorousMedia::derive: no support for more than one Aqueous component");
        }
    }
    else if (name=="Aqueous_Pressure") {

        // The pressure field is the Aqueous pressure in atm
        // (assumes nphase==1,2) 
        int ncomp = 1;
        int ngrow = mf.nGrow();
        AmrLevel::derive("pressure",time,mf,dcomp);
        if (model == model_list["richard"] || model == model_list["steady-saturated"]) {
            mf.mult(BL_ONEATM,dcomp,ncomp,ngrow);
            mf.plus(BL_ONEATM,dcomp,ncomp,ngrow);
        }
        else {
  	    BoxLib::Abort(std::string("PorousMedia::derive: Aqueous_Pressure not yet implemented for " + model).c_str());
        }
    }
    else if (name=="Aqueous_Volumetric_Flux_X" || name=="Aqueous_Volumetric_Flux_Y" || name=="Aqueous_Volumetric_Flux_Z")
    {
        int dir = ( name=="Aqueous_Volumetric_Flux_X"  ?  0  :
                    name == "Aqueous_Volumetric_Flux_Y" ? 1 : 2);

        BL_ASSERT(dir < BL_SPACEDIM);
        if (model == model_list["richard"] || model == model_list["steady-saturated"]) {
	  MultiFab tmf(grids,BL_SPACEDIM,0);
	  umac_edge_to_cen(u_mac_curr,tmf); 
	  MultiFab::Copy(mf,tmf,dir,0,1,0);
        }
        else {
	  BoxLib::Abort(std::string("PorousMedia::derive: Aqueous_Volumetric_Flux not yet implemented for "+model).c_str());
        }
    }
    else if (name=="Porosity") {
        
        const BoxArray& BA = mf.boxArray();
        BL_ASSERT(rec->deriveType() == BA[0].ixType());
        int ngrow = mf.nGrow();
        int ncomp = 1; // just porosity
        BL_ASSERT(rec->numDerive()==ncomp);
        BL_ASSERT(mf.nGrow()<=3); // rock_phi only has this many
        MultiFab::Copy(mf,*rock_phi,0,dcomp,ncomp,ngrow);

    } else {

        not_found_yet = true;
    }

    if (not_found_yet) {

        for (int n=0; n<ntracers && not_found_yet; ++n) {
            std::string tname = "Aqueous_" + tNames[n] + "_Concentration";
            if (name==tname) {
                AmrLevel::derive(tNames[n],time,mf,dcomp);
                not_found_yet = false;
            }
        }
    }

    if (not_found_yet)
    {
        AmrLevel::derive(name,time,mf,dcomp);
    }
}

void
PorousMedia::manual_tags_placement (TagBoxArray&    tags,
				    Array<IntVect>& bf_lev)
{
  //
  // Tag inflow and outflow faces for refinement
  // 
  Array<Orientation> Faces;
  const BCRec& p_bc = desc_lst[Press_Type].getBC(0);
  getDirichletFaces(Faces,Press_Type,p_bc);

  if (Faces.size()>0)
    {
      for (int j =0; j<4; ++j)
	{
	  for (int i=0; i<Faces.size(); ++i)
	    {
	      const Orientation& Face = Faces[i];
	      const int oDir = Face.coordDir();
	      const Box& crse_domain = BoxLib::coarsen(geom.Domain(),bf_lev[level]);
	      const int mult = (Face.isLow() ? +1 : -1);

	      
	      // Refine entire boundary if new boxes within grid_tol
	      // from outflow
        
	      const int grid_tol = 2;
	      Box flowBox = Box(BoxLib::adjCell(crse_domain,Face,grid_tol));
	      flowBox.shift(oDir,mult*grid_tol);
	      

	      // Only refine if there are already tagged cells in the region
	      
	      bool hasTags = false;
	      for (MFIter tbi(tags); !hasTags && tbi.isValid(); ++tbi)
		if (tags[tbi].numTags(flowBox) > 0) hasTags = true;

	      ParallelDescriptor::ReduceBoolOr(hasTags);
	      
	      // hack to make sure inlet is always refined.
	      if (hasTags)
		tags.setVal(BoxArray(&flowBox,1),TagBox::SET);
	    }
	}
    }	
}

void
PorousMedia::create_umac_grown (MultiFab* u_mac, MultiFab* u_macG)
{

  // This complicated copy handles the periodic boundary condition properly.

  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::create_umac_grown1()");
  BL_ASSERT(level==0);
	    
  for (int n = 0; n < BL_SPACEDIM; ++n)
    {
      MultiFab u_ghost(u_mac[n].boxArray(),1,1);
      u_ghost.setVal(1.e40);
      u_ghost.copy(u_mac[n]);
      u_ghost.FillBoundary();
      geom.FillPeriodicBoundary(u_ghost);
      for (MFIter mfi(u_macG[n]); mfi.isValid(); ++mfi)
	{
	  u_macG[n][mfi].copy(u_ghost[mfi]);
	}
    }
}

void
PorousMedia::create_umac_grown (MultiFab* u_mac, 
				PArray<MultiFab>& u_mac_crse, 
				MultiFab* u_macG) 
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::create_umac_grown2()");

  BL_ASSERT(level>0);

  const BoxArray& fgrids = grids;
  BoxList         bl     = BoxLib::GetBndryCells(fgrids,1);

  BoxArray f_bnd_ba(bl);

  bl.clear();

  BoxArray c_bnd_ba = BoxArray(f_bnd_ba.size());

  for (int i = 0; i < f_bnd_ba.size(); ++i)
    {
      c_bnd_ba.set(i,Box(f_bnd_ba[i]).coarsen(crse_ratio));
      f_bnd_ba.set(i,Box(c_bnd_ba[i]).refine(crse_ratio));
    }

  for (int n = 0; n < BL_SPACEDIM; ++n)
    {
      //
      // crse_src & fine_src must have same parallel distribution.
      // We'll use the KnapSack distribution for the fine_src_ba.
      // Since fine_src_ba should contain more points, this'll lead
      // to a better distribution.
      //
      BoxArray crse_src_ba(c_bnd_ba);
      BoxArray fine_src_ba(f_bnd_ba);

      crse_src_ba.surroundingNodes(n);
      fine_src_ba.surroundingNodes(n);

      std::vector<long> wgts(fine_src_ba.size());

      for (unsigned int i = 0; i < wgts.size(); i++)
	{
	  wgts[i] = fine_src_ba[i].numPts();
	}
      DistributionMapping dm;
      //
      // This call doesn't invoke the MinimizeCommCosts() stuff.
      // There's very little to gain with these types of coverings
      // of trying to use SFC or anything else.
      // This also guarantees that these DMs won't be put into the
      // cache, as it's not representative of that used for more
      // usual MultiFabs.
      //
      dm.KnapSackProcessorMap(wgts,ParallelDescriptor::NProcs());

      MultiFab crse_src,  fine_src; 

      crse_src.define(crse_src_ba, 1, 0, dm, Fab_allocate);
      fine_src.define(fine_src_ba, 1, 0, dm, Fab_allocate);
	    
      crse_src.setVal(1.e200);
      fine_src.setVal(1.e200);
	
      //
      // We want to fill crse_src from lower level u_mac including u_mac's grow cells.
      // Gotta do it in steps since parallel copy only does valid region.
      //
      const MultiFab& u_macLL = u_mac_crse[n];
	  
      BoxArray edge_grids = u_macLL.boxArray();
      edge_grids.grow(1);
      
      MultiFab u_macC(edge_grids,1,0);
      
      for (MFIter mfi(u_macLL); mfi.isValid(); ++mfi)
	u_macC[mfi].copy(u_macLL[mfi]);

      crse_src.copy(u_macC);
      
      for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
	{
	  const int  nComp = 1;
	  const Box& box   = crse_src[mfi].box();
	  const int* rat   = crse_ratio.getVect();
	  FORT_PC_CF_EDGE_INTERP(box.loVect(), box.hiVect(), &nComp, rat, &n,
				 crse_src[mfi].dataPtr(),
				 ARLIM(crse_src[mfi].loVect()),
				 ARLIM(crse_src[mfi].hiVect()),
				 fine_src[mfi].dataPtr(),
				 ARLIM(fine_src[mfi].loVect()),
				 ARLIM(fine_src[mfi].hiVect()));
	}
      crse_src.clear();
      //
      // Replace pc-interpd fine data with preferred u_mac data at
      // this level u_mac valid only on surrounding faces of valid
      // region - this op will not fill grow region.
      //
      fine_src.copy(u_mac[n]);

      for (MFIter mfi(fine_src); mfi.isValid(); ++mfi)
	{
	  //
	  // Interpolate unfilled grow cells using best data from
	  // surrounding faces of valid region, and pc-interpd data
	  // on fine edges overlaying coarse edges.
	  //
	  const int  nComp = 1;
	  const Box& fbox  = fine_src[mfi.index()].box(); 
	  const int* rat   = crse_ratio.getVect();
	  FORT_EDGE_INTERP(fbox.loVect(), fbox.hiVect(), &nComp, rat, &n,
			   fine_src[mfi].dataPtr(),
			   ARLIM(fine_src[mfi].loVect()),
			   ARLIM(fine_src[mfi].hiVect()));
	  
	}

      // This complicated copy handles the periodic boundary condition properly.
      MultiFab u_ghost(u_mac[n].boxArray(),1,1);
      u_ghost.setVal(1.e40);
      u_ghost.copy(u_mac[n]);     
      u_ghost.FillBoundary();
      geom.FillPeriodicBoundary(u_ghost);
      for (MFIter mfi(u_macG[n]); mfi.isValid(); ++mfi)
	{
	  u_macG[n][mfi].copy(u_ghost[mfi]);
	}
      u_macG[n].copy(fine_src);
    }
}

void
PorousMedia::GetCrseUmac(PArray<MultiFab>& u_mac_crse,
                         Real              time          ) const
{
  BL_ASSERT(level>0);
  BL_ASSERT(u_mac_crse.size() == BL_SPACEDIM);

  const PorousMedia* pm = dynamic_cast<const PorousMedia*>(&parent->getLevel(level-1));

  Real t_old = pm->state[State_Type].prevTime();
  Real t_new = pm->state[State_Type].curTime(); 
  Real alpha = (time - t_old)/(t_new - t_old);
  const Geometry& cgeom  = parent->Geom(level-1);
  for (int i=0; i<BL_SPACEDIM; ++i)
    {
      BL_ASSERT(!u_mac_crse.defined(i));
      const BoxArray eba = BoxArray(pm->boxArray()).surroundingNodes(i);
   
      u_mac_crse.set(i,new MultiFab(eba, 1, 1));

      // This complicated copy is to ensure we copy the boundary
      // data of the coarse grid to ensure periodic boundary
      // condition is correct.
      BoxArray edge_grids = u_mac_crse[i].boxArray();
      edge_grids.grow(1);
      MultiFab u_macC(edge_grids,1,0);
      MultiFab u_macD(edge_grids,1,0);
      MultiFab u_macE(eba,1,1);
      for (MFIter mfi(u_mac_crse[i]); mfi.isValid(); ++mfi)
	{
	  u_macC[mfi].copy(pm->u_macG_prev[i][mfi]);
	  Real omalpha = 1.0 - alpha;
	  u_macC[mfi].mult(omalpha);
 
	  u_macD[mfi].copy(pm->u_macG_curr[i][mfi]);
	  u_macD[mfi].mult(alpha);
	}
      for (MFIter mfi(u_macC); mfi.isValid(); ++mfi)
	{
	  u_mac_crse[i][mfi].copy(u_macC[mfi]);
	  u_macE[mfi].copy(u_macD[mfi]);
	}
      MultiFab::Add(u_mac_crse[i],u_macE,0,0,1,1);

      //        FArrayBox UmacCrseTemp;
      //         for (MFIter mfi(u_mac_crse[i]); mfi.isValid(); ++mfi)
      //         {
      //             UmacCrseTemp.resize(mfi.validbox(),1);

      //             UmacCrseTemp.copy(pm->u_macG_prev[i][mfi]);
      //             Real omalpha = 1.0 - alpha;
      //             UmacCrseTemp.mult(omalpha);

      //             u_mac_crse[i][mfi].copy(pm->u_macG_curr[i][mfi]);
      //             u_mac_crse[i][mfi].mult(alpha);
      //             u_mac_crse[i][mfi].plus(UmacCrseTemp);
      //         }
      u_mac_crse[i].FillBoundary();
      cgeom.FillPeriodicBoundary(u_mac_crse[i],false);
    }
}

void
PorousMedia::GetCrsePressure (MultiFab& phi_crse,
                              Real      time      ) const
{
  if (level==0) return;

  const PorousMedia* pm = dynamic_cast<const PorousMedia*>(&parent->getLevel(level-1));
    
  Real t_old = pm->state[Press_Type].prevTime();
  Real t_new = pm->state[Press_Type].curTime();
  Real alpha = (time - t_old)/(t_new - t_old);
  const Geometry& cgeom  = parent->Geom(level-1);
    
  phi_crse.clear();
  phi_crse.define(pm->boxArray(), 1, 1, Fab_allocate); 

  // BUT NOTE we don't trust phi's ghost cells.
  FArrayBox PhiCrseTemp;

  if (std::fabs(time-t_new)<1.e-10 ) {
    const MultiFab& P_crse_new = pm->get_new_data(Press_Type);
    //MultiFab::Copy(phi_crse,P_crse_new,0,0,1,1);
    for (MFIter mfi(phi_crse); mfi.isValid(); ++mfi)
      phi_crse[mfi].copy(P_crse_new[mfi]);
      
  } 
  else if (std::fabs(time- t_old)<1.e-10) 
    {
      const MultiFab& P_crse_old = pm->get_old_data(Press_Type);
      //MultiFab::Copy(phi_crse,P_crse_old,0,0,1,1);
      for (MFIter mfi(phi_crse); mfi.isValid(); ++mfi)
	phi_crse[mfi].copy(P_crse_old[mfi]);
    
    } 
  else 
    {
      const MultiFab& P_crse_old = pm->get_old_data(Press_Type);
      const MultiFab& P_crse_new = pm->get_new_data(Press_Type);
      for (MFIter mfi(phi_crse); mfi.isValid(); ++mfi)
	{
	  PhiCrseTemp.resize(phi_crse[mfi].box(),1);

	  PhiCrseTemp.copy(P_crse_old[mfi]);
	  Real omalpha = 1.0 - alpha;
	  PhiCrseTemp.mult(omalpha);

	  phi_crse[mfi].copy(P_crse_new[mfi]);
	  phi_crse[mfi].mult(alpha);
	  phi_crse[mfi].plus(PhiCrseTemp);
	 
	}
    }

  phi_crse.FillBoundary();
  cgeom.FillPeriodicBoundary(phi_crse,true);
}

// ============
// IO Functions
// ============

void
PorousMedia::fill_from_plotfile (MultiFab&          mf,
                                 int                dcomp,
                                 const std::string& pltfile,
                                 const std::string& varname)
{
  const Real strt_time = ParallelDescriptor::second();

  if (pltfile.empty())
    BoxLib::Abort("fill_from_plotfile(): pltfile not specified");

  if (varname.empty())
    BoxLib::Abort("fill_from_plotfile(): varname not specified");

  if (verbose>1 && ParallelDescriptor::IOProcessor())
    std::cout << "fill_from_plotfile(): reading data from: " << pltfile << '\n';

  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  DataServices dataServices(pltfile, fileType);

  if (!dataServices.AmrDataOk())
    //
    // This calls ParallelDescriptor::EndParallel() and exit()
    //
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
  AmrData&           amrData   = dataServices.AmrDataRef();
  Array<std::string> plotnames = amrData.PlotVarNames();

  if (amrData.FinestLevel() < level)
    BoxLib::Abort("fill_from_plotfile(): not enough levels in plotfile");

  if (amrData.ProbDomain()[level] != Domain())
    BoxLib::Abort("fill_from_plotfile(): problem domains do not match");

  int idx = -1;
  for (int i = 0; i < plotnames.size(); ++i)
    if (plotnames[i] == varname) idx = i;

  if (idx == -1)
    {
      std::string msg = "fill_from_plotfile(): could not find '";
      msg += varname;
      msg += "' in the plotfile";
      BoxLib::Abort(msg.c_str());
    }

  amrData.FillVar(mf, level, varname, dcomp);
  amrData.FlushGrids(idx);

  if (verbose>1 && ParallelDescriptor::IOProcessor())
    std::cout << "fill_from_plotfile(): finished init from plotfile" << '\n';

  if (verbose>3)
  {
    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "PorousMedia::fill_from_plotfile(): lev: "
                  << level
                  << ", time: " << run_time << '\n';
  }
}

void
PorousMedia::checkPoint (const std::string& dir,
                         std::ostream&  os,
                         VisMF::How     how,
                         bool           dump_old)
{
  AmrLevel::checkPoint(dir,os,how,dump_old);

  std::string Level = BoxLib::Concatenate("Level_", level, 1);
  std::string uxfile = "/umac_x";
  std::string uyfile = "/umac_y";
  std::string FullPath = dir;
  if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
    {
      FullPath += '/';
    }
  FullPath += Level;
  uxfile = FullPath + uxfile;
  uyfile = FullPath + uyfile;
  VisMF::Write(u_mac_curr[0], uxfile);
  VisMF::Write(u_mac_curr[1], uyfile);

  std::string utxfile = "/umact_x";
  std::string utyfile = "/umact_y";
  utxfile = FullPath + utxfile;
  utyfile = FullPath + utyfile;
  VisMF::Write(u_macG_trac[0], utxfile);
  VisMF::Write(u_macG_trac[1], utyfile);

#if (BL_SPACEDIM == 3)
  std::string uzfile = "/umac_z";
  uzfile = FullPath + uzfile;
  VisMF::Write(u_mac_curr[2], uzfile);
  std::string utzfile = "/umact_z";
  utzfile = FullPath + utzfile;
  VisMF::Write(u_macG_trac[2], utzfile);
#endif 

#ifdef MG_USE_FBOXLIB
  if (model != model_list["richard"])
    {
      std::string rxfile = "/rhs_RhoD_x";
      std::string ryfile = "/rhs_RhoD_y";
      rxfile = FullPath + rxfile;
      ryfile = FullPath + ryfile;
      VisMF::Write(rhs_RhoD[0], rxfile);
      VisMF::Write(rhs_RhoD[1], ryfile);
#if (BL_SPACEDIM == 3)
      std::string rzfile = "/rhs_RhoD_z";
      rzfile = FullPath + rzfile;
      VisMF::Write(rhs_RhoD[2], rzfile);
#endif 
    }
#endif

  os << dt_eig << '\n';
}

// =================
// Utility functions
// =================

void 
PorousMedia::check_sum()
{
  // gathering some statistics of the solutions.

  Real minmax[2] = {1,1};

  MultiFab& S_new = get_new_data(State_Type);
  FArrayBox tmp,tmp2;

  for (MFIter mfi(S_new);mfi.isValid();++mfi) 
    {
      tmp.resize(mfi.validbox(),1);
      tmp2.resize(mfi.validbox(),1);
      tmp.setVal(0);
      tmp2.setVal(0);
    
      for (int kk=0; kk < ncomps; kk++)
	{
	  if (solid.compare(pNames[pType[kk]]) != 0) {
	    tmp2.copy(S_new[mfi],mfi.validbox(),kk,mfi.validbox(),0,1);
	    tmp2.mult(1.0/density[kk]);
	    tmp.plus(tmp2,mfi.validbox(),0,0,1);
	  }
	}
      minmax[0] = std::min(minmax[0],tmp.min(mfi.validbox(),0));
      minmax[1] = std::max(minmax[1],tmp.max(mfi.validbox(),0));
    }
    
  const int IOProc = ParallelDescriptor::IOProcessorNumber();

  ParallelDescriptor::ReduceRealMax(&minmax[0],2,IOProc);

  if (verbose>3 && ParallelDescriptor::IOProcessor())
    {
      std::cout << "   SUM SATURATION MAX/MIN = " 
		<< minmax[1] << ' ' << minmax[0] << '\n';
    }
}

void 
PorousMedia::check_minmax()
{
  MultiFab* rho;
  MultiFab& S_new = get_new_data(State_Type);
  
  rho = new MultiFab(grids,1,0);
  MultiFab::Copy(*rho,S_new,0,0,1,0);

  for (int kk = 1; kk<ncomps; kk++)
    {
      if (solid.compare(pNames[pType[kk]]) != 0) 
	MultiFab::Add(*rho,S_new,kk,0,1,0);
    }
 
  Array<Real> smin(ncomps,1.e20), smax(ncomps,-1.e20);

  for (int kk = 0; kk < ncomps; kk++)
    {
      for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
	{
	  smax[kk] = std::max(smax[kk],S_new[mfi].max(mfi.validbox(),kk));
	  smin[kk] = std::min(smin[kk],S_new[mfi].min(mfi.validbox(),kk));
	}
    }
  const int IOProc = ParallelDescriptor::IOProcessorNumber();

  ParallelDescriptor::ReduceRealMax(smax.dataPtr(), ncomps, IOProc);
  ParallelDescriptor::ReduceRealMin(smin.dataPtr(), ncomps, IOProc);
  
  if (verbose > 3 && ParallelDescriptor::IOProcessor())
    {
      for (int kk = 0; kk < ncomps; kk++)
	{
	  std::cout << "   SNEW MAX/MIN OF COMP " << kk
		    << ' ' << smax[kk] << "  " << smin[kk] << '\n';
	}
    }

  Real rhomaxmin[2] = {-1.e20,+1.e20};
  for (MFIter mfi(*rho); mfi.isValid(); ++mfi)
    {
      rhomaxmin[0] = std::max(rhomaxmin[0],(*rho)[mfi].max(mfi.validbox(),0));
      rhomaxmin[1] = std::min(rhomaxmin[1],(*rho)[mfi].min(mfi.validbox(),0));
    }

  ParallelDescriptor::ReduceRealMax(&rhomaxmin[0], 1, IOProc);
  ParallelDescriptor::ReduceRealMin(&rhomaxmin[1], 1, IOProc);

  if (verbose > 3 && ParallelDescriptor::IOProcessor())
    {  
      std::cout << "   RHO MAX/MIN "
		<< ' ' << rhomaxmin[0] << "  " << rhomaxmin[1] << '\n';
    }

  delete rho;
}

void 
PorousMedia::check_minmax(int fscalar, int lscalar)
{
  MultiFab& S_new = get_new_data(State_Type);
  
  const int nscal = lscalar - fscalar + 1;

  Array<Real> smin(nscal,1.e20), smax(nscal,-1.e20);

  for (int kk = 0; kk < nscal; kk++)
    {
      for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
	{
            smax[kk] = std::max(smax[kk], S_new[mfi].max(mfi.validbox(),fscalar+kk));
            smin[kk] = std::min(smin[kk], S_new[mfi].min(mfi.validbox(),fscalar+kk));
	}
    }
  const int IOProc = ParallelDescriptor::IOProcessorNumber();
  ParallelDescriptor::ReduceRealMax(smax.dataPtr(), nscal, IOProc);
  ParallelDescriptor::ReduceRealMin(smin.dataPtr(), nscal, IOProc);
  
  if (verbose>3 && ParallelDescriptor::IOProcessor())
    {
        for (int kk = 0; kk < nscal; kk++)
	{
	  std::cout << "   SNEW MAX/MIN OF COMP "
                    << fscalar+kk
		    << ' ' << smax[kk] 
		    << ' ' << smin[kk] << '\n';
	}
    }
}

void 
PorousMedia::check_minmax(MultiFab& mf)
{
  const int ncomp = mf.nComp();
  Array<Real> smin(ncomp,1.e20), smax(ncomp,-1.e20);

  for (int kk = 0; kk < ncomp; kk++)
    {
      for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	  smax[kk] = std::max(smax[kk],mf[mfi].max(mfi.validbox(),kk));
	  smin[kk] = std::min(smin[kk],mf[mfi].min(mfi.validbox(),kk));
	}
    }
  const int IOProc = ParallelDescriptor::IOProcessorNumber();

  ParallelDescriptor::ReduceRealMax(smax.dataPtr(), ncomp, IOProc);
  ParallelDescriptor::ReduceRealMin(smin.dataPtr(), ncomp, IOProc);
  
  if (verbose>3 && ParallelDescriptor::IOProcessor())
    {
      for (int kk = 0; kk < ncomp; kk++)
	{
	  std::cout << " MAX/MIN OF MF " << kk
		    << ' ' << smax[kk] << "  " << smin[kk] << '\n';
	}
    }
}

void 
PorousMedia::check_minmax(MultiFab* u_mac)
{
  //
  // Write out the min and max of the MAC velocities.
  //
  Real umax[BL_SPACEDIM] = {D_DECL(-1.e20,-1.e20,-1.e20)};
  Real umin[BL_SPACEDIM] = {D_DECL(+1.e20,+1.e20,+1.e20)};

  for (MFIter mfi(u_mac[0]); mfi.isValid(); ++mfi)
    {
      const int i = mfi.index();

      umax[0] = std::max(umax[0],u_mac[0][i].max(u_mac[0].boxArray()[i]));
      umin[0] = std::min(umin[0],u_mac[0][i].min(u_mac[0].boxArray()[i]));
      umax[1] = std::max(umax[1],u_mac[1][i].max(u_mac[1].boxArray()[i]));
      umin[1] = std::min(umin[1],u_mac[1][i].min(u_mac[1].boxArray()[i]));
#if(BL_SPACEDIM == 3)
      umax[2] = std::max(umax[2],u_mac[2][i].max(u_mac[2].boxArray()[i]));
      umin[2] = std::min(umin[2],u_mac[2][i].min(u_mac[2].boxArray()[i]));
#endif
    }

  const int IOProc = ParallelDescriptor::IOProcessorNumber();

  ParallelDescriptor::ReduceRealMax(&umax[0], BL_SPACEDIM, IOProc);
  ParallelDescriptor::ReduceRealMin(&umin[0], BL_SPACEDIM, IOProc);

  if (verbose>3 && ParallelDescriptor::IOProcessor())
  {
      D_TERM(std::cout << "   UMAC MAX/MIN  " << umax[0] << "  " << umin[0] << '\n';,
             std::cout << "   VMAC MAX/MIN  " << umax[1] << "  " << umin[1] << '\n';,
             std::cout << "   WMAC MAX/MIN  " << umax[2] << "  " << umin[2] << '\n';);
  }
}

void
PorousMedia::umac_edge_to_cen(MultiFab* u_mac, MultiFab& U_cc)
{
  // average velocity onto cell center
  for (MFIter mfi(U_cc); mfi.isValid(); ++mfi)
    {
      const int* lo     = mfi.validbox().loVect();
      const int* hi     = mfi.validbox().hiVect();
    
      const int* u_lo   = U_cc[mfi].loVect();
      const int* u_hi   = U_cc[mfi].hiVect();
      const Real* udat  = U_cc[mfi].dataPtr();
	  
      const int* um_lo  = (u_mac[0])[mfi].loVect();
      const int* um_hi  = (u_mac[0])[mfi].hiVect();
      const Real* umdat = (u_mac[0])[mfi].dataPtr();
	
      const int* vm_lo  = (u_mac[1])[mfi].loVect();
      const int* vm_hi  = (u_mac[1])[mfi].hiVect();
      const Real* vmdat = (u_mac[1])[mfi].dataPtr();
	
#if (BL_SPACEDIM == 3)
      const int* wm_lo  = (u_mac[2])[mfi].loVect();
      const int* wm_hi  = (u_mac[2])[mfi].hiVect();
      const Real* wmdat = (u_mac[2])[mfi].dataPtr();
#endif

      FORT_AVG_UMAC(umdat,ARLIM(um_lo),ARLIM(um_hi),
		    vmdat,ARLIM(vm_lo),ARLIM(vm_hi),
#if (BL_SPACEDIM == 3)
		    wmdat,ARLIM(wm_lo),ARLIM(wm_hi),
#endif
		    udat ,ARLIM( u_lo),ARLIM( u_hi),lo,hi);       
    }
}

void
PorousMedia::umac_cpy_edge_to_cen(MultiFab* u_mac, int idx_type, int ishift)
{
  // average velocity onto cell center
  MultiFab&  U_cor  = get_new_data(idx_type);
  for (MFIter mfi(U_cor); mfi.isValid(); ++mfi)
    {
      const int* lo     = mfi.validbox().loVect();
      const int* hi     = mfi.validbox().hiVect();
    
      const int* u_lo   = U_cor[mfi].loVect();
      const int* u_hi   = U_cor[mfi].hiVect();
      const Real* udat  = U_cor[mfi].dataPtr();
	  
      const int* um_lo  = (u_mac[0])[mfi].loVect();
      const int* um_hi  = (u_mac[0])[mfi].hiVect();
      const Real* umdat = (u_mac[0])[mfi].dataPtr();
	
      const int* vm_lo  = (u_mac[1])[mfi].loVect();
      const int* vm_hi  = (u_mac[1])[mfi].hiVect();
      const Real* vmdat = (u_mac[1])[mfi].dataPtr();
	
#if (BL_SPACEDIM == 3)
      const int* wm_lo  = (u_mac[2])[mfi].loVect();
      const int* wm_hi  = (u_mac[2])[mfi].hiVect();
      const Real* wmdat = (u_mac[2])[mfi].dataPtr();
#endif

      FORT_CPY_UMAC(umdat,ARLIM(um_lo),ARLIM(um_hi),
		    vmdat,ARLIM(vm_lo),ARLIM(vm_hi),
#if (BL_SPACEDIM == 3)
		    wmdat,ARLIM(wm_lo),ARLIM(wm_hi),
#endif
		    udat ,ARLIM( u_lo),ARLIM( u_hi),lo,hi, &ishift); 
    }
}

void
PorousMedia::compute_divu (MultiFab& soln,
			   MultiFab* umac)
{
  //
  // This compute the divergence of umac
  //

  const Real* dx   = geom.CellSize();

  for (MFIter fpi(soln); fpi.isValid(); ++fpi)
    {
      const int i = fpi.index();
      const int* lo = fpi.validbox().loVect();
      const int* hi = fpi.validbox().hiVect();

      const Real* sdat = soln[i].dataPtr();
      const int* s_lo  = soln[i].loVect();
      const int* s_hi  = soln[i].hiVect();
    
      const Real* uxdat = umac[0][i].dataPtr();
      const int*  uxlo  = umac[0][i].loVect();
      const int*  uxhi  = umac[0][i].hiVect();

      const Real* uydat = umac[1][i].dataPtr();
      const int*  uylo  = umac[1][i].loVect();
      const int*  uyhi  = umac[1][i].hiVect();

#if (BL_SPACEDIM == 3)
      const Real* uzdat = umac[2][i].dataPtr();
      const int*  uzlo  = umac[2][i].loVect();
      const int*  uzhi  = umac[2][i].hiVect();
#endif

      FORT_DIV_UMAC (sdat, ARLIM(s_lo),ARLIM(s_hi),
		     uxdat,ARLIM(uxlo),ARLIM(uxhi),
		     uydat,ARLIM(uylo),ARLIM(uyhi),
#if (BL_SPACEDIM == 3)
		     uzdat,ARLIM(uzlo),ARLIM(uzhi),
#endif
		     lo,hi,dx);
    }
}

