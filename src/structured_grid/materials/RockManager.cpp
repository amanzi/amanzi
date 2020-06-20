
#include <fstream>
#include <iomanip>

#include <RockManager.H>
#include <RockManager_F.H>

#include <Utility.H>
#include <ParmParse.H>

// Note that the permeability may be specified in m^2 or in mDa
static bool Permeability_in_m2 = true;
static bool Permeability_in_mDa = false;

// Note that the capillary pressure model alpha can be in Pa^-1 or atm^-1
static bool Capillary_Pressure_alpha_in_invPa = true;
static bool Capillary_Pressure_alpha_in_invAtm = false;

static std::string CapillaryPressureName    = "capillary_pressure";
static std::string PorosityName             = "porosity";
static std::string PermeabilityName         = "permeability";
static std::string RelativePermeabilityName = "relative_permeability";
static std::string HydCondName              = "hydraulic_conductivity";

static std::string CP_model_None = "None";
static std::string CP_model_vG = "VanGenuchten";
static std::string CP_model_BC = "BrooksCorey";

static std::string Kr_model_None = "None";
static std::string Kr_model_Mualem = "Mualem";
static std::string Kr_model_Burdine = "Burdine";

static std::string Kr_model_vG_Mualem  = CP_model_vG + "_" + Kr_model_Mualem;
static std::string Kr_model_vG_Burdine = CP_model_vG + "_" + Kr_model_Burdine;
static std::string Kr_model_BC_Mualem  = CP_model_BC + "_" + Kr_model_Mualem;
static std::string Kr_model_BC_Burdine = CP_model_BC + "_" + Kr_model_Burdine;

static int MAX_CPL_PARAMS = 7; // Must be set to accommodate the model with the most parameters
static int CPL_MODEL_ID = 0;

static int VG_M                   = 1;
static int VG_ALPHA               = 2;
static int VG_SR                  = 3;
static int VG_ELL                 = 4;
static int VG_KR_MODEL_ID         = 5;
static int VG_KR_SMOOTHING_MAX_PC = 6;

static int BC_LAMBDA              = 1;
static int BC_ALPHA               = 2;
static int BC_SR                  = 3;
static int BC_ELL                 = 4;
static int BC_KR_MODEL_ID         = 5;
static int BC_KR_SMOOTHING_MAX_PC = 6;

static Real Kr_ell_vG_Mualem_DEF = 0.5;
static Real Kr_ell_vG_Burdine_DEF = 2.0;

static Real Kr_smoothing_max_pcap_DEF = -1;
static Real Kr_smoothing_min_seff_DEF = 2;

// Interpolators
static int NUM_INIT_INTERP_EVAL_PTS_DEF = 101;
static int NUM_ADD_INTERP_EVAL_PTS_DEF = 0;
static Real ERR_FOR_ADAPTIVE_INTERPS_DEF = 1.e-10;
static Real pc_at_Sr = 1.e11;

static int Rock_Mgr_ID_ctr=0;
static std::vector<RockManager*> Rock_Mgr_Ptrs;
static std::vector<std::pair<bool,Real> > Kr_smoothing_min_seff; // Bool says whether value needs to be updated

static int max_grid_size_fine_gen_DEF = 96; // Blocking size for generating GSLib datafiles
static int ngrow_fine_gen_DEF = 9; // nGrow for generating GSLib datafiles, note really a fn of correl search radii

RockManager::RockManager(const RegionManager*      _region_manager,
                         const Array<std::string>* _solute_names,
			 Real                      _liquid_density,
			 Real                      _liquid_viscosity,
			 Real                      _gravity)
  : region_manager(_region_manager),
    liquid_density(_liquid_density),
    liquid_viscosity(_liquid_viscosity),
    gravity(_gravity),
    finalized(false)
{
  BL_PROFILE("RockManager::RockManager()");

  Initialize(_solute_names);

  rock_mgr_ID = Rock_Mgr_ID_ctr++;
  Rock_Mgr_Ptrs.resize(rock_mgr_ID+1); Rock_Mgr_Ptrs[rock_mgr_ID] = this;
}

extern "C" {
  void ROCK_MANAGER_PCAP(int* rockMgrID, const Real* saturation, int* matID, Real* time, Real* capillaryPressure, int* npts) {
    Rock_Mgr_Ptrs[*rockMgrID]->CapillaryPressure(saturation,matID,*time,capillaryPressure,*npts);
  }

  void ROCK_MANAGER_INVPCAP(int* rockMgrID, const Real* capillaryPressure, int* matID, Real* time, Real* saturation, int* npts) {
    Rock_Mgr_Ptrs[*rockMgrID]->InverseCapillaryPressure(capillaryPressure,matID,*time,saturation,*npts);
  }

  void ROCK_MANAGER_RELPERM(int* rockMgrID, const Real* saturation, int* matID, Real* time, Real* kappa, int* npts) {
    Rock_Mgr_Ptrs[*rockMgrID]->RelativePermeability(saturation,matID,*time,kappa,*npts);
  }

  void ROCK_MANAGER_DSDPCAP(int* rockMgrID, const Real* saturation, int* matID, Real* time, Real* dsdpc, int* npts) {
    Rock_Mgr_Ptrs[*rockMgrID]->DInverseCapillaryPressure(saturation,matID,*time,dsdpc,*npts);
  }

  void ROCK_MANAGER_RESIDSAT(int* rockMgrID, int* matID, Real* time, Real* Sr, int* npts) {
    Rock_Mgr_Ptrs[*rockMgrID]->ResidualSaturation(matID,*time,Sr,*npts);
  }
}

static
void
EnsureFolderExists(const std::string& full_path)
{
  // Find folder name first, and ensure folder exists
  // FIXME: Will fail on Windows
  const std::vector<std::string>& tokens = BoxLib::Tokenize(full_path,"/");
  std::string dir = (full_path[0] == '/' ? "/" : "");
  for (int i=0; i<tokens.size()-1; ++i) {
    dir += tokens[i];
    if (i<tokens.size()-2) dir += "/";
  }
  
  if(!BoxLib::FileExists(dir)) {
    if ( ! BoxLib::UtilCreateDirectory(dir, 0755)) {
      BoxLib::CreateDirectoryFailed(dir);
    }
  }
}

int
RockManager::NComp(const std::string& property_name) const
{
  BL_PROFILE("RockManager::NComp()");

  int nc = 0;
  if (materialFiller->CanDerive(property_name)) {
    nc = materialFiller->NComp(property_name);
  }
  return nc;
}

int
RockManager::FindMaterialNum(const std::string& name) const
{
  BL_PROFILE("RockManager::FindMaterial()");

  bool found=false;
  int iMat = -1;
  const PArray<Material>& mats = GetMaterials();
  for (int i=0; i<mats.size() && !found; ++i) {
    const Material& mat = mats[i];
    if (name == mat.Name()) {
      found = true;
      iMat = i;
    }
  } 
  if (iMat < 0) {
    std::string m = "Named material not found " + name;
    BoxLib::Abort(m.c_str());
  }
  return iMat;
}

const Material&
RockManager::FindMaterial(const std::string& name) const
{
  return GetMaterials()[FindMaterialNum(name)];
}

const MonotCubicInterpolator&
RockManager::CPInterpolator(const std::string& name) const
{
  return CP_s_interps[FindMaterialNum(name)];
}

const MonotCubicInterpolator&
RockManager::KrInterpolator(const std::string& name) const
{
  return Kr_s_interps[FindMaterialNum(name)];
}

const Array<std::string>
RockManager::MaterialNames() const
{
  const PArray<Material>& mats = GetMaterials();
  Array<std::string> names(mats.size());
  for (int i=0; i<mats.size(); ++i) {
    names[i] = mats[i].Name();
  }
  return names;
}

const Material&
RockManager::FindMaterialInRegions(const Array<std::string>& region_names) const
{
  BL_PROFILE("RockManager::FindMaterialInRegions()");

  int iMat = -1;
  const PArray<Material>& mats = GetMaterials();
  for (int i=0; i<mats.size() && iMat<0; ++i) {
    const Material& mat = mats[i];
    const Array<const Region*>& thisRegions = mat.Regions();

    bool contains_one = false;
    bool contains_all = true;
    for (int j=0; j<thisRegions.size() && iMat<0; ++j) {
      for (int k=0; k<region_names.size(); ++k) {
        if (region_names[k] == thisRegions[j]->name) {
          contains_one = true;
        }
        else if (contains_one) {
          contains_all = false;
        }
      }
      if (contains_one) {
        if (contains_all) {
          iMat = i;
        }
        else {
          BoxLib::Abort("Single Material does contain all regions for this boundary");
        }
      }
    }
  }
  if (iMat < 0) {
    BoxLib::Abort("No single material found containing all specified regions");
  }
  return mats[iMat];
}

void
RockManager::BuildInterpolators()
{
  BL_PROFILE("RockManager::BuildInterpolators()");

  CP_s_interps.resize(rock.size(),PArrayManage);
  Kr_s_interps.resize(rock.size(),PArrayManage);

  int nComp = materialFiller->NComp(CapillaryPressureName);
  static IntVect iv(D_DECL(0,0,0));
  static Box bx(iv,iv);
  FArrayBox pc_params(bx,nComp);
  int level=0; //not really used
  int dComp=0;
  Real time = 0;

  for (int n=0; n<rock.size(); ++n) {

    Array<Real> s(NUM_INIT_INTERP_EVAL_PTS_DEF);
    int Npts = s.size();
    Array<Real> pc(Npts);
    Array<int> mat(Npts,n);

    pc[0] = pc_at_Sr;
    InverseCapillaryPressure(pc.dataPtr(),&n,time,&(s[0]),1);

    Real ds = 1 - s[0];
    for (int i=1; i<s.size(); ++i) {
      s[i] = s[0] + ds*Real(i)/(NUM_INIT_INTERP_EVAL_PTS_DEF - 1);
    }
    CapillaryPressure(s.dataPtr(),mat.dataPtr(),time,pc.dataPtr(),Npts);
    CP_s_interps.set(n, new MonotCubicInterpolator(std::vector<Real>(s),std::vector<Real>(pc)));

    Array<Real> kr(s.size());
    RelativePermeability(s.dataPtr(),mat.dataPtr(),time,kr.dataPtr(),Npts);

    int NaddMAX = NUM_ADD_INTERP_EVAL_PTS_DEF;
    Real errForAdapt = ERR_FOR_ADAPTIVE_INTERPS_DEF;
    if (NaddMAX > 0) {
      MonotCubicInterpolator* t = new MonotCubicInterpolator(std::vector<Real>(s),std::vector<Real>(kr));
      Real diff;
      for (int i=0; i<NaddMAX; ++i) {
        std::pair<double,double> newval = t->getMissingX();
	Real interpval = (*t)(newval.first);
	RelativePermeability(&newval.first,mat.dataPtr(),time,&newval.second,1);
	diff = std::abs(newval.second-interpval);
	if (diff < errForAdapt) {
	  i=NaddMAX;
        }
	t->addPair(newval.first,newval.second);
      }
      //std::cout << n << " last err " << diff << std::endl;
      Kr_s_interps.set(n, t);
    } else {
      Kr_s_interps.set(n, new MonotCubicInterpolator(std::vector<Real>(s),std::vector<Real>(kr)));
    }
  }

  for (int n=0; n<rock.size(); ++n) {
    int Npts = WRM_plot_file[n].first;
    if (ParallelDescriptor::IOProcessor() && Npts > 0) {
      Array<Real> s(Npts);
      Array<Real> pc(Npts);
      Array<int> mat(Npts,n);

      pc[0] = pc_at_Sr;
      InverseCapillaryPressure(pc.dataPtr(),&n,time,&(s[0]),1);

      Real ds = 1 - s[0];
      for (int i=1; i<s.size(); ++i) {
        s[i] = std::max(s[0], std::min(1.0, s[0] + ds*Real(i)/(Npts - 1)));
      }
      CapillaryPressure(s.dataPtr(),mat.dataPtr(),time,pc.dataPtr(),Npts);

      Array<Real> kr(s.size());
      RelativePermeability(s.dataPtr(),mat.dataPtr(),time,kr.dataPtr(),Npts);

      const std::string& file = WRM_plot_file[n].second;
      std::cout << "Writing WRM data for material \"" << rock[n].Name()
                << "\" to file \"" << file << "\"" << std::endl;

      EnsureFolderExists(file);
      std::ofstream osf; osf.open(file.c_str());
      osf << std::setprecision(15);
      for (int i=1; i<s.size(); ++i) {
        osf << s[i] << " " << pc[i] << " " << kr[i] << std::endl;
      }
      osf.close();
    }
  }
}


// FIXME: Replace with call to MatFiller::FillCellsOutsideDomain
void
RockManager::FillBoundary(Real      time,
                          int       level,
                          MultiFab& mf,
                          int       dComp,
                          int       nComp,
                          int       nGrow) const
{
  BL_PROFILE("RockManager::FillBoundary()");

  const Geometry& geom = materialFiller->Geom(level);

  if (nGrow>0) {
    BL_ASSERT(mf.nGrow() >= nGrow);
    const Box& domain = geom.Domain();
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      FArrayBox& fab = mf[mfi];
      const Box& vbox = mfi.validbox();
      Box gbox = Box(vbox).grow(nGrow);
      if ( !(domain.contains(gbox)) ) {
        for (int d=0; d<BL_SPACEDIM; ++d) {
          Box adjCellLo = BoxLib::adjCellLo(vbox,d,1);
          Box intCellLo = Box(adjCellLo).shift(d,1);
          for (int i=0; i<nGrow; ++i) {
            Box ladjCellLo = Box(adjCellLo).shift(d,-i);
            for (int dd=0; dd<d; ++dd) {
              ladjCellLo.grow(dd,nGrow);
            }
            Box lintCellLo = Box(ladjCellLo).shift(d,i+1);
            fab.copy(fab,lintCellLo,dComp,ladjCellLo,dComp,nComp);
          }

          Box adjCellHi = BoxLib::adjCellHi(vbox,d,1);
          Box intCellHi = Box(adjCellHi).shift(d,-1);
          for (int i=0; i<nGrow; ++i) {
            Box ladjCellHi = Box(adjCellHi).shift(d,i);
            for (int dd=0; dd<d; ++dd) {
              ladjCellHi.grow(dd,nGrow);
            }
            Box lintCellHi = Box(ladjCellHi).shift(d,-i-1);
            fab.copy(fab,lintCellHi,dComp,ladjCellHi,dComp,nComp);
          }
        }
      }
    }
    bool local = false;
    bool corner=true;
    mf.FillBoundary(dComp,nComp,!corner);
    materialFiller->Geom(level).FillPeriodicBoundary(mf,dComp,nComp,corner,local);
  }
}

void
RockManager::FinalizeBuild(const Array<Geometry>& geomArray,
                           const Array<IntVect>&  refRatio,
                           int                    nGrow,
			   bool                   restart)
{
  BL_PROFILE("RockManager::FinalizeBuild()");

  for (int i=0; i<rock.size(); ++i) {
    Array<std::string> propNames = rock[i].PropertyNames();
    for (int j=0; j<propNames.size(); ++j) {      
      const Property* p = rock[i].Prop(propNames[j]);
      GSLibProperty* t = dynamic_cast<GSLibProperty*>(const_cast<Property*>(p));
      if (t!=0) {
	if (!restart && ParallelDescriptor::IOProcessor()) {
	  std::cout << "WARNING: Building GSLib file with ngrow_fine_gen: " << ngrow_fine_gen << std::endl;
	  std::cout << "   It is up to you to ensure that this is consistent with the search radii!" << std::endl;
	}
        t->BuildDataFile(geomArray,refRatio,ngrow_fine_gen,max_grid_size_fine_gen,p->coarsenRule(),propNames[j],restart);
      }
    }
  }
  materialFiller = new MatFiller(geomArray,refRatio,rock);
  finalized = true;
  BuildInterpolators();
}

void
RockManager::Initialize(const Array<std::string>* solute_names)
{
  BL_PROFILE("RockManager::Initialize()");

  if (solute_names != 0) {
    int num_solutes = solute_names->size();
    known_solutes.resize(num_solutes);
    for (int i=0; i<num_solutes; ++i) {
      known_solutes[i] = (*solute_names)[i];
    }
  }

  is_saturated = false;

  static int CP_cnt = 0;
  CP_models[CP_model_None] = CP_cnt++;
  CP_models[CP_model_vG] = CP_cnt++;
  CP_models[CP_model_BC] = CP_cnt++;

  static int Kr_cnt = 0;
  Kr_models[Kr_model_None] = Kr_cnt++;
  Kr_models[Kr_model_vG_Mualem] = Kr_cnt++;
  Kr_models[Kr_model_vG_Burdine] = Kr_cnt++;
  Kr_models[Kr_model_BC_Mualem] = Kr_cnt++;
  Kr_models[Kr_model_BC_Burdine] = Kr_cnt++;

  ParmParse pp("rock");
  int nrock = pp.countval("rocks");
  if (nrock <= 0) {
    BoxLib::Abort("At least one rock type must be defined.");
  }
  Array<std::string> r_names;  pp.getarr("rocks",r_names,0,nrock);

  ParmParse ppp("prob");
  max_grid_size_fine_gen = max_grid_size_fine_gen_DEF;
  ppp.query("max_grid_size_fine_gen",max_grid_size_fine_gen);
  ngrow_fine_gen = ngrow_fine_gen_DEF;
  ppp.query("ngrow_fine_gen",ngrow_fine_gen);

  rock.clear();
  rock.resize(nrock,PArrayManage);
  Array<std::string> material_regions;

  // Scan rock for properties that must be defined for all
  //   if defined for one. 
  bool user_specified_dispersivity = false;
  bool user_specified_tortuosity = false;
  bool user_specified_specific_storage = false;
  for (int i = 0; i<nrock; i++) {
    const std::string& rname = r_names[i];
    const std::string prefix("rock." + rname);
    ParmParse ppr(prefix.c_str());
    user_specified_tortuosity = ppr.countval("tortuosity.val");
    user_specified_dispersivity = ppr.countval("dispersivity.alphaL");
    user_specified_specific_storage = ppr.countval("specific_storage.val");
  }

  if (is_saturated) {
    user_specified_specific_storage = true; // Will use default if not specified
  }

  do_diffusion = ( user_specified_tortuosity || user_specified_dispersivity);
  
  do_tensor_diffusion = do_diffusion && user_specified_dispersivity;

  // setup static database for smoothing interval
  Kr_smoothing_min_seff.resize(nrock,std::pair<bool,Real>(true,Kr_smoothing_min_seff_DEF));

  // set up static database for WRM plot files
  WRM_plot_file.resize(nrock,std::pair<int,std::string>(0,""));

  for (int i = 0; i<nrock; i++) {

    const std::string& rname = r_names[i];
    const std::string prefix("rock." + rname);
    ParmParse ppr(prefix.c_str());

    bool generate_porosity_gslib_file = false;
    bool generate_perm_gslib_file = false;
        
    static Property::CoarsenRule arith_crsn = Property::Arithmetic;
    static Property::CoarsenRule harm_crsn = Property::ComponentHarmonic;
    static Property::RefineRule pc_refine = Property::PiecewiseConstant;

    Real rdensity = -1; // ppr.get("density",rdensity); // not actually used anywhere

    Property* Dmolec_func = 0;

    Array<Real> rDispersivity(2,0);
    Property* Dispersivity_func = 0;
    if (user_specified_dispersivity) {
      ppr.query("dispersivity.alphaL",rDispersivity[0]);
      ppr.query("dispersivity.alphaT",rDispersivity[1]);
      std::string Dispersivity_str = "dispersivity";
      Dispersivity_func = new ConstantProperty(Dispersivity_str,rDispersivity,harm_crsn,pc_refine);
    }

    Array<Real> rTortuosity(BL_SPACEDIM,1);
    Property* Tortuosity_func = 0;
    if (user_specified_tortuosity) {
      ppr.query("tortuosity.val",rTortuosity[0]);
      for (int d=1; d<BL_SPACEDIM; ++d) {
        rTortuosity[d] = rTortuosity[0];
      }
      std::string Tortuosity_str = "tortuosity";
      Tortuosity_func = new ConstantProperty(Tortuosity_str,rTortuosity,harm_crsn,pc_refine);
    }

    Real rSpecificStorage = 0;
    Property* SpecificStorage_func = 0;
    if (user_specified_specific_storage) {
      ppr.query("specific_storage.val",rSpecificStorage);
      std::string SpecificStorage_str = "specific_storage";
      SpecificStorage_func = new ConstantProperty(SpecificStorage_str,rSpecificStorage,arith_crsn,pc_refine);
    }

    Property* phi_func = 0;
    Array<Real> rpvals(1), rptimes;
    Array<std::string> rpforms;
    const std::string porosity_prefix(prefix+"."+PorosityName);
    ParmParse pprp(porosity_prefix.c_str());


    std::string PorosityValsName = "vals";
    std::string PorosityValName = "val";
    std::string PorosityStdName = "std";
    std::string PorosityTimesName = "times";
    std::string PorosityFormsName = "forms";
    std::string PorosityDistName = "distribution_type";
    std::string PorosityGSParamFileName = "gslib_param_file";
    std::string PorosityGSDataFileName = "gslib_data_file";
    std::string PorosityGSFileShiftName = "gslib_file_shift";

    std::string phi_dist = "uniform"; pprp.query(PorosityDistName.c_str(),phi_dist);

    if (phi_dist != "uniform" && phi_dist!="gslib") {
      BoxLib::Abort(std::string("Unrecognized distribution_type for rock: \""+rname).c_str());
    }
    if (phi_dist == "uniform") {
      if (pprp.countval(PorosityValsName.c_str()) > 0) {
        pprp.getarr(PorosityValsName.c_str(),rpvals,0,pprp.countval(PorosityValsName.c_str()));
        int nrpvals = rpvals.size();
        if (nrpvals>1) {
          pprp.getarr(PorosityTimesName.c_str(),rptimes,0,nrpvals);
          pprp.getarr(PorosityFormsName.c_str(),rpforms,0,nrpvals-1);
          TabularFunction pft(rptimes,rpvals,rpforms);
          phi_func = new TabularInTimeProperty(PorosityName.c_str(),pft,arith_crsn,pc_refine);
        }
        else {
          phi_func = new ConstantProperty(PorosityName,rpvals[0],arith_crsn,pc_refine);
        }
      } else if (pprp.countval(PorosityValName.c_str()) == 1) {
        rpvals.resize(1); pprp.get(PorosityValName.c_str(),rpvals[0]);
        phi_func = new ConstantProperty(PorosityName,rpvals[0],arith_crsn,pc_refine);
      } else {
        BoxLib::Abort("No valid porosity values given for uniform distribution");
      }
    }
    else {
      // phi_dist == gslib
      std::string gslib_param_file, gslib_data_file;
      generate_porosity_gslib_file = (pprp.countval(PorosityGSParamFileName.c_str()) != 0);
      if (pprp.countval(PorosityGSDataFileName.c_str()) == 0) {
        pprp.get(PorosityGSParamFileName.c_str(),gslib_param_file);
        gslib_data_file="porosity.gslib";
      }
      else {
        pprp.query(PorosityGSParamFileName.c_str(),gslib_param_file);
        pprp.get(PorosityGSDataFileName.c_str(),gslib_data_file);
      }

      Array<Real> gslib_file_shift(BL_SPACEDIM,0);
      pprp.queryarr(PorosityGSFileShiftName.c_str(),gslib_file_shift,0,BL_SPACEDIM);

      Real avg; pprp.get(PorosityValName.c_str(),avg);

      phi_func = new GSLibProperty(PorosityName,avg,gslib_param_file,gslib_data_file,gslib_file_shift,arith_crsn,pc_refine);
    }

    Property* kappa_func = 0;
    std::string KValName = "val";
    std::string KValsName = "vals";
    std::string KStdName = "std";
    std::string KDistName = "distribution_type";
    std::string KGSParamFileName = "gslib_param_file";
    std::string KGSDataFileName = "gslib_data_file";
    std::string KGSFileShiftName = "gslib_file_shift";

    const std::string perm_prefix(prefix+"."+PermeabilityName);
    const std::string HC_prefix(prefix+"."+HydCondName);
    ParmParse pprk(perm_prefix.c_str());
    ParmParse pprhc(HC_prefix.c_str());

    std::string VertValName = "vertical."+KValsName;
    bool hasPerm =  pprk.countval(VertValName.c_str()) > 0;
    bool hasHC   = pprhc.countval(VertValName.c_str()) > 0;
    if (!(hasPerm ^ hasHC)) {
      BoxLib::Abort(std::string("Must specify either permeability or hydraulic conductivy for rock: \""+rname).c_str());
    }

    std::string KName = hasPerm ? PermeabilityName : HydCondName;
    ParmParse& pprK = hasPerm ? pprk : pprhc;

    if (hasHC) {
      if (liquid_density<=0 || liquid_viscosity<=0 || gravity <= 0) {
	BoxLib::Abort("hydraulic_conductivity requires non-negative values for gravity and for the density, viscosity of liquid.");
      }
    }

    std::string K_dist = "uniform"; pprK.query(KDistName.c_str(),K_dist);

    if (K_dist != "uniform" && K_dist!="gslib") {
      BoxLib::Abort(std::string("Unrecognized perm/HC distribution_type for rock: \""+rname).c_str());
    }
    if (K_dist == "uniform") {

      Array<Real> rvpvals(1), rhpvals(1), rh1pvals(1), rvptimes(1), rhptimes(1), rh1ptimes(1);
      Array<std::string> rvpforms, rhpforms, rh1pforms;

      Array<Real> rperm_in(BL_SPACEDIM);
      if (ppr.countval(PermeabilityName.c_str())) {
	ppr.getarr(PermeabilityName.c_str(),rperm_in,0,BL_SPACEDIM);
	rhpvals[0] = rperm_in[0];
	rvpvals[0] = rperm_in[BL_SPACEDIM-1];
#if BL_SPACEDIM==3
	rh1pvals[0] = rperm_in[0];
#endif
      }
      else {

	std::string KVertValName = KName+".vertical.vals";
	std::string KVertTimesName = KName+".vertical.times";
	std::string KVertFormsName = KName+".vertical.forms";

	int nrvpvals = ppr.countval(KVertValName.c_str());
	if (nrvpvals>0) {
	  ppr.getarr(KVertValName.c_str(),rvpvals,0,nrvpvals);
	  if (nrvpvals>1) {
	    ppr.getarr(KVertTimesName.c_str(),rvptimes,0,nrvpvals);
	    ppr.getarr(KVertFormsName.c_str(),rvpforms,0,nrvpvals-1);
	  }
	} else {
	  BoxLib::Abort(std::string("No vertical permeability function specified for rock: \""+rname).c_str());
	}

	std::string KHoriValName   = KName+".horizontal.vals";
	std::string KHoriTimesName = KName+".horizontal.times";
	std::string KHoriFormsName = KName+".horizontal.forms";
	int nrhpvals = ppr.countval(KHoriValName.c_str());
	if (nrhpvals>0) {
	  ppr.getarr(KHoriValName.c_str(),rhpvals,0,nrhpvals);
	  if (nrhpvals>1) {
	    ppr.getarr(KHoriTimesName.c_str(),rhptimes,0,nrhpvals);
	    ppr.getarr(KHoriFormsName.c_str(),rhpforms,0,nrhpvals-1);
	  }

	} else {
	  BoxLib::Abort(std::string("No horizontal permeability function specified for rock: \""+rname).c_str());
	}

#if BL_SPACEDIM==3
	std::string KHori1ValName   = KName+".horizontal1.vals";
	std::string KHori1TimesName = KName+".horizontal1.times";
	std::string KHori1FormsName = KName+".horizontal1.forms";
	int nrh1pvals = ppr.countval(KHori1ValName.c_str());
	if (nrh1pvals>0) {
	  ppr.getarr(KHori1ValName.c_str(),rh1pvals,0,nrh1pvals);
	  if (nrh1pvals>1) {
	    ppr.getarr(KHori1TimesName.c_str(),rh1ptimes,0,nrh1pvals);
	    ppr.getarr(KHori1FormsName.c_str(),rh1pforms,0,nrh1pvals-1);
	  }

	} else {
	  BoxLib::Abort(std::string("No horizontal1 permeability function specified for rock: \""+rname).c_str());
	}
#endif
      }

      // Convert input permeability values to mDa, if necessary
      if (Permeability_in_m2) {
	BL_ASSERT(!Permeability_in_mDa);
	for (int j=0; j<rvpvals.size(); ++j) {
	  rvpvals[j] *= 1.01325e15; // mDa/m^2
	}
	for (int j=0; j<rhpvals.size(); ++j) {
	  rhpvals[j] *= 1.01325e15; // mDa/m^2
	}
#if BL_SPACEDIM==3
	for (int j=0; j<rh1pvals.size(); ++j) {
	  rh1pvals[j] *= 1.01325e15; // mDa/m^2
	}
#endif
      }

      // Convert hydraulic conductivity to intrinsic permeability, if necessary
      // FIXME: Reflect this conversion in the GSLib file as well
      if (hasHC) {
	Real HC_to_IP_conversion_factor = liquid_viscosity / (liquid_density * gravity);

	for (int j=0; j<rvpvals.size(); ++j) {
	  rvpvals[j] *= HC_to_IP_conversion_factor;
	}
	for (int j=0; j<rhpvals.size(); ++j) {
	  rhpvals[j] *= HC_to_IP_conversion_factor;
	}
#if BL_SPACEDIM==3
	for (int j=0; j<rh1pvals.size(); ++j) {
	  rh1pvals[j] *= HC_to_IP_conversion_factor;
	}
#endif
      }

      // The permeability is now in mDa.  
      // This needs to be multiplied with 1e-10 to be consistent 
      // with the other units in the code.  What this means is that
      // we will be evaluating the darcy velocity as:
      //
      //  u_Darcy [m/s] = ( kappa [X . mD] / mu [Pa.s] ).Grad(p) [atm/m]
      //
      // where X is the factor necessary to have this formula be dimensionally
      // consistent.  X here is 1.e-10, and can be combined with kappa for the 
      // the moment because no other derived quantities depend directly on the 
      // value of kappa  (NOTE: We will have to know that this is done however
      // if kappa is used as a diagnostic or in some way for a derived quantity).
      //
      for (int j=0; j<rvpvals.size(); ++j) {
	rvpvals[j] *= 1.e-10;
      }
      for (int j=0; j<rhpvals.size(); ++j) {
	rhpvals[j] *= 1.e-10;
      }
#if BL_SPACEDIM==3
      for (int j=0; j<rh1pvals.size(); ++j) {
	rh1pvals[j] *= 1.e-10;
      }
#endif

      if (rvpvals.size()>1 || rhpvals.size()>1) {
	Array<TabularFunction> pft(BL_SPACEDIM);
	pft[0] = TabularFunction(rhptimes,rhpvals,rhpforms);
	pft[BL_SPACEDIM-1] = TabularFunction(rvptimes,rvpvals,rvpforms);
#if BL_SPACEDIM==3
	pft[1] = TabularFunction(rh1ptimes,rh1pvals,rh1pforms);
#endif
	kappa_func = new TabularInTimeProperty(PermeabilityName,pft,harm_crsn,pc_refine);
      }
      else {
	Array<Real> vals(BL_SPACEDIM);
	vals[0] = rhpvals[0]; vals[BL_SPACEDIM-1] = rvpvals[0];
#if BL_SPACEDIM==3
	vals[1] = rh1pvals[0];
#endif
	kappa_func = new ConstantProperty(PermeabilityName,vals,harm_crsn,pc_refine);
      }

      // Set old-style values
      Array<Real> rpermeability(BL_SPACEDIM,rvpvals[0]);
      for (int j=0;j<BL_SPACEDIM-1;j++) rpermeability[j] = rhpvals[0];

    }
    else {
      // K_dist == gslib
      std::string gslib_param_file, gslib_data_file;
      generate_perm_gslib_file = (pprK.countval(KGSParamFileName.c_str()) != 0);
      if (pprK.countval(KGSDataFileName.c_str()) == 0) {
        pprK.get(KGSParamFileName.c_str(),gslib_param_file);
        gslib_data_file="permeability.gslib";
      }
      else {
        pprK.query(KGSParamFileName.c_str(),gslib_param_file);
        pprK.get(KGSDataFileName.c_str(),gslib_data_file);
      }

      Array<Real> gslib_file_shift(BL_SPACEDIM,0);
      pprK.queryarr(KGSFileShiftName.c_str(),gslib_file_shift,0,BL_SPACEDIM);

      Real avg; pprK.get(KValName.c_str(),avg);

      // Scale (as above)
      if (Permeability_in_m2) {
	BL_ASSERT(!Permeability_in_mDa);
	avg *= 1.01325e15; // to mDa
      }
      avg *= 1.e-10;

      kappa_func = new GSLibProperty(PermeabilityName,avg,gslib_param_file,gslib_data_file,gslib_file_shift,harm_crsn,pc_refine);
    }

    // capillary pressure: include cpl_coef, residual_saturation, sigma
    const std::string cpl_prefix(prefix+".cpl");
    ParmParse pp_cpl(cpl_prefix.c_str());
    std::string cpl_model = CP_model_None; pp_cpl.query("type",cpl_model);
    std::map<std::string,int>::const_iterator it = CP_models.find(cpl_model);
    int rcplType = -1;
    int rKrType = -1;

    Array<Real> rcplParam(MAX_CPL_PARAMS);

    if (it != CP_models.end()) {
      rcplType = it->second;

      bool is_vG = Is_CP_model_XX(rcplType,CP_model_vG);
      bool is_BC = Is_CP_model_XX(rcplType,CP_model_BC);

      if (rcplType == CP_models[CP_model_None]) {
        rKrType = Kr_models[Kr_model_None];
      }
      else if (is_vG || is_BC)
      {
        Real m, lambda;
        if (is_vG) {
          pp_cpl.get("m",m);
          if (m <= 0) {
            if (ParallelDescriptor::IOProcessor()) {
              std::cerr << "Invalid m (= " << m << " ) for Capillary Pressure model in material: \"" << rname << "\"" << std::endl;
            } BoxLib::Abort();
          }
        } else {
          pp_cpl.get("lambda",lambda);
          if (lambda <= 0) {
            if (ParallelDescriptor::IOProcessor()) {
              std::cerr << "Invalid lambda (= " << lambda << " ) for Capillary Pressure model in material: \"" << rname << "\"" << std::endl;
            } BoxLib::Abort();
          }
        }          
        Real Sr; pp_cpl.get("Sr",Sr);
        if (Sr < 0 || Sr > 1) {
          if (ParallelDescriptor::IOProcessor()) {
            std::cerr << "Invalid Sr (= " << Sr << " ) for Capillary Pressure model in material: \"" << rname << "\"" << std::endl;
          } BoxLib::Abort();
        }

        Real alpha; pp_cpl.get("alpha",alpha);
        if (alpha < 0 ) {
          if (ParallelDescriptor::IOProcessor()) {
            std::cerr << "Invalid alpha (= " << m << " ) for Capillary Pressure model in material: \"" << rname << "\"" << std::endl;
          } BoxLib::Abort();
        }

        std::string Kr_model; ppr.get("Kr_model", Kr_model);
	if (Kr_model != Kr_model_None
	    && Kr_model != Kr_model_Mualem
	    && Kr_model != Kr_model_Burdine) {
	  BoxLib::Abort(std::string("Invalid Kr model: \""+Kr_model
				    +"\" given for material: \""+rname+"\"").c_str());
	}
        std::string Kr_full_model_name = cpl_model + "_" + Kr_model;
        std::map<std::string,int>::const_iterator itKr = Kr_models.find(Kr_full_model_name);
        if (it != Kr_models.end()) {
          rKrType = itKr->second;
        }
        else {
          if (ParallelDescriptor::IOProcessor()) {
            std::cerr << "Invalid Kr model (= \"" << Kr_model
                      << "\") for Relative Permeability with Capillary Pressure model (\"" << cpl_model
                      << "\") in material: \"" << rname << "\""<< std::endl;
          } BoxLib::Abort();
        }

        // Get the ell value (optional for vG+{Mualem,Burdine}, required for all others
        Real Kr_ell;
        if (Is_Kr_model_XX(rKrType,Kr_model_vG_Mualem)) {
          Kr_ell = Kr_ell_vG_Mualem_DEF; ppr.query("Kr_ell",Kr_ell);
        }
        else if (Is_Kr_model_XX(rKrType,Kr_model_vG_Burdine)) {
          Kr_ell = Kr_ell_vG_Burdine_DEF; ppr.query("Kr_ell",Kr_ell);
        }
        else {
          ppr.get("Kr_ell",Kr_ell);
        }

        Real Kr_smoothing_max_pcap = Kr_smoothing_max_pcap_DEF;
        ppr.query("Kr_smoothing_max_pcap",Kr_smoothing_max_pcap);

        if (ppr.countval("WRM_plot_file")>0) {
          ppr.get("WRM_plot_file",WRM_plot_file[i].second);

          WRM_plot_file[i].first = NUM_INIT_INTERP_EVAL_PTS_DEF;
          ppr.query("WRM_plot_file_num_pts",WRM_plot_file[i].first);
        }

	// Convert input alpa values to invAtm, if necessary	
	if (Capillary_Pressure_alpha_in_invPa) {
	  BL_ASSERT(!Capillary_Pressure_alpha_in_invAtm);
	  alpha *= 1.01325e5;
	}

        // Finally, load array of Real numbers for this model
        rcplParam[CPL_MODEL_ID] = (Real)rcplType;
        if (is_vG) {
          rcplParam[VG_M]     = m;
          rcplParam[VG_ALPHA] = alpha;
          rcplParam[VG_SR]    = Sr;
          rcplParam[VG_ELL]   = Kr_ell;
          rcplParam[VG_KR_MODEL_ID]  = (Real)rKrType;
          rcplParam[VG_KR_SMOOTHING_MAX_PC] = Kr_smoothing_max_pcap;
        }
        else {
          rcplParam[BC_LAMBDA] = lambda;
          rcplParam[BC_ALPHA]  = alpha;
          rcplParam[BC_SR]     = Sr;
          rcplParam[BC_ELL]    = Kr_ell;
          rcplParam[BC_KR_MODEL_ID]  = (Real)rKrType;
          rcplParam[BC_KR_SMOOTHING_MAX_PC] = Kr_smoothing_max_pcap;
        }
      }
      else {
        if (ParallelDescriptor::IOProcessor()) {
          std::cerr << "Unknown capillary pressure (" << cpl_model << ") model for " << rname << std::endl;
        }
      }
    }

    Property* cpl_func = new ConstantProperty(CapillaryPressureName,rcplParam,arith_crsn,pc_refine);

    Array<std::string> region_names;
    ppr.getarr("regions",region_names,0,ppr.countval("regions"));
    Array<const Region*> rregions = region_manager->RegionPtrArray(region_names);
    for (int j=0; j<region_names.size(); ++j) {
      material_regions.push_back(region_names[j]);
    }

    std::vector<Property*> properties;
    if (phi_func) {
      properties.push_back(phi_func);
    }
    properties.push_back(kappa_func);
    if (Dmolec_func) {
      properties.push_back(Dmolec_func);
    }
    if (Dispersivity_func) {
      properties.push_back(Dispersivity_func);
    }
    if (Tortuosity_func) {
      properties.push_back(Tortuosity_func);
    }
    if (SpecificStorage_func) {
      properties.push_back(SpecificStorage_func);
    }
    properties.push_back(cpl_func);
    rock.set(i,new Material(rname,rregions,properties));

    delete phi_func;
    delete kappa_func;
    delete Dmolec_func;
    delete Tortuosity_func;
    delete Dispersivity_func;
    delete SpecificStorage_func;
    delete cpl_func;
  }

  nsorption_isotherms = 0;
  ncation_exchange = 0;
  if (known_solutes.size() > 0) {

    // Read rock parameters associated with chemistry
    using_sorption = false; // TRUE if CationExchangeCapacity, SorptionIsotherms, or SorptionSites specified
    aux_chem_variables.clear();
    
    ParmParse ppm("mineral");
    nminerals = ppm.countval("minerals");
    minerals.resize(nminerals);
    if (nminerals>0) {
      ppm.getarr("minerals",minerals,0,nminerals);
    }

    ParmParse pps("sorption_site");
    nsorption_sites = pps.countval("sorption_sites");
    sorption_sites.resize(nsorption_sites);
    if (nsorption_sites>0) {
      using_sorption = true;
      pps.getarr("sorption_sites",sorption_sites,0,nsorption_sites);
    }

    ICParmPair sorption_isotherm_options;
    sorption_isotherm_options[          "Kd"] = -1;
    sorption_isotherm_options[  "Langmuir_b"] = 1;
    sorption_isotherm_options["Freundlich_n"] = 1;
	
    for (int k=0; k<known_solutes.size(); ++k) {
      for (ICParmPair::const_iterator it=sorption_isotherm_options.begin();
           it!=sorption_isotherm_options.end(); ++it) {
        const std::string& str = it->first;
        bool found = false;
        for (int i=0; i<nrock; ++i) {
          const std::string prefix("rock."+r_names[i]+".sorption_isotherms."+known_solutes[k]);
          ParmParse pprs(prefix.c_str());
          if (pprs.countval(str.c_str())) {
            pprs.get(str.c_str(),sorption_isotherm_ics[r_names[i]][known_solutes[k]][str]);
            found = true;
          }
        }

        if (found) {
          using_sorption = true;
          nsorption_isotherms = known_solutes.size();
          for (int i=0; i<nrock; ++i) {
            if (sorption_isotherm_ics[r_names[i]][known_solutes[k]].count(str) == 0) {
              sorption_isotherm_ics[r_names[i]][known_solutes[k]][str] = it->second; // set to default value
            }
          }
          const std::string label = str+"_"+known_solutes[k];
          if (aux_chem_variables.find(label) == aux_chem_variables.end()) {
            sorption_isotherm_label_map[known_solutes[k]][str] = aux_chem_variables.size();
            aux_chem_variables[label]=aux_chem_variables.size()-1;
          }
        }
      }
    }

    // Make a final pass to be sure that if any isotherm parameters are set for a material, they are defaulted for all remaining materials
    if (nsorption_isotherms > 0) {
      for (int k=0; k<known_solutes.size(); ++k) {
	for (ICParmPair::const_iterator it=sorption_isotherm_options.begin();
	     it!=sorption_isotherm_options.end(); ++it) {
	  const std::string& str = it->first;
	  nsorption_isotherms = known_solutes.size();
	  for (int i=0; i<nrock; ++i) {
	    if (sorption_isotherm_ics[r_names[i]][known_solutes[k]].count(str) == 0) {
              sorption_isotherm_ics[r_names[i]][known_solutes[k]][str] = it->second; // set to default value
            }
          }
	}
      }
    }

    ICParmPair cation_exchange_options;
    cation_exchange_options["Cation_Exchange_Capacity"] = 0;
    cation_exchange_options["Ion_Exchange_Reference_Cation_Concentration_0"] = 1.e-9;

    bool ceopt_found = false;
    for (ICParmPair::const_iterator it=cation_exchange_options.begin(); it!=cation_exchange_options.end(); ++it) {
      //const std::string& parm_name = it->first;
      const std::string& in_parm_name = it->first;
      const std::string& out_parm_name = 
        (it->first == "Cation_Exchange_Capacity" ? "Ion_Exchange_Site_Density_0": it->first);
      for (int i=0; i<nrock; ++i) {
        const std::string prefix("rock."+r_names[i]);
        ParmParse pprs(prefix.c_str());
        if (pprs.countval(in_parm_name.c_str())) {
          pprs.get(in_parm_name.c_str(),cation_exchange_ics[r_names[i]][out_parm_name]);
          ceopt_found = true;
        }
      }
    }

    // If any of these found, set for all of them
    if (ceopt_found) {
      for (ICParmPair::const_iterator it=cation_exchange_options.begin(); it!=cation_exchange_options.end(); ++it) {
        const std::string& parm_name =
          (it->first == "Cation_Exchange_Capacity" ? "Ion_Exchange_Site_Density_0": it->first);
        using_sorption = true;
        ncation_exchange = 1;
        for (int i=0; i<nrock; ++i) {
          if (cation_exchange_ics[r_names[i]].count(parm_name) == 0) {
            cation_exchange_ics[r_names[i]][parm_name] = it->second; // set to default value
          }
        }
        
        // Add a aux_chem variable slot for this quantity
        if (aux_chem_variables.find(parm_name) == aux_chem_variables.end())  {
          cation_exchange_label_map[parm_name] = aux_chem_variables.size();
          aux_chem_variables[parm_name]=aux_chem_variables.size()-1;
        }
      }
    }

    ICParmPair mineralogy_options;
    mineralogy_options[      "Volume_Fraction"] = 0;
    mineralogy_options["Specific_Surface_Area"] = 0;
    for (int k=0; k<minerals.size(); ++k) {
      for (ICParmPair::const_iterator it=mineralogy_options.begin(); it!=mineralogy_options.end(); ++it) {
        const std::string& str = it->first;
        bool found = false;
        for (int i=0; i<nrock; ++i) {
          const std::string prefix("rock."+r_names[i]+".Mineralogy."+minerals[k]);
          ParmParse pprs(prefix.c_str());
          if (pprs.countval(str.c_str())) {
            pprs.get(str.c_str(),mineralogy_ics[r_names[i]][minerals[k]][str]);
            found = true;
          }
        }
	    
        if (found) {
          for (int i=0; i<nrock; ++i) {
            if (mineralogy_ics[r_names[i]][minerals[k]].count(str) == 0) {
              mineralogy_ics[r_names[i]][minerals[k]][str] = it->second; // set to default value
            }
          }
          //std::cout << "****************** mineralogy_ics[" << r_names[i] << "][" << minerals[k] 
          //	  << "][" << str << "] = " << mineralogy_ics[r_names[i]][minerals[k]][str] 
          //	  << std::endl;

          const std::string label = minerals[k]+"_"+str;
          if (aux_chem_variables.find(label) == aux_chem_variables.end()) {
            mineralogy_label_map[minerals[k]][str] = aux_chem_variables.size();
            aux_chem_variables[label]=aux_chem_variables.size()-1;
          }
        }
      }
    }

    ICParmPair complexation_options;
    complexation_options["Site_Density"] = 0;
    for (int k=0; k<sorption_sites.size(); ++k) {
      for (ICParmPair::const_iterator it=complexation_options.begin(); it!=complexation_options.end(); ++it) {
        const std::string& str = it->first;
        bool found = false;
        for (int i=0; i<nrock; ++i) {
          const std::string prefix("rock."+r_names[i]+".Surface_Complexation_Sites."+sorption_sites[k]);
          ParmParse pprs(prefix.c_str());
          if (pprs.countval(str.c_str())) {
            pprs.get(str.c_str(),surface_complexation_ics[r_names[i]][sorption_sites[k]][str]);
            found = true;
          }
        }
	    
        if (found) {
          using_sorption = true;
          for (int i=0; i<nrock; ++i) {
            if (surface_complexation_ics[r_names[i]][sorption_sites[k]].count(str) == 0) {
              surface_complexation_ics[r_names[i]][sorption_sites[k]][str] = it->second; // set to default value
            }
          }
          //std::cout << "****************** surface_complexation_ics[" << r_names[i] << "][" << sorption_sites[k] 
          //	  << "][" << str << "] = " << sorption_isotherm_ics[r_names[i]][sorption_sites[k]][str] 
          //	  << std::endl;
	      
          const std::string label = str+"_"+sorption_sites[k];
          if (aux_chem_variables.find(label) == aux_chem_variables.end()) {
            surface_complexation_label_map[sorption_sites[k]][str] = aux_chem_variables.size();
            aux_chem_variables[label]=aux_chem_variables.size()-1;
          }
        }
      }
    }
  }
}

RockManager::~RockManager()
{
  BL_PROFILE("RockManager::~RockManager()");

  delete materialFiller;
}

bool
RockManager::CanDerive(const std::string& property_name) const
{
  BL_PROFILE("RockManager::CanDerive()");

  return materialFiller->CanDerive(property_name);
}

bool
RockManager::GetProperty(Real               time,
                         int                level,
                         MultiFab&          mf,
                         const std::string& pname,
                         int                dComp,
                         int                nGrow,
                         void*              ctx,
                         bool               ignore_mixed) const
{
  BL_PROFILE("RockManager::GetProperty()");

  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }

  if (materialFiller->NComp(pname) <= 0) {
    return false;
  }

  BL_ASSERT(mf.nGrow()>=nGrow);

  bool ret = materialFiller->SetProperty(time,level,mf,pname,dComp,nGrow,0,ignore_mixed);

  // Fill any GSLib-bsaed quantities
  int nComp = NComp(pname);
  Array<Real> dummy(nComp);
  bool do_touchup = false;
  for (int i=0; i<rock.size(); ++i) {
    const Property* prop = rock[i].Prop(pname);
    const GSLibProperty* gslib_prop = dynamic_cast<const GSLibProperty*>(prop);
    if (gslib_prop != 0) {
      const AmrData* this_const_amrData = gslib_prop->GetAmrData();
      AmrData* this_amrData = const_cast<AmrData*>(this_const_amrData);
      const Geometry& geom = materialFiller->Geom(level);

      Array<int> destFillComps(nComp);
      for (int n=0; n<nComp; ++n) {
	destFillComps[n] = n;
      }

      // Build a boxarray that includes grow cells, except where they extend out the domain
      BoxArray bavals = BoxArray(mf.boxArray()).grow(nGrow);
      const Box& this_domain = this_amrData->ProbDomain()[level];
      for (int j=0; j<bavals.size(); ++j) {
	bavals.set(j,Box(bavals[j]) & this_domain);
      }
      MultiFab valstmp(bavals,nComp,0);
      this_amrData->FillVar(valstmp,level,gslib_prop->PlotfileVars(),destFillComps);

      MultiFab vals(mf.boxArray(),nComp,nGrow);
      vals.setVal(-1); // Put something computable outside domain, touchup later

      // Copy filled property over default
      for (MFIter mfi(vals); mfi.isValid(); ++mfi) {
	vals[mfi].copy(valstmp[mfi]);
      }
      valstmp.clear(); // No longer needed

      // Set id mask in order to pull out the cells of vals in this material
      iMultiFab id(mf.boxArray(),1,nGrow);
      materialFiller->SetMaterialID(level,id,nGrow,ignore_mixed);

      for (MFIter mfi(vals); mfi.isValid(); ++mfi) {
	Box bx = BoxLib::grow(mfi.validbox(),nGrow);
        const FArrayBox& vfab = vals[mfi];
        const IArrayBox& idfab = id[mfi];
        FArrayBox& mfab = mf[mfi];

	BL_ASSERT(mfab.nComp() >= nComp+dComp);
	BL_ASSERT(vfab.nComp() >= nComp);

        FORT_FILLPMAT (mfab.dataPtr(dComp), ARLIM( mfab.loVect()), ARLIM( mfab.hiVect()),
                       idfab.dataPtr(),     ARLIM(idfab.loVect()), ARLIM(idfab.hiVect()),
                        vfab.dataPtr(),     ARLIM( vfab.loVect()), ARLIM( vfab.hiVect()),
                       &i, bx.loVect(), bx.hiVect(), &nComp);
	do_touchup = true;
      }
    }
  }
  ret = true;
  ParallelDescriptor::ReduceBoolOr(do_touchup);
  if (do_touchup) {
    FillBoundary(time,level,mf,dComp,nComp,nGrow);
  }

  return ret;
}

void
RockManager::GetMaterialID(int level, iMultiFab& mf, int nGrow, bool ignore_mixed) const
{
  BL_PROFILE("RockManager::GetMaterialID()");
  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }
  return materialFiller->SetMaterialID(level,mf,nGrow,ignore_mixed);
}


static Real vgMKr(Real seff, Real m, Real mI, Real ell) {
  return std::pow(seff, ell) * std::pow(1-std::pow(1-std::pow(seff,mI),m),2);
}

static Real vgBKr(Real seff, Real m, Real mI, Real ell) {
  return std::pow(seff, ell) * ( 1 - std::pow(1-std::pow(seff,mI),m) );
}

static Real vgPc(Real seff, Real mI, Real nI, Real alphaI) {
  return alphaI * std::pow( (seff==0 ? 1 : std::pow(seff,-mI)) - 1, nI);
}

static Real vgPcInv(Real pc, Real m, Real n, Real alpha) {
  return std::pow(1 + std::pow(alpha*pc,n),-m);
}

static Real bcMKr(Real seff, Real lambda, Real ell) {
  return std::pow(seff,ell+2+2/lambda);
}

static Real bcBKr(Real seff, Real lambda, Real ell) {
  return std::pow(seff,ell+1+2/lambda);
}

static Real bcPc(Real seff, Real lambdaI, Real alphaI) {
  return alphaI * std::pow(seff,-lambdaI);
}

static Real bcPcInv(Real pc, Real lambda, Real alpha) {
  return std::pow(alpha*pc,-lambda);
}

// Capillary Pressure (given Saturation)
void
RockManager::CapillaryPressure(const Real* saturation, int* matID, Real time, Real* capillaryPressure, int npts) const
{
  BL_PROFILE("RockManager::CapillaryPressure()");
  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }
  Array<Array<int> > mat_pts = SortPtsByMaterial(matID,npts);

  // Make temp structure to interact with Property interface
  int nComp = materialFiller->NComp(CapillaryPressureName);
  Array<Real> pc_params(nComp);

  for (int j=0; j<rock.size(); ++j) {

    int N = mat_pts[j].size();

    if (N>0) {
      const Property* p = rock[j].Prop(CapillaryPressureName); BL_ASSERT(p!=0);
      p->Evaluate(time,pc_params);
      bool is_vG = Is_CP_model_XX(pc_params[CPL_MODEL_ID],CP_model_vG);
      bool is_BC = Is_CP_model_XX(pc_params[CPL_MODEL_ID],CP_model_BC);

      if (is_vG) {
        Real m      = pc_params[VG_M];
        Real alpha  = pc_params[VG_ALPHA];
        Real Sr     = pc_params[VG_SR];

        Real mI     = 1/m;
        Real omSrI  = 1/(1-Sr);
        Real alphaI = 1/alpha;
        Real n      = 1/(1-m);
        Real nI     = 1-m;

        for (int i=0; i<N; ++i) {
          int idx = mat_pts[j][i];
          Real s = std::min(1.0, std::max(0.0, saturation[idx]));        
          Real seff = std::max(0.,std::min(1.,(s - Sr)*omSrI));
          capillaryPressure[idx] = vgPc(seff,mI,nI,alphaI);
        }
      }
      else if (is_BC) {
        Real lambda   = pc_params[BC_LAMBDA];
        Real alpha    = pc_params[BC_ALPHA];
        Real Sr       = pc_params[BC_SR];

        Real lambdaI  = 1/lambda;
        Real alphaI   = 1/alpha;
        Real omSrI    = 1/(1-Sr);

        for (int i=0; i<N; ++i) {
          int idx = mat_pts[j][i];
          Real s = std::min(1.0, std::max(0.0, saturation[idx]));        
          Real seff = std::max(0.,std::min(1.,(s - Sr)*omSrI));
          capillaryPressure[idx] = bcPc(seff,lambdaI,alphaI);
        }
      }
      else {
        for (int i=0; i<N; ++i) {
          int idx = mat_pts[j][i];
	  capillaryPressure[idx] = 0;
	}
      }
    }
  }
}
 
void
RockManager::CapillaryPressure(const MultiFab&  saturation,
                               const iMultiFab& matID,
                               Real             time,
                               MultiFab&        capillaryPressure,
                               int              sComp,
                               int              dComp,
                               int              nGrow) const
{
  BL_PROFILE("RockManager::CapillaryPressure1()");
  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }
  BL_ASSERT(saturation.boxArray() == capillaryPressure.boxArray());
  BL_ASSERT(sComp < saturation.nComp()  && dComp < capillaryPressure.nComp());
  BL_ASSERT(nGrow <= saturation.nGrow() && nGrow <= capillaryPressure.nGrow());
  FArrayBox s, pc;
  IArrayBox id;
  for (MFIter mfi(saturation); mfi.isValid(); ++mfi) {
    Box box = Box(mfi.validbox()).grow(nGrow);
    s.resize(box,1); s.copy(saturation[mfi],box,sComp,box,0,1);
    id.resize(box,1); id.copy(matID[mfi],box,0,box,0,1);
    pc.resize(box,1);
    CapillaryPressure(s.dataPtr(),id.dataPtr(),time,pc.dataPtr(),box.numPts());
    capillaryPressure[mfi].copy(pc,box,0,box,dComp,1);
  }
}

Array<Array<int> >
RockManager::SortPtsByMaterial(int* matID, int npts) const
{
  BL_PROFILE("RockManager::SortPtsByMaterial()");

  Array<Array<int> > mat_pts(rock.size());
  for (int i=0; i<npts; ++i) {
    mat_pts[matID[i]].push_back(i);
  }
  return mat_pts;
}

// Inverse Capillary Pressure (Saturation given Capillary Pressure)
void
RockManager::InverseCapillaryPressure(const Real* capillaryPressure, int* matID, Real time, Real* saturation, int npts) const
{
  BL_PROFILE("RockManager::InverseCapillaryPressure()");
  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }
  Array<Array<int> > mat_pts = SortPtsByMaterial(matID,npts);

  // Make temp structure to interact with Property interface
  int nComp = materialFiller->NComp(CapillaryPressureName);
  Array<Real> pc_params(nComp);

  for (int j=0; j<rock.size(); ++j) {
    const Property* p = rock[j].Prop(CapillaryPressureName); BL_ASSERT(p!=0);
    p->Evaluate(time,pc_params);
    bool is_vG = Is_CP_model_XX(pc_params[CPL_MODEL_ID],CP_model_vG);
    bool is_BC = Is_CP_model_XX(pc_params[CPL_MODEL_ID],CP_model_BC);

    if (is_vG) {
      Real m     = pc_params[VG_M];
      Real alpha = pc_params[VG_ALPHA];
      Real Sr    = pc_params[VG_SR];
      Real n     = 1./(1-m);
      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Real seff = (capillaryPressure[idx] <= 0  ? 1 : 
                     std::pow( std::pow(alpha*capillaryPressure[idx],n) + 1, -m));
        seff = std::max(0.,std::min(1.,seff));
        saturation[idx] = seff*(1 - Sr) + Sr;
      }
    }
    else if (is_BC) {
      Real mLambda = -pc_params[BC_LAMBDA];
      Real alpha   =  pc_params[BC_ALPHA];
      Real Sr      =  pc_params[BC_SR];
      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Real seff = (capillaryPressure[idx] <= 0  ? 1 : 
                     std::pow(alpha*capillaryPressure[idx],mLambda));
        seff = std::max(0.,std::min(1.,seff));
        saturation[idx] = seff*(1 - Sr) + Sr;
      }
    }
    else {
      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
	saturation[idx] = 1.0;
      }
    }
  }
}

void
RockManager::InverseCapillaryPressure(const MultiFab&  capillaryPressure,
                                      const iMultiFab& matID,
                                      Real             time,
                                      MultiFab&        saturation,
                                      int              sComp,
                                      int              dComp,
                                      int              nGrow) const
{
  BL_PROFILE("RockManager::InverseCapillaryPressure1()");
  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }
  BL_ASSERT(saturation.boxArray() == capillaryPressure.boxArray());
  BL_ASSERT(dComp < saturation.nComp()  && sComp < capillaryPressure.nComp());
  BL_ASSERT(nGrow <= saturation.nGrow() && nGrow <= capillaryPressure.nGrow());
  FArrayBox s, pc;
  IArrayBox id;
  for (MFIter mfi(saturation); mfi.isValid(); ++mfi) {
    Box box = Box(mfi.validbox()).grow(nGrow);
    pc.resize(box,1); pc.copy(capillaryPressure[mfi],box,sComp,box,0,1);
    id.resize(box,1); id.copy(matID[mfi],box,0,box,0,1);
    s.resize(box,1);
    InverseCapillaryPressure(pc.dataPtr(),id.dataPtr(),time,s.dataPtr(),box.numPts());
    saturation[mfi].copy(s,box,0,box,dComp,1);
  }
}

// D (Saturation) / D(CapillaryPressure)
void
RockManager::DInverseCapillaryPressure(const Real* saturation, int* matID, Real time, Real* DsaturationDpressure, int npts) const
{
  BL_PROFILE("RockManager::DInverseCapillaryPressure()");
  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }
  Array<Array<int> > mat_pts = SortPtsByMaterial(matID,npts);

  // Make temp structure to interact with Property interface
  int nComp = materialFiller->NComp(CapillaryPressureName);
  Array<Real> pc_params(nComp);

  for (int j=0; j<rock.size(); ++j) {
    const Property* p = rock[j].Prop(CapillaryPressureName); BL_ASSERT(p!=0);
    p->Evaluate(time,pc_params);
    bool is_vG = Is_CP_model_XX(pc_params[CPL_MODEL_ID],CP_model_vG);
    bool is_BC = Is_CP_model_XX(pc_params[CPL_MODEL_ID],CP_model_BC);

    if (is_vG) {
      Real m     = pc_params[VG_M];
      Real alpha = pc_params[VG_ALPHA];
      Real Sr    = pc_params[VG_SR];
      Real n     = 1./(1-m);
      Real b     = -1/m;
      Real fac   = - (1 - Sr)*alpha*m*n;
      Real omSrI = 1/(1 - Sr);

      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Real seff = std::max(0.,std::min(1.,(saturation[idx] - Sr)*omSrI));
        Real sb = std::pow(seff,b);
        DsaturationDpressure[idx] = fac*std::pow(sb-1,m)*seff/sb;
      }
    }
    else if (is_BC) {
      Real lambda = pc_params[BC_LAMBDA];
      Real alpha  = pc_params[BC_ALPHA];
      Real Sr     = pc_params[BC_SR];
      Real fac    = -alpha*lambda;
      Real oplI   = 1+1/lambda;
      Real omSrI = 1/(1 - Sr);

      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Real seff = std::max(0.,std::min(1.,(saturation[idx] - Sr)*omSrI));
        DsaturationDpressure[idx] = fac*std::pow(seff,oplI);
      }
    }
  }
}

void
RockManager::DInverseCapillaryPressure(const MultiFab&  saturation,
                                       const iMultiFab& matID,
                                       Real             time,
                                       MultiFab&        DsaturationDpressure,
                                       int              sComp,
                                       int              dComp,
                                       int              nGrow) const
{
  BL_PROFILE("RockManager::DInverseCapillaryPressure1()");
  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }
  BL_ASSERT(saturation.boxArray() == DsaturationDpressure.boxArray());
  BL_ASSERT(dComp < saturation.nComp()  && sComp < DsaturationDpressure.nComp());
  BL_ASSERT(nGrow <= saturation.nGrow() && nGrow <= DsaturationDpressure.nGrow());
  FArrayBox s, dsdp;
  IArrayBox id;
  for (MFIter mfi(saturation); mfi.isValid(); ++mfi) {
    Box box = Box(mfi.validbox()).grow(nGrow);
    s.resize(box,1); s.copy(saturation[mfi],box,sComp,box,0,1);
    id.resize(box,1); id.copy(matID[mfi],box,0,box,0,1);
    dsdp.resize(box,1);
    DInverseCapillaryPressure(s.dataPtr(),id.dataPtr(),time,dsdp.dataPtr(),box.numPts());
    DsaturationDpressure[mfi].copy(dsdp,box,0,box,dComp,1);
  }
}

static Real KrInterp(Real se, Real seth, Real m_th, Real m_int) {
  Real dels  = 1 - se;
  Real dels1 = 1 - seth;
  return 1 + dels*dels * m_int/dels1 +  dels*dels * (dels-dels1) * (m_th - 2*m_int) / (dels1*dels1);
}

void
RockManager::RelativePermeability(const Real* saturation, int* matID, Real time, Real* kappa, int npts) const
{
  BL_PROFILE("RockManager::RelativePermeability()");
  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }
  Array<Array<int> > mat_pts = SortPtsByMaterial(matID,npts);

  // Make temp structure to interact with Property interface
  int nComp = materialFiller->NComp(CapillaryPressureName);
  Array<Real> pc_params(nComp);

  for (int j=0; j<rock.size(); ++j) {
    const Property* p = rock[j].Prop(CapillaryPressureName); BL_ASSERT(p!=0);
    p->Evaluate(time,pc_params);
    bool is_vG = Is_CP_model_XX(pc_params[CPL_MODEL_ID],CP_model_vG);
    bool is_BC = Is_CP_model_XX(pc_params[CPL_MODEL_ID],CP_model_BC);

    if (is_vG) {

      Real m                     = pc_params[VG_M];
      Real alpha                 = pc_params[VG_ALPHA];
      Real Sr                    = pc_params[VG_SR];
      Real ell                   = pc_params[VG_ELL];
      int  Kr_model              = pc_params[VG_KR_MODEL_ID];
      Real Kr_smoothing_max_pcap = pc_params[BC_KR_SMOOTHING_MAX_PC];

      Real omSrI = 1/(1 - Sr);
      Real mI = 1/m;
      Real n = 1/(1-m);
      bool is_Mualem  = Is_Kr_model_XX(Kr_model,Kr_model_vG_Mualem);
      bool is_Burdine = Is_Kr_model_XX(Kr_model,Kr_model_vG_Burdine);

      // Find smoothing interval in terms of seff (need to evaluate here since properties may be time-dependent)
      std::pair<bool,Real>& Kr_min_seff = Kr_smoothing_min_seff[j];
      if ( (Kr_smoothing_max_pcap > 0)
           && (p->isTimeDependent() || Kr_min_seff.first) )
      {
        Kr_min_seff.second = vgPcInv(Kr_smoothing_max_pcap,m,n,alpha);
        Kr_min_seff.first = false;

        if (ParallelDescriptor::IOProcessor()) {
          std::cout << "For material \"" << rock[j].Name() << "\" seff thresh is " << Kr_min_seff.second << std::endl;
        }
      }

      if (is_Mualem) {
        Real Kr_slope_thresh, Kr_slope_interval;
        if (Kr_smoothing_max_pcap > 0) {
          Real mI = 1/m;
          Real ds = (1 - Kr_min_seff.second)*0.001;
          Real Kr_thresh1 = vgMKr(Kr_min_seff.second   ,m,mI,ell);
          Real Kr_thresh2 = vgMKr(Kr_min_seff.second+ds,m,mI,ell);
          Kr_slope_thresh = -(Kr_thresh2 - Kr_thresh1)/ds;
          Kr_slope_interval = (Kr_thresh1 - 1)/(1 - Kr_min_seff.second);
        }

        for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
          int idx = mat_pts[j][i];
          Real seff = std::max(0.,std::min(1.,(saturation[idx] - Sr)*omSrI));
          if (seff > Kr_min_seff.second && seff < 1) {
            kappa[idx] = KrInterp(seff, Kr_min_seff.second, Kr_slope_thresh, Kr_slope_interval);
          }
          else {
            kappa[idx] = vgMKr(seff,m,mI,ell);
          }
        }
      }
      else if (is_Burdine) {
        Real Kr_slope_thresh, Kr_slope_interval;
        if (Kr_smoothing_max_pcap > 0) {
          Real mI = 1/m;
          Real ds = (1 - Kr_min_seff.second)*0.001;
          Real Kr_thresh1 = vgBKr(Kr_min_seff.second   ,m,mI,ell);
          Real Kr_thresh2 = vgBKr(Kr_min_seff.second+ds,m,mI,ell);
          Kr_slope_thresh = -(Kr_thresh2 - Kr_thresh1)/ds;
          Kr_slope_interval = (Kr_thresh1 - 1)/(1 - Kr_min_seff.second);
        }

        for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
          int idx = mat_pts[j][i];
          Real seff = std::max(0.,std::min(1.,(saturation[idx] - Sr)*omSrI));
          if (seff > Kr_min_seff.second && seff < 1) {
            kappa[idx] = KrInterp(seff, Kr_min_seff.second, Kr_slope_thresh, Kr_slope_interval);
          }
          else {
            kappa[idx] = vgBKr(seff,m,mI,ell);
          }
        }
      }

    } else if (is_BC) {

      Real lambda = pc_params[BC_LAMBDA];
      Real alpha  = pc_params[BC_ALPHA];
      Real Sr     = pc_params[BC_SR];
      Real ell    = pc_params[BC_ELL];
      Real omSrI  = 1/(1 - Sr);

      bool is_Mualem  = Is_Kr_model_XX(pc_params[VG_KR_MODEL_ID],Kr_model_BC_Mualem);
      bool is_Burdine = Is_Kr_model_XX(pc_params[VG_KR_MODEL_ID],Kr_model_BC_Burdine);

      BL_ASSERT(is_Mualem || is_Burdine);
      Real f = (is_Mualem ? ell + 2 + 2/lambda : ell + 1 + 2/lambda);

      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Real seff = std::max(0.,std::min(1.,(saturation[idx] - Sr)*omSrI));
        kappa[idx] = std::pow(seff,f);
      }

    }
  }
}

void
RockManager::RelativePermeability(const MultiFab&  saturation,
                                  const iMultiFab& matID,
                                  Real             time,
                                  MultiFab&        kappa,
                                  int              sComp,
                                  int              dComp,
                                  int              nGrow) const
{
  BL_PROFILE("RockManager::RelativePermeability1()");
  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }
  BL_ASSERT(saturation.boxArray() == kappa.boxArray());
  BL_ASSERT(dComp < saturation.nComp()  && sComp < kappa.nComp());
  BL_ASSERT(nGrow <= saturation.nGrow() && nGrow <= kappa.nGrow());
  FArrayBox s, k;
  IArrayBox id;
  for (MFIter mfi(saturation); mfi.isValid(); ++mfi) {
    Box box = Box(mfi.validbox()).grow(nGrow);
    s.resize(box,1); s.copy(saturation[mfi],box,sComp,box,0,1);
    id.resize(box,1); id.copy(matID[mfi],box,0,box,0,1);
    k.resize(box,1);
    RelativePermeability(s.dataPtr(),id.dataPtr(),time,k.dataPtr(),box.numPts());
    kappa[mfi].copy(k,box,0,box,dComp,1);
  }
}


void 
RockManager::ResidualSaturation(int* matID, Real time, Real* Sr, int npts) const
{
  BL_PROFILE("RockManager::ResidualSaturation()");
  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }
  Array<Array<int> > mat_pts = SortPtsByMaterial(matID,npts);

  // Make temp structure to interact with Property interface
  int nComp = materialFiller->NComp(CapillaryPressureName);
  Array<Real> pc_params(nComp);

  for (int j=0; j<rock.size(); ++j) {
    const Property* p = rock[j].Prop(CapillaryPressureName); BL_ASSERT(p!=0);
    p->Evaluate(time,pc_params);
    bool is_vG = Is_CP_model_XX(pc_params[CPL_MODEL_ID],CP_model_vG);
    bool is_BC = Is_CP_model_XX(pc_params[CPL_MODEL_ID],CP_model_BC);

    if (is_vG) {
      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Sr[idx] = pc_params[VG_SR];
      }
    } else if (is_BC) {
      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Sr[idx] = pc_params[BC_SR];
      }
    }
  }
}

void
RockManager::ResidualSaturation(const iMultiFab& matID,
                                Real             time,
                                MultiFab&        Sr,
                                int              dComp,
                                int              nGrow) const
{
  BL_PROFILE("RockManager::ResidualSaturation1()");
  if (!finalized) {
    BoxLib::Abort("RockManager not yet functional, must call RockManager::FinalizeBuild");
  }
  BL_ASSERT(dComp < Sr.nComp());
  BL_ASSERT(nGrow <= Sr.nGrow());
  FArrayBox sr;
  IArrayBox id;
  for (MFIter mfi(Sr); mfi.isValid(); ++mfi) {
    Box box = Box(mfi.validbox()).grow(nGrow);
    id.resize(box,1); id.copy(matID[mfi],box,0,box,0,1);
    sr.resize(box,1);
    ResidualSaturation(id.dataPtr(),time,sr.dataPtr(),box.numPts());
    Sr[mfi].copy(sr,box,0,box,dComp,1);
  }
}

// Annoying set of functions necessary make [] operator of maps const
bool
RockManager::Is_CP_model_XX(int model_id, const std::string& str) const
{
  std::map<std::string,int>::const_iterator it = CP_models.find(str);
  if (it != CP_models.end()) {
    return model_id == it->second;
  }
  return false;
}

bool
RockManager::Is_Kr_model_XX(int model_id, const std::string& str) const
{
  std::map<std::string,int>::const_iterator it = Kr_models.find(str);
  if (it != Kr_models.end()) {
    return model_id == it->second;
  }
  return false;
}

void
RockManager::RockChemistryProperties(FArrayBox&  fab,
                                     const Real* dx,
                                     const std::map<std::string,int>& aux_chem_variables_map)
{
  BL_PROFILE("RockManager::RockChemistryProperties()");

  int nSolutes = known_solutes.size();
  if (nSolutes>0) {

    // Note: sorption_isotherm_ics[rockname][solute][property] = val
    for (ChemICMap::const_iterator it=sorption_isotherm_ics.begin(); it!=sorption_isotherm_ics.end(); ++it) {
      const std::string& material_name = it->first;
      const ICLabelParmPair& solute_to_pp = it->second; 
      for (ICLabelParmPair::const_iterator it1=solute_to_pp.begin(); it1!=solute_to_pp.end(); ++it1) {
        const std::string& solute_name = it1->first;
        const ICParmPair& parm_pairs = it1->second;
        for (ICParmPair::const_iterator it2=parm_pairs.begin(); it2!=parm_pairs.end(); ++it2) {
          const std::string& parameter = it2->first;
          std::string key = ChemistryHelper_Structured::BuildPropertyParameterName(solute_name,"Isotherm",parameter);
          std::map<std::string,int>::const_iterator it3 = aux_chem_variables_map.find(key);
          if (it3 == aux_chem_variables_map.end() && ParallelDescriptor::IOProcessor()) {
            std::cout << "RockManager::RockChemistryProperties Unable to locate parameter in aux_data (material,solute,parameter): "
                      << material_name << "   "
                      << solute_name << "    "
                      << parameter << "  key: " << key << std::endl;
            BoxLib::Abort();
          }
          int comp = it3->second;
          Real value = it2->second;
          const Array<const Region*>& rock_regions = FindMaterial(material_name).Regions();
          for (int j=0; j<rock_regions.size(); ++j) {
            rock_regions[j]->setVal(fab,value,comp,dx,0);
          }
        }
      }
    }
  }

  if (nminerals>0) {
    // Note: mineralogy_ics[rockname][mineralname][property] = val
    for (ChemICMap::const_iterator it=mineralogy_ics.begin(); it!=mineralogy_ics.end(); ++it) {
      const std::string& material_name = it->first;
      const ICLabelParmPair& mineral_to_pp = it->second; 
      for (ICLabelParmPair::const_iterator it1=mineral_to_pp.begin(); it1!=mineral_to_pp.end(); ++it1) {
        const std::string& mineral_name = it1->first;
        const ICParmPair& parm_pairs = it1->second;
        for (ICParmPair::const_iterator it2=parm_pairs.begin(); it2!=parm_pairs.end(); ++it2) {
          const std::string& parameter = it2->first;
          std::string key = ChemistryHelper_Structured::BuildPropertyParameterName(mineral_name,parameter);
          std::map<std::string,int>::const_iterator it3 = aux_chem_variables_map.find(key);
          if (it3 == aux_chem_variables_map.end() && ParallelDescriptor::IOProcessor()) {
            std::cout << "RockManager::RockChemistryProperties  Unable to locate parameter in aux_data (material,mineral,parameter): "
                      << material_name << "   "
                      << mineral_name << "    "
                      << parameter << "  key: " << key << std::endl;
            BoxLib::Abort();
          }
          int comp = it3->second;
          Real value = it2->second;
          const Array<const Region*>& rock_regions = FindMaterial(material_name).Regions();
          for (int j=0; j<rock_regions.size(); ++j) {
            rock_regions[j]->setVal(fab,value,comp,dx,0);
          }
        }
      }
    }
  }     

  if (nsorption_sites>0) {
    // Note: surface_complexation_ics[rockname][sorptionsitename][property]) = val
    for (ChemICMap::const_iterator it=surface_complexation_ics.begin(); it!=surface_complexation_ics.end(); ++it) {
      const std::string& material_name = it->first;
      const ICLabelParmPair& sorption_site_to_pp = it->second; 
      for (ICLabelParmPair::const_iterator it1=sorption_site_to_pp.begin(); it1!=sorption_site_to_pp.end(); ++it1) {
        const std::string& sorption_site_name = it1->first;
        const ICParmPair& parm_pairs = it1->second;
        for (ICParmPair::const_iterator it2=parm_pairs.begin(); it2!=parm_pairs.end(); ++it2) {
          const std::string& parameter = it2->first;
          std::string key = ChemistryHelper_Structured::BuildPropertyParameterName(sorption_site_name,"Surface",parameter);
          std::map<std::string,int>::const_iterator it3 = aux_chem_variables_map.find(key);
          if (it3 == aux_chem_variables_map.end() && ParallelDescriptor::IOProcessor()) {
            std::cout << "RockManager::RockChemistryProperties  Unable to locate parameter in aux_data (material,sorption_site,parameter): "
                      << material_name << "   "
                      << sorption_site_name << "    "
                      << parameter << "  key: " << key << std::endl;
            BoxLib::Abort();
          }
          int comp = it3->second;
          Real value = it2->second;
          const Array<const Region*>& rock_regions = FindMaterial(material_name).Regions();
          for (int j=0; j<rock_regions.size(); ++j) {
            rock_regions[j]->setVal(fab,value,comp,dx,0);
          }
        }
      }
    }
  }
 
  if (ncation_exchange>0) {
    // Note: cation_exchange_ics[rockname][property] = val
    for (ICLabelParmPair::const_iterator it=cation_exchange_ics.begin(); it!=cation_exchange_ics.end(); ++it) {
      const std::string& material_name = it->first;
      const ICParmPair& parm_pairs = it->second;
      for (ICParmPair::const_iterator it1=parm_pairs.begin(); it1!=parm_pairs.end(); ++it1) {
        const std::string& key = it1->first;
        std::map<std::string,int>::const_iterator it3 = aux_chem_variables_map.find(key);
        if (it3 == aux_chem_variables_map.end() && ParallelDescriptor::IOProcessor()) {
          std::cout << "RockManager::RockChemistryProperties Unable to locate parameter in aux_data: " << key << std::endl;
          BoxLib::Abort();
        }
        int comp = it3->second;
        Real value = it1->second;
        const Array<const Region*>& rock_regions = FindMaterial(material_name).Regions();
        for (int j=0; j<rock_regions.size(); ++j) {
          rock_regions[j]->setVal(fab,value,comp,dx,0);
        }
      }
    }
  }
}

void
RockManager::RockChemistryProperties(FArrayBox&         fab,
                                     const std::string& material_name,
                                     const std::map<std::string,int>& aux_chem_variables_map)
{
  BL_PROFILE("RockManager::RockChemistryProperties1()");

  int nSolutes = known_solutes.size();
  if (nSolutes>0) {

    // Note: sorption_isotherm_ics[rockname][solute][property] = val
    ChemICMap::const_iterator it = sorption_isotherm_ics.find(material_name);
    if (it != sorption_isotherm_ics.end()) {
      const ICLabelParmPair& solute_to_pp = it->second; 
      for (ICLabelParmPair::const_iterator it1=solute_to_pp.begin(); it1!=solute_to_pp.end(); ++it1) {
        const std::string& solute_name = it1->first;
        const ICParmPair& parm_pairs = it1->second;
        for (ICParmPair::const_iterator it2=parm_pairs.begin(); it2!=parm_pairs.end(); ++it2) {
          const std::string& parameter = it2->first;
          std::string key = ChemistryHelper_Structured::BuildPropertyParameterName(solute_name,"Isotherm",parameter);
          std::map<std::string,int>::const_iterator it3 = aux_chem_variables_map.find(key);
          if (it3 == aux_chem_variables_map.end() && ParallelDescriptor::IOProcessor()) {
            std::cout << "RockManager::RockChemistryProperties Unable to locate parameter in aux_data (material,solute,parameter): "
                      << material_name << "   "
                      << solute_name << "    "
                      << parameter << "  key: " << key << std::endl;
            BoxLib::Abort();
          }
          int comp = it3->second;
          Real value = it2->second;
          fab.setVal(value,comp);
        }
      }
    }
  }

  if (nminerals>0) {
    // Note: mineralogy_ics[rockname][mineralname][property] = val
    ChemICMap::const_iterator it=mineralogy_ics.find(material_name);
    if (it != mineralogy_ics.end())
    {
      const ICLabelParmPair& mineral_to_pp = it->second; 
      for (ICLabelParmPair::const_iterator it1=mineral_to_pp.begin(); it1!=mineral_to_pp.end(); ++it1) {
        const std::string& mineral_name = it1->first;
        const ICParmPair& parm_pairs = it1->second;
        for (ICParmPair::const_iterator it2=parm_pairs.begin(); it2!=parm_pairs.end(); ++it2) {
          const std::string& parameter = it2->first;
          std::string key = ChemistryHelper_Structured::BuildPropertyParameterName(mineral_name,parameter);
          std::map<std::string,int>::const_iterator it3 = aux_chem_variables_map.find(key);
          if (it3 == aux_chem_variables_map.end() && ParallelDescriptor::IOProcessor()) {
            std::cout << "RockManager::RockChemistryProperties  Unable to locate parameter in aux_data (material,mineral,parameter): "
                      << material_name << "   "
                      << mineral_name << "    "
                      << parameter << "  key: " << key << std::endl;
            BoxLib::Abort();
          }
          int comp = it3->second;
          Real value = it2->second;
          fab.setVal(value,comp);
        }
      }
    }
  }     

  if (nsorption_sites>0) {
    // Note: surface_complexation_ics[rockname][sorptionsitename][property]) = val
    ChemICMap::const_iterator it=surface_complexation_ics.find(material_name);
    if (it != surface_complexation_ics.end()) {
      const ICLabelParmPair& sorption_site_to_pp = it->second; 
      for (ICLabelParmPair::const_iterator it1=sorption_site_to_pp.begin(); it1!=sorption_site_to_pp.end(); ++it1) {
        const std::string& sorption_site_name = it1->first;
        const ICParmPair& parm_pairs = it1->second;
        for (ICParmPair::const_iterator it2=parm_pairs.begin(); it2!=parm_pairs.end(); ++it2) {
          const std::string& parameter = it2->first;
          std::string key = ChemistryHelper_Structured::BuildPropertyParameterName(sorption_site_name,"Surface",parameter);
          std::map<std::string,int>::const_iterator it3 = aux_chem_variables_map.find(key);
          if (it3 == aux_chem_variables_map.end() && ParallelDescriptor::IOProcessor()) {
            std::cout << "RockManager::RockChemistryProperties  Unable to locate parameter in aux_data (material,sorption_site,parameter): "
                      << material_name << "   "
                      << sorption_site_name << "    "
                      << parameter << "  key: " << key << std::endl;
            BoxLib::Abort();
          }
          int comp = it3->second;
          Real value = it2->second;
          fab.setVal(value,comp);
        }
      }
    }
  }
 
  if (ncation_exchange>0) {
    // Note: cation_exchange_ics[rockname][property] = val
    ICLabelParmPair::const_iterator it=cation_exchange_ics.find(material_name);
    if (it != cation_exchange_ics.end()) {
      const ICParmPair& parm_pairs = it->second;
      for (ICParmPair::const_iterator it1=parm_pairs.begin(); it1!=parm_pairs.end(); ++it1) {
        const std::string& key = it1->first;
        std::map<std::string,int>::const_iterator it3 = aux_chem_variables_map.find(key);
        if (it3 == aux_chem_variables_map.end() && ParallelDescriptor::IOProcessor()) {
          std::cout << "RockManager::RockChemistryProperties Unable to locate parameter aux_data: " << key << std::endl;
          BoxLib::Abort();
        }
        int comp = it3->second;
        Real value = it1->second;
        fab.setVal(value,comp);
      }
    }
  }
}
