/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _RockManager_H_
#define _RockManager_H_

#include <RegionManager.H>
#include <MatFiller.H>
#include <MonotCubicInterpolator.H>
#include <ChemistryHelper_Structured.H>
#include <RegionData.H>

#include <PArray.H>
#include <DataServices.H>

/*
  RockManager

  A rock manager is a wrapper around the functions used to compute properties of the
  porous media. The most significant data member is a MatFiller object, which is used
  to store user input data (region definition, parameters of material properties, etc).
  A MatFiller however was only really designed to manage properties with parameters
  with a specified time variation.  This class can additionally evaluate a select
  set of properties that are a function of the "state".  Construction of a RockManager
  object requires a description of the computational domain (including AMR refinement
  layout, and a description of the regions that define material boundaries (via
  a RegionManager object). In the constructor, "materials" are defined over subsets of
  the regions, via domain-specific parameterizations documented elsewhere.  The
  parameters of each thusly-defined property are available via the "GetProperty"
  memebers.  The set of derived properties (ones that require additonal data for
  evaluation) currently include:

  1. Capillary Pressure, as a function of saturation
  2. Saturation, as a function of capillary pressure
  3. The derivative of 2, as a function of saturation
  4. Mobility (relative permeability divided by fluid velocity)
  5. The derivative of mobility, as a function of saturation

  Material properties are available on sets of points, given a corresponding
  set of ID values that select the associated material ID.  Alternatively
  properties can be computed in the context of the AMR hierarchy, where an integer
  refinement level is used to disambiguate the index space of the box arrays requested
  to be filled.

  Also, the class returns the material ID in structures associated with the AMR
  hierarchy.  The intention is that the work to assign a material index to a large
  set of cells can be done infrequently and reused.  It also provides a mechanism for
  use of this structure in multi-language (e.g. C++/Fortran) situations.
 */

class RockManager
{
public:
  typedef std::map<std::string,Real> ICParmPair; // ic parameter and value
  typedef std::map<std::string, ICParmPair > ICLabelParmPair; // ic parameter and value associated with a labeled entity
  typedef std::map<std::string, ICLabelParmPair> ChemICMap; //
  typedef std::map<std::string, std::map<std::string,int> > LabelIdx;

  RockManager(const RegionManager*      _region_manager,
	      const Array<std::string>* _solute_names = 0,
	      Real                      _liquid_density = 0,
	      Real                      _liquid_viscosity = 0,
	      Real                      _gravity = 0);

  ~RockManager();

  void FinalizeBuild(const Array<Geometry>& geomArray,
                     const Array<IntVect>&  refRatio,
                     int                    nGrow,
		     bool                   restart);

  bool CanDerive(const std::string& property_name) const;

  bool GetProperty(Real               time,
                   int                level,
                   MultiFab&          mf,
                   const std::string& pname,
                   int                dComp,
                   int                nGrow,
                   void*              ctx = 0,
                   bool               ignore_mixed = false) const;

  /* Fill a MultiFab with the ID of the material at that location, given a
     refinement level (as defined in the ctr), including nGrow cells. The ID
     returned is suitable to be passed as an identifier for evaluating
     state-dependent properties (such as relative permeability) */
  void GetMaterialID(int level, iMultiFab& mf, int nGrow, bool ignore_mixed = false) const;

  void RockChemistryProperties(FArrayBox&  fab,
                               const Real* dx,
                               const std::map<std::string,int>& aux_chem_variables_map);

  void RockChemistryProperties(FArrayBox&         fab,
                               const std::string& material_name,
                               const std::map<std::string,int>& aux_chem_variables_map);

  // Capillary Pressure (given Saturation)
  void CapillaryPressure(const Real* saturation, int* matID, Real time, Real* capillaryPressure, int npts) const;
  void CapillaryPressure(const MultiFab&  saturation,
                         const iMultiFab& matID,
                         Real             time,
                         MultiFab&        capillaryPressure,
                         int              sComp,
                         int              dComp,
                         int              nGrow) const;

  // Inverse Capillary Pressure (Saturation given Capillary Pressure)
  void InverseCapillaryPressure(const Real* capillaryPressure, int* matID, Real time, Real* saturation, int npts) const;
  void InverseCapillaryPressure(const MultiFab&  capillaryPressure,
                                const iMultiFab& matID,
                                Real             time,
                                MultiFab&        saturation,
                                int              sComp,
                                int              dComp,
                                int              nGrow) const;

  // D (Saturation) / D(CapillaryPressure)
  void DInverseCapillaryPressure(const Real* saturation, int* matID, Real time, Real* DsaturationDpressure, int npts) const;
  void DInverseCapillaryPressure(const MultiFab&  saturation,
                                 const iMultiFab& matID,
                                 Real             time,
                                 MultiFab&        DsaturationDpressure,
                                 int              sComp,
                                 int              dComp,
                                 int              nGrow) const;

  // Relative Permeability (given saturation)
  void RelativePermeability(const Real* saturation, int* matID, Real time, Real* kappa, int npts) const;
  void RelativePermeability(const MultiFab&  saturation,
                            const iMultiFab& matID,
                            Real             time,
                            MultiFab&        kappa,
                            int              sComp,
                            int              dComp,
                            int              nGrow) const;

  // Residual saturation
  void ResidualSaturation(int* matID, Real time, Real* Sr, int npts) const;
  void ResidualSaturation(const iMultiFab& matID,
                          Real             time,
                          MultiFab&        Sr,
                          int              dComp,
                          int              nGrow) const;

  const MatFiller* GetMatFiller() const {return materialFiller;}
  int ID() const {return rock_mgr_ID;}
  int NComp(const std::string& property_name) const;
  int FinestLevel() const {return materialFiller->NumLevels()-1;}
  const PArray<Material>& GetMaterials() const {return rock;}
  const Material& GetMaterial(int i) const {return rock[i];}
  const Material& FindMaterial(const std::string& name) const;
  const Material& FindMaterialInRegions(const Array<std::string>& region_names) const;

  bool Is_CP_model_XX(int model_id, const std::string& str) const;
  bool Is_Kr_model_XX(int model_id, const std::string& str) const;

  int NumMinerals() const {return nminerals;}
  const Array<std::string>& MineralNames() const {return minerals;}
  int NumSorptionSites() const {return nsorption_sites;}
  const Array<std::string>& SorptionSiteNames() const {return sorption_sites;}
  int NumCationExchange() const {return ncation_exchange;}
  int NumSorptionIsotherms() const {return nsorption_isotherms;}
  const std::map<std::string,int> AuxChemVariablesMap() const {return aux_chem_variables;}
  bool UsingSorption() const {return using_sorption;}

  const ChemICMap& SorptionIsothermICs() const {return sorption_isotherm_ics;}  // isotherm[rockname][solute][property] = val
  const ChemICMap& MineralogyICs() const {return mineralogy_ics;}           // mineralogy[rockname][mineralsname][property] = val
  const ChemICMap& SurfaceComplexationICs() const {return surface_complexation_ics;} // sc[rockname][sorptionsitename][property]) = val
  const ICLabelParmPair& CationExchangeCapacityICs() const {return cation_exchange_ics;} // ceics[rockname][property] = val

  const LabelIdx& MineralogyLabelMap() const {return mineralogy_label_map;}
  const LabelIdx& SorptionIsothermLabelMap() const {return sorption_isotherm_label_map;}
  const LabelIdx& SurfaceComplexationLabelMap() const {return surface_complexation_label_map;}
  const LabelIdx& SorptionChemistryLabelMap() const {return sorption_chem_label_map;} // sc[solute][parameter] = val
  const std::map<std::string,int>& CationExchangeLabelMap() const {return cation_exchange_label_map;}

  bool DoDiffusion() const {return do_diffusion;}
  bool DoTensorDiffusion() const {return do_tensor_diffusion;}

  const Array<std::string> MaterialNames() const;
  const MonotCubicInterpolator& CPInterpolator(const std::string& name) const;
  const MonotCubicInterpolator& KrInterpolator(const std::string& name) const;

protected:
  void Initialize(const Array<std::string>* soluteNames);
  int FindMaterialNum(const std::string& name) const;

  Array<Array<int> > SortPtsByMaterial(int* matID, int npts) const;
  void BuildInterpolators();
  void FillBoundary(Real      time,
                    int       level,
                    MultiFab& mf,
                    int       dComp,
                    int       nComp,
                    int       nGrow) const;

  const RegionManager* region_manager;
  Real liquid_density, liquid_viscosity, gravity;
  MatFiller* materialFiller;
  PArray<Material> rock;
  bool is_saturated;
  bool do_diffusion, do_tensor_diffusion;
  std::map<std::string,int> CP_models;
  std::map<std::string,int> Kr_models;
  PArray<MonotCubicInterpolator> CP_s_interps;
  PArray<MonotCubicInterpolator> Kr_s_interps;
  Array<std::pair<int,std::string> > WRM_plot_file;
  int max_grid_size_fine_gen, ngrow_fine_gen;
  bool finalized;

  int rock_mgr_ID;

  // Minerals/chemistry
  int nminerals;
  Array<std::string> minerals;
  int nsorption_sites;
  Array<std::string> sorption_sites;
  int ncation_exchange;
  int nsorption_isotherms;
  std::map<std::string,int> aux_chem_variables;
  bool using_sorption;

  ChemICMap sorption_isotherm_ics;    // isotherm[rockname][solute][property] = val
  ChemICMap mineralogy_ics;           // mineralogy[rockname][mineralsname][property] = val
  ChemICMap surface_complexation_ics; // sc[rockname][sorptionsitename][property]) = val
  ICLabelParmPair cation_exchange_ics; // ceics[rockname][parameter] = val
  Array<Real> tInflow;

  LabelIdx mineralogy_label_map;
  LabelIdx sorption_isotherm_label_map;
  LabelIdx surface_complexation_label_map;
  LabelIdx sorption_chem_label_map; // sc[solute][parameter] = val
  std::map<std::string,int> cation_exchange_label_map;

  Array<std::string> known_solutes;
};

#endif
