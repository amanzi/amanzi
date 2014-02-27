#ifndef _RockManager_H_
#define _RockManager_H_

#include <RegionManager.H>
#include <MatFiller.H>
#include <MonotCubicInterpolator.H>

#include <PArray.H>

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
  RockManager(const RegionManager*   _region_manager,
              const Array<Geometry>& geomArray,
              const Array<IntVect>&  refRatio);
  ~RockManager();

  bool CanDerive(const std::string& property_name) const;

  bool GetProperty(Real               time,
                   int                level,
                   MultiFab&          mf,
                   const std::string& pname,
                   int                dComp,
                   int                nGrow,
                   void*              ctx = 0,
                   bool               ignore_mixed = false) const;

  bool GetPropertyDirect(Real               t,
                         int                level,
                         FArrayBox&         fab,
                         const Box&         box,
                         const std::string& pname,
                         int                dComp,
                         void*              ctx = 0) const;

  /* Fill a MultiFab with the ID of the material at that location, given a
     refinement level (as defined in the ctr), including nGrow cells. The ID
     returned is suitable to be passed as an identifier for evaluating
     state-dependent properties (such as relative permeability) */
  void GetMaterialID(int level, iMultiFab& mf, int nGrow, bool ignore_mixed = false) const;


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
  /*
    Number of components for this property
   */
  int NComp(const std::string& property_name) const;
  const PArray<Material>& GetMaterials() const {return rock;}

  bool Is_CP_model_XX(int model_id, const std::string& str) const;
  bool Is_Kr_model_XX(int model_id, const std::string& str) const;

protected:
  void Initialize(const Array<Geometry>& geomArray,
                  const Array<IntVect>&  refRatio);

  Array<Array<int> > SortPtsByMaterial(int* matID, int npts) const;
  void BuildInterpolators();

  const RegionManager* region_manager;
  MatFiller* materialFiller;
  PArray<Material> rock;
  bool is_saturated;
  bool is_diffusive;
  bool tensor_diffusion;
  std::map<std::string,int> CP_models;
  std::map<std::string,int> Kr_models;
  PArray<MonotCubicInterpolator> CP_s_interps;
  PArray<MonotCubicInterpolator> Kr_s_interps;
  Array<std::pair<int,std::string> > WRM_plot_file;
  bool interps_built;

  int rock_mgr_ID;
};

#endif
