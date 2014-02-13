#ifndef _RockManager_H_
#define _RockManager_H_

#include <RegionManager.H>
#include <MatFiller.H>

#include <PArray.H>

class RockManager
{
public:
  RockManager(const RegionManager*   _region_manager,
              const Array<Geometry>& geomArray,
              const Array<IntVect>&  refRatio);

  bool SetProperty(Real               time,
                   int                level,
                   MultiFab&          mf,
                   const std::string& pname,
                   int                dComp,
                   int                nGrow,
                   void*              ctx = 0,
                   bool               ignore_mixed = false) {return false;}

  bool SetPropertyDirect(Real               t,
                         int                level,
                         FArrayBox&         fab,
                         const Box&         box,
                         const std::string& pname,
                         int                dComp,
                         void*              ctx = 0) {return false;}

  /* Fill a MultiFab with the ID of the material at that location, given a
     refinement level (as defined in the ctr), including nGrow cells. The ID
     returned is suitable to be passed as an identifier for evaluating
     state-dependent properties (such as relative permeability) */
  void SetMaterialID(int level, iMultiFab& mf, int nGrow, bool ignore_mixed = false) {}

protected:
  void Initialize(const Array<Geometry>& geomArray,
                  const Array<IntVect>&  refRatio);

  const RegionManager* region_manager;
  MatFiller* materialFiller;
  PArray<Material> materials;
  bool is_saturated;
  bool is_diffusive;
  bool tensor_diffusion;
  bool use_shifted_Kr_eval;
  bool saturation_threshold_for_vg_Kr;
};

#endif
