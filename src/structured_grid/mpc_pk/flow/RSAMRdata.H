/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _RSAMRdata_H_
#define _RSAMRdata_H_

#include <PorousMedia.H>
#include <RSdata.H>

struct RSAMRdata
  : public RSdata
{
  /*
    An RSdata that manages memory and data evals via a PMAmr structure
   */
  RSAMRdata(int slev, int nlevs, Layout& layout, PMAmr* amrp, NLScontrol& nlsc, const RockManager* rm);
  virtual ~RSAMRdata();
  virtual void SetUpMemory(NLScontrol& nlsc);
  virtual void ResetRhoSat();
  virtual void SetInflowVelocity(PArray<MFTower>& velocity,
				 Real             time);
  virtual void FillStateBndry (MFTower& press,
                               Real time);
  virtual void calcInvPressure (MFTower&       N,
                                const MFTower& P,
                                Real           time,
                                int            sComp,
                                int            dComp,
                                int            nGrow) const;

  virtual void calcLambda (MFTower&       Lambda,
                           const MFTower& N,
                           Real           time,
                           int            sComp,
                           int            dComp,
                           int            nGrow) const;

  virtual void calcRichardAlpha (MFTower&       Alpha,
                                 const MFTower& N,
                                 Real           time,
                                 int            sComp,
                                 int            dComp,
                                 int            nGrow) const;

  virtual Array<int>& rinflowBCLo();
  virtual Array<int>& rinflowBCHi();
  virtual void SetDensity();
  virtual void SetGravity();
  virtual void SetViscosity();
  virtual void SetIsSaturated();

  enum RSdata_Property {
    RSdata_KappaCCdir = 0,
    RSdata_Num_Properties = 1
  };

  struct PropertyManagerCtx {
    const RockManager* rockMgr;
    std::string property_name;
  };

  struct PropertyManager
  {
    PropertyManager(MFTower                           *pdat_ptr,
                    const std::set<PropertyManager*>& dependents,
                    PropertyManagerCtx                ctx)
      : property_dataPtr(pdat_ptr),
        property_dependents(dependents),
        property_ctx(ctx),
        eval_time(-666) {}
    virtual const MFTower* UpdateProperty(Real t);
    virtual void UpdateDependents(Real t)
      {
        for (std::set<PropertyManager*>::const_iterator it=property_dependents.begin(),
               End=property_dependents.end(); it!=End; ++it) {
          (*it)->UpdateProperty(t);
        }
      }
    MFTower* property_dataPtr;
    std::set<PropertyManager*> property_dependents;
    PropertyManagerCtx property_ctx;
    Real eval_time;
  };

  virtual const MFTower* GetKappaCCdir(Real t);
  virtual const PArray<MFTower>& GetKappaEC(Real t);
  virtual const MFTower* GetSource(Real t);

  PorousMedia& GetPMlevel(int level) {return pm[level];}
  PMAmr* pm_amr;
  PArray<PorousMedia> pm;
  std::map<RSdata_Property,PropertyManager*> managed_properties;
  Real eval_time_for_source;
};

#endif
