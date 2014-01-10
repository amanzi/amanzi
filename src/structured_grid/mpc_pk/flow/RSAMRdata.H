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
  RSAMRdata(int slev, int nlevs, Layout& layout, PMAmr* amrp, NLScontrol& nlsc);
  virtual ~RSAMRdata();
  virtual void SetUpMemory(NLScontrol& nlsc);
  virtual void ResetRhoSat();
  virtual void SetInflowVelocity(PArray<MFTower>& velocity,
				 Real             time);
  virtual void FillStateBndry (MFTower& press,
                               Real time);
  virtual void calcInvPressure (MFTower&       N,
				const MFTower& P) const;
  virtual void calcLambda (MFTower&       lbd,
			   const MFTower& N);
  virtual void calcRichardAlpha(MFTower&       alpha,
                                const MFTower& rhoSat,
                                Real           time);
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
    MatFiller* matFiller;
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
        eval_time(-1) {}
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
