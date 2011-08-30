#include <winstd.H>

#include "PorousMedia.H"

Amr* Observation::amrp = 0;

Observation::Observation()
{
  times_idx = 0;
  val_old = 0.;

  obs_type_list["average"]          = 0;
  obs_type_list["integral"]         = 1;
  obs_type_list["time_integral"]    = 2;
  obs_type_list["squared_integral"] = 3;
  obs_type_list["flux"]             = 4;
}  

void 
Observation::process(Real t_old, 
		     Real t_new)
{
  // Must set the amr pointer prior to use via Observation::setAmrPtr
  BL_ASSERT(amrp); 

  int nObs = vals.size();

  // determine observation at time t_new
  switch (obs_type_list[obs_type] )
    {
    case 0:
      val_new  = average(t_new);

      break;
      
    case 1:
      val_new = volume_integral(t_new);
      break;
      
    case 2:
      // times[nObs] = 0.5*(t_old+t_new);
      val_new = volume_time_integral(t_old,t_new);
      break;

    default:
      ;// Do nothing
    }

  // determine with the observation is requested at this time step
  switch (obs_type_list[obs_type] )
    {
    case 2:
      if (vals.size() < 1) 
	{
	  vals.resize(1);
	  vals[0] = 0.;
	}
      if (times[0] <= t_old && times[1] >= t_new)
	vals[0] = vals[0] + val_new;   
      break;

    default:
      for (int i=0; i<times.size(); i++)
	{
	  if (times[i] >= t_old && times[i] <= t_new)
	    {
	      vals.resize(nObs+1);
	      vals[nObs] = val_old + (times[i]-t_old)/(t_new-t_old)*(val_new-val_old);
	    }
	}
    }
  val_old = val_new;
}

std::pair<Real,Real>
Observation::integral_and_volume (Real time)
{
  const int finest_level = amrp->finestLevel();
  const Array<IntVect>& refRatio = amrp->refRatio();

  Real int_inside = 0.0;
  Real vol_inside = 0;  
  
  Real vol_scale_lev = 1;

  for (int lev = 0; lev <= finest_level; lev++)
    {
      if (lev>0)
        {
          for (int d=0; d<BL_SPACEDIM; ++d)
            {
              vol_scale_lev *= 1./refRatio[lev][d];
            }
        }

      const MultiFab& S_new = amrp->getLevel(lev).get_new_data(State_Type);
      BoxArray baf;
      if (lev < finest_level)
        {
          BoxArray baf = amrp->boxArray(lev+1);
          baf.coarsen(refRatio[lev]);
        }

      FArrayBox fabVOL, fabINT;
      const Real* dx = amrp->Geom(lev).CellSize();
      Real vol = 1.;
      for (int i=0;i<BL_SPACEDIM;i++) vol *= dx[i];

      for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
          const Box& cbox = mfi.validbox();
          fabVOL.resize(cbox,1);

          // Initialize to zero everywhere, then set to 1 inside region
          fabVOL.setVal(0);
          PorousMedia::region_array[region]->setVal(fabVOL,vol,0,dx,0);

          // Zero where covered with better data
          if (lev < finest_level)
            {
              std::vector< std::pair<int,Box> > isects = baf.intersections(cbox);
              
              for (int ii = 0, N = isects.size(); ii < N; ii++)
                {
                  fabVOL.setVal(0,isects[ii].second,0,fabVOL.nComp());
                }
            }

          // Now increment volume and integral inside
          vol_inside += fabVOL.sum(0) * vol_scale_lev;
          
          fabINT.resize(cbox,1);
          fabINT.copy(S_new[mfi],id,0,1);
          fabINT.mult(fabVOL);
          int_inside += fabINT.sum(0);
        }
    }
  ParallelDescriptor::ReduceRealSum(vol_inside);
  ParallelDescriptor::ReduceRealSum(int_inside);
  return std::pair<Real,Real>(int_inside,vol_inside);
}

Real
Observation::average (Real time)
{
  std::pair<Real,Real> int_vol = integral_and_volume(time);
  return (int_vol.second == 0 ? 0 : int_vol.first / int_vol.second);
}

Real
Observation::volume_integral (Real time)
{
  return integral_and_volume(time).first;
}

Real
Observation::volume_time_integral (Real t_old, Real t_new)
{
  // Use trapezoid rule in time on volume integrals.
  std::pair<Real,Real> int_vol_0 = integral_and_volume(t_old);
  std::pair<Real,Real> int_vol_1 = integral_and_volume(t_new);

  return 0.5*(t_new-t_old)*(int_vol_0.first + int_vol_1.first);
}
