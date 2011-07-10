#ifndef _OBSERVATION_H_
#define _OBSERVATION_H_

#include <string>

#include <Array.H>
#include <FArrayBox.H>
#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <LevelBld.H>
#include <AmrLevel.H>
#include <Region.H>

/*
  @ManDoc: Class defining all possible types of observations
*/
class Observation
{
public: 
  /* This is the identifier for a particular observation
  */
  std::string name;

  /* These specify which variable the source is contributing to. 
       var_type: comp or tracer.
       var_id  : name listed in the comp or tracer. 
                 ALL: all of the components in var_type.
		 USER: user specified routine in PROB_$D.F
       id: based on var_type, var_id, this is converted to index in State_Type.
  */
  std::string var_type;
  std::string var_id;
  Array<int> id;

  /* This specifies the region in which the observation is applied. 
     The user give the name but this is internally change 
     to corresponding index in region_array.
   */
  int region;

  /* This specifies the observation type
     0: production
     1: mass fraction
     2: mole fraction
     3: saturation
   */
  int obs_type;

  /* This specifies the operation type
     0: average
     1: integral over space
     2: integral over space and time
   */
  int op_type;

  Array<Real> times;
  Array<Real> vals;

  void process(Real t_old, Real t_new);
  static void setAmrPtr(Amr* amr) {amrp = amr;}

protected:
  std::pair<Real,Real> integral_and_volume (Real time);
  static Amr* amrp; // Should be const, but cannot be
  Real average(Real time);
  Real volume_integral(Real time);  
  Real volume_time_integral(Real t_old, Real t_new);
};
#endif /*_OBSERVATION_H_*/
