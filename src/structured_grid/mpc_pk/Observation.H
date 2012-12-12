#ifndef _OBSERVATION_H_
#define _OBSERVATION_H_

#include <string>
#include <map>

#include <Array.H>
#include <FArrayBox.H>
#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <LevelBld.H>
#include <AmrLevel.H>
#include <Region.H>

#include <cmath>


/*
  @ManDoc: Class defining all possible types of observations
*/
class Observation
{
public: 
  /* This is the identifier for a particular observation
  */
  std::string name;

  /* This is the string that is sent to the derive function of AmrLevel
  */
  std::string field;

  /* This specifies the region in which the observation is applied. 
   */
  const Region& region;

  /* This specifies the observation type
   */
  std::string obs_type;
  static std::map<std::string, int> obs_type_list;

  Array<Real> times;
  std::map<int,Real> vals;

  std::string event_label;

  Observation(const std::string& name,
              const std::string& field,
              const Region&      region,
              const std::string& obs_type,
              const std::string& eventLabel);

  void process(Real t_old, Real t_new, int iter);
  static void setAmrPtr(Amr* amr) {amrp = amr;}
  static void Cleanup ();
  static void Initialize ();

protected:
  std::pair<Real,Real> integral_and_volume (Real time);
  static Amr* amrp; // Should be const, but cannot be
  Real average(Real time);
  Real volume_integral(Real time);  
  Real volume_time_integral(Real t_old, Real t_new);
  Real point_sample (Real time);

  Real val_old, val_new;  
  bool obs_data_initialized;
};
#endif /*_OBSERVATION_H_*/
