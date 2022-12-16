/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

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
#include <Amr.H>
#include <Region.H>

#include <cmath>

class PMAmr;

/*
  @ManDoc: Class defining all possible types of observations
*/
class Observation
{
public:
  typedef std::pair<std::ostream&,std::ostream&> ostream; // first:text, second:binary
  typedef std::pair<std::istream&,std::istream&> istream;

  /* This is the identifier for a particular observation
  */
  std::string name;

  /* This is the string that is sent to the derive function of AmrLevel
  */
  std::string field;

  /* This specifies the name of the region in which the observation is applied.
   */
  std::string region_name;

  /* This specifies the observation type
   */
  std::string obs_type;
  static std::map<std::string, int> obs_type_list;

  Array<Real> times;
  std::map<int,Real> vals;

  std::string event_label;

  Observation(const std::string& name,
              const std::string& field,
              const std::string& region_name,
              const std::string& obs_type,
              const std::string& eventLabel);

  Observation(Observation::istream& is);
  void CheckPoint(Observation::ostream& os) const;

  void process(Real t_old, Real t_new, int iter, int verbose);
  static void setPMAmrPtr(PMAmr* amr) {amrp = amr;}
  static PMAmr* PMAmrPtr() {return amrp;}

  static void Cleanup ();
  static void Initialize ();

protected:
  const Region& GetRegion(const std::string& name) const;

  std::pair<Real,Real> integral_and_volume (Real time);
  static PMAmr* amrp; // Should be const, but cannot be
  Real average(Real time);
  Real volume_integral(Real time);
  Real volume_time_integral(Real t_old, Real t_new);
  Real point_sample (Real time);
  Real peak_sample (Real time);

  Real val_old, val_new;
  bool obs_data_initialized;
};
#endif /*_OBSERVATION_H_*/
