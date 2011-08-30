#ifndef _ROCK_H_
#define _ROCK_H_

#include <string>

#include <Array.H>
#include <FArrayBox.H>
#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <LevelBld.H>
#include <AmrLevel.H>
#include <Region.H>

//
//  Rock class stores rock properties and related function.
//
class Rock
{
public: 

  std::string name;
  Real density;
  Real porosity;
  Array<Real> permeability;
  
  /* The relative permeability types are
     0: none
     2: quadratic 
     3: van Genuchten
     krParam gives the parameters needed to evaluate
     relative permeability.
   */
  int krType;
  Array<Real> krParam;

  /* The capillary pressure types are
     0: none
     3: van Genuchten
     cplParam gives the parameters needed to evaluate
     capillary pressure.
   */
  int cplType;
  Array<Real> cplParam;

  // regions in which we have this rock type and
  // the distributions in each region.
  Array<int> region;

  /* The dist_type are as follows:
     0: read from file
     1: constant
     2: random
     3: sequential gaussian simulation
     *_dist_param gives the parameters needed to generate 
     the above distribtuion.
  */
  Array<int> porosity_dist_type;
  Array<int> permeability_dist_type;
  Array< Array<Real> > porosity_dist_param;
  Array< Array<Real> > permeability_dist_param;

  // parameters for constructing maps
  static int twoexp;
  static int max_level;
  static Array<int>  n_cell, fratio;
  static Array<Real> problo, probhi;
  
  void build_pmap(MultiFab&      mfdata, 
		  Array<Region*>& region_array, 
		  std::string&    gsfile);

  void build_kmap(MultiFab&      mfdata, 
		  Array<Region*>& region_array, 
		  std::string&    gsfile);

  //
  // @ManDoc Set mfdata to constant values given by val.
  //
  void set_constant_val(MultiFab&    mfdata, 
			Region*      region_local, 
			Array<Real>& val);

  //
  // @ManDoc Set mfdata to constant permeability values.
  //
  void set_constant_kval(MultiFab& mfdata, 
			 Region*   region_local);

  //
  // @ManDoc Set mfdata to constant porosity value.
  //
  void set_constant_pval(MultiFab& mfdata, 
			 Region*   region_local);

  //
  // @ManDoc Set fab to parameters relevant to relative permeability.
  //
  void set_constant_krval(FArrayBox& fab,
			  Array<Region*>& region_array, 
			  const Real* dx);
  //
  // @ManDoc Set fab to parameters relevant to capillary pressure.
  //
  void set_constant_cplval(FArrayBox& fab, 
			   Array<Region*>& region_array,
			   const Real* dx);

  static BoxArray build_finest_data(int        max_level, 
				    Array<int>& n_cell,
				    Array<int>& fratio,
				    int        maxBaseGrid = 32);

private:
  
  // ideally these are specific to the regions
  int r_max_level;
  Array<int>  r_n_cell, r_fratio;
  Array<Real> r_problo, r_probhi;
 
};
#endif /*_ROCK_H_*/
