#ifndef _MatFiller_H_
#define _MatFiller_H_

#include <Geometry.H>
#include <BoxArray.H>
#include <PArray.H>
#include <VisMF.H>
#include <Utility.H>

#include <Region.H>
#include <Material.H>

/*
  A data structure used to manage the evaluation of material properties for a block-structured
  adaptive grid application.  A problem is defined by the physical domain and extent of its
  consituent regions, and its discretization, given by the grid spacing (constant at each
  "level"), and a set of ratios that define the refinement relative to each successive level.
  Finally, a rule is required to assign material properties to a cell--here, implementation of
  that rule is deferred to the innerds of the Material class, through its implementation of the
  "setVal" function.

  The primary members of this class are "SetProperty" (and "SetMaterialID", if that information
  will be used by other procedures to set material-based information).  The SetProperty
  function simply evaluates a property over the cells requested.  However, if material boundaries
  do not coincide with the interfaces of the discrete cells at the coarsest level, there is 
  some ambiguity for how to do this.  This class implements the rule that a the material is
  defined by the material ID determined at the finest grid resolution desired.  If, when a 
  cell on a coarser level is refined to that level, its constituent cells share the same 
  material ID, we store that ID in a data structure at the coarse level. The SetMaterialID 
  function is designed to return data according to this rule, with the understanding that
  a value of -1 indicates that the cell is "mixed" at that level.  Any of the properties
  can be evaluated directly in non-mixed cells.  For mixed cells, the properties are evaluated
  at the level required to uniquify the material.  The material property is then averaged
  to the desired level of refinement and returned.

  The purpose of this class is to attempt to manage this process in an "optimal" way.  While
  correct definition requires an initial evaluation of the material ID on every fine cell in
  the domain, this is only required once.  After this, that information is cached in an
  efficient manner consistent with the block-structured AMR approach to data compression.

  Finally, a few of the members of this class are virtual so that alternative rules may
  be implemented for the various key pieces of the algorithm (evualate-then-avaerage, etc),
  although default behavior is provided, as discussed in the function descriptions below).  
  Also, the Property class provides for the passing of a "user-defined" context in order to
  generalize the notion of a "Property".  This class respects that generality, but defaults
  the context to a null object pointer by default.
 */
class MatFiller
{
public:

  /* Define the problem of interest */
  MatFiller(const Array<Geometry>&  geomArray,
            const Array<IntVect>&   refRatio,
            const PArray<Material>& materials);

  /* Fill a MultiFab with the value of a property, given a time, 
     and refinement level (as defined in the ctr), including
     nGrow cells.  ctx is a user-defined context which is passed
     through to the Property "eval" function. */
  void SetProperty(Real               time,
                   int                level,
                   MultiFab&          mf,
                   const std::string& pname,
                   int                dComp,
                   int                nGrow,
                   void*              ctx = 0);

  /* Fill a MultiFab with the ID (a Real-casted integer) of the 
     material at that location, given a refinement level (as defined 
     in the ctr), including nGrow cells. The ID is consistent 
     with MatIdx()*/
  void SetMaterialID(int level, MultiFab& mf, int nGrow);


  /* Access functions follow */

  /* The coarsen rule registered for a given property */
  Property::CoarsenRule coarsenRule (const std::string& property_name) const;

  /* The refine rule registered for a given property */
  Property::RefineRule refineRule (const std::string& property_name) const;

  /* The number of components this property will fill */
  int nComp(const std::string& property_name) const;

  /* The Geometry object at this level  */
  const Geometry& Geom(int level) const {return geomArray[level];}

  /* The refinement ratio between this level and the next finer. 
     Note, level<NumLevels()-1 */
  const IntVect& RefRatio(int crse_level) {return refRatio[crse_level];}

  /* The set of cells that are "mixed" at this level....i.e., when cells in this set are
     refined to the next finer level, they will be composed of finer cells which do not
     all share the same material ID */
  const BoxArray& Mixed(int level) const;

  /* The internally cached (minimal) data structure at this level containing the 
     (Real-casted integer) ID of the material at this location.  The ID is consistent 
     with MatIdx(). */
  MultiFab& MaterialID(int level) {return materialID[level];}

  /* A map of material name to integer used to identify Materials in this class */
  const std::map<std::string,int>& MatIdx() const {return matIdx;}

  /* The number of AMR levels used to uniquely determine the material ID of the cells. */
  const int NumLevels() const {return nLevs;}

protected:

  /* Currently, not really critical, since there is no default constructor */
  bool Initialized() const {return initialized;}

  /* Manage all actions necessary to define this structure */
  void Initialize();

  /* Manage the cell-ID parts of the construction process */
  void define(const Array<Geometry>& _geomArray,
              const Array<IntVect>&  _refRatio,
              const PArray<Material>& materials);

  /* Given fine property values, generate coarse property values
     per the coarsen rule registered for this property. */
  virtual void CoarsenData(const FArrayBox&   fineFab,
                           int                sComp,
                           FArrayBox&         crseFab,
                           const Box&         crseBox,
                           int                dComp,
                           int                nComp,
                           const IntVect&     ref,
                           const std::string& pname) const;

  /* Given crse property values, generate fine property values 
     per the refine rule registered for this property */
  virtual void RefineData(const FArrayBox&   crseFab,
                          int                sComp,
                          FArrayBox&         fineFab,
                          const Box&         crseBox,
                          int                dComp,
                          int                nComp,
                          const IntVect&     ref,
                          const std::string& pname) const;

  /* Manage the generation of property value based on evaluations
     at finer levels.  This is factored out to allow permuting the
     order of evaluation and averaging.  The default implementation
     evaluates the property at the finer level, then calls 
     CoarsenData to generate coarsened values. */
  virtual void FillCoarseCells(Real               t,
			       int                level,
			       MultiFab&          mfc,
			       const std::string& pname,
			       int                dComp,
			       int                nComp,
			       void*              ctx = 0);

  /* Fill all cells that are outside the "domain" as defined by
     the Geometry objects passed into the ctr.  By default, this
     is done by piecewise-constant extrapolation normal to the 
     boundary interface. */
  virtual void FillCellsOutsideDomain(Real               t,
				      int                level,
				      MultiFab&          mf,
				      const std::string& pname,
				      int                dComp,
				      int                nComp,
				      int                nGrow);

  /* Ensure that the properties passed in for all materials are consistent,
     e.g. that they return a consistent number of components, and have 
     consistent inter-level transfer instructions. */
  void VerifyProperties();

  /* For all levels necessary, find the set of cells that are "mixed",
     in that if refined to the finest level, the material ID of 
     the consituent cells will not all be the same */
  Array<BoxArray> FindMixedCells();

  /* The array of mixed cells, per refinement level */
  Array<BoxArray> ba_mixed;

  /* The material ID of cells at each level.  This is the minimal set
     necessary to represent each level.  The ID is consistent with 
     MatIdx(), except that mixed cells at this level will have the value
     of -1.
  */
  PArray<MultiFab> materialID;

  /* Complete copies of the objects passed in on construction */
  Array<Geometry> geomArray;
  Array<IntVect> refRatio;
  PArray<Material> materials;

  /* The number of levels required to uniquely determine cell IDs */
  int nLevs;

  /* The map used to identify materials */
  std::map<std::string,int> matIdx;

  /* A list of registered material names */
  Array<std::string> matNames;

  /* Has this object been initialized */
  bool initialized;

  /* Number of components returned by each property */
  std::map<std::string,int> property_nComps;

  /* The coarsen rule for each property */
  std::map<std::string,Property::CoarsenRule> property_cRules;

  /* The refine rule for each property */
  std::map<std::string,Property::RefineRule> property_rRules;
};

#endif

