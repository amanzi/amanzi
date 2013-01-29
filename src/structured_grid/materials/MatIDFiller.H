#ifndef _MatIDFiller_H_
#define _MatIDFiller_H_

#include <Geometry.H>
#include <BoxArray.H>
#include <PArray.H>
#include <VisMF.H>
#include <Utility.H>

#include <Region.H>
#include <Material.H>

class MatIDFiller
{
public:
  MatIDFiller();

  MatIDFiller(const Array<Geometry>& _geomArray,
              const Array<IntVect>&  _refRatio,
              const PArray<Material>& materials);

  void define(const Array<Geometry>& _geomArray,
              const Array<IntVect>&  _refRatio,
              const PArray<Material>& materials);

  bool Initialized() const {return initialized;}

  const Geometry& Geom(int level) const {return geomArray[level];}
  const IntVect& RefRatio(int crse_level) {return refRatio[crse_level];}
  void SetMaterialID(int level, MultiFab& mf, int nGrow);

  const BoxArray& Mixed(int level) const;
  MultiFab& MaterialID(int level) {return materialID[level];}
  const std::map<std::string,int>& MatIdx() const {return matIdx;}
  const int NumLevels() const {return nLevs;}

  void Initialize();

protected:
  Array<BoxArray> FindMixedCells();

  Array<BoxArray> ba_mixed;
  PArray<MultiFab> materialID;

  Array<Geometry> geomArray;
  Array<IntVect> refRatio;
  int nLevs;
  PArray<Material> materials;
  std::map<std::string,int> matIdx;
  Array<std::string> matNames;
  bool initialized;
};

#endif
