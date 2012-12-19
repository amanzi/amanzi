#ifndef _SATURATIONCURVE_HH_
#define _SATURATIONCURVE_HH_

namespace Amanzi {
class SaturationCurve {

public:
  // S_*, functional form for a multi-phase model, assuming
  // s_1/(s_1-s_2) = S_*(p_c21) = S_*(p_2 - p_1)
  virtual double s_star(double pc) = 0;

  // access to mesh block on which the model is defined
  const int mesh_block() { return meshblock; };

protected:
  void set_mesh_block(const int meshblock_) { meshblock = meshblock_; };
  int meshblock;

};
}
#endif

