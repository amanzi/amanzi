#include "vanGenuchtenModel.hh"
#include "math.h"
#include <iostream>

namespace Amanzi {
vanGenuchtenModel::vanGenuchtenModel(int meshblock, double lambda, double P0) :
  lambda_(lambda), P0_(P0) {
  nu_ = 1.0/(1.0-lambda_);
  set_mesh_block(meshblock_);
}

double s_star(double pc) {
  return pow(1.0 + pow(pc/P0_, nu_), -lambda_);
};
}
