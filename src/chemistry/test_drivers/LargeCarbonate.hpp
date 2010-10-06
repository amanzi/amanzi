/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __LargeCarbonate_hpp__
#define __LargeCarbonate_hpp__

#include <vector>

#include "Beaker.hpp"

// Driver class for evalating geochemical related processes at a
// single computational node

class LargeCarbonate : public Beaker {
 public:
  LargeCarbonate();
  ~LargeCarbonate();

  void setup(std::vector<double> &total, const std::string mineral_kinetics_file);


 private:
};

#endif  // __LargeCarbonate_hpp__
