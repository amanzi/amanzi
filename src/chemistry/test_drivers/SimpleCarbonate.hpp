/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __SimpleCarbonate_hpp__
#define __SimpleCarbonate_hpp__

#include <vector>

#include "Beaker.hpp"

// Driver class for evalating geochemical related processes at a
// single computational node

class SimpleCarbonate : public Beaker {
 public:
  SimpleCarbonate();
  ~SimpleCarbonate();

  void setup(std::vector<double> &total);

 private:
};

#endif  // __SimpleCarbonate_hpp__
