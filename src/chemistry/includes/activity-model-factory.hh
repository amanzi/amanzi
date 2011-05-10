/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __ACTIVITY_MODEL_FACTORY_HH__
#define __ACTIVITY_MODEL_FACTORY_HH__

#include <string>

#include "activity-model.hh"

class ActivityModelFactory 
{
 public:
  ActivityModelFactory();
  ~ActivityModelFactory();

  ActivityModel* Create(const std::string& model);

  static const std::string debye_huckel;
  static const std::string unit;

 protected:

 private:
};

#endif  // __ACTIVITY_MODEL_FACTORY_HPP__

