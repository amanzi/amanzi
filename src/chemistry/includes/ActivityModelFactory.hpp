/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __ACTIVITY_MODEL_FACTORY_HPP__
#define __ACTIVITY_MODEL_FACTORY_HPP__

#include <string>

#include "ActivityModel.hpp"

class ActivityModelFactory 
{
 public:
  ActivityModelFactory();  
  ~ActivityModelFactory();

  ActivityModel* create(std::string model);

  static const string debye_huckel;
  static const string unit;

 protected:

 private:

};

#endif  // __ACTIVITY_MODEL_FACTORY_HPP__

