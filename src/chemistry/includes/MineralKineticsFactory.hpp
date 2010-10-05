/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __MINERAL_KINETICS_FACTORY_HPP__

#define __MINERAL_KINETICS_FACTORY_HPP__

/*******************************************************************************
**
**  File Name: MineralKineticsFactory.h
**
**  Description: factory class for reading mineral rates from a file
**  and creating a kinetic rate object.
**
*******************************************************************************/

#include <string>

#include "StringTokenizer.hpp"
#include "Verbosity.hpp"

class KineticRate;

class MineralKineticsFactory
{
 public:
  MineralKineticsFactory(void);
  ~MineralKineticsFactory(void);
  std::vector<KineticRate*> Create(std::string file_name);
 
  void verbosity(const Verbosity s_verbosity) { this->verbosity_ = s_verbosity; };
  Verbosity verbosity(void) const { return this->verbosity_; };

 protected:

 private:
  Verbosity verbosity_;
  static const std::string kTST;

  std::vector<KineticRate*> rates;

  void ReadFile(const std::string file_name);
  void ParseRate(StringTokenizer rate, KineticRate** kinetic_rate);
};

#endif     /* __MINERAL_KINETICS_FACTORY_HPP__ */

