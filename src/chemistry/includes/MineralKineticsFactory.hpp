/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __MINERAL_KINETICS_FACTORY_HPP__

#define __MINERAL_KINETICS_FACTORY_HPP__

/*******************************************************************************
**
**  File Name: MineralKineticsFactory.h
**
**  Description: factory class for building a mineral kinetic rate object
**
*******************************************************************************/
#include <vector>
#include <string>

#include "Species.hpp"
#include "Mineral.hpp"
#include "verbosity.hh"
#include "string-tokenizer.hh"

class KineticRate;

class MineralKineticsFactory
{
 public:
  MineralKineticsFactory(void);
  ~MineralKineticsFactory(void);
 
  KineticRate* Create(const std::string& rate_type, 
                      const StringTokenizer& rate_data,
                      const Mineral& mineral,
                      const SpeciesArray& primary_species);

  SpeciesId VerifyMineralName(const std::string mineral_name,
                              const std::vector<Mineral>& minerals) const;


  void set_verbosity(const Verbosity s_verbosity) { this->verbosity_ = s_verbosity; };
  Verbosity verbosity(void) const { return this->verbosity_; };

 protected:

 private:
  Verbosity verbosity_;
  static const std::string kTST;

};

#endif     /* __MINERAL_KINETICS_FACTORY_HPP__ */

