/* -*-  mode: c++; c-default-style: "google-c-style"; indent-tabs-mode: nil -*- */
#ifndef __MINERAL_KINETICS_CREATOR_HPP__

#define __MINERAL_KINETICS_CREATOR_HPP__

/*******************************************************************************
**
**  File Name: MineralKineticsCreator.h
**
**  Description: factory class for reading mineral rates from a file
**  and creating a kinetic rate object.
**
*******************************************************************************/

#include <string>

#include "StringTokenizer.hpp"

class MineralKineticsCreator
{

public:

  MineralKineticsCreator(void);
  ~MineralKineticsCreator(void);

  void readFile(const std::string file_name);

  void verbosity(const int s_verbosity) { this->verbosity_ = s_verbosity; };
  int verbosity(void) const { return this->verbosity_; };

protected:

private:
  int verbosity_;

  void ParseReaction(const::std::string rxn_string);
  void ParseRate(StringTokenizer rate);
  void ParseTstParameters(StringTokenizer rate);

};


#endif     /* __MINERAL_KINETICS_CREATOR_HPP__ */

