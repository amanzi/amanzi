/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Description: Base class for reactions (aqueous equilibrium
** complexes, minerals)
*/

#include <iostream>
#include <iomanip>

#include "Species.hpp"
#include "Reaction.hpp"
#include "Block.hpp"

Reaction::Reaction()
    : name_("BaseReactionClass")
{
} // end Reaction() constructor

Reaction::Reaction(const std::string in_name)
    : name_(in_name)
{
} // end Reaction() constructor


Reaction::~Reaction() 
{
} // end Reaction() destructor

/*
**
**  required interface functions
**
*/
void Reaction::Update(const std::vector<Species> primary_species) 
{
  static_cast<void>(primary_species);
} // end Update()

void Reaction::AddContributionToTotal(std::vector<double>& total) 
{
  static_cast<void>(total);
} // end AddContributionToTotal()

void Reaction::AddContributionToDTotal(const std::vector<Species> primary_species,
                                       Block* dtotal) 
{
  static_cast<void>(primary_species);
  static_cast<void>(dtotal);
} // end AddContributionToDTotal()

void Reaction::Display(void) const
{
  std::cout << "  Reaction Type: " << name() << std::endl;
} // end Display()
