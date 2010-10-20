/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cmath>

#include <iostream>
#include <iomanip>

#include "IonExchangeSite.hpp"

IonExchangeSite::IonExchangeSite() 
    : Species(),
      cation_exchange_capacity_(0.0)
{
} // end IonExchangeSite constructor


IonExchangeSite::IonExchangeSite(SpeciesName exchanger_name, 
                                 SpeciesId exchanger_id, 
                                 double exchanger_charge, 
                                 double mol_wt, double size)
    : Species(exchanger_id, exchanger_name, exchanger_charge, mol_wt, size),
      cation_exchange_capacity_(0.0)
{
} // end IonExchangeSite constructor

IonExchangeSite::~IonExchangeSite() 
{  
} // end IonExchangeSite destructor


void IonExchangeSite::Update(void)
{
} // end Update()


void IonExchangeSite::Display(void) const
{
  std::cout << std::setw(15) << name()
            << std::setw(10) << charge()
            << std::setw(10) << cation_exchange_capacity() 
            << std::endl;
} // end Display()
