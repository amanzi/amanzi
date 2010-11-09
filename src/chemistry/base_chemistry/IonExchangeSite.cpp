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


IonExchangeSite::IonExchangeSite(const SpeciesName exchanger_name, 
                                 const SpeciesId exchanger_id, 
                                 const double exchanger_charge, 
                                 const std::string exchanger_location,
                                 const double mol_wt, const double size)
    : Species(exchanger_id, exchanger_name, exchanger_charge, mol_wt, size),
      cation_exchange_capacity_(0.0),
      location_(exchanger_location)
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
            << std::setw(15) << location()
            << std::setw(10) << charge()
            << std::setw(10) << cation_exchange_capacity() 
            << std::endl;
} // end Display()

void IonExchangeSite::DisplayResultsHeader(void) const
{
  std::cout << std::setw(15) << "Name" 
            << std::setw(15) << "CEC"
            << std::endl;
} // end DisplayResultsHeader()

void IonExchangeSite::DisplayResults(void) const
{
  std::cout << std::setw(15) << name()
            << std::scientific << std::setprecision(5)
            << std::setw(15) << cation_exchange_capacity()
            << std::endl;
} // end DisplayResults()
