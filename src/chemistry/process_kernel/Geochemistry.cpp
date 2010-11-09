/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "Geochemistry.hpp"
#include "BeakerFactory.hpp"
#include "Beaker.hpp"
#include "Verbosity.hpp"

Geochemistry::Geochemistry()
    : name_("DefaultGeochemistry"),
      verbosity_(kSilent),
      beaker(NULL),
{
} // end Geochemistry() constructor

Geochemistry::~Geochemistry() 
{
  if (beaker != NULL) {
    delete beaker;
  }
} // end Geochemistry destructor

void Geochemistry::Setup(void) 
{
  // what format is the input data....? just in input file, a struct...?
  

  // use a beaker factory to create the appropriate type of beaker
  // depending on the thermo database file format?

  beaker = new ThermoDatabase();
  beaker->verbosity(verbosity);

  // extract the parameters from the input and set them here
  this->parameters = beaker->GetDefaultParameters();
  this->parameters.XXX = YYY;

  // size component vectors here, assign an initial value?
  this->components.primaries.clear();
  this->components.minerals.clear();
  this->components.ion_exchange_sites.clear();

  this->components.primaries.resize(XXX);
  this->components.minerals.resize(XXX);
  this->components.ion_exchange_sites.resize(XXX);

  // now setup the beaker
  beaker->Setup(components, parameters);

  // initial speciation to make sure everything is working
  beaker->Speciate(components, parameters);
} // end Setup()


void Geochemistry::Advance(void)
{
  // unpack the data here
  // some loop to extract current data from the entire domain data
  // assign the data to the appropriate location in parameters and components
  for (point in strang_data_format) {
    parameters.porosity = point.XXX;
    for (p in strange_data_format) {
      components.primaries[p] = XXX;
    }
    for (mineral in strange_data_format) {
      components.minerals[m] = XXX;
    }
  }
  
  beaker->ReactionStep(components, parameters, delta_time);

  // package the results appropriately

}  // end Advance()


void Geochemistry::CommitState(void);
{
  // save the input values as the internal state of the system....

}  // end CommitState()
