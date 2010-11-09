/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __Geochemistry_hpp__
#define __Geochemistry_hpp__

/*

What does this family of classes need to do:

This is basically the chemistry process kernel that was planned for amanzi. 

Geochemistry is where the Beaker class will live and all setup, etc
takes place. The interface inherits from Geochemistry and translates
from the external data structures to the internal ones.

Setup() : create a beaker object, initialized with the desired
chemistry, and perform the initial speciation to verify that
everything is working correctly.

Advance() : receive the new state from the driver routine, call reaction step etc.

Notes: 

- my gut feeling is that a template class is best for this, but need
to see how it goes.

- Do we still need the Chemistry_State object that was origionally planned?

- what do we do about memory for the (I assume) massive arrays we are getting. we can't copy them....

- would it be easier if every calling program wrote it's own wrapper 

*/

#include "Beaker.hpp"

#include <vector>

// chemistry process kernel

class Geochemistry {
 public:
  ~Geochemistry();

  virtual void Setup(void);

  // can't make these abstract unless it is a template because we don't know what the input vector type is?
  virtual void Advance() = 0;
  virtual void CommitState() = 0;

  // having something like this is going to be annoying:
  /*
  virtual void Advance(Petsc data);
  virtual void Advance(std::vector data);
  virtual void Advance(Teuchos::RPC<Epetra_MultiVector>);
  virtual void Advance(YourFavoriteFlavorVector data);
  */

  virtual void set_verbosity(const Verbosity s_verbosity) { this->verbosity_ = s_verbosity; };
  virtual Verbosity verbosity(void) const { return this->verbosity_; };

  void set_name(const std::string& s) { this->name_ = s; };
  std::string name(void) const { return this->name_; };

 protected:
  Geochemistry();

 private:
  std::string name_;
  Verbosity verbosity_;
  Beaker* beaker;
  Beaker::BeakerComponents components;
  Beaker::BeakerParameters parameters;

};

#endif // __Geochemistry_hpp__
