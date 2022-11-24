/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! Base factory for self-registering classes of a given base type.
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*!
  
  In many cases in Amanzi/ATS, we may have multiple options that inherit a
  common (likely purely) virtual class.  For instance, many implementations
  of the equations of state class will provide a basic method for rho(T,p),
  including both real "fits" to data, analytic expressions, and fake EOS
  classes for testing which may be constant.  We would like to be able to:

   1. choose the implementation at run time
   2. easily add new implementations

  To do #1, we use a factory design pattern.  Like most factories, an
  implementation must be "registered" with the factory.  To do #2, this
  registration must NOT be done in the factory's source code itself.

  This is made a little easier by the fact that nearly all of these things
  will be constructed using a single interface for the constructor, which
  (explicitly) takes a single argument -- a Teuchos::ParameterList -- and
  parses that list for its actual parameters.  While it is usually a good
  idea to have a factory take the input list, do the parsing, and call the
  model's constructor with the parameters, that would require every model
  implementation to have its own factory.  To simply things for scientists
  writing these models, we choose to do the parsing within the
  constructor/initialization.

  The obvious exception to this is the model type parameter, which must get
  read by a factory and mapped to an implementation's constructor.

  This factory is a general purpose factory which is templated to take a
  single base class.  Implementations of that base class then "register"
  themselves with the factory instance (which is stored statically since we
  cannot correctly manage the cleanup).  This factory assumes all
  constructors for all implementations of all base classes take a single
  Teuchos::ParameterList as an argument.

  Simplest usage (for our EOS example):

   // eos_factory.cc  (no .hh file necessary)
   #include "eos.hh" // header for class EOS, a purely virtual base class
   #include "factory.hh" // this file
   template <> Factory<EOS>::map_type* Factory<EOS>::map_; // explicity
                                                           // instantiate the
                                                           // static registry

   // eos_implementation.hh
   #include "eos.hh"
   class DerivedEOS : public EOS {
     DerivedEOS(Teuchos::ParameterList& plist);
     ...
   private:
     static RegisteredFactory<EOS,DerivedEOS> factory_; // my factory
     ...
   };


   // pk_using_an_eos.cc
   #include "eos.hh"

   void init(...) {
     ...
     Factory<EOS> eos_factory;
     my_eos_ = eos_factory.CreateInstance("my_eos_type", eos_plist);
     ...

  You're not supposed to understand this, just find another example that uses
  it and copy it.
*/

#ifndef AMANZI_FACTORY_HH_
#define AMANZI_FACTORY_HH_

#include <iostream>
#include <map>
#include <string>
#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Utils {

template <typename TBase>
class Factory {
 public:
  typedef std::map<std::string, TBase* (*)(Teuchos::ParameterList&)> map_type;

  static TBase* CreateInstance(const std::string& s, Teuchos::ParameterList& plist)
  {
    typename map_type::iterator iter = GetMap()->find(s);
    if (iter == GetMap()->end()) {
      std::cout << "Factory: cannot get item of type: " << s << std::endl;

      for (typename map_type::iterator piter = GetMap()->begin(); piter != GetMap()->end();
           ++piter) {
        std::cout << "  option: " << piter->first << std::endl;
      }
      return 0;
    }
    return iter->second(plist);
  }

 protected:
  static map_type* GetMap()
  {
    static map_type* map_;
    if (!map_) { map_ = new map_type; }
    return map_;
  }

 private:
  static map_type* map_;
};

template <typename TBase>
typename Factory<TBase>::map_type* Factory<TBase>::map_;

template <typename TBase, typename TDerived>
TBase*
CreateT(Teuchos::ParameterList& plist)
{
  return new TDerived(plist);
}


template <typename TBase, typename TDerived>
class RegisteredFactory : public Factory<TBase> {
 public:
  // Constructor for the registered factory.  Needs some error checking in
  // case a name s is already in the map? (i.e. two implementations trying to
  // call themselves the same thing) --etc
  RegisteredFactory(const std::string& s)
  {
    for (typename Factory<TBase>::map_type::iterator iter = Factory<TBase>::GetMap()->begin();
         iter != Factory<TBase>::GetMap()->end();
         ++iter) {}
    Factory<TBase>::GetMap()->insert(
      std::pair<std::string, TBase* (*)(Teuchos::ParameterList&)>(s, &CreateT<TBase, TDerived>));
    for (typename Factory<TBase>::map_type::iterator iter = Factory<TBase>::GetMap()->begin();
         iter != Factory<TBase>::GetMap()->end();
         ++iter) {}
  }
};

} // namespace Utils
} // namespace Amanzi

#endif
