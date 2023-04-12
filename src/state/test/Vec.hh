/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#pragma once
#include "Data_Helpers.hh"

int g_constructor_calls_default = 0;
int g_constructor_calls_copy = 0;
int g_constructor_calls_main = 0;

struct Vec {
  Vec() { g_constructor_calls_default++; }

  Vec(const Vec& other) : v(other.v) { g_constructor_calls_copy++; }

  Vec(int size) : v(size, 0.0) { g_constructor_calls_main++; }

  std::vector<double> v;
};

class VecFactory {
 public:
  VecFactory() : size_(-1) {}

  void set_size(int size) { size_ = size; }

  Teuchos::RCP<Vec> Create()
  {
    assert(size_ >= 0);
    return Teuchos::rcp(new Vec(size_));
  }

 private:
  int size_;
};


template<>
bool
Amanzi::Helpers::Initialize<Vec>(Teuchos::ParameterList& plist,
           Vec& t)
{
  std::cout << "found it!" << std::endl;
  return true;
}

template<>
void
Amanzi::Helpers::WriteVis<Vec>(const Amanzi::Visualization& vis,
         Teuchos::ParameterList& attrs,
         const Vec& vec)
{}

template<>
void
Amanzi::Helpers::WriteCheckpoint<Vec>(const Amanzi::Checkpoint& vis,
         Teuchos::ParameterList& attrs,
         const Vec& vec)
{}

template<>
void
Amanzi::Helpers::ReadCheckpoint<Vec>(const Amanzi::Checkpoint& vis,
         Teuchos::ParameterList& attrs,
         Vec& vec)
{}

