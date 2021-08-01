/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

// Helper functions for matrix traits.

#pragma once

namespace Amanzi {
namespace AmanziSolvers {
namespace Impl {



//
// trait: is_assembled
// ------------------------------
//
// Default is false
//
template<typename Operator, typename Dummy=void>
struct is_assembled
{
  const static bool value = false;
};

//
// is_assembled::value == true if Operator has method FillComplete()
//
template<typename Operator>
struct is_assembled<Operator,
                    typename std::enable_if<std::is_member_function_pointer<
                      decltype(static_cast<int (Operator::*)(bool)>(&Operator::FillComplete))>::value>::type>
{
  const static bool value = true;
};



//
// trait: is_assembling
// ------------------------------
//
// Default is false
//
template<typename Operator, typename Dummy=void>
struct is_assembling
{
  const static bool value = false;
};

//
// is_assembling::value == true if Operator has method AssembleMatrix()
//
template<typename Operator>
struct is_assembling<Operator,
                     typename std::enable_if<std::is_member_function_pointer<
                                               decltype(static_cast<Teuchos::RCP<Matrix_type> (Operator::*)(void)>(&Operator::A))>::value>::type>
{
  const static bool value = true;
};


} // namespace Impl
} // namespace AmanziSolvers
} // namespace Amanzi
