//
// If this compiles, the test passes...
//
#include "UnitTest++.h"

//#include "InverseFactory.hh"
#include "AmanziTypes.hh"
using namespace Amanzi; 

// class Tpetra::CsrMatrix {
//  public:
//  void 	fillComplete (const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) 
//  void 	fillComplete (const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)
// };


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
// template<typename Operator>
// struct is_assembled<Operator,
//                     typename std::enable_if<!std::is_void<Operator>::value>::type>
// {
//   const static bool value = true;
// };
template<typename Operator>
struct is_assembled<Operator,
                    typename std::enable_if<std::is_member_function_pointer<
                      decltype(
                        static_cast<void (Operator::*) (const Teuchos::RCP<Teuchos::ParameterList> &)>
                        (&Operator::fillComplete))>::value>::type>
{
  const static bool value = true;
};



SUITE(SOLVERS) {
  TEST(static_asserts) {
    static_assert(is_assembled<Matrix_type>::value, "Matrix_type is assembled");
  }
}
