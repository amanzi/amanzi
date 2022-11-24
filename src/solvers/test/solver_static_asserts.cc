//
// If this compiles, the test passes...
//
#include "UnitTest++.h"
#include "Epetra_CrsMatrix.h"

//#include "InverseFactory.hh"


// class Epetra_CrsMatrix {
//  public:
//   int FillComplete(bool OptimizeDataStorage) { return 1;}
//   int FillComplete(std::string& other1, std::string& other2) { return 1;}
// };


//
// trait: is_assembled
// ------------------------------
//
// Default is false
//
template <typename Operator, typename Dummy = void>
struct is_assembled {
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
template <typename Operator>
struct is_assembled<
  Operator,
  typename std::enable_if<std::is_member_function_pointer<
    decltype(static_cast<int (Operator::*)(bool)>(&Operator::FillComplete))>::value>::type> {
  const static bool value = true;
};


SUITE(SOLVERS)
{
  TEST(static_asserts)
  {
    //    static_assert(!is_assembling<Epetra_CrsMatrix>::value, "Epetra_CrsMatrix is not assembling");
    static_assert(is_assembled<Epetra_CrsMatrix>::value, "Epetra_CrsMatrix is assembled");
  }
}
