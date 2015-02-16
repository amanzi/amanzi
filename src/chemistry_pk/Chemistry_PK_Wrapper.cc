/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors:: Daniil Svyatskiy

  Temporary wrapper converting the Chemistry_PK, which inherits from 
  BDFFnBase<CompositeVector>, to use TreeVectors.

*/


#include "chemistry_pk.hh"
#include "Chemistry_PK_Wrapper.hh"

namespace Amanzi {
namespace AmanziChemistry{


Chemistry_PK_Wrapper::Chemistry_PK_Wrapper(Teuchos::ParameterList& pk_tree,
					   const Teuchos::RCP<Teuchos::ParameterList>& global_list,
					   const Teuchos::RCP<State>& S,
					   const Teuchos::RCP<TreeVector>& soln) :
    S_(S),
    soln_(soln),
    glist_(global_list)
{


  std::string pk_name = pk_tree.name();

  const char* result = pk_name.data();
  while ((result = std::strstr(result, "->")) != NULL) {
    result += 2;
    pk_name = result;   
  }


// grab the component names

  if (glist_->isSublist("Cycle Driver")){
    if (glist_->sublist("Cycle Driver").isParameter("component names")){   
      comp_names_ = glist_->sublist("Cycle Driver").get<Teuchos::Array<std::string> >("component names").toVector();
    }
    else{
      Errors::Message msg("Chemistry_PK_Wrapper: Cycle Driver has no input parameter component names.");
      Exceptions::amanzi_throw(msg);
    }
  }


  Teuchos::ParameterList chemistry_parameter_list = glist_->sublist("PKs").sublist("Chemistry");
  

  CS = Teuchos::rcp(new AmanziChemistry::Chemistry_State(chemistry_parameter_list, comp_names_, S));

  //CS->Initialize();
 
// construct
  pk_ = Teuchos::rcp(new Chemistry_PK(chemistry_parameter_list, CS));

}


}//namespace
}//namespace
