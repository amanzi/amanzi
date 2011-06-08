
#include <string>
#include <sstream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "test_input.hh"

typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;


const std::string TestInput::num_sublist_keyword_  = "Number of sublists";
const std::string TestInput::sublist_root_name_    = "Sublist"; 

/* Default constructor */
TestInput::TestInput()
  : InputBaseClass("Test Parameters")
{

  /* Default parameter values */
  Int_  = 0;
  Dbl_  = 0.0;
  Bool_ = false; 

  num_sublist_ = 0;   /* Number of block regions */

  /* Set the private parameter list */
  Teuchos::RCP<Teuchos::ParameterList> plist =
      defineValidTestParameterList();

  setMyParamList(plist);

}

/* Constructor with all parameters set */
TestInput::TestInput(int x0, double x1, bool x2) 
  : InputBaseClass("Test Parameters")
{
  Int_  = x0;
  Dbl_  = x1;
  Bool_ = x2;

  num_sublist_ = 0;

  /* Set the private parameter list */
  Teuchos::RCP<Teuchos::ParameterList> plist =
      defineValidTestParameterList();

  setMyParamList(plist);

}

#if 0
/* Override setParameterList from ParameterListAcceptor */
void TestInput::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const
                                    &plist)
{
  //TEST_FOR_EXPECT(is_null(plist));

  /* Validate the incoming list */
  Teuchos::RCP<const Teuchos::ParameterList> plist_valid = getValidParameters();
  plist->validateParametersAndSetDefaults(*plist_valid,0);

  /* Set the private parameter list from AcceptorDefaultBase */
  setMyParamList(plist);

  /* Set private variables */
  nx_ = plist->get<int>("Number of cells in X");
  ny_ = plist->get<int>("Number of cells in Y");
  nz_ = plist->get<int>("Number of cells in Z");

  x0_ = plist->get<double>("X_Min");
  x1_ = plist->get<double>("X_Max");
  y0_ = plist->get<double>("Y_Min");
  y1_ = plist->get<double>("Y_Max");
  z0_ = plist->get<double>("Z_Min");
  z1_ = plist->get<double>("Z_Max");

  num_blocks_ = plist->get<int>(num_blocks_keyword_);


}

Teuchos::RCP<const Teuchos::ParameterList>
TestInput::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validParams;

  if ( is_null(validParams) ) {

    validParams = defineValidTestParameterList();

  }

  return validParams;
}

std::string TestInput::generate_mesh_block_name(int label) const
{
  std::stringstream ss;
  ss << "Mesh block " << label;
  return ss.str();
}

void TestInput::add_block(double z0, double z1)
{
  /* Will use Validators that compare zmin and zmax when available */
  Teuchos::ParameterEntry z0_entry;
  z0_entry.setValue<double>(z0,false,"",this->DoubleOnlyValidator_);

  Teuchos::ParameterEntry z1_entry;
  z1_entry.setValue<double>(z1,false,"",this->DoubleOnlyValidator_);


  /* Update the number of blocks */
  Teuchos::RCP<Teuchos::ParameterList> plist = getNonconstParameterList();
  int num_blocks = plist->get<int>(num_blocks_keyword_); 
  num_blocks++;
  plist->set<int>(num_blocks_keyword_,num_blocks);

  /* Now add the new block sublist */
  std::string block_name = generate_mesh_block_name(num_blocks);
  Teuchos::ParameterList& new_block =
      plist->sublist(blocks_sublist_name_,true).sublist(block_name);
  new_block.setEntry("Z0",z0_entry);
  new_block.setEntry("Z1",z1_entry);

  setParameterList(plist);

}
#endif

/* Private */

Teuchos::RCP<Teuchos::ParameterList> 
TestInput::defineValidTestParameterList() const
{
  Teuchos::RCP<Teuchos::ParameterList> 
      plist = Teuchos::rcp(new Teuchos::ParameterList(paramListName_));
  Teuchos::ParameterEntry entry;

  /* Int parameter */
  entry.setValue(Int_,
                 true,
                 "",
                 this->IntOnlyValidator_);
  plist->setEntry("Integer Parameter",entry);

  entry.setValue(Dbl_,
                 true,
                 "",
                 this->DoubleOnlyValidator_);
  plist->setEntry("Double",entry);

  /* bool Parameter */
  entry.setValue(Bool_,
                 true,
                 "",
                 Teuchos::null);
  plist->setEntry("Boolean Parameter",entry);

  /* Number of sublists */
  entry.setValue(num_sublist_,
                 true,
                 "Number of sublists",
                 this->IntOnlyValidator_);
  plist->setEntry(num_sublist_keyword_,entry);
  plist->sublist("Sublists");


  return plist;
}
