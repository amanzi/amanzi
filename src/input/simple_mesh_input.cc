
#include <string>
#include <sstream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "simple_mesh_input.hh"

typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;


const std::string SimpleMeshInput::num_blocks_keyword_   = "Number of blocks";
const std::string SimpleMeshInput::blocks_sublist_name_  = "Blocks"; 

/* Default constructor */
SimpleMeshInput::SimpleMeshInput()
  : InputBaseClass("Simple Mesh Parameters")
{

  /* Default parameter values */
  nx_=0;   ny_=0;   nz_=0;     /* Number of cells*/
  x0_=0.0; y0_=0.0; z0_=0.0;   /* Lower mesh bounds */
  x1_=0.0; y1_=0.0; z1_=0.0;   /* Upper mesh bounds */

  num_blocks_ = 0;   /* Number of block regions */

  /* Set the private parameter list */
  Teuchos::RCP<Teuchos::ParameterList> plist =
      defineValidSimpleMeshParameterList();

  setMyParamList(plist);

}

/* Constructor with all parameters set */
SimpleMeshInput::SimpleMeshInput(double x0, double y0, double z0, 
                                 double  x1, double y1, double z1,
                                 int nx, int ny, int nz)
  : InputBaseClass("Simple Mesh Parameters")
{
  nx_=nx; ny_=ny; nz_=nz;
  x0_=x0; y0_=y0; z0_=z0;
  x1_=x1; y1_=y1; z1_=z1;

  num_blocks_ = 0;

  /* Set the private parameter list */
  Teuchos::RCP<Teuchos::ParameterList> plist =
      defineValidSimpleMeshParameterList();

  setMyParamList(plist);

}

/* Override setParameterList from ParameterListAcceptor */
void SimpleMeshInput::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const
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
SimpleMeshInput::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validParams;

  if ( is_null(validParams) ) {

    validParams = defineValidSimpleMeshParameterList();

  }

  return validParams;
}

std::string SimpleMeshInput::generate_mesh_block_name(int label) const
{
  std::stringstream ss;
  ss << "Mesh block " << label;
  return ss.str();
}

void SimpleMeshInput::add_block(double z0, double z1)
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


/* Private */

Teuchos::RCP<Teuchos::ParameterList> 
SimpleMeshInput::defineValidSimpleMeshParameterList() const
{
  Teuchos::RCP<Teuchos::ParameterList> 
      plist = Teuchos::rcp(new Teuchos::ParameterList(paramListName_));
  Teuchos::ParameterEntry entry;

  /* X direction */
  entry.setValue(x0_,
                 true,
                 "",
                 this->DoubleOnlyValidator_);
  plist->setEntry("X_Min",entry);

  entry.setValue(x1_,
                 true,
                 "",
                 this->DoubleOnlyValidator_);
  plist->setEntry("X_Max",entry);

  /* Y direction */
  entry.setValue(y0_,
                 true,
                 "",
                 this->DoubleOnlyValidator_);
  plist->setEntry("Y_Min",entry);

  entry.setValue(y1_,
                 true,
                 "",
                 this->DoubleOnlyValidator_);
  plist->setEntry("Y_Max",entry);

  /* Z direction */
  entry.setValue(z0_,
                 true,
                 "",
                 this->DoubleOnlyValidator_);
  plist->setEntry("Z_Min",entry);

  entry.setValue(z1_,
                 true,
                 "",
                 this->DoubleOnlyValidator_);
  plist->setEntry("Z_Max",entry);

  /* X direction */
  entry.setValue(nx_,
                 true,
                 "",
                 this->IntOnlyValidator_);
  plist->setEntry("Number of cells in X",entry);

  /* Y direction */
  entry.setValue(ny_,
                 true,
                 "",
                 this->IntOnlyValidator_);
  plist->setEntry("Number of cells in Y",entry);

  /* Z direction */
  entry.setValue(nz_,
                 true,
                 "",
                 this->IntOnlyValidator_);
  plist->setEntry("Number of cells in Z",entry);

  /* Number of blocks */
  entry.setValue(num_blocks_,
                 true,
                 "Number of blocks (regions)",
                 this->IntOnlyValidator_);
  plist->setEntry(num_blocks_keyword_,entry);
  plist->sublist("Blocks");


  return plist;
}











