
#include "NKADirFactory.H"


NKADirFactory::NKADirFactory(const Teuchos::RCP<NOX::GlobalData> &gd, 
			     Teuchos::ParameterList &params, 
			     const NOX::Abstract::Vector &initvec)

{
  my_dir = Teuchos::rcp(new NKADirection(gd, params, initvec));
};


NKADirFactory::~NKADirFactory() 
{
};


Teuchos::RCP<NOX::Direction::Generic> 
NKADirFactory::buildDirection(const Teuchos::RCP<NOX::GlobalData> &gd, 
			      Teuchos::ParameterList &params) const
{
  return my_dir;
};
