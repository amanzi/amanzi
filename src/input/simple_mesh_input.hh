#ifndef SIMPLE_MESH_INPUT_HPP
#define SIMPLE_MESH_INPUT_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

#include "input_base_class.hh"

class SimpleMeshInput : public InputBaseClass
{
  public:

    /* Constructors */
    SimpleMeshInput();
    SimpleMeshInput(double x0, double y0, double z0,
                    double x1, double y1, double z1,
                    int nx, int ny, int nz);

    /* Methods */
    void add_block(double z0,double z1);


    /* Overridden methods from ParameterListAcceptorDefaultBase */
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const & );
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;


  private:
  
    /* Keywords */
    static const std::string num_blocks_keyword_;
    static const std::string blocks_sublist_name_;
    int num_blocks_;

    int nx_, ny_, nz_;
    double x0_, x1_, y0_, y1_, z0_, z1_ ;
    
    Teuchos::RCP<Teuchos::ParameterList>
        defineValidSimpleMeshParameterList() const; 

    std::string generate_mesh_block_name(int) const;
    Teuchos::RCP<Teuchos::ParameterList> create_block_list() const;


};

        
#endif
