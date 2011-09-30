#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "FlowBC.hpp"

#include "float.h"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"


namespace Amanzi
{

  FlowBC::FlowBC(Teuchos::ParameterList &list, 
		 const Teuchos::RCP<AmanziMesh::Mesh> &mesh) : mesh_(mesh)
  {
    
    // set the line prefix for output
    this->setLinePrefix("Amanzi::FlowBC      ");
    // make sure that the line prefix is printed
    this->getOStream()->setShowLinePrefix(true);

    // Read the sublist for verbosity settings.
    Teuchos::readVerboseObjectSublist(&list,this);    

    // set up the verbose output stream
    using Teuchos::OSTab;
    Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab       
    

    // first figure out the number of sublists, this number 
    // equals the number of boundary conditions
    int nbc = 0;
    for (Teuchos::ParameterList::ConstIterator it = list.begin(); it != list.end(); it++)
      {
	if (list.isSublist(list.name(it))) 
	  {
	    if ( list.name(it) != "VerboseObject" )
	      {
		nbc++;
	      }
	  }
      }

    if (nbc == 0) 
      {
	Errors::Message m;
	m << "Error in Flow BC list.";
	m << "  No flow boundary conditions have been defined.";
	Exceptions::amanzi_throw(m);	
      }


    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))
      {
	*out << "found " << nbc << " flow boundary condition parameter lists." << std::endl;
      }    


    // create the array that holds all the data we will read
    bc_.resize(nbc);
    
    
    // map to keep track of side set ids
    std::map<int,int> ss_id_map;

    int i = 0;
    for (Teuchos::ParameterList::ConstIterator it = list.begin(); it != list.end(); it++)
      {
	if (list.isSublist(list.name(it))) 
	  { 
	    if (list.name(it) != "VerboseObject") 
	      {
		
		Teuchos::ParameterList &bc_param = list.sublist( list.name(it) );
		
		// Get the Set ID and verify that it is valid.
		if ( !bc_param.isParameter("Side set ID") )
		  {
		    Errors::Message m;
		    m << "Syntax error in Flow BC list.";
		    m << "  Sublist " << list.name(it).c_str();
		    m << " does not contain required parameter: Side set ID";
		    Exceptions::amanzi_throw(m);	
		  }
		
		bc_[i].SetID = bc_param.get<int>("Side set ID");
		
		if (!mesh->valid_set_id(bc_[i].SetID, AmanziMesh::FACE)) 
		  {
		    std::stringstream ss;
		    ss << bc_[i].SetID;
		    
		    Errors::Message m;
		    m << "Error in Flow BC list.";
		    m << " Side set ID " <<  ss.str().c_str();
		    m << " refered to in sublist " << list.name(it).c_str();
		    m << " does not exist in the mesh.";
		    Exceptions::amanzi_throw(m);	
		  }
		
		// Get the corresponding list of (local) face IDs.
		bc_[i].Faces.resize(mesh->get_set_size(bc_[i].SetID, AmanziMesh::FACE, AmanziMesh::USED));
		mesh->get_set(bc_[i].SetID, AmanziMesh::FACE, AmanziMesh::USED, bc_[i].Faces.begin(), bc_[i].Faces.end());
		
		// Get the BC type and check it against the list of defined types.
		
		if ( !bc_param.isParameter("Type") )
		  {
		    Errors::Message m;
		    m << "Syntax error in Flow BC list.";
		    m << " Sublist " << list.name(it).c_str();
		    m << " does not contain required parameter: Type";
		    Exceptions::amanzi_throw(m);
		  }
		
		std::string type = bc_param.get<string>("Type");
		
		bool read_value, need_aux;
		if (type == "Pressure Constant") 
		  {
		    bc_[i].Type = PRESSURE_CONSTANT;
		    need_aux = true;
		    read_value = true;
		  }
		else if (type == "No Flow") 
		  {
		    bc_[i].Type = NO_FLOW;
		    need_aux = false;
		    read_value = false;
		  }
		else if (type == "Darcy Constant") 
		  {
		    bc_[i].Type = DARCY_CONSTANT;
		    need_aux = false;
		    read_value = true;
		  }
		else if (type == "Static Head") 
		  {
		    bc_[i].Type = STATIC_HEAD;
		    need_aux = true;
		    read_value = true;
		  }
		else if (type == "Time Dependent Pressure Constant" ) 
		  {
		    bc_[i].Type = TIME_DEPENDENT_PRESSURE_CONSTANT;
		    need_aux = true;
		    read_value = true;
		  }
		else 
		  { 
		    Errors::Message m;
		    m << "Error in Flow BC list, in ";
		    m << " sublist " << list.name(it).c_str() << ".";
		    m << " Type must be one of Pressure Constant, No Flow, Darcy Constant,";
		    m << " Static Head, or Time Dependent Pressure Constant";
		    Exceptions::amanzi_throw(m);
		  }
		
		
		// Temp storage needed for Dirichlet-type conditions.
		if (need_aux) 
		  {
		    bc_[i].Aux.resize(bc_[i].Faces.size());
		  }
		
		// Get the BC data value if required.
		if (read_value) 
		  {
		    if ( !bc_param.isParameter("BC value") )
		      {
			Errors::Message m;
			m << "Syntax error in Flow BC list.";
			m << " Sublist " << list.name(it).c_str();
			m << " does not contain required parameter: BC value";
			Exceptions::amanzi_throw(m);
		      }
		    
		    bc_[i].Value = bc_param.get<double>("BC value");
		  }
		
		
		// now for a time dependent pressure constant we need to 
		// read three additional parameters
		
		if (type == "Time Dependent Pressure Constant") 
		  {
		    
		    if ( !bc_param.isParameter("Initial BC value") )
		      {
			Errors::Message m;
			m << "Syntax error in Flow BC list.";
			m << " Sublist " << list.name(it).c_str();
			m << " does not contain required parameter: Initial BC value";
			Exceptions::amanzi_throw(m);
		      }
		    
		    bc_[i].InitialValue = bc_param.get<double>("Initial BC value");
		    
		    if ( !bc_param.isParameter("Initial Time") )
		      {
			Errors::Message m;
			m << "Syntax error in Flow BC list.";
			m << " Sublist " << list.name(it).c_str();
			m << " does not contain required parameter: Initial Time";
			Exceptions::amanzi_throw(m);
		      }		
		    
		    bc_[i].InitialTime = bc_param.get<double>("Initial Time");
		    
		    if ( !bc_param.isParameter("Final Time") )
		      {
			Errors::Message m;
			m << "Syntax error in Flow BC list.";
			m << " Sublist " << list.name(it).c_str();
			m << " does not contain required parameter: Final Time";
			Exceptions::amanzi_throw(m);
		      }				
		    
		    bc_[i].FinalTime = bc_param.get<double>("Final Time");
		  }

		// only proceed if there is not already a flow boundary condition
		// on the current side set
		if (ss_id_map.count(bc_[i].SetID) == 0)
		  {
		    ss_id_map.insert( std::pair<int,int>(bc_[i].SetID,i) );
		  }
		else
		  {
		    std::stringstream ss;
		    ss << bc_[i].SetID;
		    
		    Errors::Message m;
		    m << "Error in Flow BC list,";
		    m << " Sublist " << list.name(it).c_str() << ".";
		    m << " Side set ID " << ss.str().c_str();
		    m << " is refered to in another sublist of the Flow BC list.";
		    Exceptions::amanzi_throw(m);
		  }
		
		// increment the boundary condition counter
		i++;

	      }
	  }

      }

    //dump the data structure
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))
      {
    	*out << "Dumping the BC data structure..." << std::endl;
    	for (std::vector<bc_spec>::iterator it = bc_.begin(); 
    	     it != bc_.end();
    	     it++)
    	  {
    	    *out << "Type          = " << it->Type << std::endl;
    	    *out << "SetID         = " << it->SetID << std::endl;
    	    *out << "#faces        = " << it->Faces.size() << std::endl;
    	    *out << "Value         = " << it->Value << std::endl;
    	    *out << "Initial Value = " << it->InitialValue << std::endl;
    	    *out << "Final Value   = " << it->InitialValue << std::endl;
    	    *out << "Initial Time  = " << it->InitialValue << std::endl << std::endl;
    	  }
      }        



  }
} // close namespace Amanzi
