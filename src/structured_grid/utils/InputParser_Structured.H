#include "Teuchos_ParameterList.hpp"
#include "Teuchos_MPISession.hpp"

namespace Amanzi {
  namespace AmanziInput {
    
    static int ndim;
    static double simulation_time;
    static Teuchos::Array<double> domlo;
    static Teuchos::Array<double> domhi;      

    typedef struct {
	std::string phase;
	std::string component;
	std::string solute;
      } tridata;

    std::string underscore(const std::string& instring);

    Teuchos::ParameterList convert_to_structured(const Teuchos::ParameterList& parameter_list);

    Teuchos::ParameterList setup_structured();
    
    void convert_to_structured_mesh    (const Teuchos::ParameterList& parameter_list, 
					Teuchos::ParameterList&       struc_list);
    void convert_to_structured_control (const Teuchos::ParameterList& parameter_list, 
					Teuchos::ParameterList&       struc_list);
    void convert_to_structured_region  (const Teuchos::ParameterList& parameter_list, 
					Teuchos::ParameterList&       struc_list);
    void convert_to_structured_material(const Teuchos::ParameterList& parameter_list, 
					Teuchos::ParameterList&       struc_list);
    void convert_to_structured_comp    (const Teuchos::ParameterList& parameter_list, 
					Teuchos::ParameterList&       struc_list);   
    void convert_to_structured_tracer  (const Teuchos::ParameterList& parameter_list, 
					Teuchos::ParameterList&       struc_list);
    void convert_to_structured_output  (const Teuchos::ParameterList& parameter_list, 
					Teuchos::ParameterList&       struc_list);
  }
}
