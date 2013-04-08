#include <algorithm>
#include "errors.hh"

#include "domain-function.hh"

namespace Amanzi {
namespace Functions {

/* ******************************************************************
* Calculate pairs <list of cells, function>
****************************************************************** */
void DomainFunction::Define(const std::vector<std::string>& regions,
                            const Teuchos::RCP<const VectorFunction>& f) 
{
 
  // Create the domain
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  
  // add the spec
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
}


void DomainFunction::Define(std::string region,
                            const Teuchos::RCP<const VectorFunction>& f) 
{
  RegionList regions(1,region);
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
}



/* ******************************************************************
* Compute and normalize the result, so far by volume
****************************************************************** */
void DomainFunction::Compute(double time)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }
  
  if (specs_and_ids_.size() == 0) return;

  // create the input tuple
  int dim = mesh_->space_dimension();
  double *args = new double[1+dim];
  double *xargs = args+1;
  args[0] = time;
 

  for (SpecAndIDsList::const_iterator
           spec_and_ids=specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_and_ids!=specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_and_ids) {
    
    Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;
    for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
      AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
      for (int i=0; i!=dim; ++i) xargs[i] = xc[i];
      // Careful tracing of the typedefs is required here: spec_and_ids->first
      //  is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*id] = (*(*spec_and_ids)->first->second)(args)[0];
    }
  }

  delete [] args;
}

  
} // namespace Functions
}  // namespace Amanzi
