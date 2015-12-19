/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   GeometricModel.cc
 * @author Rao V. Garimella
 * @date Mon Aug  1 10:05:25 2011
 * 
 * @brief  
 * 
 * 
 */

#include "GeometricModel.hh"
#include "dbc.hh"
#include "errors.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class GeometricModel
// -------------------------------------------------------------

// -------------------------------------------------------------
// GeometricModel:: constructors / destructor
// -------------------------------------------------------------

// Constructor

GeometricModel::GeometricModel(const unsigned int dim, 
                               const VerboseObject *verbobj) : 
    topo_dimension_(dim), verbosity_obj_(verbobj)
{
  if (dim != 2 && dim != 3) {
    std::cerr << "Only 2D and 3D domains are supported" << std::endl;
    throw std::exception();
  }
  Regions.clear();
}

// Copy constructor

GeometricModel::GeometricModel(const GeometricModel& old) :
  topo_dimension_(old.dimension()), verbosity_obj_(old.verbosity_obj())
{
  int i, nr;

  nr = old.Num_Regions();

  for (i = 0; i < nr; i++) {
    RegionPtr r = old.Region_i(i);
    Regions.push_back(r);
  }
}


GeometricModel::GeometricModel(const unsigned int dim,
                               Teuchos::ParameterList gm_params,
                               const Epetra_MpiComm *comm,
                               const VerboseObject *verbobj) :
  topo_dimension_(dim), verbosity_obj_(verbobj)
{
  
  if (dim != 2 && dim != 3) {
    Errors::Message mesg("Only 2D and 3D domains are supported");
    amanzi_throw(mesg);
  }

  const int region_id_offset = 59049; // arbitrary number to avoid clashing
                                      // with IDs of LabeledSet regions
  int ngregions = 0; // Number of regions 

  // Go through the parameter list and populate the geometric model with regions

  for (Teuchos::ParameterList::ConstIterator i = gm_params.begin(); i != gm_params.end(); i++)
    {
      if (gm_params.isSublist(gm_params.name(i))) 
        {

          // Region name - get that from parameter list

          std::string region_name = gm_params.name(i);

          // Region ID - our internal numerical identifier

          unsigned int region_id = ++ngregions+region_id_offset;


          // Extract sublist specifying region

          const Teuchos::ParameterList &reg_spec = gm_params.sublist(gm_params.name(i));

          // While the XML file does not prevent it, there should only be one
          // specification of what the region looks like

          unsigned int k = 0;
          for (Teuchos::ParameterList::ConstIterator j = reg_spec.begin(); j != reg_spec.end(); j++, k++) 
            {
          
              if (k > 1) 
                {
                  std::stringstream sstream;
                  sstream << "ERROR: Region " << region_name << 
                    " described in multiple ways";
                  Errors::Message mesg(sstream.str());
                  amanzi_throw(mesg);
                }



              // Shape of the region

              std::string shape = reg_spec.name(j);


              // Create the region

              Amanzi::AmanziGeometry::RegionPtr regptr = 
                RegionFactory(region_name, region_id, reg_spec, dim, comm,
                              verbosity_obj());
              
              
              // Add it to the geometric model
              
              Regions.push_back(regptr);

            }
          
        }
      else 
        {
          Errors::Message mesg("Error: Improper region specification");
          amanzi_throw(mesg);
        }
    }
}


// Destructor

GeometricModel::~GeometricModel(void)
{

  // If a geometric model is deleted, we will not delete all the 
  // the regions added to it because someone else may be holding 
  // on to a pointer to the regions. For now, the top level routine
  // deleting the geometric model has to delete the regions first
  // to prevent a memory leak

  // Once we can get RegionFactory to work with Reference Counted
  // Pointers we can remove this comment

  Regions.clear();
}


// Constructor with Region List

GeometricModel::GeometricModel(const unsigned int dim, 
                               const std::vector<RegionPtr>& in_Regions,
                               const VerboseObject *verbobj) : 
  topo_dimension_(dim), verbosity_obj_(verbobj)
{
  Regions.clear();
  Regions = in_Regions;
}


// Add a Region

void GeometricModel::Add_Region(const RegionPtr& regptr)
{
  if (topo_dimension_ < regptr->dimension()) {
    Errors::Message mesg("Topological dimension of geometric model less than that of the region");
    amanzi_throw(mesg);
  }

  Regions.push_back(regptr);
}


// Number of Regions

int GeometricModel::Num_Regions(void) const
{
  return Regions.size();
}


// Get the i'th Region

RegionPtr GeometricModel::Region_i(const int i) const 
{
  return Regions[i];
}


// Get a region by its ID

RegionPtr GeometricModel::FindRegion(const int id) const
{
  for (std::vector<RegionPtr>::const_iterator r = Regions.begin();
       r != Regions.end(); ++r) {
    if ((*r)->id() == id) return *r;
  }
  return NULL;
}


// Get a region by its name

RegionPtr GeometricModel::FindRegion(const std::string name) const
{
  for (std::vector<RegionPtr>::const_iterator r = Regions.begin();
       r != Regions.end(); ++r) {
    if ((*r)->name() == name) return *r;
  }
  return NULL;
}


// Check if regions cover the domain extents.  
// This will work perfectly for domains with rectangular regions
// but not so for other types of regions

bool GeometricModel::Rough_Check_Tiling(void) const
{
  return true;
}


} // namespace AmanziGeometry
} // namespace Amanzi
