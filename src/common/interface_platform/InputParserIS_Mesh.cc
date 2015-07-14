#include <sstream>
#include <string>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputParserIS.hh"
#include "InputParserIS_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateMeshList_(Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  Teuchos::ParameterList msh_list;
  Teuchos::RCP<Teuchos::ParameterList> mlist = Teuchos::sublist(plist, "Mesh", true);
  Teuchos::RCP<Teuchos::ParameterList> ulist = Teuchos::sublist(mlist, "Unstructured", true);

  if (ulist->isSublist("Generate Mesh")) {
    Teuchos::ParameterList& generate = ulist->sublist("Generate Mesh").sublist("Uniform Structured");
    Teuchos::Array<int> ncells = generate.get<Teuchos::Array<int> >("Number of Cells");
    Teuchos::Array<double> low = generate.get<Teuchos::Array<double> >("Domain Low Coordinate");
    Teuchos::Array<double> high = generate.get<Teuchos::Array<double> >("Domain High Coordinate");

    Teuchos::ParameterList& msh_gen = msh_list.sublist("Unstructured").sublist("Generate Mesh");

    msh_gen.set< Teuchos::Array<int> >("Number of Cells",ncells);
    msh_gen.set< Teuchos::Array<double> >("Domain Low Coordinate",low);
    msh_gen.set< Teuchos::Array<double> >("Domain High Coordinate",high);

  } else if (ulist->isSublist("Read Mesh File")) {
    std::string format = ulist->sublist("Read Mesh File").get<std::string>("Format");
        
    if (format == "Exodus II") {
      // process the file name to replace .exo with .par in the case of a parallel run
      Teuchos::ParameterList& fn_list = msh_list.sublist("Unstructured").sublist("Read Mesh File");
      fn_list.set<std::string>("Format", "Exodus II");
      std::string file = ulist->sublist("Read Mesh File").get<std::string>("File");
      std::string suffix(file.substr(file.size()-4, 4));

      if (suffix != ".exo") {
        Exceptions::amanzi_throw(Errors::Message("Exodus II was specified as a mesh file format but the suffix of the file that was specified is not .exo"));
      }

      // figure out if this is a parallel run
      std::string framework("Unspecified");
      if (ulist->isSublist("Expert")) {
        framework = ulist->sublist("Expert").get<std::string>("Framework");
      }

      // Assume that if the framework is unspecified then stk::mesh is used
      // This is obviously a kludge but I don't know how to get around it
      //
      // We also have to be able to tell if we have prepartitioned
      // files or if we have one file that we want to partition

      if (numproc_ > 1) {
        std::string par_file(file);
        par_file.replace(file.size()-4,4,std::string(".par"));

        // attach the right extensions as required by Nemesis file naming conventions
        // in which files are named as mymesh.par.N.r where N = numproc and r is rank

        int rank = numproc_-1;
        int ndigits = (int)floor(log10(numproc_)) + 1;
        std::string fmt = boost::str(boost::format("%%s.%%d.%%0%dd") % ndigits);
        std::string par_file_w_ext = boost::str(boost::format(fmt) % par_file % numproc_ % rank);
        boost::filesystem::path p(par_file_w_ext);

        if (boost::filesystem::exists(p))
          fn_list.set<std::string>("File",par_file); // Nemesis file exists. Use the .par extension
        else
          fn_list.set<std::string>("File",file); // Use original .exo file extension

      } else {
        // don't translate the suffix if this is a serial run
        fn_list.set<std::string>("File",file);
      }

      msh_list.sublist("Unstructured").sublist("Read Mesh File") = fn_list;

    } else {
      msh_list.sublist("Unstructured").sublist("Read Mesh File") = ulist->sublist("Read Mesh File");
    }
  }

  if (ulist->isSublist("Expert")) {
    msh_list.sublist("Unstructured").sublist("Expert") = ulist->sublist("Expert");
  }

  return msh_list;
}


/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputParserIS::CopyDomainList_(Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  Teuchos::ParameterList dom_list = *Teuchos::sublist(plist, "Domain", false);
  return dom_list;
}


/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputParserIS::CopyRegionsList_(Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  Teuchos::ParameterList reg_list;

  // in the Mesh list
  // first find the mesh file name specified in the Mesh list
  std::string meshfile;
  if (plist->sublist("Mesh").sublist("Unstructured").isSublist("Read Mesh File")) {
    meshfile = plist->sublist("Mesh").sublist("Unstructured").sublist("Read Mesh File").get<std::string>("File");
  }

  Teuchos::ParameterList rlist = *Teuchos::sublist(plist, "Regions", true);
  // loop over all regions and find possible labeled set definitions
  for (Teuchos::ParameterList::ConstIterator i = rlist.begin(); i != rlist.end(); i++) {
    // only count sublists
    if (rlist.isSublist(rlist.name(i))) {
      if (rlist.sublist((rlist.name(i))).isSublist("Region: Labeled Set")) {
        std::string file = rlist.sublist((rlist.name(i))).sublist("Region: Labeled Set").get<std::string>("File");
        boost::filesystem::path meshfile_path(meshfile);
        boost::filesystem::path labeled_set_meshfile_path(file);
        if (meshfile != file) {
          Errors::Message msg("There is a labeled set region that refers to a mesh file that is different from the mesh file that is defined in the Mesh list: " + file);
          Exceptions::amanzi_throw(msg);
        }
      }
    }
  }

  // all is well, return the Regions list for insertion into the native list
  reg_list = plist->sublist("Regions");

  return reg_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
