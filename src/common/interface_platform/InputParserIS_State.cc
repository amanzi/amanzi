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
* STATE sublist
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateStateList_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList stt_list;

  // first we write initial conditions for scalars and vectors, not region-specific

  Teuchos::ParameterList& stt_ic = stt_list.sublist("initial conditions");
  //
  // --- gravity
  Teuchos::Array<double> gravity(spatial_dimension_);
  for (int i = 0; i != spatial_dimension_-1; ++i) gravity[i] = 0.0;
  gravity[spatial_dimension_-1] = -GRAVITY_MAGNITUDE;
  stt_ic.sublist("gravity").set<Teuchos::Array<double> >("value", gravity);
  //
  // --- viscosity
  //
  Teuchos::ParameterList& phase_list = plist->sublist("Phase Definitions");
  double viscosity = phase_list.sublist(phase_name).sublist("Phase Properties").sublist("Viscosity: Uniform").get<double>("Viscosity");
  stt_ic.sublist("fluid_viscosity").set<double>("value", viscosity);
  //
  // --- density
  //
  double density = phase_list.sublist(phase_name).sublist("Phase Properties").sublist("Density: Uniform").get<double>("Density");
  stt_ic.sublist("fluid_density").set<double>("value", density);
  // this is stupid, but for some reason we also have an array for water density, so here it goes...
  stt_ic.sublist("water_density").sublist("function").sublist("All")
      .set<std::string>("region","All")
      .set<std::string>("component","cell")
      .sublist("function").sublist("function-constant")
      .set<double>("value", density);

  constant_density = density;  // save it for PKs
  //
  // --- region specific initial conditions from material properties
  //
  std::map<std::string,int> region_to_matid;
  std::map<int,std::string> matid_to_material;
  int matid_ctr = 0;
  // loop over the material properties
  Teuchos::ParameterList& matprop_list = plist->sublist("Material Properties");
  for (Teuchos::ParameterList::ConstIterator i = matprop_list.begin(); i != matprop_list.end(); i++) {
    // get the regions
    Teuchos::Array<std::string> regions = matprop_list.sublist(matprop_list.name(i)).get<Teuchos::Array<std::string> >("Assigned Regions");

    // record the material ID for each region that this material occupies
    matid_ctr++;
    for (int ii = 0; ii < regions.size(); ii++) {
      if (region_to_matid.find(regions[ii]) == region_to_matid.end()) {
        region_to_matid[regions[ii]] = matid_ctr;
        matid_to_material[matid_ctr] = matprop_list.name(i);
      } else {
        std::stringstream ss;
        ss << "There is more than one material assinged to region " << regions[ii] << ".";
        Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
      }
    }

    // read the parameters that need to be written for each region in the Assigned regions list

    // create regions string
    std::string reg_str;
    for (Teuchos::Array<std::string>::const_iterator ireg=regions.begin(); ireg!=regions.end(); ++ireg) {
      reg_str = reg_str + *ireg;
    }

    // porosity...
    double porosity;
    if (matprop_list.sublist(matprop_list.name(i)).isSublist("Porosity: Uniform")) {
      porosity = matprop_list.sublist(matprop_list.name(i)).sublist("Porosity: Uniform").get<double>("Value");
    } else {
      Exceptions::amanzi_throw(Errors::Message("Porosity must be specified as Intrinsic Porosity: Uniform, for every region."));
    }
    Teuchos::ParameterList &porosity_ic = stt_ic.sublist("porosity");
    porosity_ic.sublist("function").sublist(reg_str)
        .set<Teuchos::Array<std::string> >("regions",regions)
        .set<std::string>("component","cell")
        .sublist("function").sublist("function-constant")
        .set<double>("value", porosity);

    // permeability...
    double perm_x, perm_y, perm_z;
    std::string perm_file, perm_attribute, perm_format;
    bool perm_init_from_file = false;
    Teuchos::ParameterList& mplist = matprop_list.sublist(matprop_list.name(i));
    if ((mplist.isSublist("Intrinsic Permeability: Uniform") || 
	 mplist.isSublist("Intrinsic Permeability: Anisotropic Uniform")) && 
	(mplist.isSublist("Hydraulic Conductivity: Uniform") || 
	 mplist.isSublist("Hydraulic Conductivity: Anisotropic Uniform"))) {
      Exceptions::amanzi_throw(Errors::Message("Permeability can only be specified either Intrinsic Permeability or Hydraulic Conductivity, but not both."));
    }

    if (mplist.isSublist("Intrinsic Permeability: Uniform")) {
      perm_x = perm_y = perm_z = mplist.sublist("Intrinsic Permeability: Uniform").get<double>("Value");
    } else if (mplist.isSublist("Intrinsic Permeability: Anisotropic Uniform")) {
      perm_x = mplist.sublist("Intrinsic Permeability: Anisotropic Uniform").get<double>("x");
      perm_y = mplist.sublist("Intrinsic Permeability: Anisotropic Uniform").get<double>("y");
      if (mplist.sublist("Intrinsic Permeability: Anisotropic Uniform").isParameter("z")) {
        if (spatial_dimension_ == 3) {
          perm_z = mplist.sublist("Intrinsic Permeability: Anisotropic Uniform").get<double>("z");
        } else {
          Exceptions::amanzi_throw(Errors::Message("Intrinsic Permeability: Anisotropic Uniform defines a value for z, while the spatial dimension of the problem not 3."));
        }
      }
    } else if (mplist.isSublist("Hydraulic Conductivity: Uniform")) {
      perm_x = perm_y = perm_z = mplist.sublist("Hydraulic Conductivity: Uniform").get<double>("Value");
      // now scale with rho, g and mu to get correct permeability values
      perm_x *= viscosity/(density*GRAVITY_MAGNITUDE);
      perm_y *= viscosity/(density*GRAVITY_MAGNITUDE);
      perm_z *= viscosity/(density*GRAVITY_MAGNITUDE);
    } else if (mplist.isSublist("Hydraulic Conductivity: Anisotropic Uniform")) {
      perm_x = mplist.sublist("Hydraulic Conductivity: Anisotropic Uniform").get<double>("x");
      perm_y = mplist.sublist("Hydraulic Conductivity: Anisotropic Uniform").get<double>("y");
      if (mplist.sublist("Hydraulic Conductivity: Anisotropic Uniform").isParameter("z")) {
        if (spatial_dimension_ == 3) {
          perm_z = mplist.sublist("Hydraulic Conductivity: Anisotropic Uniform").get<double>("z");
        } else {
          Exceptions::amanzi_throw(Errors::Message("Hydraulic Conductivity: Anisotropic Uniform defines a value for z, while the spatial dimension of the problem not 3."));
        }
      }
      // now scale with rho, g and mu to get correct permeablity values
      perm_x *= viscosity/(density*GRAVITY_MAGNITUDE);
      perm_y *= viscosity/(density*GRAVITY_MAGNITUDE);
      perm_z *= viscosity/(density*GRAVITY_MAGNITUDE);
    } else if (mplist.isSublist("Intrinsic Permeability: File")) {
      Teuchos::ParameterList &aux_list = mplist.sublist("Intrinsic Permeability: File");
      bool input_error = false;
      if (aux_list.isParameter("File")) {
        perm_file = aux_list.get<std::string>("File");
      } else {
        input_error = true;
      }
      if (aux_list.isParameter("Attribute")) {
        perm_attribute = aux_list.get<std::string>("Attribute");
      } else {
        input_error = true;
      }
      if (aux_list.isParameter("Format")) {
        perm_format = mplist.sublist("Intrinsic Permeability: File").get<std::string>("Format");
      } else {
        input_error = true;
      }
      if (input_error) {
        Exceptions::amanzi_throw(Errors::Message("The list 'Input Permeability: File' could not be parsed, a required parameter is missing. Check the input specification."));
      }
      perm_init_from_file = true;
    } else {
      Exceptions::amanzi_throw(Errors::Message("Permeability can only be specified as 'Intrinsic Permeability: Uniform', 'Intrinsic Permeability: Anisotropic Uniform', 'Hydraulic Conductivity: Uniform', 'Hydraulic Conductivity: Anisotropic Uniform', or 'Intrinsic Permeability: File'"));
    }

    Teuchos::ParameterList &permeability_ic = stt_ic.sublist("permeability");
    if (perm_init_from_file) {
      if (perm_format == std::string("exodus")) {
        // first make sure the file actually exists

        boost::filesystem::path p;
        if (numproc_>1) {
          // attach the right extensions as required by Nemesis file naming conventions
          // in which files are named as mymesh.par.N.r where N = numproc and r is rank
          // and check if files exist
          std::string perm_file_par = perm_file.substr(0,perm_file.size()-4) + std::string(".par");
          int rank = numproc_-1;
          int ndigits = static_cast<int>(floor(log10(numproc_))) + 1;
          std::string fmt = boost::str(boost::format("%%s.%%d.%%0%dd") % ndigits);
          std::string par_file_w_ext = boost::str(boost::format(fmt) % perm_file_par % numproc_ % rank);
          p = par_file_w_ext;
          if (!boost::filesystem::exists(p)) {
            Exceptions::amanzi_throw(Errors::Message("Permeability initialization from file: not all the partitioned files '" + perm_file_par + ".<numranks>.<rank>' exist, and possibly none of them exist."));
          }
          perm_file = perm_file_par;
        } else { // numproc_ == 1
          std::string suffix = perm_file.substr(perm_file.size()-4);
          if (suffix != std::string(".exo")) {          
            Exceptions::amanzi_throw(Errors::Message("Permeability initialization from file: in serial the exodus mesh file must have the suffix .exo"));
          }
          p = perm_file;
          if (!boost::filesystem::exists(p)) {
            Exceptions::amanzi_throw(Errors::Message("Permeability initialization from file: the file '" + perm_file + "' does not exist."));
          }
        }
        permeability_ic.sublist("exodus file initialization")
            .set<std::string>("file",perm_file)
            .set<std::string>("attribute",perm_attribute);
      } else {
        Exceptions::amanzi_throw(Errors::Message("Permeabily initialization from file, incompatible format specified: '" + perm_format + "', only 'exodus' is supported."));
      }
     
    } else {
      Teuchos::ParameterList& aux_list =
        permeability_ic.sublist("function").sublist(reg_str)
        .set<Teuchos::Array<std::string> >("regions",regions)
        .set<std::string>("component","cell")
        .sublist("function");
      aux_list.set<int>("Number of DoFs",spatial_dimension_)
        .set<std::string>("Function type","composite function");
      aux_list.sublist("DoF 1 Function").sublist("function-constant")
        .set<double>("value", perm_x);
      if (spatial_dimension_ >= 2) {
        aux_list.sublist("DoF 2 Function").sublist("function-constant")
          .set<double>("value", perm_y);
      }
      if (spatial_dimension_ == 3) {
        aux_list.sublist("DoF 3 Function").sublist("function-constant")
          .set<double>("value", perm_z);
      }
    }

    // specific_yield...
    if (matprop_list.sublist(matprop_list.name(i)).isSublist("Specific Yield: Uniform")) {
      Teuchos::ParameterList& spec_yield_ic = stt_ic.sublist("specific_yield");
      double specific_yield = matprop_list.sublist(matprop_list.name(i)).sublist("Specific Yield: Uniform").get<double>("Value");
      spec_yield_ic.sublist("function").sublist(reg_str)
          .set<Teuchos::Array<std::string> >("regions",regions)
          .set<std::string>("component","cell")
          .sublist("function").sublist("function-constant")
          .set<double>("value", specific_yield);
    }

    // specific_storage...
    if (matprop_list.sublist(matprop_list.name(i)).isSublist("Specific Storage: Uniform")) {
      Teuchos::ParameterList& spec_stor_ic = stt_ic.sublist("specific_storage");
      double specific_storage = matprop_list.sublist(matprop_list.name(i)).sublist("Specific Storage: Uniform").get<double>("Value");
      spec_stor_ic.sublist("function").sublist(reg_str)
          .set<Teuchos::Array<std::string> >("regions",regions)
          .set<std::string>("component","cell")
          .sublist("function").sublist("function-constant")
          .set<double>("value", specific_storage);
    }

    // particle_density...
    if (matprop_list.sublist(matprop_list.name(i)).isSublist("Particle Density: Uniform")) {
      Teuchos::ParameterList& part_dens_ic = stt_ic.sublist("particle_density");
      double particle_density = matprop_list.sublist(matprop_list.name(i)).sublist("Particle Density: Uniform").get<double>("Value");
      part_dens_ic.sublist("function").sublist(reg_str)
          .set<Teuchos::Array<std::string> >("regions",regions)
          .set<std::string>("component","cell")
          .sublist("function").sublist("function-constant")
          .set<double>("value", particle_density);
    }
  }

  // and now the initialization of fields that are initialized via an Initial Condition in the Akuna spec

  // loop over the initial conditions
  Teuchos::ParameterList& ic_list = plist->sublist("Initial Conditions");

  if (! ic_list.isParameter("Init from Checkpoint File")) {
    // only process initial conditions if we are not initializing from
    // a checkpoint file
    for (Teuchos::ParameterList::ConstIterator iic = ic_list.begin(); iic != ic_list.end(); ++iic) {
      // get the regions
      Teuchos::Array<std::string> regions = ic_list.sublist(ic_list.name(iic)).get<Teuchos::Array<std::string> >("Assigned Regions");

      Teuchos::ParameterList* ic_for_region = &ic_list.sublist(ic_list.name(iic));

      // create regions string
      std::string reg_str;
      for (Teuchos::Array<std::string>::const_iterator ireg=regions.begin(); ireg!=regions.end(); ++ireg) {
        reg_str = reg_str + *ireg;
      }

      // pressure...
      if (ic_for_region->isSublist("IC: Uniform Pressure") || ic_for_region->isSublist("IC: Linear Pressure")) {

        Teuchos::ParameterList& pressure_ic = stt_ic.sublist("pressure");
        if (ic_for_region->isSublist("IC: Uniform Pressure")) {
          double p = ic_for_region->sublist("IC: Uniform Pressure").get<double>("Value");

          pressure_ic.sublist("function").sublist(reg_str)
              .set<Teuchos::Array<std::string> >("regions",regions)
              .set<std::string>("component","cell")
              .sublist("function").sublist("function-constant")
              .set<double>("value", p);

        } else if (ic_for_region->isSublist("IC: Linear Pressure")) {
          Teuchos::Array<double> grad = ic_for_region->sublist("IC: Linear Pressure").get<Teuchos::Array<double> >("Gradient Value");
          Teuchos::Array<double> refcoord = ic_for_region->sublist("IC: Linear Pressure").get<Teuchos::Array<double> >("Reference Point");
          double refval = ic_for_region->sublist("IC: Linear Pressure").get<double>("Reference Value");

          Teuchos::Array<double> grad_with_time(grad.size()+1);
          grad_with_time[0] = 0.0;
          for (int j = 0; j != grad.size(); ++j) {
            grad_with_time[j + 1] = grad[j];
          }

          Teuchos::Array<double> refcoord_with_time(refcoord.size()+1);
          refcoord_with_time[0] = 0.0;
          for (int j = 0; j != refcoord.size(); ++j) {
            refcoord_with_time[j + 1] = refcoord[j];
          }

          pressure_ic.sublist("function").sublist(reg_str)
              .set<Teuchos::Array<std::string> >("regions",regions)
              .set<std::string>("component","cell")
              .sublist("function").sublist("function-linear")
              .set<double>("y0", refval)
              .set<Teuchos::Array<double> >("x0",refcoord_with_time)
              .set<Teuchos::Array<double> >("gradient",grad_with_time);

        } else if (ic_for_region->isSublist("IC: File Pressure")) {
          Exceptions::amanzi_throw(Errors::Message("IC: File Pressure cannot currently be used to initialize pressure in a region."));
        }
      }

      // saturation...
      if (ic_for_region->isSublist("IC: Uniform Saturation") || ic_for_region->isSublist("IC: Linear Saturation")) {
        Teuchos::ParameterList& saturation_ic = stt_ic.sublist("water_saturation");
        if (ic_for_region->isSublist("IC: Uniform Saturation")) {
          double s = ic_for_region->sublist("IC: Uniform Saturation").get<double>("Value");
          saturation_ic.sublist("function").sublist(reg_str)
              .set<Teuchos::Array<std::string> >("regions",regions)
              .set<std::string>("component","cell")
              .sublist("function").sublist("function-constant")
              .set<double>("value", s);
        } else if (ic_for_region->isSublist("IC: Linear Saturation")) {
          Teuchos::ParameterList& saturation_ic = ic_for_region->sublist("IC: Linear Saturation");
          Teuchos::Array<double> grad = saturation_ic.get<Teuchos::Array<double> >("Gradient Value");
          Teuchos::Array<double> refcoord = saturation_ic.get<Teuchos::Array<double> >("Reference Point");
          double refval = saturation_ic.get<double>("Reference Value");

          Teuchos::Array<double> grad_with_time(grad.size() + 1);
          grad_with_time[0] = 0.0;
          for (int j = 0; j != grad.size(); ++j) {
            grad_with_time[j + 1] = grad[j];
          }

          Teuchos::Array<double> refcoord_with_time(refcoord.size()+1);
          refcoord_with_time[0] = 0.0;
          for (int j = 0; j != refcoord.size(); ++j) {
            refcoord_with_time[j + 1] = refcoord[j];
          }

          saturation_ic.sublist("function").sublist(reg_str)
              .set<Teuchos::Array<std::string> >("regions",regions)
              .set<std::string>("component","cell")
              .sublist("function").sublist("function-linear")
              .set<double>("y0", refval)
              .set<Teuchos::Array<double> >("x0",refcoord_with_time)
              .set<Teuchos::Array<double> >("gradient",grad_with_time);
        } else if (ic_for_region->isSublist("IC: File Saturation")) {
          Exceptions::amanzi_throw(Errors::Message("IC: File Saturation cannot currently be used to initialize saturation in a region."));
        }
      }

      // darcy_flux...
      if (ic_for_region->isSublist("IC: Uniform Velocity")) {
        Teuchos::ParameterList &darcy_flux_ic = stt_ic.sublist("darcy_flux");
        Teuchos::Array<double> vel_vec = ic_for_region->sublist("IC: Uniform Velocity").get<Teuchos::Array<double> >("Velocity Vector");

        if (vel_vec.size() != spatial_dimension_) {
          Exceptions::amanzi_throw(Errors::Message("The velocity vector defined in in the IC: Uniform Velocity list does not match the spatial dimension of the problem."));
        }

        Teuchos::ParameterList& d_list =
            darcy_flux_ic.set<bool>("dot with normal", true)
            .sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions",regions)
            .set<std::string>("component","face")
            .sublist("function")
            .set<int>("Number of DoFs", vel_vec.size())
            .set<std::string>("Function type", "composite function");

        for (int ii = 0; ii != vel_vec.size(); ++ii) {
          std::stringstream dof_str;
          dof_str << "DoF " << ii+1 << " Function";

          d_list.sublist(dof_str.str())
              .sublist("function-constant")
              .set<double>("value", vel_vec[ii]);
        }
      }

      // total_component_concentration...
      if (ic_for_region->isSublist("Solute IC")) {

        Teuchos::ParameterList &concentration_ic = stt_ic.sublist("total_component_concentration");
        if (plist->sublist("Execution Control").get<std::string>("Transport Model") != std::string("Off")  ||
            plist->sublist("Execution Control").get<std::string>("Chemistry Model") != std::string("Off")) {
          // write the initial conditions for the solutes, note that we hardcode for there only being one phase, with one phase component

          Teuchos::ParameterList& dof_list = concentration_ic.sublist("function").sublist(reg_str)
              .set<Teuchos::Array<std::string> >("regions",regions)
              .set<std::string>("component","cell")
              .sublist("function")
              .set<int>("Number of DoFs", comp_names.size())
              .set<std::string>("Function type", "composite function");


          for (int ii = 0; ii < comp_names.size(); ii++) {
            if (ic_for_region->sublist("Solute IC").sublist(phase_name).sublist(phase_comp_name).isSublist(comp_names[ii])) {
              Teuchos::ParameterList& conc_ic = ic_for_region->sublist("Solute IC").sublist(phase_name).sublist(phase_comp_name).sublist(comp_names[ii]).sublist("IC: Uniform Concentration");

              if (conc_ic.isParameter("Geochemical Condition")) { // Geochemical condition?
                // Add an entry to State->initial conditions->geochemical conditions.
                Teuchos::ParameterList& geochem_cond_list = stt_ic.sublist("geochemical conditions");
                std::string geochem_cond_name = conc_ic.get<std::string>("Geochemical Condition");
                Teuchos::ParameterList& geochem_cond = geochem_cond_list.sublist(geochem_cond_name);
                geochem_cond.set<Teuchos::Array<std::string> >("regions", regions);

                // Now add fake initial concentrations, since the State requires entries 
                // of this sort to initialize fields. The chemistry PK will simply 
                // overwrite these values when it initializes its data.
                std::stringstream dof_str;
                dof_str << "DoF " << ii+1 << " Function";
                dof_list.sublist(dof_str.str())
                  .sublist("function-constant")
                  .set<double>("value",0.0);
              }
              else { // ordinary initial concentration value.
                double conc = conc_ic.get<double>("Value");

                std::stringstream dof_str;
                dof_str << "DoF " << ii+1 << " Function";

                dof_list.sublist(dof_str.str())
                  .sublist("function-constant")
                  .set<double>("value",conc);
              }
            }
          }
        }
      }
    }
  }

  // add mesh partitions to the state list
  stt_list.sublist("mesh partitions") = CreatePartitionList_(plist);

  return stt_list;
}


/* ******************************************************************
* Mesh patition sublist based on materials
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreatePartitionList_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList part_list;
  Teuchos::ParameterList& tmp_list = part_list.sublist("materials");

  Teuchos::Array<std::string> regions;
  if (plist->isSublist("Material Properties")) {
    for (Teuchos::ParameterList::ConstIterator it = plist->sublist("Material Properties").begin();
         it != plist->sublist("Material Properties").end(); ++it) {
      if ((it->second).isList()) {
        Teuchos::ParameterList& mat_sublist = plist->sublist("Material Properties").sublist(it->first);
        Teuchos::Array<std::string> names = mat_sublist.get<Teuchos::Array<std::string> >("Assigned Regions");

        for (int i = 0; i < names.length(); i++) {
          regions.append(names[i]);
        } 
      }
    }
  }
  tmp_list.set<Teuchos::Array<std::string> >("region list", regions);
  
  return part_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
