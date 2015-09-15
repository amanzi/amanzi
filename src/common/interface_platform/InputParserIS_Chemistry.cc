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
* Chemisty sublist
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateChemistryList_(Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  Teuchos::ParameterList chem_list;
  if (plist->isSublist("Chemistry")) {
    chem_list = plist->sublist("Chemistry");
    chem_list.set<std::string>("chemistry model", chemistry_model_);

    chem_list.set<std::string>("activity model", chem_list.get<std::string>("Activity Model", "unit"));
    chem_list.set<int>("maximum Newton iterations", chem_list.get<int>("Maximum Newton Iterations", 100));
    chem_list.set<double>("tolerance", chem_list.get<double>("Tolerance", 1.0e-12));
    chem_list.set<double>("max time step (s)", chem_list.get<double>("Max Time Step (s)", 1.0e+10));
    chem_list.set<double>("min time step (s)", chem_list.get<double>("Min Time Step (s)", 1.0e-10));
    chem_list.set<double>("initial time step (s)", chem_list.get<double>("Initial Time Step (s)", 1.0e-10));
    chem_list.set<std::string>("time step control method",
        chem_list.get<std::string>("Time Step Control Method", "fixed"));
    chem_list.set<int>("time step cut threshold", chem_list.get<int>("Time Step Cut Threshold", 8));
    chem_list.set<double>("time step cut factor", chem_list.get<double>("Time Step Cut Factor", 2.0));
    chem_list.set<int>("time step increase threshold", chem_list.get<int>("Time Step Increase Threshold", 4));
    chem_list.set<double>("time step increase factor", chem_list.get<double>("Time Step Increase Factor", 1.2));
    if (chem_list.isParameter("Auxiliary Data"))
        chem_list.set<Teuchos::Array<std::string> >("auxiliary data",
        chem_list.get<Teuchos::Array<std::string> >("Auxiliary Data"));
   
    Teuchos::ParameterList& chem_ic = chem_list.sublist("initial conditions");

    //
    // read the minerals
    //

    // Teuchos::Array<std::string> minerals;
    // if (plist->sublist("Phase Definitions").isSublist("Solid")) {
    //minerals = plist->sublist("Phase Definitions").sublist("Solid")
    // .get<Teuchos::Array<std::string> >("Minerals");
    if (mineral_names_.size() > 0) {
      chem_list.set<Teuchos::Array<std::string> >("Minerals", mineral_names_);
    }
    if (sorption_site_names_.size() > 0) {
      chem_list.set<Teuchos::Array<std::string> >("Sorption Sites", sorption_site_names_);
    }

    chem_list.set<int>("number of component concentrations", comp_names_all_.size());

    //
    // --- region specific initial conditions
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

      if (mineral_names_.size() > 0) {
        double mvf(0.0), msa(0.0);

        // mineral volume fractions
        Teuchos::ParameterList &mineral_volfrac_ic = chem_ic.sublist("mineral_volume_fractions");
        // mineral specific surface area
        Teuchos::ParameterList &mineral_surfarea_ic = chem_ic.sublist("mineral_specific_surface_area");

        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {

          // mineral volume fractions and specific surface area

          Teuchos::ParameterList& aux1_list =
              mineral_volfrac_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          Teuchos::ParameterList& aux2_list =
              mineral_surfarea_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          aux1_list.set<int>("Number of DoFs", mineral_names_.size())
              .set("Function type", "composite function");

          aux2_list.set<int>("Number of DoFs", mineral_names_.size())
              .set("Function type", "composite function");

          for (int j = 0; j<mineral_names_.size(); ++j) {
            std::stringstream ss;
            ss << "DoF " << j + 1 << " Function";

            mvf = matprop_list.sublist(matprop_list.name(i)).sublist("Mineralogy").sublist(mineral_names_[j])
                .get<double>("Volume Fraction");

            msa = matprop_list.sublist(matprop_list.name(i)).sublist("Mineralogy").sublist(mineral_names_[j])
                .get<double>("Specific Surface Area");

            aux1_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", mvf);

            aux2_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", msa);
          }
        }
      }


      if (matprop_list.sublist(matprop_list.name(i)).isParameter("Cation Exchange Capacity")) {

        double cec = matprop_list.sublist(matprop_list.name(i)).get<double>("Cation Exchange Capacity");

        Teuchos::ParameterList &ion_exchange_sites_ic = chem_ic.sublist("ion_exchange_sites");
        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {
          Teuchos::ParameterList& aux1_list =
              ion_exchange_sites_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          // this needs to be more general... for now we initialize with one DoF
          aux1_list.set<int>("Number of DoFs", 1)
              .set("Function type", "composite function");

          aux1_list.sublist("DoF 1 Function").sublist("function-constant")
              .set<double>("value", cec);
        }

        Teuchos::ParameterList &ion_exchange_ref_cation_conc_ic = chem_ic.sublist("ion_exchange_ref_cation_conc");
        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {
          Teuchos::ParameterList& aux1_list =
              ion_exchange_ref_cation_conc_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          // this needs to be more general... for now we initialize with one DoF
          aux1_list.set<int>("Number of DoFs", 1)
              .set("Function type", "composite function");

          aux1_list.sublist("DoF 1 Function").sublist("function-constant")
              .set<double>("value", 1.0);  // this should be read from the input file??... TODO
        }
      }

      if (matprop_list.sublist(matprop_list.name(i)).isSublist("Sorption Isotherms")) {

        Teuchos::ParameterList& isotherm_kd_ic = chem_ic.sublist("isotherm_kd");
        Teuchos::ParameterList& isotherm_langmuir_b_ic = chem_ic.sublist("isotherm_langmuir_b");
        Teuchos::ParameterList& isotherm_freundlich_n_ic = chem_ic.sublist("isotherm_freundlich_n");

        Teuchos::ParameterList& sorption_isotherms_list = matprop_list.sublist(matprop_list.name(i)).sublist("Sorption Isotherms");

        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {

          // Kd

          {
            Teuchos::ParameterList& aux1_list =
                isotherm_kd_ic.sublist("function").sublist(*ir)
                .set<std::string>("region",*ir)
                .set<std::string>("component","cell")
                .sublist("function");

            aux1_list.set<int>("Number of DoFs", comp_names_.size())
                .set("Function type", "composite function");

            for ( int ic = 0; ic != comp_names_.size(); ++ic) {

              std::stringstream ss;
              ss << "DoF " << ic + 1 << " Function";

              double kd(0.0);
              if (sorption_isotherms_list.isSublist(comp_names_[ic])) {
                kd = sorption_isotherms_list.sublist(comp_names_[ic]).get<double>("Kd",0.0);
              }

              aux1_list.sublist(ss.str()).sublist("function-constant")
                  .set<double>("value", kd);

            }
          }

          // Langmuir
          {
            Teuchos::ParameterList& aux2_list =
                isotherm_langmuir_b_ic.sublist("function").sublist(*ir)
                .set<std::string>("region",*ir)
                .set<std::string>("component","cell")
                .sublist("function");

            aux2_list.set<int>("Number of DoFs", comp_names_.size())
                .set("Function type", "composite function");

            for ( int ic = 0; ic != comp_names_.size(); ++ic) {

              std::stringstream ss;
              ss << "DoF " << ic + 1 << " Function";

              double langmuir_b(1.0);
              if (sorption_isotherms_list.isSublist(comp_names_[ic])) {
                langmuir_b = sorption_isotherms_list.sublist(comp_names_[ic]).get<double>("Langmuir b",1.0);
              }

              aux2_list.sublist(ss.str()).sublist("function-constant")
                  .set<double>("value", langmuir_b);

            }
          }

          // Freundlich
          {
            Teuchos::ParameterList& aux3_list =
                isotherm_freundlich_n_ic.sublist("function").sublist(*ir)
                .set<std::string>("region",*ir)
                .set<std::string>("component","cell")
                .sublist("function");

            aux3_list.set<int>("Number of DoFs", comp_names_.size())
                .set("Function type", "composite function");

            for ( int ic = 0; ic != comp_names_.size(); ++ic) {

              std::stringstream ss;
              ss << "DoF " << ic + 1 << " Function";

              double freundlich_n(1.0);
              if (sorption_isotherms_list.isSublist(comp_names_[ic])) {
                freundlich_n = sorption_isotherms_list.sublist(comp_names_[ic]).get<double>("Freundlich n",1.0);
              }

              aux3_list.sublist(ss.str()).sublist("function-constant")
                  .set<double>("value", freundlich_n);
            }
          }
        }
      }

      if (sorption_site_names_.size() > 0) {
        Teuchos::ParameterList& sorption_sites_ic = chem_ic.sublist("sorption_sites");

        Teuchos::ParameterList& sorption_sites_list = matprop_list.sublist(matprop_list.name(i)).sublist("Surface Complexation Sites");

        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {

          Teuchos::ParameterList& aux1_list =
              sorption_sites_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          aux1_list.set<int>("Number of DoFs", sorption_site_names_.size())
              .set("Function type", "composite function");

          for ( int ic = 0; ic != sorption_site_names_.size(); ++ic) {

            std::stringstream ss;
            ss << "DoF " << ic + 1 << " Function";

            double value(0.0);
            if (sorption_sites_list.isSublist(sorption_site_names_[ic])) {
              value = sorption_sites_list.sublist(sorption_site_names_[ic]).get<double>("Site Density",0.0);
            }

            aux1_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", value);
          }
        }
      }
    }

    Teuchos::ParameterList& ic_list = plist->sublist("Initial Conditions");

    for (Teuchos::ParameterList::ConstIterator ic = ic_list.begin(); ic != ic_list.end(); ++ic) {
      if (ic_list.isSublist(ic->first)) {

        Teuchos::ParameterList& ics = ic_list.sublist(ic->first);

        Teuchos::Array<std::string> ass_regions = ics.get<Teuchos::Array<std::string> >("Assigned Regions");

        Teuchos::ParameterList &free_ion_species_ic = chem_ic.sublist("free_ion_species");

        for (Teuchos::Array<std::string>::const_iterator ir=ass_regions.begin(); ir!=ass_regions.end(); ir++) {

          Teuchos::ParameterList& aux1_list =
              free_ion_species_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          aux1_list.set<int>("Number of DoFs", comp_names_.size())
              .set("Function type", "composite function");

          for (int j = 0; j<comp_names_.size(); ++j) {
            std::stringstream ss;
            ss << "DoF " << j + 1 << " Function";

            double value(1.0e-9);
            value = ics.sublist("Solute IC").sublist(phases_[0].name)
                .sublist(phases_[0].solute_name).sublist(comp_names_[j])
                .sublist("IC: Uniform Concentration").get<double>("Free Ion Guess",1.0e-9);

            aux1_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", value);
          }
        }
      }
    }
  }

  chem_list.sublist("VerboseObject") = CreateVerbosityList_(verbosity_level);

  return chem_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
