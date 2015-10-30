#include <sstream>
#include <string>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputAnalysis.hh"

namespace Amanzi {

/* ******************************************************************
* Initilization.
****************************************************************** */
void InputAnalysis::Init(Teuchos::ParameterList& plist) 
{
  plist_ = &plist;

  if (plist.isSublist("Analysis")) {
    Teuchos::ParameterList vo_list = plist.sublist("Analysis");
    vo_ = new VerboseObject("InputAnalysis", vo_list); 
  } 
}


/* ******************************************************************
* Analysis of collected regions
****************************************************************** */
void InputAnalysis::RegionAnalysis() 
{
  if (!plist_->isSublist("Analysis")) return;
  Teuchos::ParameterList alist = plist_->sublist("Analysis");

  Errors::Message msg;
  Teuchos::OSTab tab = vo_->getOSTab();

  if (alist.isParameter("used source regions")) {
    std::vector<std::string> regions = alist.get<Teuchos::Array<std::string> >("used source regions").toVector();

    for (int i = 0; i < regions.size(); i++) {
      AmanziMesh::Entity_ID_List block;
      mesh_->get_set_entities(regions[i], AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      int nblock = block.size();

      double volume(0.0);
      for (int n = 0; n < nblock; n++) {
        volume += mesh_->cell_volume(block[n]);
      }
#ifdef HAVE_MPI
      int nblock_tmp = nblock;
      double volume_tmp = volume;
      mesh_->get_comm()->SumAll(&nblock_tmp, &nblock, 1);
      mesh_->get_comm()->SumAll(&volume_tmp, &volume, 1);
#endif

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
        std::string name(regions[i]);
        name.resize(std::min(40, (int)name.size()));
        *vo_->os() << "src: \"" << name << "\" has " << nblock << " cells" 
                   << " of " << volume << " [m^3]" << std::endl;
      }

      if (nblock == 0) {
        msg << "Used source region is empty.";
        Exceptions::amanzi_throw(msg);
      } 
      // AmanziMesh::Entity_ID_List::iterator i;
      // for (i = block.begin(); i != block.end(); i++) {};
    }
  }

  if (alist.isParameter("used boundary condition regions")) {
    std::vector<std::string> regions = alist.get<Teuchos::Array<std::string> >("used boundary condition regions").toVector();

    for (int i = 0; i < regions.size(); i++) {
      AmanziMesh::Entity_ID_List block;
      mesh_->get_set_entities(regions[i], AmanziMesh::FACE, AmanziMesh::OWNED, &block);
      int nblock = block.size();

      double area(0.0);
      for (int n = 0; n < nblock; n++) {
        area += mesh_->face_area(block[n]);
      }
#ifdef HAVE_MPI
      int nblock_tmp = nblock;
      double area_tmp = area;
      mesh_->get_comm()->SumAll(&nblock_tmp, &nblock, 1);
      mesh_->get_comm()->SumAll(&area_tmp, &area, 1);
#endif
      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
        std::string name(regions[i]);
        name.resize(std::min(40, (int)name.size()));
        *vo_->os() << "bc: \"" << name << "\" has " << nblock << " faces"
                   << " of " << area << " [m^2]" << std::endl;
      }

      if (nblock == 0) {
        msg << "Used boundary condition region is empty.";
        Exceptions::amanzi_throw(msg);
      } 
    }
  }

  if (alist.isParameter("used observation regions")) {
    std::vector<std::string> regions = alist.get<Teuchos::Array<std::string> >("used observation regions").toVector();

    int nblock;
    for (int i = 0; i < regions.size(); i++) {
      double volume(0.0);
      AmanziMesh::Entity_ID_List block;
      std::string type;

      if (mesh_->valid_set_name(regions[i], AmanziMesh::CELL)) {
        mesh_->get_set_entities(regions[i], AmanziMesh::CELL, AmanziMesh::OWNED, &block);
        nblock = block.size();
        type = "cells";
        for (int n = 0; n < nblock; n++) 
            volume += mesh_->cell_volume(block[n]);
      }
      else if (mesh_->valid_set_name(regions[i], AmanziMesh::FACE)) {
        mesh_->get_set_entities(regions[i], AmanziMesh::FACE, AmanziMesh::OWNED, &block);
        nblock = block.size();
        type = "faces";
        for (int n = 0; n < nblock; n++) 
            volume += mesh_->face_area(block[n]);
      } 
      else {
        nblock = 0;
        std::string name(regions[i]);
        name.resize(std::min(40, (int)name.size()));
        *vo_->os() << "Observation region: \"" << name << "\" has unknown type." << std::endl;
      }

#ifdef HAVE_MPI
      int nblock_tmp = nblock;
      double volume_tmp = volume;
      mesh_->get_comm()->SumAll(&nblock_tmp, &nblock, 1);
      mesh_->get_comm()->SumAll(&volume_tmp, &volume, 1);
#endif

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
        std::string name(regions[i]);
        name.resize(std::min(40, (int)name.size()));
        *vo_->os() << "obs: \"" << name << "\" has " << nblock << " " << type 
                   << ", size: " << volume << std::endl;
      }

      if (nblock == 0) {
        msg << "Used observation region is empty.";
        Exceptions::amanzi_throw(msg);
      } 
    }
  }
}


/* ******************************************************************
* DEBUG output: boundary conditions.
****************************************************************** */
void InputAnalysis::OutputBCs() 
{
  if (!plist_->isSublist("Analysis")) return;
  if (vo_->getVerbLevel() < Teuchos::VERB_EXTREME) return;

  int bc_counter = 0;

  if (plist_->isSublist("Flow")) {
    Teuchos::ParameterList& flow_list = plist_->sublist("Flow");

    Teuchos::ParameterList richards_list;
    if (flow_list.isSublist("Richards problem")) {
      richards_list = flow_list.sublist("Richards problem");

      Teuchos::ParameterList bc_list;
      if (richards_list.isSublist("boundary conditions")) {
        bc_list = richards_list.sublist("boundary conditions");
      }

      Teuchos::ParameterList mass_flux_list, pressure_list, seepage_list, head_list;
      if (bc_list.isSublist("mass flux")) {
        mass_flux_list = bc_list.sublist("mass flux");
        for (Teuchos::ParameterList::ConstIterator i = mass_flux_list.begin(); i != mass_flux_list.end(); i++) {

          if (mass_flux_list.isSublist(mass_flux_list.name(i))) {
            Teuchos::ParameterList& bc = mass_flux_list.sublist(mass_flux_list.name(i));

            if ((bc.sublist("outward mass flux")).isSublist("function-tabular")) {
              std::stringstream ss;
              ss << "BCmassflux" << bc_counter++;

              Teuchos::ParameterList& f_tab = (bc.sublist("outward mass flux")).sublist("function-tabular");

              Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double> >("x values");
              Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double> >("y values");
              Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string> >("forms");

              int np = times.size()*2 - 1;
              Teuchos::Array<double> times_plot(np);
              Teuchos::Array<double> values_plot(np);

              for (int i = 0; i < times.size() - 1; i++) {
                times_plot[2*i] = times[i];
                values_plot[2*i] = values[i];
                times_plot[2*i + 1] = 0.5*(times[i] + times[i+1]);
              }
              times_plot[np - 1] = times[times.size() - 1];
              values_plot[np - 1] = values[times.size() - 1];

              for (int i = 0; i < time_fns.size(); i++) {
                if (time_fns[i] == "linear") {
                  values_plot[2*i + 1] = 0.5 *( values[i] + values[i+1]);
                } else if (time_fns[i] == "constant") {
                  values_plot[2*i + 1] = values[i];
                  times_plot[2*i + 1] = times[i+1];
                } else {
                  Exceptions::amanzi_throw(Errors::Message("In the definition of BCs: tabular function can only be Linear or Constant"));
                }
              }

              std::string filename = ss.str() + ".dat";
              std::ofstream ofile(filename.c_str());

              ofile << "# "<<"time "<< "flux"<<std::endl;
              for (int i = 0; i < np; i++) {
                ofile <<times_plot[i] << " " << values_plot[i] << std::endl;
              }

              ofile.close();
            }
          }
        }
      }
      if (bc_list.isSublist("pressure")) {
        pressure_list = bc_list.sublist("pressure");
        for (Teuchos::ParameterList::ConstIterator i = pressure_list.begin(); i != pressure_list.end(); i++) {

          if (pressure_list.isSublist(pressure_list.name(i))) {
            Teuchos::ParameterList& bc = pressure_list.sublist(pressure_list.name(i));
            if ((bc.sublist("boundary pressure")).isSublist("function-tabular")) {
              std::stringstream ss;
              ss << "BCpressure" << bc_counter++;


              Teuchos::ParameterList& f_tab = (bc.sublist("boundary pressure")).sublist("function-tabular");

              Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double> >("x values");
              Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double> >("y values");
              Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string> >("forms");

              int np = times.size()*2 - 1;
              Teuchos::Array<double> times_plot(np);
              Teuchos::Array<double> values_plot(np);

              for (int i = 0; i < times.size() - 1; i++) {
                times_plot[2*i] = times[i];
                values_plot[2*i] = values[i];
                times_plot[2*i + 1] = 0.5*(times[i] + times[i+1]);
              }
              times_plot[np - 1] = times[times.size() - 1];
              values_plot[np - 1] = values[times.size() - 1];

              for (int i = 0; i<time_fns.size(); i++) {
                if (time_fns[i] == "linear") {
                  values_plot[2*i + 1] = 0.5 *( values[i] + values[i+1]);
                } else if (time_fns[i] == "constant") {
                  values_plot[2*i + 1] = values[i];
                  times_plot[2*i + 1] = times[i+1];
                } else {
                  Exceptions::amanzi_throw(Errors::Message("In the definition of BCs: tabular function can only be Linear or Constant"));
                }
              }

              std::string filename = ss.str() + ".dat";
              std::ofstream ofile(filename.c_str());

              ofile << "# time "<<"pressure"<<std::endl;
              for (int i = 0; i < np; i++) {
                ofile << times_plot[i] << " " << values_plot[i] << std::endl;
              }

              ofile.close();
            }
          }
        }
      }

      if (bc_list.isSublist("seepage face")) {
        seepage_list = bc_list.sublist("seepage face");
        for (Teuchos::ParameterList::ConstIterator i = seepage_list.begin(); i != seepage_list.end(); i++) {

          if (seepage_list.isSublist(seepage_list.name(i))) {
            Teuchos::ParameterList& bc = seepage_list.sublist(seepage_list.name(i));
            if ((bc.sublist("outward mass flux")).isSublist("function-tabular")) {
              std::stringstream ss;
              ss << "BCseepage" << bc_counter++;

              Teuchos::ParameterList& f_tab = (bc.sublist("outward mass flux")).sublist("function-tabular");

              Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double> >("x values");
              Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double> >("y values");
              Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string> >("forms");

              int np = times.size()*2 - 1;
              Teuchos::Array<double> times_plot(np);
              Teuchos::Array<double> values_plot(np);

              for (int i = 0; i < times.size() - 1; i++) {
                times_plot[2*i] = times[i];
                values_plot[2*i] = values[i];
                times_plot[2*i + 1] = 0.5*(times[i] + times[i+1]);
              }
              times_plot[np - 1] = times[times.size() - 1];
              values_plot[np - 1] = values[times.size() - 1];

              for (int i = 0; i<time_fns.size(); i++) {
                if (time_fns[i] == "linear") {
                  values_plot[2*i + 1] = 0.5 *( values[i] + values[i+1]);
                } else if (time_fns[i] == "constant") {
                  values_plot[2*i + 1] = values[i];
                  times_plot[2*i + 1] = times[i+1];
                } else {
                  Exceptions::amanzi_throw(Errors::Message("In the definition of BCs: tabular function can only be Linear or Constant"));
                }
              }

              std::string filename = ss.str() + ".dat";
              std::ofstream ofile(filename.c_str());

              ofile << "# time " << "flux" << std::endl;
              for (int i = 0; i < np; i++) {
                ofile << times_plot[i] << " " << values_plot[i] << std::endl;
              }

              ofile.close();
            }
          }
        }
      }

      if (bc_list.isSublist("static head")) {
        head_list = bc_list.sublist("static head");
        for (Teuchos::ParameterList::ConstIterator i = head_list.begin(); i != head_list.end(); i++) {

          if (head_list.isSublist(head_list.name(i))) {
            Teuchos::ParameterList& bc = head_list.sublist(head_list.name(i));
            if ((bc.sublist("water table elevation")).isSublist("function-tabular")) {
              std::stringstream ss;
              ss << "BChead" << bc_counter++;

              Teuchos::ParameterList& f_tab = (bc.sublist("water table elevation")).sublist("function-tabular");

              Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double> >("x values");
              Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double> >("y values");
              Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string> >("forms");

              int np = times.size()*2 - 1;
              Teuchos::Array<double> times_plot(np);
              Teuchos::Array<double> values_plot(np);

              for (int i = 0; i < times.size() - 1; i++) {
                times_plot[2*i] = times[i];
                values_plot[2*i] = values[i];
                times_plot[2*i + 1] = 0.5*(times[i] + times[i+1]);
              }
              times_plot[np - 1] = times[times.size() - 1];
              values_plot[np - 1] = values[times.size() - 1];

              for (int i = 0; i < time_fns.size(); i++) {
                if (time_fns[i] == "linear") {
                  values_plot[2*i + 1] = 0.5 *( values[i] + values[i+1]);
                } else if (time_fns[i] == "constant") {
                  values_plot[2*i + 1] = values[i];
                  times_plot[2*i + 1] = times[i+1];
                } else {
                  Exceptions::amanzi_throw(Errors::Message("In the definition of BCs: tabular function can only be Linear or Constant"));
                }
              }

              std::string filename = ss.str() + ".dat";
              std::ofstream ofile(filename.c_str());

              ofile << "# time " << "head" << std::endl;
              for (int i = 0; i < np; i++) {
                ofile << times_plot[i] << " " << values_plot[i] << std::endl;
              }

              ofile.close();
            }
          }
        }
      }
    }
  }
}

}  // namespace Amanzi
