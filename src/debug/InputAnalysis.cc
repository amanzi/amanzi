#include <fstream>
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

  if (plist.isSublist("analysis")) {
    Teuchos::ParameterList vo_list = plist.sublist("analysis");
    vo_ = new VerboseObject("InputAnalysis:" + domain_, vo_list); 
  } 
}


/* ******************************************************************
* Analysis of collected regions
****************************************************************** */
void InputAnalysis::RegionAnalysis() 
{
  if (!plist_->isSublist("analysis")) return;
  Teuchos::ParameterList alist = plist_->sublist("analysis").sublist(domain_);

  Errors::Message msg;
  Teuchos::OSTab tab = vo_->getOSTab();

  if (alist.isParameter("used source regions")) {
    std::vector<std::string> regions = alist.get<Teuchos::Array<std::string> >("used source regions").toVector();
    regions.erase(SelectUniqueEntries(regions.begin(), regions.end()), regions.end());

    for (int i = 0; i < regions.size(); i++) {
      AmanziMesh::Entity_ID_List block;
      std::vector<double> vofs;
      mesh_->get_set_entities_and_vofs(regions[i], AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block, &vofs);
      int nblock = block.size();
      int nvofs = vofs.size();

      double volume(0.0), frac;
      for (int n = 0; n < nblock; n++) {
        frac = (nvofs == 0) ? 1.0 : vofs[n];
        volume += mesh_->cell_volume(block[n]) * frac;
      }

      double vofs_min(1.0), vofs_max(0.0);
      for (int n = 0; n < nvofs; ++n) {
        vofs_min = std::min(vofs_min, vofs[n]);
        vofs_max = std::max(vofs_max, vofs[n]);
      }
      if (nvofs == 0) vofs_max = 1.0;
#ifdef HAVE_MPI
      int nblock_tmp(nblock), nvofs_tmp(nvofs);
      double volume_tmp(volume), vofs_min_tmp(vofs_min), vofs_max_tmp(vofs_max);

      mesh_->get_comm()->SumAll(&nblock_tmp, &nblock, 1);
      mesh_->get_comm()->SumAll(&nvofs_tmp, &nvofs, 1);
      mesh_->get_comm()->SumAll(&volume_tmp, &volume, 1);
      mesh_->get_comm()->MinAll(&vofs_min_tmp, &vofs_min, 1);
      mesh_->get_comm()->MaxAll(&vofs_max_tmp, &vofs_max, 1);
#endif

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
        std::string name(regions[i]);
        name.resize(std::min(40, (int)name.size()));
        *vo_->os() << "src: \"" << name << "\" has " << nblock << " cells" 
                   << " of " << volume << " [m^3]";
        if (nvofs > 0) *vo_->os() << ", vol.fractions: " << vofs_min << "/" << vofs_max;
        *vo_->os() << std::endl;
      }

      if (nblock == 0) {
        msg << "Used source region is empty.";
        Exceptions::amanzi_throw(msg);
      } 
    }
  }

  if (alist.isParameter("used boundary condition regions")) {
    std::vector<std::string> regions = alist.get<Teuchos::Array<std::string> >("used boundary condition regions").toVector();
    regions.erase(SelectUniqueEntries(regions.begin(), regions.end()), regions.end());

    for (int i = 0; i < regions.size(); i++) {
      AmanziMesh::Entity_ID_List block;
      std::vector<double> vofs;
      mesh_->get_set_entities_and_vofs(regions[i], AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED, &block, &vofs);
      int nblock = block.size();
      int nvofs = vofs.size();

      double frac, area(0.0);
      for (int n = 0; n < nblock; n++) {
        frac = (nvofs == 0) ? 1.0 : vofs[n];
        area += mesh_->face_area(block[n]) * frac;
      }

      double vofs_min(1.0), vofs_max(0.0);
      for (int n = 0; n < nvofs; ++n) {
        vofs_min = std::min(vofs_min, vofs[n]);
        vofs_max = std::max(vofs_max, vofs[n]);
      }

      // verify that all faces are boundary faces
      int bc_flag(1);
      AmanziMesh::Entity_ID_List cells;

      for (int n = 0; n < nblock; ++n) {
        area += mesh_->face_area(block[n]) * frac;
        mesh_->face_get_cells(block[n], AmanziMesh::Parallel_type::ALL, &cells);
        if (cells.size() != 1) bc_flag = 0;
      }

#ifdef HAVE_MPI
      int nblock_tmp(nblock), nvofs_tmp(nvofs), bc_flag_tmp(bc_flag);
      double area_tmp(area), vofs_min_tmp(vofs_min), vofs_max_tmp(vofs_max);

      mesh_->get_comm()->SumAll(&nblock_tmp, &nblock, 1);
      mesh_->get_comm()->SumAll(&nvofs_tmp, &nvofs, 1);
      mesh_->get_comm()->SumAll(&area_tmp, &area, 1);
      mesh_->get_comm()->MinAll(&vofs_min_tmp, &vofs_min, 1);
      mesh_->get_comm()->MaxAll(&vofs_max_tmp, &vofs_max, 1);
      mesh_->get_comm()->MinAll(&bc_flag_tmp, &bc_flag, 1);
#endif
      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
        std::string name(regions[i]);
        name.resize(std::min(40, (int)name.size()));
        *vo_->os() << "bc: \"" << name << "\" has " << nblock << " faces"
                   << " of " << area << " [m^2]";
        if (nvofs > 0) *vo_->os() << ", vol.fractions: " << vofs_min << "/" << vofs_max;
        *vo_->os() << std::endl;
      }

      if (nblock == 0) {
        msg << "Used boundary region is empty.";
        Exceptions::amanzi_throw(msg);
      } 
      if (bc_flag == 0) {
        msg << "Used boundary region has non-boundary entries.";
        Exceptions::amanzi_throw(msg);
      }
    }
  }

  if (alist.isParameter("used observation regions")) {
    std::vector<std::string> regions = alist.get<Teuchos::Array<std::string> >("used observation regions").toVector();

    int nblock, nblock_tmp;
    for (int i = 0; i < regions.size(); i++) {
      double volume(0.0), volume_tmp;
      std::string type;
      AmanziMesh::Entity_ID_List block;
      std::vector<double> vofs;

      // observation region may use either cells of faces
      if (! mesh_->valid_set_name(regions[i], AmanziMesh::CELL) &&
          ! mesh_->valid_set_name(regions[i], AmanziMesh::FACE)) {
        std::string name(regions[i]);
        name.resize(std::min(40, (int)name.size()));
        *vo_->os() << "Observation region: \"" << name << "\" has unknown type." << std::endl;
      }

      mesh_->get_set_entities_and_vofs(regions[i], AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block, &vofs);
      nblock_tmp = nblock = block.size();
      type = "cells";
      for (int n = 0; n < nblock; n++) 
          volume += mesh_->cell_volume(block[n]);

#ifdef HAVE_MPI
      volume_tmp = volume;
      mesh_->get_comm()->SumAll(&nblock_tmp, &nblock, 1);
      mesh_->get_comm()->SumAll(&volume_tmp, &volume, 1);
#endif

      if (nblock == 0) {
        mesh_->get_set_entities_and_vofs(regions[i], AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED, &block, &vofs);
        nblock_tmp = nblock = block.size();
        type = "faces";
        for (int n = 0; n < nblock; n++) 
            volume += mesh_->face_area(block[n]);

#ifdef HAVE_MPI
        volume_tmp = volume;
        mesh_->get_comm()->SumAll(&nblock_tmp, &nblock, 1);
        mesh_->get_comm()->SumAll(&volume_tmp, &volume, 1);
#endif
      } 

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
  if (!plist_->isSublist("analysis")) return;
  if (vo_->getVerbLevel() < Teuchos::VERB_EXTREME) return;

  int bc_counter = 0;

  if (plist_->isSublist("flow")) {
    Teuchos::ParameterList& flow_list = plist_->sublist("flow");

    Teuchos::ParameterList bc_list;
    if (flow_list.isSublist("boundary conditions")) {
      bc_list = flow_list.sublist("boundary conditions");
    }

    Teuchos::ParameterList mass_flux_list, pressure_list, seepage_list, head_list;
    if (bc_list.isSublist("mass flux")) {
      mass_flux_list = bc_list.sublist("mass flux");
      for (auto it = mass_flux_list.begin(); it != mass_flux_list.end(); ++it) {

        if (mass_flux_list.isSublist(mass_flux_list.name(it))) {
          Teuchos::ParameterList& bc = mass_flux_list.sublist(mass_flux_list.name(it));

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
      for (auto it = pressure_list.begin(); it != pressure_list.end(); ++it) {

        if (pressure_list.isSublist(pressure_list.name(it))) {
          Teuchos::ParameterList& bc = pressure_list.sublist(pressure_list.name(it));
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
      for (auto it = seepage_list.begin(); it != seepage_list.end(); ++it) {

        if (seepage_list.isSublist(seepage_list.name(it))) {
          Teuchos::ParameterList& bc = seepage_list.sublist(seepage_list.name(it));
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
      for (auto it = head_list.begin(); it != head_list.end(); ++it) {

        if (head_list.isSublist(head_list.name(it))) {
          Teuchos::ParameterList& bc = head_list.sublist(head_list.name(it));
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


/* ******************************************************************
* Selects unique entries and places them in [first, last)
****************************************************************** */
template<class Iterator>
Iterator InputAnalysis::SelectUniqueEntries(Iterator first, Iterator last)
{
  while (first != last) {
    Iterator next(first);
    last = std::remove(++next, last, *first);
    first = next;
  }
  return last;
}

}  // namespace Amanzi
