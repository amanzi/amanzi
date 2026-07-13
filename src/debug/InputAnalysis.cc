/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

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
void
InputAnalysis::Init(Teuchos::ParameterList& plist)
{
  plist_ = &plist;
  vo_ = new VerboseObject("InputAnalysis:" + domain_, plist);
}


/* ******************************************************************
* Analysis of collected regions
****************************************************************** */
void
InputAnalysis::RegionAnalysis()
{
  if (vo_->getVerbLevel() < Teuchos::VERB_MEDIUM) return;

  Teuchos::OSTab tab = vo_->getOSTab();

  // Helper: print cell-region stats (count + volume)
  auto analyzeCellRegion = [&](const std::string& region) {
    int nblock(0), nblock_tmp, nvofs(0);
    double volume(0.0), frac;
    AmanziMesh::cEntity_ID_View block;
    AmanziMesh::cDouble_View vofs;

    try {
      Kokkos::tie(block, vofs) = mesh_->getSetEntitiesAndVolumeFractions(
        region, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
      nblock = block.size();
      nvofs = vofs.size();
      for (int n = 0; n < nblock; n++) {
        frac = (nvofs == 0) ? 1.0 : vofs[n];
        volume += mesh_->getCellVolume(block[n]) * frac;
      }
    } catch (...) {
      nblock = -1;
    }

    double vofs_min(1.0), vofs_max(0.0);
    for (int n = 0; n < nvofs; ++n) {
      vofs_min = std::min(vofs_min, vofs[n]);
      vofs_max = std::max(vofs_max, vofs[n]);
    }
    if (nvofs == 0) vofs_max = 1.0;

    int nvofs_tmp(nvofs);
    double volume_tmp(volume), vofs_min_tmp(vofs_min), vofs_max_tmp(vofs_max);
    mesh_->getComm()->MinAll(&nblock, &nblock_tmp, 1);
    if (nblock_tmp < 0) { nblock = 0; nvofs = 0; volume = 0.0; }
    else { mesh_->getComm()->SumAll(&nblock, &nblock_tmp, 1);  nblock = nblock_tmp; }
    mesh_->getComm()->SumAll(&nvofs_tmp, &nvofs, 1);
    mesh_->getComm()->SumAll(&volume_tmp, &volume, 1);
    mesh_->getComm()->MinAll(&vofs_min_tmp, &vofs_min, 1);
    mesh_->getComm()->MaxAll(&vofs_max_tmp, &vofs_max, 1);

    std::string name(region);
    name.resize(std::min(40, (int)name.size()));
    *vo_->os() << "cell: \"" << name << "\" has " << nblock << " cells, vol: " << volume << " [m^3]";
    if (nvofs > 0) *vo_->os() << ", vol.fractions: " << vofs_min << "/" << vofs_max;
    *vo_->os() << std::endl;
  };

  // Helper: print face-region stats (count + area)
  auto analyzeFaceRegion = [&](const std::string& region) {
    AmanziMesh::cEntity_ID_View block;
    AmanziMesh::cDouble_View vofs;
    Kokkos::tie(block, vofs) = mesh_->getSetEntitiesAndVolumeFractions(
      region, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    int nblock = block.size();
    int nvofs = vofs.size();

    double frac, area(0.0);
    for (int n = 0; n < nblock; n++) {
      frac = (nvofs == 0) ? 1.0 : vofs[n];
      area += mesh_->getFaceArea(block[n]) * frac;
    }

    double vofs_min(1.0), vofs_max(0.0);
    for (int n = 0; n < nvofs; ++n) {
      vofs_min = std::min(vofs_min, vofs[n]);
      vofs_max = std::max(vofs_max, vofs[n]);
    }

    int nblock_tmp(nblock), nvofs_tmp(nvofs);
    double area_tmp(area), vofs_min_tmp(vofs_min), vofs_max_tmp(vofs_max);
    mesh_->getComm()->SumAll(&nblock_tmp, &nblock, 1);
    mesh_->getComm()->SumAll(&nvofs_tmp, &nvofs, 1);
    mesh_->getComm()->SumAll(&area_tmp, &area, 1);
    mesh_->getComm()->MinAll(&vofs_min_tmp, &vofs_min, 1);
    mesh_->getComm()->MaxAll(&vofs_max_tmp, &vofs_max, 1);

    std::string name(region);
    name.resize(std::min(40, (int)name.size()));
    *vo_->os() << "face: \"" << name << "\" has " << nblock << " faces, area: " << area << " [m^2]";
    if (nvofs > 0) *vo_->os() << ", vol.fractions: " << vofs_min << "/" << vofs_max;
    *vo_->os() << std::endl;
  };

  // Helper: print node-region stats (count only)
  auto analyzeNodeRegion = [&](const std::string& region) {
    AmanziMesh::cEntity_ID_View block;
    AmanziMesh::cDouble_View vofs;
    Kokkos::tie(block, vofs) = mesh_->getSetEntitiesAndVolumeFractions(
      region, AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
    int nblock = block.size(), nblock_tmp;
    mesh_->getComm()->SumAll(&nblock, &nblock_tmp, 1);  nblock = nblock_tmp;

    std::string name(region);
    name.resize(std::min(40, (int)name.size()));
    *vo_->os() << "node: \"" << name << "\" has " << nblock << " nodes" << std::endl;
  };

  // cell regions
  if (plist_->isParameter("cell regions")) {
    std::vector<std::string> regions =
      plist_->get<Teuchos::Array<std::string>>("cell regions").toVector();
    regions.erase(SelectUniqueEntries(regions.begin(), regions.end()), regions.end());
    if (regions.size() == 1 && regions[0] == "{*}")
      regions = mesh_->getResolvedSetNames(AmanziMesh::Entity_kind::CELL);
    for (const auto& r : regions) analyzeCellRegion(r);
  }

  // face regions
  if (plist_->isParameter("face regions")) {
    std::vector<std::string> regions =
      plist_->get<Teuchos::Array<std::string>>("face regions").toVector();
    regions.erase(SelectUniqueEntries(regions.begin(), regions.end()), regions.end());
    if (regions.size() == 1 && regions[0] == "{*}")
      regions = mesh_->getResolvedSetNames(AmanziMesh::Entity_kind::FACE);
    for (const auto& r : regions) analyzeFaceRegion(r);
  }

  // node regions
  if (plist_->isParameter("node regions")) {
    std::vector<std::string> regions =
      plist_->get<Teuchos::Array<std::string>>("node regions").toVector();
    regions.erase(SelectUniqueEntries(regions.begin(), regions.end()), regions.end());
    if (regions.size() == 1 && regions[0] == "{*}")
      regions = mesh_->getResolvedSetNames(AmanziMesh::Entity_kind::NODE);
    for (const auto& r : regions) analyzeNodeRegion(r);
  }

  // "used source regions" -> cell analysis
  if (plist_->isParameter("used source regions")) {
    std::vector<std::string> regions =
      plist_->get<Teuchos::Array<std::string>>("used source regions").toVector();
    regions.erase(SelectUniqueEntries(regions.begin(), regions.end()), regions.end());
    for (const auto& r : regions) analyzeCellRegion(r);
  }

  // "used boundary condition regions" -> face analysis
  if (plist_->isParameter("used boundary condition regions")) {
    std::vector<std::string> regions =
      plist_->get<Teuchos::Array<std::string>>("used boundary condition regions").toVector();
    regions.erase(SelectUniqueEntries(regions.begin(), regions.end()), regions.end());
    for (const auto& r : regions) analyzeFaceRegion(r);
  }

  // "used observation regions" -> cell-then-face analysis
  if (plist_->isParameter("used observation regions")) {
    std::vector<std::string> regions =
      plist_->get<Teuchos::Array<std::string>>("used observation regions").toVector();
    for (const auto& r : regions) {
      if (mesh_->isValidSetName(r, AmanziMesh::Entity_kind::CELL))
        analyzeCellRegion(r);
      else if (mesh_->isValidSetName(r, AmanziMesh::Entity_kind::FACE))
        analyzeFaceRegion(r);
      else
        *vo_->os() << "obs: \"" << r << "\" has unknown type." << std::endl;
    }
  }
}


/* ******************************************************************
* DEBUG output: boundary conditions.
****************************************************************** */
void
InputAnalysis::OutputBCs()
{
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

            Teuchos::ParameterList& f_tab =
              (bc.sublist("outward mass flux")).sublist("function-tabular");

            Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double>>("x values");
            Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double>>("y values");
            Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string>>("forms");

            int np = times.size() * 2 - 1;
            Teuchos::Array<double> times_plot(np);
            Teuchos::Array<double> values_plot(np);

            for (int i = 0; i < times.size() - 1; i++) {
              times_plot[2 * i] = times[i];
              values_plot[2 * i] = values[i];
              times_plot[2 * i + 1] = 0.5 * (times[i] + times[i + 1]);
            }
            times_plot[np - 1] = times[times.size() - 1];
            values_plot[np - 1] = values[times.size() - 1];

            for (int i = 0; i < time_fns.size(); i++) {
              if (time_fns[i] == "linear") {
                values_plot[2 * i + 1] = 0.5 * (values[i] + values[i + 1]);
              } else if (time_fns[i] == "constant") {
                values_plot[2 * i + 1] = values[i];
                times_plot[2 * i + 1] = times[i + 1];
              } else {
                Exceptions::amanzi_throw(Errors::Message(
                  "In the definition of BCs: tabular function can only be Linear or Constant"));
              }
            }

            std::string filename = ss.str() + ".dat";
            std::ofstream ofile(filename.c_str());

            ofile << "# "
                  << "time "
                  << "flux" << std::endl;
            for (int i = 0; i < np; i++) {
              ofile << times_plot[i] << " " << values_plot[i] << std::endl;
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


            Teuchos::ParameterList& f_tab =
              (bc.sublist("boundary pressure")).sublist("function-tabular");

            Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double>>("x values");
            Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double>>("y values");
            Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string>>("forms");

            int np = times.size() * 2 - 1;
            Teuchos::Array<double> times_plot(np);
            Teuchos::Array<double> values_plot(np);

            for (int i = 0; i < times.size() - 1; i++) {
              times_plot[2 * i] = times[i];
              values_plot[2 * i] = values[i];
              times_plot[2 * i + 1] = 0.5 * (times[i] + times[i + 1]);
            }
            times_plot[np - 1] = times[times.size() - 1];
            values_plot[np - 1] = values[times.size() - 1];

            for (int i = 0; i < time_fns.size(); i++) {
              if (time_fns[i] == "linear") {
                values_plot[2 * i + 1] = 0.5 * (values[i] + values[i + 1]);
              } else if (time_fns[i] == "constant") {
                values_plot[2 * i + 1] = values[i];
                times_plot[2 * i + 1] = times[i + 1];
              } else {
                Exceptions::amanzi_throw(Errors::Message(
                  "In the definition of BCs: tabular function can only be Linear or Constant"));
              }
            }

            std::string filename = ss.str() + ".dat";
            std::ofstream ofile(filename.c_str());

            ofile << "# time "
                  << "pressure" << std::endl;
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

            Teuchos::ParameterList& f_tab =
              (bc.sublist("outward mass flux")).sublist("function-tabular");

            Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double>>("x values");
            Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double>>("y values");
            Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string>>("forms");

            int np = times.size() * 2 - 1;
            Teuchos::Array<double> times_plot(np);
            Teuchos::Array<double> values_plot(np);

            for (int i = 0; i < times.size() - 1; i++) {
              times_plot[2 * i] = times[i];
              values_plot[2 * i] = values[i];
              times_plot[2 * i + 1] = 0.5 * (times[i] + times[i + 1]);
            }
            times_plot[np - 1] = times[times.size() - 1];
            values_plot[np - 1] = values[times.size() - 1];

            for (int i = 0; i < time_fns.size(); i++) {
              if (time_fns[i] == "linear") {
                values_plot[2 * i + 1] = 0.5 * (values[i] + values[i + 1]);
              } else if (time_fns[i] == "constant") {
                values_plot[2 * i + 1] = values[i];
                times_plot[2 * i + 1] = times[i + 1];
              } else {
                Exceptions::amanzi_throw(Errors::Message(
                  "In the definition of BCs: tabular function can only be Linear or Constant"));
              }
            }

            std::string filename = ss.str() + ".dat";
            std::ofstream ofile(filename.c_str());

            ofile << "# time "
                  << "flux" << std::endl;
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

            Teuchos::ParameterList& f_tab =
              (bc.sublist("water table elevation")).sublist("function-tabular");

            Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double>>("x values");
            Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double>>("y values");
            Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string>>("forms");

            int np = times.size() * 2 - 1;
            Teuchos::Array<double> times_plot(np);
            Teuchos::Array<double> values_plot(np);

            for (int i = 0; i < times.size() - 1; i++) {
              times_plot[2 * i] = times[i];
              values_plot[2 * i] = values[i];
              times_plot[2 * i + 1] = 0.5 * (times[i] + times[i + 1]);
            }
            times_plot[np - 1] = times[times.size() - 1];
            values_plot[np - 1] = values[times.size() - 1];

            for (int i = 0; i < time_fns.size(); i++) {
              if (time_fns[i] == "linear") {
                values_plot[2 * i + 1] = 0.5 * (values[i] + values[i + 1]);
              } else if (time_fns[i] == "constant") {
                values_plot[2 * i + 1] = values[i];
                times_plot[2 * i + 1] = times[i + 1];
              } else {
                Exceptions::amanzi_throw(Errors::Message(
                  "In the definition of BCs: tabular function can only be Linear or Constant"));
              }
            }

            std::string filename = ss.str() + ".dat";
            std::ofstream ofile(filename.c_str());

            ofile << "# time "
                  << "head" << std::endl;
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
Iterator
InputAnalysis::SelectUniqueEntries(Iterator first, Iterator last)
{
  while (first != last) {
    Iterator next(first);
    last = std::remove(++next, last, *first);
    first = next;
  }
  return last;
}

} // namespace Amanzi
