/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
*/


#ifndef AMANZI_MESH_DATA_HH_
#define AMANZI_MESH_DATA_HH_

#include <iostream>
#include <memory>
#include <vector>

#include "Element_block.hh"
#include "Side_set.hh"
#include "Node_set.hh"
#include "Coordinates.hh"
#include "Parameters.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace Data {


//
// This is bad memory management and should be refactored. --etc
//
template <typename P>
void reclaim(std::vector<P*>& v)
{
  for (typename std::vector<P*>::iterator it = v.begin(); it != v.end(); ++it)
  delete *it;
}


class Data {
 private:

  // Mixing raw pointers and unique_ptr?  Really?  Why?
  // fixme --etc
  std::unique_ptr<Parameters> params_;
  std::unique_ptr<Coordinates<double> > coords_;
  std::vector<Element_block*> element_blocks_;
  std::vector<Side_set*> side_sets_;
  std::vector<Node_set*> node_sets_;

 public:
  int element_blocks() const { return element_blocks_.size(); }
  int side_sets() const      { return side_sets_.size(); }
  int node_sets() const      { return node_sets_.size(); }

  Element_block& element_block(int id);
  const Element_block& element_block(int id) const;
  const Side_set& side_set(int id) const;
  const Node_set& node_set(int id) const;

  const Parameters& parameters() const { return *params_; }

  const Coordinates<double>& coordinates() const { return *coords_; }

  void take_data_from(Parameters* params,
                      Coordinates<double>* coords,
                      std::vector<Element_block*> blocks,
                      std::vector<Side_set*> side_sets,
                      std::vector<Node_set*> node_sets);

  static Data* build_from(Parameters* params,
                          Coordinates<double>* coords, 
                          std::vector<Element_block*> blocks,
                          std::vector<Side_set*> side_sets,
                          std::vector<Node_set*> node_sets);
    
  void to_stream (std::ostream& stream, const bool& verbose = false) const;

  virtual ~Data();
};


inline std::ostream& operator<<(std::ostream& stream, const Data& m)
{
  m.to_stream (stream);
  return stream;
}

}  // namespace Data
}  // namespace AmanziMesh
}  // namespace Amanzi

#endif


