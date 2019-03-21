/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
*/

#include <algorithm>
#include <boost/format.hpp>
#include "exodusII.h"

#include "AmanziComm.hh"
#include "Parallel_Exodus_file.hh"
#include "Exodus_readers.hh"


namespace Amanzi {
namespace Exodus {

// -------------------------------------------------------------
//  class Parallel_Exodus_file
// -------------------------------------------------------------

// -------------------------------------------------------------
// Parallel_Exodus_file:: constructors / destructor
// -------------------------------------------------------------
/** 
 * 
 * The files are name as @c nem_spread names them, namely @e
 * basename.N.n, where @e N is the total number of processors and @e n
 * is the local processor rank.
 * 
 * @param comm_ parallel environment
 * @param basename ExodusII file set base name
 */
Parallel_Exodus_file::Parallel_Exodus_file(const Comm_ptr_type& comm,
                                           const std::string& basename)
  : comm_(comm),
    basename_(basename)
{
  const int np(comm_->getSize());
  const int me(comm_->getRank());

  std::string s(basename_);

  // create the format string that pads me with the correct number of zeros
  int ndigits = (int)floor(log10(np)) + 1;
  std::stringstream frmt;
  frmt << ".%d.%0" << ndigits << "d";

  s += boost::str(boost::format( frmt.str().c_str() ) % np % me);
  //s += boost::str(boost::format(".%d.%d") % np % me);

  std::cerr << "Process " << me << " of " << np << ": trying file " << s << std::endl;

  int ierr(0);
  ExodusError ex_all;

  try {
    file_.reset(new Exodus_file(s.c_str()));
  } catch (const ExodusError& e) {
    std::string msg;
    msg = boost::str(boost::format("Process %d: error opening \"%s\": %s") %
                     me % s % e.what());
    std::cerr << msg << std::endl;
    ierr += 1;
    ex_all.add_data(msg.c_str());
  }

  int gerr(0);
  Teuchos::reduceAll(comm_, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  if (gerr > 0) {
    Exceptions::amanzi_throw( ex_all );
  }
}

// -------------------------------------------------------------
// Parallel_Exodus_file::read_mesh
// -------------------------------------------------------------
Teuchos::RCP<AmanziMesh::Data::Data> 
Parallel_Exodus_file::read_mesh(void)
{

  mesh_.reset(read_exodus_file(*file_));

  // Even though element blocks are defined in all local files, if the
  // local block is empty, it will not have a cell type specified, so
  // we need get that from the other processes

  comm_->Barrier();
  const int np(comm_->getSize());
  const int me(comm_->getRank());

  std::vector<int> byproc(np);
  int nblk(mesh_->element_blocks());

  // check the number of blocks; should be the same on all processes

  int ierr(0);
  comm_->GatherAll(&nblk, &byproc[0], 1);
  for (int p = 1; p < np; p++) {
    if (byproc[p] != byproc[p-1]) ierr++;
  }
  int aerr(0);
  Teuchos::reduceAll(comm_, Teuchos::REDUCE_SUM, 1, &ierr, &aerr);

  if (aerr) {
    std::string msg(basename_);
    msg += ": mismatched numbers of element blocks";
    Exceptions::amanzi_throw( ExodusError(msg.c_str()) );
  }

  for (int b = 0; b < nblk; b++) {
    int mytype(mesh_->element_block(b).element_type());
    std::vector<int> alltype(np, AmanziMesh::CELLTYPE_UNKNOWN);
    comm_->GatherAll(&mytype, &alltype[0], 1);

    std::vector<int>::iterator junk;
    junk = std::remove(alltype.begin(), alltype.end(),
                       AmanziMesh::CELLTYPE_UNKNOWN);
    alltype.erase(junk, alltype.end());
    junk = std::unique(alltype.begin(), alltype.end());
    alltype.erase(junk, alltype.end());
    
    if (alltype.empty()) {

      // this means that all processes reported the AmanziMesh::Data::UNKNOWN
      // type for this element block; this is OK as long as it's empty
      // on all processes.

      if (mesh_->element_block(b).num_elements() > 0) {
        std::string msg = 
          boost::str(boost::format("Process %d: %s: element block %d: block element type unknown") %
                     me % basename_ %  mesh_->element_block(b).id());
        std::cerr << msg << std::endl;
        ierr++;
      }
    } else if (alltype.size() > 1) {
      
      // this means that at least two processes reported different,
      // not unknown, element types for this block; this is an error
      // and needs to be reported.

      std::string msg = 
        boost::str(boost::format("Process %d: %s: element block %d: block element type mismatch") %
                   me % basename_ %  mesh_->element_block(b).id());
      std::cerr << msg << std::endl;
      ierr++;

    } else {

      AmanziMesh::Cell_type
        thetype(static_cast<AmanziMesh::Cell_type>(alltype.front()));

      if (mesh_->element_block(b).element_type() == AmanziMesh::CELLTYPE_UNKNOWN) {
        mesh_->element_block(b).element_type(thetype);
      }
    }
  }

  Teuchos::reduceAll(comm_, Teuchos::REDUCE_SUM, 1, &ierr, &aerr);
  if (aerr) {
    std::string msg = 
      boost::str(boost::format("%s: element block type errors") % basename_);
    Exceptions::amanzi_throw( ExodusError(msg.c_str()) );
  }

  return mesh_;
}

// -------------------------------------------------------------
// Parallel_Exodus_file::cellmap
// -------------------------------------------------------------
Map_ptr_type
Parallel_Exodus_file::cellmap(void)
{
  if (mesh_.is_null()) {
    read_mesh();
  }
  
  int myelem(mesh_->parameters().num_elements_);

  std::vector<int> gids(myelem);

  int ret_val = ex_get_id_map(file_->id, EX_ELEM_MAP, &gids[0]);
  if (ret_val < 0) {
    std::string msg = 
      boost::str(boost::format("%s: error: cannot read element number map (%d)") %
                 file_->filename % ret_val);
    Exceptions::amanzi_throw( ExodusError (msg.c_str()) );
  }

  comm_->Barrier();

  return Teuchos::rcp(new Map_type(-1, myelem, &gids[0], 1, *comm_));
}


// -------------------------------------------------------------
// Parallel_Exodus_file::vertexmap
// -------------------------------------------------------------
Map_ptr_type
Parallel_Exodus_file::vertexmap(void)
{
  if (mesh_.is_null()) {
    read_mesh();
  }
  
  int myvert(mesh_->parameters().num_nodes_);

  std::vector<int> gids(myvert);

  int ret_val = ex_get_id_map(file_->id, EX_NODE_MAP, &gids[0]);
  if (ret_val < 0) {
    std::string msg = 
      boost::str(boost::format("%s: error: cannot read vertex number map (%d)") %
                 file_->filename % ret_val);
    Exceptions::amanzi_throw( ExodusError (msg.c_str()) );
  }

  comm_->Barrier();

  return Teuchos::rcp(new Map_type(-1, myvert, &gids[0], 1, *comm_));
}

} // namespace Exodus
} // namespace Amanzi


