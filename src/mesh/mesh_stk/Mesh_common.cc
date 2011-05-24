// -------------------------------------------------------------
/**
 * @file   Mesh_common.cc
 * @author William A. Perkins
 * @date Wed Dec 29 12:12:43 2010
 * 
 * @brief  Implementations of some stk::mesh utility routines
 * 
 * 
 */
// -------------------------------------------------------------
// Created December 29, 2010 by William A. Perkins
// Last Change: Wed Dec 29 12:12:43 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>

#include "dbc.hh"
#include "Mesh_common.hh"

namespace stk {
  namespace mesh {


    // -------------------------------------------------------------
    // check_face_ownership
    // -------------------------------------------------------------
    /** 
     * Collective
     * This makes sure that each face is owned by a process that owns
     * one of the connecting cells.  The convention is to pick the
     * process with the lowest rank.
     * 
     * @param bulk_data 
     */
    void 
    check_face_ownership(BulkData& bulk_data)
    {
      const int me(bulk_data.parallel_rank());
      const MetaData& meta_data = bulk_data.mesh_meta_data();
      EntityVector faces;
      Selector owned(meta_data.locally_owned_part());
      get_selected_entities(owned, bulk_data.buckets(Face), faces);

      std::vector<EntityProc> eproc;

      for (EntityVector::iterator f = faces.begin(); 
           f != faces.end(); f++) {

        EntityVector theface, nodes, cells;
        theface.push_back(*f);
        get_entities_through_relations(theface, Node, nodes);
        get_entities_through_relations(nodes, Element, cells);

        // there better be at least one cell related to this face

        ASSERT(!cells.empty());

        // There should be at most 2 cells related to this face

        ASSERT(cells.size() <= 2);

        unsigned int cowner = 
          std::min<unsigned int>(cells.front()->owner_rank(), 
                                 cells.back()->owner_rank());
        

        // the face should be owned by the same process as the cell,
        // if not, change the owner; only the face owner can request
        // this

        if (cowner != me) {
          eproc.push_back(EntityProc(*f, cowner));
          // std::string msg = 
          //     boost::str(boost::format("%d: face %d: changing owner from %d to %d") %
          //                me % (*f)->identifier() % 
          //                (*f)->owner_rank() % cowner);
          // std::cerr << msg << std::endl;
        } else {
          // std::string msg = 
          //     boost::str(boost::format("%d: face %d, owner %d: cell %d, owner %d") %
          //                me % (*f)->identifier() % (*f)->owner_rank() % 
          //                cells.front()->identifier() % cells.front()->owner_rank());
          // std::cerr << msg << std::endl;
        }            
      }
      bulk_data.modification_begin();
      bulk_data.change_entity_owner(eproc);
      bulk_data.modification_end();
    }

    // -------------------------------------------------------------
    // check_node_ownership
    // -------------------------------------------------------------
    /** 
     * Collective
     * This makes sure that each node is owned by a process that owns
     * one of the connecting cells.  The convention is to pick the
     * process with the lowest rank.
     * 
     * @param bulk_data 
     */
    void 
    check_node_ownership(BulkData& bulk_data)
    {
      const int me(bulk_data.parallel_rank());
      const MetaData& meta_data = bulk_data.mesh_meta_data();
      EntityVector nodes;
      Selector owned(meta_data.locally_owned_part());
      get_selected_entities(owned, bulk_data.buckets(Node), nodes);

      std::vector<EntityProc> eproc;

      for (EntityVector::iterator n = nodes.begin(); n != nodes.end(); n++) {

        EntityVector thenode, cells;
        thenode.push_back(*n);
        get_entities_through_relations(nodes, Element, cells);

        unsigned int new_owner(10000000);

        if (!cells.empty()) {

          // find the lowest owner rank amongst the connected cells
          // the node should be owned by the same process as the cell,
          // if not, change the owner; only the node owner can request
          // this

          for (EntityVector::iterator c = cells.begin(); c != cells.end(); c++) {
            if ((*c)->owner_rank() < new_owner) new_owner = (*c)->owner_rank();
          }

        } else {

          // If there are a no cells related to this node, then we
          // definitely need to get it off this processor.  Look
          // through the list of processes that share it and choose
          // the lowest ranked one.

          const PairIterEntityComm sharing = (*n)->sharing();
          
          
          if (sharing.empty()) continue;

          for ( size_t j = 0 ; j < sharing.size() ; ++j ) {
            if (sharing[j].proc != me && sharing[j].proc < new_owner) {
              new_owner = sharing[j].proc;
            }
          }          

        }

        if (new_owner != me) {
          eproc.push_back(EntityProc(*n, new_owner));
        }            
      }
      bulk_data.modification_begin();
      bulk_data.change_entity_owner(eproc);
      bulk_data.modification_end();
    }


  } // close namespace mesh
} // close namespace stk
