#ifndef STK1DDISC_HPP
#define STK1DDISC_HPP

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"

#include "STKMesh1D.hpp"

class STK1DDisc
{
  
public:
  
  STK1DDisc (Teuchos::RCP<STKMesh1D>& stkMesh_,
		  const Teuchos::RCP<const Epetra_Comm>& comm);
  
  ~STK1DDisc ();

  int getnodeDOF (stk::mesh::Entity& node) const;
    
  Teuchos::RCP<const Epetra_Map> getMap() const
  { return map; };

  Teuchos::RCP<const Epetra_Map> getElementMap() const
  { return element_map; };

  Teuchos::RCP<const Epetra_Map> getOverlapMap() const
  { return overlap_map; };

  Teuchos::ArrayRCP<double>&  getCoordinates() const;

  Teuchos::RCP<const Epetra_CrsGraph>  getJacobianGraph() const
  { return graph; };

  Teuchos::RCP<const Epetra_CrsGraph>  getOverlapJacobianGraph() const
  { return overlap_graph; };
  
  const int getnumMyElements () const
  { return numMyElements; };

  const Teuchos::RCP<STKMesh1D>& getMesh() const
  { return stkMesh; };

  const Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> >& getElNodeID() const
  { return elNodeID; };

protected:
  //! Epetra communicator
  Teuchos::RCP<const Epetra_Comm> comm;

  //! Element map
  Teuchos::RCP<Epetra_Map> elem_map;
  
  //! Node map
  Teuchos::RCP<Epetra_Map> node_map;

  //! Element map
  Teuchos::RCP<Epetra_Map> element_map;

  //! Unknown Map
  Teuchos::RCP<Epetra_Map> map;
  
  //! Overlapped unknown map
  Teuchos::RCP<Epetra_Map> overlap_map;
  
  //! Jacobian matrix graph
  Teuchos::RCP<Epetra_CrsGraph> graph;
  
  //! Overlapped Jacobian matrix graph
  Teuchos::RCP<Epetra_CrsGraph> overlap_graph;
  
  //! Processor ID
  unsigned int myPID;

  //! Number of elements on this processor
  unsigned int numMyElements;

  //! Number of nodes per element
  unsigned int nodes_per_element;

  Teuchos::ArrayRCP< Teuchos::ArrayRCP<int> > elNodeID;
  
  mutable Teuchos::ArrayRCP<double> coordinates;

  Teuchos::RCP<STKMesh1D> stkMesh;

  //! list of all overlap nodes, saved for getting coordinates for mesh motion
  std::vector< stk::mesh::Entity * > overlapnodes ;

};


#endif
