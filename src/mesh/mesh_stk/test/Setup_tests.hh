#ifndef _SETUP_TESTS_H_
#define _SETUP_TESTS_H_

#include <stk_util/parallel/Parallel.hpp>
#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>

#include "../Mesh_STK_Impl.hh"
#include "../Mesh_STK_factory.hh"

#include "Element_block.hh"
#include "Coordinates.hh"
#include "Element_types.hh"
#include "Field_data.hh"

#include "Example_Mesh.hh"

extern stk::ParallelMachine parallel_machine;

typedef Teuchos::RCP<Amanzi::AmanziMesh::STK::Mesh_STK_Impl> Mesh_p;
typedef Teuchos::RCP<Amanzi::AmanziMesh::STK::Mesh_STK_factory> Factory_p;
typedef Teuchos::RCP<Epetra_MpiComm> Comm_p;


static double real_coordinates [3];
static void set_real_coordinates (double x, double y, double z)
{
    real_coordinates [0] = x;
    real_coordinates [1] = y;
    real_coordinates [2] = z;
}


struct Mesh_setup : public Test_mesh
{

    Factory_p factory;
    Mesh_p mesh;
    Comm_p communicator;
    int my_pid;
    

    Mesh_setup ()
    {
        Amanzi::AmanziMesh::Data::Fields fields;

        communicator = Comm_p (new Epetra_MpiComm(parallel_machine));
        my_pid = communicator->MyPID ();

        const int bucket_size = 20;
        
        factory = Factory_p(new Amanzi::AmanziMesh::STK::Mesh_STK_factory(parallel_machine, bucket_size));
        mesh    = Mesh_p (factory->build_mesh (*data, fields, NULL));

    }

    ~Mesh_setup ()
    {
        communicator->Barrier ();
    }

};


#endif /* _SETUP_TESTS_H_ */
