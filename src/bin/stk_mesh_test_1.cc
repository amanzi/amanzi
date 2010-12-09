#include <iostream>
#include "mpi.h"

// Trilinos header files
#include "Teuchos_GlobalMPISession.hpp"

// Amanzi header files
#include "Exodus_readers.hh"
#include "Mesh_factory.hh"
#include "Mesh_maps_stk.hh"
//#include "Data.hh"


Teuchos::RCP<STK_mesh::Mesh_maps_stk> import_exodus_mesh (const char *file)
{
    const int bucket_size = 20; // no clue what this means
    STK_mesh::Mesh_factory factory(MPI_COMM_WORLD, bucket_size);
    // Need the comm as an argument

    // Read the mesh data from the Exodus mesh file.
    Mesh_data::Data *data = ExodusII::read_exodus_file(file);
    //std::cout << *data << std::endl;

    // Create the low-level STK mesh.
    Mesh_data::Fields fields; // required but not currently used
    Teuchos::RCP<STK_mesh::Mesh> stkmesh(factory.build_mesh(*data, fields));

    // Create the associated mesh maps, which is what Amanzi regards as the mesh.
    Teuchos::RCP<STK_mesh::Mesh_maps_stk> mesh(new STK_mesh::Mesh_maps_stk(stkmesh));
    return mesh;
}


void dump_mesh_topology (STK_mesh::Mesh_maps_stk &mesh)
{
    // Dump the numbers of the various mesh features.
    int nnode_own = mesh.count_entities(Mesh_data::NODE, OWNED);
    int nnode_use = mesh.count_entities(Mesh_data::NODE, USED);
    int nface_own = mesh.count_entities(Mesh_data::FACE, OWNED);
    int nface_use = mesh.count_entities(Mesh_data::FACE, USED);
    int ncell_own = mesh.count_entities(Mesh_data::CELL, OWNED);
    int ncell_use = mesh.count_entities(Mesh_data::CELL, USED);

    std::cout << "Number of nodes (owned/used) = "
              << nnode_own << "/" << nnode_use << std::endl;
    std::cout << "Number of faces (owned/used) = "
              << nface_own << "/" << nface_use << std::endl;
    std::cout << "Number of cells (owned/used) = "
              << ncell_own << "/" << ncell_use << std::endl;

    // Dump the connectivity between mesh features.
    unsigned int cnode[8], cface[6], fnode[4];

    for (int j = 0; j < ncell_use; ++j) 
    {

        mesh.cell_to_nodes(j, cnode, cnode+8);
        std::cout << std::endl << "cell " << j << " has nodes:";
        for (int i = 0; i < 8; ++i) std::cout << " " << cnode[i];
        std::cout << std::endl;

        mesh.cell_to_faces(j, cface, cface+6);

        for (int k = 0; k < 6; ++k) 
        {
            std::cout << "  side " << k << " is face " << cface[k] << " and has nodes:";
            mesh.face_to_nodes(cface[k], fnode, fnode+4);
            for (int i = 0; i < 4; ++i) std::cout << " " << fnode[i];
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}


void dump_mesh_geometry (STK_mesh::Mesh_maps_stk &mesh)
{
    std::cout << "Node coordinates:" << std::endl;int nnode_use = mesh.count_entities(Mesh_data::NODE, USED);
    double x[3];
    for (int j = 0; j < nnode_use; ++j) 
    {
        std::cout << "node " << j << ":";
        mesh.node_to_coordinates(j, x, x+3);
        for (int i = 0; i < 3; ++i)
            std::cout << " " << x[i];
        std::cout << std::endl;
    }
}


int main (int argc, char* argv[])
{
    // Initialize MPI.
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    // Create a mesh that corresponds to the Exodus II mesh file.
    Teuchos::RCP<STK_mesh::Mesh_maps_stk> mesh = import_exodus_mesh (argv[1]);

    // Dump what we've created to stdout.
    dump_mesh_topology (*mesh);
    dump_mesh_geometry (*mesh);

    return 0;
}
