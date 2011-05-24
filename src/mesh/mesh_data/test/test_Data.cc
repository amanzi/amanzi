#include <UnitTest++.h>

#include "../Data.hh"

SUITE (Data)
{

    TEST (Construct)
    {
        Amanzi::AmanziMesh::Data::Data mesh ();
    }

}
