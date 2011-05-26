#include <UnitTest++.h>

#include "../Field_data.hh"

#include <string>

SUITE (Field_data)
{

    TEST (Construct)
    {

        std::string name = "Rumplestiltskin";
        Amanzi::AmanziMesh::Entity_kind location = Amanzi::AmanziMesh::CELL;
        Amanzi::AmanziMesh::Data::FIELD_TYPE type = Amanzi::AmanziMesh::Data::SCALAR;

        Amanzi::AmanziMesh::Data::Field f (name, type, location);

        CHECK_EQUAL (f.name (), name);
        CHECK_EQUAL (f.type (), type);
        CHECK_EQUAL (f.location (), location);

    }

    TEST (Default)
    {
        std::string name = "Bob's crab shack";
        Amanzi::AmanziMesh::Data::Field f (name);

        CHECK_EQUAL (f.name (), name);
        CHECK_EQUAL (f.location (), Amanzi::AmanziMesh::CELL);
        CHECK_EQUAL (f.type (), Amanzi::AmanziMesh::Data::SCALAR);
    }

}

SUITE (Fields)
{

    TEST (Empty)
    {
        Amanzi::AmanziMesh::Data::Fields fields;
        CHECK (fields.empty ());
    }

    TEST (Add)
    {
        Amanzi::AmanziMesh::Data::Fields fields;

        Amanzi::AmanziMesh::Data::Field f ("Smith", Amanzi::AmanziMesh::Data::SCALAR, 
                                     Amanzi::AmanziMesh::NODE);
        fields.add_field (f);

        Amanzi::AmanziMesh::Data::Field g ("Jones", Amanzi::AmanziMesh::Data::VECTOR, 
                                     Amanzi::AmanziMesh::EDGE);
        fields.add_field (g);

        Amanzi::AmanziMesh::Data::Field h ("Brown", Amanzi::AmanziMesh::Data::SCALAR, 
                                     Amanzi::AmanziMesh::FACE);
        fields.add_field (h);

        CHECK (!fields.empty ());
        int count = 0;
        for (Amanzi::AmanziMesh::Data::Fields::const_iterator it = fields.begin ();
             it != fields.end ();
             ++it)
        {
            CHECK_EQUAL (it->location (), count);
            ++count;
        }
        CHECK_EQUAL (3, count);

     }
}
