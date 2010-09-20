#include <UnitTest++.h>

#include "../Field_data.hh"
#include "../Entity_kind.hh"

#include <string>

SUITE (Field_data)
{

    TEST (Construct)
    {

        std::string name = "Rumplestiltskin";
        Mesh_data::Entity_kind location = Mesh_data::CELL;
        Mesh_data::FIELD_TYPE type = Mesh_data::SCALAR;

        Mesh_data::Field f (name, type, location);

        CHECK_EQUAL (f.name (), name);
        CHECK_EQUAL (f.type (), type);
        CHECK_EQUAL (f.location (), location);

    }

    TEST (Default)
    {
        std::string name = "Bob's crab shack";
        Mesh_data::Field f (name);

        CHECK_EQUAL (f.name (), name);
        CHECK_EQUAL (f.location (), Mesh_data::CELL);
        CHECK_EQUAL (f.type (), Mesh_data::SCALAR);
    }

}

SUITE (Fields)
{

    TEST (Empty)
    {
        Mesh_data::Fields fields;
        CHECK (fields.empty ());
    }

    TEST (Add)
    {
        Mesh_data::Fields fields;

        Mesh_data::Field f ("Smith", Mesh_data::SCALAR, Mesh_data::NODE);
        fields.add_field (f);

        Mesh_data::Field g ("Jones", Mesh_data::VECTOR, Mesh_data::EDGE);
        fields.add_field (g);

        Mesh_data::Field h ("Brown", Mesh_data::SCALAR, Mesh_data::FACE);
        fields.add_field (h);

        CHECK (!fields.empty ());
        int count = 0;
        for (Mesh_data::Fields::const_iterator it = fields.begin ();
             it != fields.end ();
             ++it)
        {
            CHECK_EQUAL (it->location (), count);
            ++count;
        }
        CHECK_EQUAL (3, count);

     }
}
