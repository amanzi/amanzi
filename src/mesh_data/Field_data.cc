#include "Field_data.hh"
#include "dbc.hh"


namespace Mesh_data
{

Field::Field (std::string name, FIELD_TYPE type, Entity_kind location) : 
    name_ (name), type_ (type), location_ (location)
{
    ASSERT (ok_field_type (type));
    ASSERT (valid_entity_kind (location));
}

}
