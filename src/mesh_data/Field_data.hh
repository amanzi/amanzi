#ifndef _FIELD_DATA_HH_
#define _FIELD_DATA_HH_

#include <list>
#include <string>

#include "Entity_kind.hh"
#include "Element_field_types.hh"

namespace Mesh_data
{

class Field
{

    Entity_kind location_;
    FIELD_TYPE type_;
    std::string name_;

public:

    FIELD_TYPE type () const { return type_; }
    Entity_kind location () const { return location_; }
    std::string name () const { return name_; }

    Field (std::string name, FIELD_TYPE t = SCALAR, Entity_kind l = CELL);

};


class Fields
{

    typedef std::list<Field> Fields_;
    Fields_ fields_;

public:

    typedef Fields_::const_iterator const_iterator;

    Fields () : fields_ () { }

    void add_field (Field& field) { fields_.push_back (field); }

    const_iterator begin () const { return fields_.begin (); }
    const_iterator end () const { return fields_.end (); }
    
    bool empty () const { return fields_.empty (); }

};

}

#endif
