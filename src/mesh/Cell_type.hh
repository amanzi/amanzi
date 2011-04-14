#ifndef _AMZ_CELL_TYPE_H_
#define _AMZ_CELL_TYPE_H_


// Names of standard element types

namespace Amz_Mesh
{


enum Cell_type
{
    TRI,
    QUAD,
    POLYGON,
    TET,
    PRISM,
    PYRAMID,
    HEX,
    POLYHED
};

bool valid_cell_type (Cell_type type);
}

#endif /* _AMZ_CELL_TYPE_H_ */
