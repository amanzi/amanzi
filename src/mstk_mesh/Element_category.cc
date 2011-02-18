#include "Element_category.hh"

namespace MSTK_mesh
{

bool valid_category (Element_Category category)
{
    return (category >= OWNED && category <= USED);
}

}
