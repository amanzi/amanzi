#ifndef _ELEMENT_CATEGORY_H_
#define _ELEMENT_CATEGORY_H_

namespace STK_mesh
{

enum Element_Category
{
    OWNED = 1,
    GHOST = 2,
    USED  = 3
};

bool valid_category (Element_Category category);

}

#endif /* _ELEMENT_CATEGORY_H_ */
