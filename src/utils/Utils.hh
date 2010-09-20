#ifndef _UTILS_HH_
#define _UTILS_HH_

#include <vector>

namespace Utils
{


template <typename P>
void reclaim (std::vector<P*>& v)
{

    for (typename std::vector<P*>::iterator it = v.begin (); it != v.end (); ++it)
        delete *it;
}

}


#endif
