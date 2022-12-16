/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _PorousMedia_Labels_H_
#define _PorousMedia_Labels_H_

#include <map>
#include <string>

namespace PMAMR {
    static std::string RpurposeDEF[7] = {"xlobc", "ylobc", "zlobc", "xhibc", "yhibc", "zhibc", "all"};
    static std::string RlabelDEF[7] = {"XLOBC", "YLOBC", "ZLOBC", "XHIBC", "YHIBC", "ZHIBC", "All"};
}

namespace Amanzi {
    namespace AmanziInput {

        struct GlobalData {
            static std::map<std::string,std::string> AMR_to_Amanzi_label_map;
        };
        inline static std::map<std::string,std::string>& AMRToAmanziLabelMap() {
            return GlobalData::AMR_to_Amanzi_label_map;
        }
    }
}

#endif
