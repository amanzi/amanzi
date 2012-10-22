#ifndef _PorousMedia_Labels_H_
#define _PorousMedia_Labels_H_

#include <map>
#include <string>

namespace PMAMR {
    static std::string RpurposeDEF[7] = {"xlobc", "ylobc", "zlobc", "xhibc", "yhibc", "zhibc", "all"};
}

namespace Amanzi {
    namespace AmanziInput {

        struct GlobalData {
            static std::map<std::string,std::string> AMR_to_Amanzi_label_map;
        };
        static std::map<std::string,std::string>& AMRToAmanziLabelMap() {
            return GlobalData::AMR_to_Amanzi_label_map;
        }
    }
}

#endif
