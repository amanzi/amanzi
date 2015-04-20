#ifndef AMANZI_UTILS_LEAST_SQUARE_HH_
#define AMANZI_UTILS_LEAST_SQUARE_HH_

#include <vector>

namespace Amanzi {
namespace Utils {

double bestLSfit(const std::vector<double>& h, const std::vector<double>& error);

}  // namespace Utils
}  // namespace Amanzi

#endif
