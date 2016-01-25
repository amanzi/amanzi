#ifndef _UTILS_MESSAGE_H_
#define _UTILS_MESSAGE_H_

#include <algorithm>
#include <string>

namespace Errors {

void encode_string(std::string& s, int n, int m, int* out) {
  for (int i = 0; i < n; ++i) out[i] = m + 200 + (int)s[i];
  out[n] = m + 200;
}

void decode_string(int* in, int n, std::string& s) {
  int m = in[n];
  for (int i = 0; i < n; ++i) s[i] = in[i] - m;
}

}  // namespace Errors
#endif
