#ifndef _EXODUS_RECORDS_HH_
#define _EXODUS_RECORDS_HH_

#include <ostream>
#include <vector>

#include "Exodus_file.hh"

namespace Amanzi {
namespace Exodus {

struct Info_records
{
  Info_records (Exodus_file file);

  int num_records;
  std::vector<std::string> records;

  void to_stream (std::ostream& stream) const;

 private:

  void read_record_count_ (int id);

};

} // namespace Exodus
} // namespace Amanzi

#endif
