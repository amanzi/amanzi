#include "Exodus_records.hh"

#include <string.h>

#include "Exodus_error.hh"
#include "exodusII.h"

namespace ExodusII
{


Info_records::Info_records (Exodus_file file)
{
    read_record_count_ (file.id);

    std::vector<char*> record_data (num_records);
    for (int record = 0; record < num_records; ++record)
        record_data [record] = new char [MAX_LINE_LENGTH];

    int ret_val = ex_get_info(file.id, &record_data[0]);

    records.resize (num_records);
    for (int record = 0; record < num_records; ++record)
    {
        char* record_ptr = record_data [record];
        records [record] = std::string (record_ptr, record_ptr + strlen (record_ptr));
        delete [] record_ptr;
    }

}

void Info_records::read_record_count_ (int id)
{
    float f_dummy;
    char c_dummy;

    int ret_val = ex_inquire (id, EX_INQ_INFO, &num_records, &f_dummy, &c_dummy);
    if (ret_val < 0) throw ExodusError (ret_val);
}


void Info_records::to_stream (std::ostream& stream) const
{
    stream << "Information Records:\n";
    stream << "  Number of records: " << num_records << "\n";
    if (num_records > 0)
    {
        stream << "  Records:\n";
        for (int record = 0; record < num_records; ++record)
        {
            stream << record << ". " << records [record];
        }
        stream << "\n";
    }
}

}
