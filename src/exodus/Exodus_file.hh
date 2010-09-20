#ifndef _EXODUS_FILE_HH_
#define _EXODUS_FILE_HH_

#include <iostream>


namespace ExodusII
{

struct Exodus_file
{
    int id;
    char const * filename;
    int exodus_word_size;
    int system_word_size;

    float version;
    
    Exodus_file (char const * filename);

    void to_stream (std::ostream& stream) const;
};





} // namespace ExodusII

#endif
