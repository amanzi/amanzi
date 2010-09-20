#include <iostream>

#include "../Exodus_readers.hh"
#include "Data.hh"

int main (int argc, const char* argv [])
{
    
    Mesh_data::Data* mesh = ExodusII::read_exodus_file ("htc_rad_test-random.exo");
    std::cout << *mesh << std::endl << std::endl;;

    mesh = ExodusII::read_exodus_file ("cubit.e");
    std::cout << *mesh << std::endl << std::endl;;

    mesh = ExodusII::read_exodus_file ("freefem.e");
    std::cout << *mesh << std::endl;

}
