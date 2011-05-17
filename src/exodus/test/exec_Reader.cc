#include <iostream>

#include "../Exodus_readers.hh"
#include "Data.hh"

void read_file (const char* filename)
{

  Amanzi::AmanziMesh::Data::Data *mesh = Amanzi::Exodus::read_exodus_file (filename);
  std::cout << *mesh << std::endl << std::endl;

  delete mesh;

}


int main (int argc, const char* argv [])
{
    
  if (argc > 1)
  {
        
    for (int i = 1; i < argc; ++i)
      read_file (argv [i]);

  }

}
