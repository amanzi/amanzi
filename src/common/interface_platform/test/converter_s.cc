/*
  This is the input component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Jeffrey Johnson (jnjohnson@lbl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "ParmParse.H"
#include "UnitTest++.h"

// Amanzi
#include "InputConverterS.hh"


void*
getLevelBld()
{
  return NULL;
}


/* **************************************************************** */
TEST(CONVERTER_S)
{
  using namespace Amanzi;

  int rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  rank = 0;
  if (rank == 0) std::cout << "Test: convert v2.x -> ParmParse test" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/converter_s_input.xml";

  Amanzi::AmanziInput::InputConverterS converter(xmlFileName);
  try {
    // Translate the input. This produces a singleton instance of ParmParse that is
    // populated with data.
    converter.Translate(rank);

    // Dump the guy to stdout.
    if (rank == 0) {
      ParmParse pp;
      std::ofstream dumpy("converter_s_test_output.dump");
      pp.dumpTable(dumpy);

      std::cout << "Successful translation. Validating the result...\n\n";
    }
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
