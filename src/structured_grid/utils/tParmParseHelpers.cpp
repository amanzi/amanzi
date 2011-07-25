#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"	
#include "Teuchos_Version.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Utility.H"
#include "ParmParseHelpers.H"
#include "ccse-mpi.H"
#include "Array.H"

#ifdef BL_USE_OMP
#include "omp.h"
#endif

static double epsilon = 1.e-10;

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    std::string xmlFileName = "sample.xml";
    std::cout << "\n# Reading list from XML file ...\n\n";
    Teuchos::ParameterList params = *Teuchos::getParametersFromXmlFile(xmlFileName);

    // Initialize the ParmParse table using a ParameterList object, where sublists of the 
    // latter correspond to Records of the former.
    BoxLib::Initialize_ParmParse(params);


    // To test, we read in a hierarchical sample file, and real some double arrays with 
    // both type of parameter parsing tools, and then assert that they get the same answer.
    // This test requires that the tool be able to recognize and properly translate strings,
    // string arrays, doubles and double arrays, so its fairly exhaustive, but can be extended.

    ParmParse pp;


    std::string regions_label("regions");
    int num_regions = pp.countval(regions_label.c_str());
    Array<std::string> regions(num_regions);
    pp.getarr(regions_label.c_str(),regions,0,num_regions);

    Teuchos::Array<std::string> tregions = 
      Teuchos::getParameter<Teuchos::Array<std::string> >(params, regions_label);
    BL_ASSERT(tregions.size()==num_regions);

    for (int i=0; i<num_regions; ++i)
      {
        BL_ASSERT(tregions[i] == regions[i]);
        ParmParse::Record ppr = pp.getRecord(regions[i]);
        
        Teuchos::ParameterEntry tregion = params.getEntry(regions[i]);
        BL_ASSERT(tregion.isList());
        Teuchos::ParameterList& sublist = params.sublist(regions[i]);
        
        std::string box_label("box");
        ParmParse::Record pprr = ppr->getRecord(box_label);
        const Teuchos::ParameterEntry& tbox = sublist.getEntry(box_label);
        BL_ASSERT(tbox.isList());
        
        std::string lo_label("lo");
        int num_lo = pprr->countval(lo_label.c_str());
        Array<double> lo(num_lo);
        pprr->getarr(lo_label.c_str(),lo,0,num_lo);
        const Teuchos::Array<double>& tlo = sublist.sublist(box_label).get<Teuchos::Array<double> >(lo_label);
        
        for (int j=0; j<num_lo; ++j)
          {
            BL_ASSERT(std::abs(lo[j] - tlo[j])/(epsilon + lo[j] + tlo[j]) < epsilon);
          }
        
        std::string hi_label("hi");
        int num_hi = pprr->countval(hi_label.c_str());
        Array<double> hi(num_hi);
        pprr->getarr(hi_label.c_str(),hi,0,num_hi);
        const Teuchos::Array<double>& thi = sublist.sublist(box_label).get<Teuchos::Array<double> >(hi_label);
        
        for (int j=0; j<num_hi; ++j)
          {
            BL_ASSERT(std::abs(hi[j] - thi[j])/(epsilon + hi[j] + thi[j]) < epsilon);
          }
        
      }

    ParmParse::Finalize();

    BoxLib::Finalize();
    
    return 0;
}
