#include <iostream>
#include <fstream>
#include <iomanip>
using std::cout;
using std::endl;
#include <ParmParse.H>

#include <MatFillerPCarithAvg.H>
#include <Region.H>

static std::map<std::string, Region*> regions;
static PArray<Material> materials;

std::ostream& operator<<(std::ostream& os, const Material& mat) {
  std::cout << "Material: " << mat.Name() << std::endl;
  std::cout << "  Regions: ";
  const PArray<Region>& regions = mat.Regions();
  for (int i=0; i<regions.size(); ++i) {
    std::cout << regions[i].name << " ";
  }
  std::cout << std::endl;
  return os;
}

PArray<Region>
build_region_PArray(const Array<std::string>& region_names)
{
  PArray<Region> ret(region_names.size(), PArrayNoManage);

  for (int i=0; i<region_names.size(); ++i)
  {
    const std::string& name = region_names[i];
    std::map<std::string, Region*>::const_iterator it = regions.find(name);
    if (it != regions.end()) {
      ret.set(i,it->second);
    }
    else {
      std::string m = "Named region not found: " + name;
      BoxLib::Error(m.c_str());
    }
  }
  return ret;
}

void SetRegionsTANK()
{
  Array<Real> lo(2), hi(2);
  std::string r_name, r_purpose, r_type;
  r_name = "SoilLower";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0;
  lo[1] = 0;
  hi[0] = 40;
  hi[1] = 10;
  regions["SoilLower"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "SoilRight";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 12;
  lo[1] = 10;
  hi[0] = 40;
  hi[1] = 18;
  regions["SoilRight"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "SoilUpper";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0;
  lo[1] = 18;
  hi[0] = 40;
  hi[1] = 24;
  regions["SoilUpper"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "TankConcFloor";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0.5;
  lo[1] = 10;
  hi[0] = 12;
  hi[1] = 10.3;
  regions["TankConcFloor"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "TankConcRoof1";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0;
  lo[1] = 17.7;
  hi[0] = 11.67;
  hi[1] = 18;
  regions["TankConcRoof1"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "TankConcRoof2";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 11.69;
  lo[1] = 17.7;
  hi[0] = 12;
  hi[1] = 18;
  regions["TankConcRoof2"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "TankConcWall";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 11.7;
  lo[1] = 10.3;
  hi[0] = 12;
  hi[1] = 17.7;
  regions["TankConcWall"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "TankFFfloor";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0;
  lo[1] = 10;
  hi[0] = 0.5;
  hi[1] = 10.31;
  regions["TankFFfloor"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "TankFFwall";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 11.67;
  lo[1] = 10.33;
  hi[0] = 11.69;
  hi[1] = 18;
  regions["TankFFwall"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "TankGrout";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0;
  lo[1] = 10.33;
  hi[0] = 11.67;
  hi[1] = 17.69;
  regions["TankGrout"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "TankLinerFloor";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0.5;
  lo[1] = 10.3;
  hi[0] = 11.7;
  hi[1] = 10.31;
  regions["TankLinerFloor"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "TankLinerRoof";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0;
  lo[1] = 17.69;
  hi[0] = 11.67;
  hi[1] = 17.7;
  regions["TankLinerRoof"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "TankLinerWall";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 11.69;
  lo[1] = 10.31;
  hi[0] = 11.7;
  hi[1] = 17.7;
  regions["TankLinerWall"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "TankWaste";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0;
  lo[1] = 10.31;
  hi[0] = 11.69;
  hi[1] = 10.33;
  regions["TankWaste"] = new boxRegion(r_name,r_purpose,r_type,lo,hi);
}

void
SetMaterialsTANK()
{
  std::string phi_str = "porosity";
  std::map<std::string,Property*> property_map; property_map[phi_str] = 0;
  materials.resize(5, PArrayManage);
  Array<std::string> region_names;
  PArray<Region> regionset;
  region_names.push_back("SoilLower");
  region_names.push_back("SoilRight");
  region_names.push_back("SoilUpper");
  regionset = build_region_PArray(region_names);
  delete property_map[phi_str]; property_map[phi_str] = new ConstantProperty(1);
  materials.set(0,new Material("Soil",regionset,property_map));

  region_names.clear();
  region_names.push_back("TankConcFloor");
  region_names.push_back("TankConcRoof1");
  region_names.push_back("TankConcRoof2");
  region_names.push_back("TankConcWall");
  regionset = build_region_PArray(region_names);
  delete property_map[phi_str]; property_map[phi_str] = new ConstantProperty(2);
  materials.set(1,new Material("TankConc",regionset,property_map));

  region_names.clear();
  region_names.push_back("TankFFfloor");
  region_names.push_back("TankFFwall");
  region_names.push_back("TankWaste");
  regionset = build_region_PArray(region_names);
  delete property_map[phi_str]; property_map[phi_str] = new ConstantProperty(3);
  materials.set(2,new Material("TankFF",regionset,property_map));

  region_names.clear();
  region_names.push_back("TankGrout");
  regionset = build_region_PArray(region_names);
  delete property_map[phi_str]; property_map[phi_str] = new ConstantProperty(4);
  materials.set(3,new Material("TankGrout",regionset,property_map));

  region_names.clear();
  region_names.push_back("TankLinerFloor");
  region_names.push_back("TankLinerRoof");
  region_names.push_back("TankLinerWall");
  regionset = build_region_PArray(region_names);
  //delete property_map[phi_str]; property_map[phi_str] = new ConstantProperty(5);
  int nvals = 2;
  Array<double> p_values(nvals); p_values[0] = 5; p_values[1] = 6;
  Array<double> p_times(nvals); p_times[0] = 1; p_times[1] = 2;
  Array<std::string> p_forms(nvals-1); p_forms[0] = "Linear";
  TabularFunction pft(p_times,p_values,p_forms);
  delete property_map[phi_str]; property_map[phi_str] = new TabularInTimeProperty(pft);

  materials.set(4,new Material("TankLiner",regionset,property_map));

  region_names.clear();
}

void SetRegionsSMALL()
{
  Array<Real> lo(2), hi(2);
  std::string r_name, r_purpose, r_type;
  r_name = "Blue";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 5;
  lo[1] = 1;
  hi[0] = 7;
  hi[1] = 3;
  regions[r_name] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "Red1";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0;
  lo[1] = 3;
  hi[0] = 8;
  hi[1] = 8;
  regions[r_name] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "Red2";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0;
  lo[1] = 1;
  hi[0] = 5;
  hi[1] = 3;
  regions[r_name] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "Red3";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 0;
  lo[1] = 0;
  hi[0] = 8;
  hi[1] = 1;
  regions[r_name] = new boxRegion(r_name,r_purpose,r_type,lo,hi);

  r_name = "Red4";
  r_purpose = "all";
  r_type = "box";
  lo[0] = 7;
  lo[1] = 1;
  hi[0] = 8;
  hi[1] = 3;
  regions[r_name] = new boxRegion(r_name,r_purpose,r_type,lo,hi);
}

void
SetMaterialsSMALL()
{
  std::string phi_str = "porosity";
  std::map<std::string,Property*> property_map; property_map[phi_str] = 0;
  materials.resize(2, PArrayManage);
  Array<std::string> region_names;
  PArray<Region> regionset;
  region_names.push_back("Blue");
  regionset = build_region_PArray(region_names);
  delete property_map[phi_str]; property_map[phi_str] = new ConstantProperty(1);
  materials.set(0,new Material("Blue",regionset,property_map));

  region_names.clear();
  region_names.push_back("Red1");
  region_names.push_back("Red2");
  region_names.push_back("Red3");
  region_names.push_back("Red4");
  regionset = build_region_PArray(region_names);
  delete property_map[phi_str]; property_map[phi_str] = new ConstantProperty(2);
  materials.set(1,new Material("Red",regionset,property_map));

  region_names.clear();
}



void
DestroyRegions()
{
  for (std::map<std::string, Region*>::iterator it=regions.begin(), End=regions.end(); it!=End; ++it) {
    delete it->second;
  }
}

void
DestroyMaterials()
{
  materials.clear();
}

void WritePlotfile(const std::string         &pfversion,
                   const PArray<MultiFab>    &data,
                   const Real                 time,
                   const Real                *probLo,
                   const Real                *probHi,
                   const Array<int>          &refRatio,
                   const Array<Box>          &probDomain,
                   const Array<Array<Real> > &dxLevel,
                   const int                  coordSys,
                   const std::string         &oFile,
                   const Array<std::string>  &names,
                   const bool                 verbose,
		   const bool                 isCartGrid,
		   const Real                *vfeps,
		   const int                 *levelSteps);

int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);

  ParmParse pp;
  std::string case_size="small"; pp.query("case_size",case_size);

  Array<int> n_cells(BL_SPACEDIM); 
  int nLevs;
  int coord = 0;
  Array<int> is_per(BL_SPACEDIM,0);
  Array<Real> prob_lo(BL_SPACEDIM,0);
  Array<Real> prob_hi(BL_SPACEDIM);
  if (case_size=="small") {
    nLevs = 2;
    D_EXPR(n_cells[0] = 2, n_cells[1] = 2, n_cells[2] = 1);
    D_EXPR(prob_hi[0] = 8, prob_hi[1] = 8, prob_hi[2] = 8);
  }
  else {
    if (case_size=="medium") {
      nLevs = 2;
    } else if (case_size=="large"){
      nLevs = 3;
    } else if (case_size=="xlarge"){
      nLevs = 4;
    }
    else {
      BoxLib::Abort("No plan for this case_size");
    }
    D_EXPR(n_cells[0] = 40, n_cells[1] = 24, n_cells[2] = 1);
    D_EXPR(prob_hi[0] = 40, prob_hi[1] = 24, prob_hi[2] = 1);
  }

  Array<int> rRatio(nLevs-1,4);
  Array<IntVect> refRatio(nLevs-1);
  for (int lev=0; lev<nLevs-1; ++lev) {
    refRatio[lev] = rRatio[lev] * IntVect::TheUnitVector();
  }

  const RealBox rb(prob_lo.dataPtr(),prob_hi.dataPtr());

  Array<Geometry> geomArray(nLevs);
  for (int lev=0; lev<nLevs; ++lev) {
    Box domain;
    if (lev==0) {
      IntVect be;
      for (int d=0; d<BL_SPACEDIM; ++d) {
        be[d] = n_cells[d] - 1;
      }
      domain = Box(IntVect(D_DECL(0,0,0)),be);
    }
    else {
      domain = Box(geomArray[lev-1].Domain()).refine(refRatio[lev-1]);
    }
    geomArray[lev] = Geometry(domain,&rb,coord,is_per.dataPtr());
  }

  Region::geometry_eps = 1.e-6;
  Region::domlo.resize(BL_SPACEDIM);
  Region::domhi.resize(BL_SPACEDIM);
  for (int d=0; d<BL_SPACEDIM; ++d) {
    Region::domlo[d] = Geometry::ProbLo()[d];
    Region::domhi[d] = Geometry::ProbHi()[d];
  }

  if (case_size=="small") {
    SetRegionsSMALL();
    SetMaterialsSMALL();
  }
  else {
    SetRegionsTANK();
    SetMaterialsTANK();
  }

  bool verbose = false; pp.query("verbose",verbose);
  if (verbose) {
    // Echo regions and materials
    if (ParallelDescriptor::IOProcessor()) {
      for (std::map<std::string,Region*>::const_iterator it=regions.begin(); it!=regions.end(); ++it) {
        std::cout << *(it->second) << std::endl;
      }
      
      for (int i=0; i<materials.size(); ++i) {
        std::cout << materials[i] << std::endl;
      }
    }
  }

  PArray<Property> properties(materials.size());
  for (int i=0, N=materials.size(); i<N; ++i) {
    Real m_porosity = (Real)(i+1)/N;
    properties.set(i, new ConstantProperty(m_porosity));
  }

  MatFillerPCarithAvg matFiller(geomArray,refRatio,materials);
  const std::map<std::string,int>& mat_map = matFiller.MatIdx();

  bool fail = false;

  int nout = 1;
  Array<std::string> propNames(nout);
  propNames[0] = "porosity";

  Real t = 1.5;
  PArray<MultiFab> data(nLevs,PArrayManage);
  for (int lev=nLevs-1; lev>=0; --lev) {
    BoxArray ba(matFiller.MaterialID(lev).boxArray());
    data.set(lev, new MultiFab(ba,nout,0));
    for (int n=0; n<nout; ++n) {
      matFiller.SetProperty(t,lev,data[lev],propNames[n],n,0);
    }
  }

  Array<Real> bins(materials.size(),0);
  for (int lev=0; lev<nLevs; ++lev) {
    for (MFIter mfi(data[lev]); mfi.isValid(); ++mfi) {
      const Box& vbox = mfi.validbox();
      const FArrayBox& fab = data[lev][mfi];
      const FArrayBox& id = matFiller.MaterialID(lev)[mfi];
      for (IntVect iv=vbox.smallEnd(), BIG=vbox.bigEnd(); iv<=BIG; vbox.next(iv)) {
	int bin = (int) id(iv,0);
	int val = fab(iv,0);
	if (bin>=0) {
	  bins[bin] += val;
	}
      }
    }
  }

  ParallelDescriptor::ReduceRealSum(bins.dataPtr(),bins.size());

  if (case_size=="small") {
    fail &= bins[0] == 8;
    fail &= bins[1] == 30;
  }
  else {
    fail &= bins[0] == 1464;
    fail &= bins[1] == 15726;
    fail &= bins[2] == 3837;
    fail &= bins[3] == 29028;
    fail &= bins[4] == 9685;
  }

  t = 2.5;
  for (int lev=nLevs-1; lev>=0; --lev) {
    for (int n=0; n<nout; ++n) {
      matFiller.SetProperty(t,lev,data[lev],propNames[n],n,0);
    }
  }

  for (int lev=0; lev<nLevs; ++lev) {
    for (MFIter mfi(data[lev]); mfi.isValid(); ++mfi) {
      const Box& vbox = mfi.validbox();
      const FArrayBox& fab = data[lev][mfi];
      const FArrayBox& id = matFiller.MaterialID(lev)[mfi];
      for (IntVect iv=vbox.smallEnd(), BIG=vbox.bigEnd(); iv<=BIG; vbox.next(iv)) {
	int bin = (int) id(iv,0);
	int val = fab(iv,0);
	if (bin>=0) {
	  bins[bin] += val;
	}
      }
    }
  }

  ParallelDescriptor::ReduceRealSum(bins.dataPtr(),bins.size());
  if (case_size=="small") {
    fail &= bins[0] == 8;
    fail &= bins[1] == 60;
  }
  else {
    fail &= bins[0] == 1464;
    fail &= bins[1] == 15726;
    fail &= bins[2] == 3837;
    fail &= bins[3] == 29028;
    fail &= bins[4] == 11622;
  }

  if (verbose) {
    // Write out result to pltfile
    std::string pfversion = "MaterialData-0.2";
    Array<Box> pDomain(nLevs);
    Array<Array<Real> > dxLevel(nLevs,Array<Real>(BL_SPACEDIM));
    for (int i=0; i<nLevs; ++i) {
      pDomain[i] = geomArray[i].Domain();
      for (int d=0; d<BL_SPACEDIM; ++d) {
        dxLevel[i][d] = geomArray[i].CellSize()[d];
      }
    }
    int coordSys = (int)CoordSys::Coord();
    std::string outFileName = "pltfile";
    bool verbose = false;
    bool isCartGrid = false;
    Array<Real> vfeps(nLevs,1.e-10);
    Array<int> levelSteps(nLevs,0);
    
    WritePlotfile(pfversion,data,t,Geometry::ProbLo(),Geometry::ProbHi(),
                  rRatio,pDomain,dxLevel,coordSys,outFileName,propNames,
                  verbose,isCartGrid,vfeps.dataPtr(),levelSteps.dataPtr());
  }

  DestroyMaterials();
  DestroyRegions();
  FabArrayBase::verbose = false;
  BoxLib::Finalize();
  if (fail) {
    BoxLib::Abort();
  }
  return 0;
}

void WritePlotfile(const std::string         &pfversion,
                   const PArray<MultiFab>    &data,
                   const Real                 time,
                   const Real                *probLo,
                   const Real                *probHi,
                   const Array<int>          &refRatio,
                   const Array<Box>          &probDomain,
                   const Array<Array<Real> > &dxLevel,
                   const int                  coordSys,
                   const std::string         &oFile,
                   const Array<std::string>  &names,
                   const bool                 verbose,
		   const bool                 isCartGrid,
		   const Real                *vfeps,
		   const int                 *levelSteps)
{
    if(ParallelDescriptor::IOProcessor()) {
      if( ! BoxLib::UtilCreateDirectory(oFile,0755)) {
         BoxLib::CreateDirectoryFailed(oFile);
      }
    }
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
    
    std::string oFileHeader(oFile);
    oFileHeader += "/Header";
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    
    std::ofstream os;
    
    //os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    
    if(verbose && ParallelDescriptor::IOProcessor()) {
      std::cout << "Opening file = " << oFileHeader << '\n';
    }
    
    os.open(oFileHeader.c_str(), std::ios::out|std::ios::binary);
    
    if(os.fail()) {
      BoxLib::FileOpenFailed(oFileHeader);
    }
    //
    // Start writing plotfile.
    //
    os << pfversion << '\n';
    int n_var = data[0].nComp();
    os << n_var << '\n';
    for (int n = 0; n < n_var; n++) os << names[n] << '\n';
    os << BL_SPACEDIM << '\n';
    os << std::setprecision(30) << time << '\n';
    const int finestLevel = data.size() - 1;
    os << finestLevel << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++) os << probLo[i] << ' ';
    os << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++) os << probHi[i] << ' ';
    os << '\n';
    for (int i = 0; i < finestLevel; i++) os << refRatio[i] << ' ';
    os << '\n';
    for (int i = 0; i <= finestLevel; i++) os << probDomain[i] << ' ';
    os << '\n';
    if(levelSteps != 0) {
      for (int i = 0; i <= finestLevel; i++) os << levelSteps[i] << ' ';
    } else {
      for (int i = 0; i <= finestLevel; i++) os << 0 << ' ';
    }
    os << '\n';
    for(int i = 0; i <= finestLevel; i++) {
      for(int k = 0; k < BL_SPACEDIM; k++) {
            os << dxLevel[i][k] << ' ';
      }
      os << '\n';
    }
    if(isCartGrid) {
      for(int i(0); i <= finestLevel; i++) {
        os << vfeps[i] << ' ';
      }
      os << '\n';
    }
    os << coordSys << '\n';
    os << 0 << '\n';                  // --------------- The bndry data width.
    //
    // Write out level by level.
    //
    for(int iLevel(0); iLevel <= finestLevel; ++iLevel) {
        //
        // Write state data.
        //
        const BoxArray& ba = data[iLevel].boxArray();
        int nGrids = ba.size();
        char buf[64];
        sprintf(buf, "Level_%d", iLevel);
        
        if(ParallelDescriptor::IOProcessor()) {
            os << iLevel << ' ' << nGrids << ' ' << time << '\n';
            if(levelSteps != 0) {
              os << levelSteps[iLevel] << '\n';
	    } else {
              os << 0 << '\n';
	    }
            
            for(int i(0); i < nGrids; ++i) {
              const Box &b = ba[i];
              for(int n(0); n < BL_SPACEDIM; ++n) {
                Real glo = b.smallEnd()[n] * dxLevel[iLevel][n];
                Real ghi = (b.bigEnd()[n]+1) * dxLevel[iLevel][n];
                os << glo << ' ' << ghi << '\n';
              }
            }
            //
            // Build the directory to hold the MultiFabs at this level.
            //
            std::string Level(oFile);
            Level += '/';
            Level += buf;
            
            if( ! BoxLib::UtilCreateDirectory(Level, 0755)) {
              BoxLib::CreateDirectoryFailed(Level);
	    }
        }
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        //
        // Now build the full relative pathname of the MultiFab.
        //
        static const std::string MultiFabBaseName("MultiFab");
        
        std::string PathName(oFile);
        PathName += '/';
        PathName += buf;
        PathName += '/';
        PathName += MultiFabBaseName;
        
        if(ParallelDescriptor::IOProcessor()) {
            //
            // The full name relative to the Header file.
            //
            std::string RelativePathName(buf);
            RelativePathName += '/';
            RelativePathName += MultiFabBaseName;
            os << RelativePathName << '\n';
        }
        VisMF::Write(data[iLevel], PathName);
    }
    
    os.close();
}
