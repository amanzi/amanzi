
#include <Property.H>
#include <GSLibInt.H>
#include <MatFiller.H>

#include <Utility.H>
#include <WritePlotfile.H>

static int crse_init_factor = 32;
static int max_grid_fine_gen = 32;

static std::string MaterialPlotFileVersion = "MaterialPlotFile-V1.0";

Property*
ConstantProperty::clone() const
{
  const ConstantProperty* t = dynamic_cast<const ConstantProperty*>(this);
  BL_ASSERT(t!=0);
  ConstantProperty* ret = new ConstantProperty(t->Name(), t->values, t->coarsen_rule, t->refine_rule);
  return ret;
}

bool
ConstantProperty::Evaluate(Real t, Array<Real>& result) const
{
  int N = values.size();
  result.resize(N);
  for (int i=0; i<N; ++i) {
    result[i] = values[i];
  }
  return true;
}

Property*
GSLibProperty::clone() const
{
  const GSLibProperty* t = dynamic_cast<const GSLibProperty*>(this);
  BL_ASSERT(t!=0);
  Array<Real> shift(BL_SPACEDIM,0);
  GSLibProperty* ret = new GSLibProperty(t->Name(), t->avg, t->std, t->param_file, t->data_file, shift, t->coarsen_rule, t->refine_rule);
  if (t->dataServices == 0) {
    ret->dataServices = 0;
  }
  else {
    ret->dataServices = new DataServices(t->dataServices->GetFileName(), t->dataServices->GetFileType());
  }
  ret->num_comps = t->num_comps;
  ret->varnames = t->PlotfileVars();
  return ret;
}

static
void
EnsureFolderExists(const std::string& full_path)
{
  // Find folder name first, and ensure folder exists
  // FIXME: Will fail on Windows
  const std::vector<std::string>& tokens = BoxLib::Tokenize(full_path,"/");
  std::string dir = (full_path[0] == '/' ? "/" : "");
  for (int i=0; i<tokens.size()-1; ++i) {
    dir += tokens[i];
    if (i<tokens.size()-2) dir += "/";
  }
  
  if(!BoxLib::FileExists(dir)) {
    if ( ! BoxLib::UtilCreateDirectory(dir, 0755)) {
      BoxLib::CreateDirectoryFailed(dir);
    }
  }
}

void
GSLibProperty::BuildGSLibFile(Real                   avg,
                              Real                   std,
                              const std::string&     gslib_param_file,
                              const std::string&     gslib_data_file,
                              const Array<Geometry>& geom_array,
                              const Array<IntVect>&  ref_ratio,
                              int                    num_grow,
                              int                    max_grid_size_fine_gen,
                              Property::CoarsenRule  crule,
                              const std::string&     varname)
{
  int nLev = geom_array.size();
  int finest_level = nLev - 1;
  const Geometry& geom = geom_array[finest_level];
  Box stat_box = Box(geom.Domain());

  if(!BoxLib::FileExists(gslib_param_file)) {
    const std::string& str = "GSLib parameter file: \"" + gslib_param_file + "\" does not exist";
    BoxLib::Abort(str.c_str());
  }

  // Original interface supports layered structure, we disable that for now
  Array<Real> avgVals(1,avg);
  Array<int> n_cell(BL_SPACEDIM);
  const Geometry& geom0 = geom_array[0];
  for (int d=0; d<BL_SPACEDIM; ++d) {
    n_cell[d] = geom0.Domain().length(d);
  }

  // Find cummulative refinement ratio
  int twoexp = 1;
  for (int i = 1; i<nLev; i++) {
    twoexp *= ref_ratio[i-1][0];  // FIXME: Assumes uniform refinement
  }

  PArray<MultiFab> stat(nLev,PArrayManage);
  BoxArray stat_ba(stat_box);
  stat_ba.maxSize(max_grid_size_fine_gen);
  int ng_cum = num_grow * twoexp;
  num_comps = (crule == ComponentHarmonic  ?  BL_SPACEDIM : 1);
  stat.set(finest_level, new MultiFab(stat_ba,num_comps,ng_cum));

  const Array<Real> prob_lo(geom0.ProbLo(),BL_SPACEDIM);
  const Array<Real> prob_hi(geom0.ProbHi(),BL_SPACEDIM);
  
  GSLibInt::rdpGaussianSim(avgVals,std,n_cell,prob_lo,prob_hi,twoexp,stat[finest_level],
                           crse_init_factor,max_grid_size_fine_gen,ng_cum,gslib_param_file);

  for (int lev=finest_level-1; lev>=0; --lev) {
    int ltwoexp = 1;
    for (int i = 1; i<lev; i++) {
      ltwoexp *= ref_ratio[i-1][0];  // FIXME: Assumes uniform refinement
    }

    const Box& domain = geom_array[lev].Domain();
    BoxArray ba(domain);
    ba.maxSize(max_grid_size_fine_gen / ref_ratio[lev][0]); // FIXME: Assumes uniform refinement
    stat.set(lev, new MultiFab(ba,num_comps,num_grow*ltwoexp));

    BoxArray baf = BoxArray(ba).refine(ref_ratio[lev]);
    MultiFab fine(baf,1,stat[lev].nGrow()*ref_ratio[lev][0]);// FIXME: Assumes uniform refinement
    BoxArray bafg = BoxArray(baf).grow(fine.nGrow());
    MultiFab fineg(bafg,1,0);
    fineg.copy(stat[lev+1]); // parallel copy
    for (MFIter mfi(fine); mfi.isValid(); ++mfi) {
      for (int n=0; n<num_comps; ++n) {
	const Box& bx = fineg[mfi].box();
	fine[mfi].copy(fineg[mfi],bx,0,bx,n,1);
      }
    }
    fineg.clear();

    for (MFIter mfi(fine); mfi.isValid(); ++mfi) {
      const FArrayBox& finefab = fine[mfi];
      FArrayBox& crsefab = stat[lev][mfi];
      const Box& cbox = crsefab.box();
      if ( !(finefab.box().contains(Box(cbox).refine(ref_ratio[lev])) ) ) {
        std::cout << "c,f: " << cbox << " " << finefab.box() << std::endl;
        BoxLib::Abort();
      }
      MatFiller::CoarsenData(fine[mfi],0,stat[lev][mfi],cbox,0,1,ref_ratio[lev],crule);
    }
  }

  EnsureFolderExists(gslib_data_file);
  
  Array<MultiFab*> data(nLev);
  Array<Box> prob_domain(nLev);
  Array<Array<Real> > dx_level(nLev,Array<Real>(BL_SPACEDIM));
  Array<int> int_ref(nLev-1);
  for (int lev=0; lev<nLev; ++lev) {
    data[lev] = &(stat[lev]);
    prob_domain[lev] = geom_array[lev].Domain();
    for (int d=0; d<BL_SPACEDIM; ++d) {
      dx_level[lev][d] = geom_array[lev].CellSize(d);
    }
    if (lev<finest_level) {
      int_ref[lev] = ref_ratio[lev][0];
    }
  }
  Real time=0;
  bool verbose=false;
  Array<Real> vfeps(BL_SPACEDIM,0);
  Array<int> level_steps(nLev,0);
  varnames.resize(num_comps);
  for (int n=0; n<num_comps; ++n) {
    varnames[n] = BoxLib::Concatenate(varname+"_",n,1);
  }
  bool is_cart_grid = false;
  WritePlotfile(MaterialPlotFileVersion,data,time,geom0.ProbLo(),geom0.ProbHi(),int_ref,prob_domain,
                dx_level,geom0.Coord(),gslib_data_file,varnames,verbose,is_cart_grid,vfeps.dataPtr(),
                level_steps.dataPtr());
  ParallelDescriptor::Barrier(); // Wait until everyone finished to avoid reading before completely written
}

void
GSLibProperty::BuildDataFile(const Array<Geometry>& geom_array,
                             const Array<IntVect>&  ref_ratio,
                             int                    num_grow,
                             int                    max_grid_size_fine_gen,
                             Property::CoarsenRule  crule,
                             const std::string&     varname)
{
  BuildGSLibFile(avg, std, param_file, data_file,
                 geom_array, ref_ratio, num_grow,
                 max_grid_fine_gen, crule, varname);

  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  delete dataServices;
  dataServices = new DataServices(data_file, fileType);
  if (!dataServices->AmrDataOk()) {
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
  }
}

const AmrData*
GSLibProperty::GetAmrData() const
{
  if (dataServices == 0) {
    BoxLib::Abort("GSLib file not initialized");
  }
  return &(dataServices->AmrDataRef());
}


bool
GSLibProperty::Evaluate(Real t, Array<Real>& result) const
{
  result.resize(num_comps);
  for (int i=0; i<num_comps; ++i) {
    result[i] = values[0];
  }
  return false;
}


Property*
TabularInTimeProperty::clone() const
{
  const TabularInTimeProperty* t = dynamic_cast<const TabularInTimeProperty*>(this);
  BL_ASSERT(t!=0);
  TabularInTimeProperty* ret = new TabularInTimeProperty(t->Name(), t->Functions(), t->coarsen_rule, t->refine_rule);
  return ret;
}

bool
TabularInTimeProperty::Evaluate(Real t, Array<Real>& result) const
{
  int N = funcs.size();
  result.resize(N);
  for (int i=0; i<N; ++i) {
    result[i] = funcs[i](t);
  }
  return true;
}
