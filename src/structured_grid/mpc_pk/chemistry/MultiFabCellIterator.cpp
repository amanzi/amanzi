#include <MultiFabCellIterator.H>

MultiFabCellIterator::MultiFabCellIterator(const DistributionMapping& dm,
					   const BoxArray&            ba,
					   int                        ng)
  : dmap(dm), grids(ba), num_grow(ng)
{
  initialize();
}

MultiFabCellIterator::MultiFabCellIterator(const MultiFabCellIterator& rhs)
{
  *this = rhs;
}

IntVect
MultiFabCellIterator::idx_to_intvect(const Box& box, long idx)
{
  IntVect retIdx;
  IntVect size = box.size();
  long newVal = idx;
#if BL_SPACEDIM==3
  retIdx[2] = (int)(newVal/(size[0]*size[1]));
  newVal -= retVal[2] * size[0]*size[1];
#endif
  retIdx[1] = (int)(newVal/size[0]);
  retIdx[0] = newVal - retIdx[1] * size[0];
  retIdx += box.smallEnd();
  return retIdx;
}

void 
MultiFabCellIterator::set_cell_index_and_boxID(long idx)
{
  if (idx >= end_index) {
    current_box_id = -1;
    current_index = end_index;
    return;
  }
  current_index = idx;
  if (current_box_id<0 || !(current_index>=offsets[current_box_id]
			    && (current_index == index_map.size()-1 || current_index<offsets[current_box_id+1]) ) )
    {
      current_box_id=-1;
      bool done = false;
      for (int i=0, End=index_map.size(); i<End && !done; ++i) {
        if (current_index >= offsets[i]) {
          current_box_id = i;
        }
        else { 
          done = true;
        }
      }
      current_box = Box(grids[index_map[current_box_id]]).grow(num_grow);
    }
  current_cell_index = idx_to_intvect(current_box,idx-offsets[current_box_id]);
}

void
MultiFabCellIterator::initialize()
{
  int myproc = ParallelDescriptor::MyProc();
  for (int i=0; i<grids.size(); ++i) {
    if (dmap[i] == myproc) {
      index_map.push_back(i);
      if (index_map.size()==1) {
	offsets.push_back(0);
      }
      else {
	Box tmpbox = Box(grids[index_map[index_map.size()-2]]).grow(num_grow);
	long numPts = tmpbox.numPts();
	offsets.push_back(offsets[offsets.size()-1] + numPts);
      }
    }
  }
  BL_ASSERT(index_map.size()==offsets.size());
  end_index = begin_index = 0;
  current_box_id = -1;    
  if (!index_map.empty()) { 
    Box tmpbox = Box(grids[index_map[index_map.size()-1]]).grow(num_grow);
    long numPts = tmpbox.numPts();
    end_index = offsets[offsets.size()-1] + numPts;
  }
}

MultiFabCellIterator
MultiFabCellIterator::begin() const
{
  MultiFabCellIterator ret(*this);
  ret.set_cell_index_and_boxID(begin_index);
  return ret;
}

MultiFabCellIterator
MultiFabCellIterator::end() const
{
  MultiFabCellIterator ret(*this);
  ret.set_cell_index_and_boxID(end_index);
  return ret;
}

void
MultiFabCellIterator::operator++()
{
  current_index++;
  set_cell_index_and_boxID(current_index);
  if (current_index>end_index) {
    BoxLib::Abort("Stepped beyond end of iterator");
  }
}

std::ostream& operator<<(std::ostream& os, const MultiFabCellIterator& it) {
  os << "current_box_id, current_cell_index: " << it.idx() << ", " << it.index();
  return os;
}

