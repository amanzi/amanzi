#ifndef MultiFabCellIterator_H
#define MultiFabCellIterator_H

#include <MultiFab.H>

class MultiFabCellIterator
{
public:
  MultiFabCellIterator(const DistributionMapping& dm,
		       const BoxArray&            ba,
		       int                        ng);

  MultiFabCellIterator(const MultiFabCellIterator& rhs);

  MultiFabCellIterator begin() const;
  MultiFabCellIterator end() const;
  virtual void operator++();

  bool operator<(const MultiFabCellIterator& rhs) const {
    return current_index < rhs.idx();
  }
  bool operator!=(const MultiFabCellIterator& rhs) const {
    return current_index != rhs.idx();
  }

  int boxID() const {return index_map[current_box_id];}
  const IntVect& index() const {return current_cell_index;}
  const Box& box() const {return current_box;}
  int nGrow() const {return num_grow;}
  const BoxArray& boxArray() const {return grids;}
  long idx() const {return current_index;}

protected:
  void initialize();
  void set_cell_index_and_boxID(long idx);
  static IntVect idx_to_intvect(const Box& box, long idx);

  DistributionMapping dmap;
  BoxArray grids;
  int num_grow;
  int begin_index, end_index;
  long current_index;
  int current_box_id;
  IntVect current_cell_index;
  Box current_box;
  std::vector<int> index_map;
  std::vector<long> offsets;
  long set_size;
};

std::ostream& operator<<(std::ostream& os, const MultiFabCellIterator& it);

#endif // MultiFabCellIterator_H
