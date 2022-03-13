/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
   Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
   Amanzi is released under the three-clause BSD License.
   The terms of use and "as is" disclaimer for this license are
   provided in the top-level COPYRIGHT file.
   See $ATS_DIR/COPYRIGHT

   Author: Ethan Coon
*/
//! Simple utility for doing intersections across a communicator.

/*

This could be made much much more general, but for now it does two jobs:

1. It performs reductions on strings, forming the "max" across a communicator
(where max means alphabetically last according to std::string::operator> ).

2. It forms the intersection of a sorted list of strings across a communicator,
finding all strings that are on all ranks.

This gets used by domain sets in visualization, to make sure that IO is done
only for vectors that live on all ranks (otherwise the IO method hangs).

*/

#pragma once

#include <vector>
#include <string>
#include <set>

#include "dbc.hh"
#include "AmanziTypes.hh"

namespace Amanzi {
namespace Utils {

namespace Impl {

#ifdef HAVE_MPI

template<int MAXLEN>
void MyMaxString(char* in, char* out, int *count, MPI_Datatype *dptr) {
  AMANZI_ASSERT(*count == 1);
  for (int i=0; i!=*count; ++i) {
    std::string instr(&in[MAXLEN*i], MAXLEN);
    std::string outstr(&out[MAXLEN*i], MAXLEN);
    if (instr > outstr) {
      for (int j=0; j!=MAXLEN; ++j) out[MAXLEN*i+j] = in[MAXLEN*i+j];
    }
  }
}

template<int MAXLEN>
void MyMaxStringWrapper(void* in, void *out, int *count, MPI_Datatype *dptr) {
  MyMaxString<MAXLEN>((char*)in, (char*)out, count, dptr);
}

#endif

} // namespace Impl


//
// A little helper class that, along with the Impl namespace, allows the
// determination of the alphabetically last string across a communicator.
//
template <int maxlen>
class StringReducer {
  static const int MAXLEN = maxlen;

 public:
  StringReducer(const Comm_ptr_type& comm)
    : comm_ptr_(comm) {

#ifdef HAVE_MPI
    // get the raw comm
    auto comm_raw = Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm);
    AMANZI_ASSERT(comm_raw);
    comm_ = comm_raw->Comm();

    // construct the data type
    MPI_Type_contiguous(MAXLEN, MPI_CHAR, &dtype_);
    MPI_Type_commit(&dtype_);

    // construct the op
    MPI_Op_create(Impl::MyMaxStringWrapper<MAXLEN>, 1, &op_);

    int ierr = MPI_Comm_rank(comm_, &rank_);
    ierr |= MPI_Comm_size(comm_, &size_);
    AMANZI_ASSERT(!ierr);
#endif
  }

  // reduceAll on strings -- returns the alphabetically last string across comm
  std::string reduceMaxString(const std::string& str) const {
#ifdef HAVE_MPI
    // buffer strings
    std::string string_g(MAXLEN, ' ');
    std::string string_l(str);
    string_l.insert(str.size(), MAXLEN-1, ' ');

    // reduce
    int ierr = MPI_Allreduce((void*) string_l.data(), (void*) string_g.data(), 1, dtype_, op_, comm_);
    AMANZI_ASSERT(!ierr);

    // erase the trailing spaces
    int i;
    for (i=MAXLEN-1; i!=-1; --i) {
      if (string_g[i] != ' ') break;
    }
    return string_g.substr(0,i+1);
#else
    return str;
#endif
  }


  int reduceMinInt(int size) const {
#ifdef HAVE_MPI
    int size_g = 0;
    int ierr = MPI_Allreduce(&size, &size_g, 1, MPI_INT, MPI_MIN, comm_);
    AMANZI_ASSERT(!ierr);
    return size_g;
#else
    return size;
#endif
  }


  // intersectAll -- given a SORTED list of strings, returns the SORTED list of
  // keys that exist on all ranks of comm.
  //
  // NOTE: if uncertain, call checkValidInput()!
  std::vector<std::string> intersectAll(const std::vector<std::string>& in) const {

#ifdef HAVE_MPI
    std::vector<std::string> out;
    int nleft = reduceMinInt(in.size());
    int i = 0;
    int keep = 0;
    while (nleft > 0) {
      std::string l_next = in[i];
      auto g_next = reduceMaxString(l_next);

      if (g_next == l_next) {
        // l_next is the largest, maybe all one that all have?
        keep = 1;
      } else {
        // l_next is not the largest, so some other process does not have l_next
        keep = 0;
        i++;
      }
      if (reduceMinInt(keep) == 1) {
        // all have l_next
        out.emplace_back(l_next);
        i++;
      }

      nleft = reduceMinInt(in.size() - i);
    }
    return out;
#else
    return in;
#endif
  }


  void checkValidInput(const std::vector<std::string>& in) const {
    std::set<std::string> set(in.begin(), in.end());
    AMANZI_ASSERT(set.size() == in.size());
    if (in.size() > 1) {
      for (int i=0; i!=(in.size()-1); ++i) {
        AMANZI_ASSERT(in[i+1] > in[i]);
      }
    }
  }


 private:
#ifdef HAVE_MPI
  MPI_Datatype dtype_;
  MPI_Op op_;
  MPI_Comm comm_;
  int rank_, size_;
#endif
  Comm_ptr_type comm_ptr_;
};

} // namespace Utils
} // namespace Amanzi
