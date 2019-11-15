#include "Epetra_SerialDenseMatrix.h"

#include "errors.hh"
#include "exceptions.hh"
#include "UniqueHelpers.hh"
#include "HDF5Reader.hh"
#include "FunctionBilinearAndTime.hh"

namespace Amanzi {

FunctionBilinearAndTime::FunctionBilinearAndTime(const std::string& filename,
        const std::string& time_header,
        const std::string& x_header,
        const std::string& y_header,
        const std::string& val_header)
    : filename_(filename),
      x_header_(x_header),
      y_header_(y_header),
      val_header_(val_header),
      t_before_(-1.0),
      t_after_(-1.0),
      current_interval_(-2)
{
  HDF5Reader reader(filename_);
  reader.ReadData(time_header, times_);
}



FunctionBilinearAndTime::FunctionBilinearAndTime(const FunctionBilinearAndTime& other)
    : t_before_(other.t_before_),
      t_after_(other.t_after_),
      current_interval_(other.current_interval_),
      times_(other.times_),
      filename_(other.filename_),
      val_before_(std::unique_ptr<FunctionBilinear>(other.val_before_->Clone())),
      val_after_(std::unique_ptr<FunctionBilinear>(other.val_after_->Clone()))
{}


double FunctionBilinearAndTime::operator()(const std::vector<double>& x) const
{
  // get the interval of the current time
  int interval = std::lower_bound(times_.begin(), times_.end(), x[0]) - times_.begin() - 1;
  if (interval != current_interval_) {
    if ((current_interval_ == -2) || (interval != current_interval_+1)) {
      // completely uninitialized, or not an increment, load both
      if (interval == -1) {
        val_before_ = nullptr;
        val_after_ = Load_(0);
      } else if (interval == times_.size() - 1) {
        val_before_ = Load_(times_.size()-1);
        val_after_ = nullptr;
      } else {
        val_before_ = Load_(interval);
        val_after_ = Load_(interval + 1);
      }

    } else {
      // an increment, we can move after --> before
      val_before_ = std::move(val_after_);
      val_after_ = Load_(interval + 1);
    }
    current_interval_ = interval;
  }

  // interpolate
  if (interval == -1) {
    return (*val_after_)(x);
  } else if (interval == times_.size()-1) {
    return (*val_before_)(x);
  } else {
    double s = (x[0] - times_[interval]) / (times_[interval+1] - times_[interval]);
    return (*val_after_)(x) * s + (*val_before_)(x) * (1.0-s);
  }
}

std::unique_ptr<FunctionBilinear>
FunctionBilinearAndTime::Load_(const int time_index) const {
  HDF5Reader reader(filename_);

  std::vector<double> x, y;
  reader.ReadData(x_header_, x);
  reader.ReadData(y_header_, y);

  Epetra_SerialDenseMatrix values;
  reader.ReadMatData(val_header_+"/"+std::to_string(time_index), values);
  return std::make_unique<FunctionBilinear>(x, y, values, 1, 2);
}

} // namespace Amanzi
