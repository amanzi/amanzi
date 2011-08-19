/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_EXAMPLE_HH_
#define AMANZI_CHEMISTRY_EXAMPLE_HH_

//
// Simple class template
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class Example {
 public:
  Example();
  virtual ~Example();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_EXAMPLE_HH_
