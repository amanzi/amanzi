/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

int g_constructor_calls_default = 0;
int g_constructor_calls_copy = 0;
int g_constructor_calls_main = 0;

struct Vec {
  Vec() { g_constructor_calls_default++; }

  Vec(const Vec& other) : v(other.v) { g_constructor_calls_copy++; }

  Vec(int size) : v(size, 0.0) { g_constructor_calls_main++; }

  std::vector<double> v;
};

class VecFactory {
 public:
  VecFactory() : size_(-1) {}

  void set_size(int size) { size_ = size; }

  Teuchos::RCP<Vec> Create()
  {
    assert(size_ >= 0);
    return Teuchos::rcp(new Vec(size_));
  }

 private:
  int size_;
};

bool
UserInitialize(Teuchos::ParameterList& plist,
               const Teuchos::ParameterList& attrs, Vec& t)
{
  std::cout << "Successfully initialized a Vec!" << std::endl;
  return true;
}

void
UserWriteVis(const Amanzi::Visualization& vis,
             const Teuchos::ParameterList& attrs, const Vec& vec)
{}

void
UserWriteCheckpoint(const Amanzi::Checkpoint& chkp,
                    const Teuchos::ParameterList& attrs, const Vec& vec)
{}
void
UserReadCheckpoint(const Amanzi::Checkpoint& chkp,
                   const Teuchos::ParameterList& attrs, Vec& vec)
{}
