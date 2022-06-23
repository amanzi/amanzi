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

  Teuchos::RCP<Vec> Create() {
    assert(size_ >= 0);
    return Teuchos::rcp(new Vec(size_));
  }

private:
  int size_;
};

bool UserInitialize(Teuchos::ParameterList& plist, Vec& t,
                    const Amanzi::Key& fieldname,
                    const std::vector<std::string>* subfieldnames) {
  std::cout << "found it!" << std::endl;
  return true;
}

void UserWriteVis(const Amanzi::Visualization& vis,
                  const Amanzi::Key& fieldname,
                  const std::vector<std::string>* subfieldnames,
                  const Vec& vec) {}

void UserWriteCheckpoint(const Amanzi::Checkpoint& chkp,
                         const Amanzi::Key& fieldname,
                         const std::vector<std::string>* subfieldnames,
                         const Vec& vec) {}
void UserReadCheckpoint(const Amanzi::Checkpoint& chkp,
                        const Amanzi::Key& fieldname,
                        const std::vector<std::string>* subfieldnames,
                        Vec& vec) {}
