// a test PK that stores its data in a struct and implements Advance.  This is
// like what Veg PFTs might do.


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "State.hh"

using namespace Amanzi;

struct PFT {
  PFT() : Bleaf(0.), Broot(0.), Bstem(0.) {}
  PFT(double leaf, double root, double stem) : Bleaf(leaf), Broot(root), Bstem(stem) {}
  double Bleaf;
  double Broot;
  double Bstem;
};


// the object to be stored/computed/etc
using PFTList = std::vector<PFT>;

// the factory 
class PFTListSpace {
 public:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;
  Teuchos::RCP<PFTList> Create() {
    return Teuchos::rcp(new PFTList(mesh->num_entities(AmanziMesh::CELL,
            AmanziMesh::OWNED)));
  }
};


bool
inline
UserInitialize(Teuchos::ParameterList& plist, PFTList& t,
               const Amanzi::Key& fieldname,
               const std::vector<std::string>& subfieldnames) { return true; }

void
UserWriteVis(const Amanzi::Visualization& vis, const Amanzi::Key& fieldname,
             const std::vector<std::string>& subfieldnames,
             const PFTList& vec) {}

void
UserWriteCheckpoint(const Amanzi::Checkpoint& chkp, const Amanzi::Key& fieldname,
             const PFTList& vec) {}
void
UserReadCheckpoint(const Amanzi::Checkpoint& chkp, const Amanzi::Key& fieldname,
                   PFTList& vec) {}



class Grow {
 public:
  Grow(double C) : C_(C) {}
  PFT operator()(const PFT& in) {
    return PFT{in.Bleaf + C_, in.Broot + C_, in.Bstem + C_};
  }

 private:
  double C_;
};

// the pk
template<class Base_t>
class PK_Veg : public Base_t {
 public:

  PK_Veg(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
         const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
         const Teuchos::RCP<State>& S) :
      Base_t(pk_tree, global_plist, S) {}

  void Setup() {
    Base_t::Setup();

    tag_new_ = "next";
    tag_old_ = "";

    // require at new and old times
    auto& fac = this->S_->template Require<PFTList, PFTListSpace>(this->key_, tag_old_, this->key_);
    fac.mesh = this->mesh_;

    // require data at the new and old times
    this->SolutionToState(tag_new_, "");
    this->SolutionToState(tag_old_, "");
  }

  bool 
  AdvanceStep(const Key& tag_old, const Key& tag_new) {
    // my local tags, used in physics PKs?  Can we get rid of these? --etc
    tag_old_ = tag_old;
    tag_new_ = tag_new;
  
    // times associated with those tags
    double t_new = S_->time(tag_new);
    double t_old = S_->time(tag_old);
    double dt = t_new - t_old;  
  
    // logging
    Teuchos::OSTab out = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "----------------------------------------------------------------" << std::endl
                 << "GROWING! t0 = " << t_old
                 << " t1 = " << t_new << " h = " << dt << std::endl
                 << "----------------------------------------------------------------" << std::endl;

    // manually take the step
    // grow a bit
    auto& data_new = this->S_->template GetW<PFTList>(this->key_, tag_new_, this->key_);
    const auto& data_old = this->S_->template Get<PFTList>(this->key_, tag_old_);

    Grow g(0.1 * dt);
    std::transform(data_old.begin(), data_old.end(), data_new.begin(), g);
    return false;
  }

  // these PKs typically assume a fixed, specified timestep
  double get_dt() {
    return 0.1;
  }


  void
  StateToSolution(TreeVector& soln, const Key& tag, const Key& suffix)
  {
    // ASSERTS for now?
    ASSERT(false);
  }

  // THOU SHALL NOT PASS!
  void FailStep(const Key& tag_old, const Key& tag_new) {
    Errors::Message message("THOU SHALL NOT PASS!");
    throw(message);
  }
  
  
 protected:
  std::string tag_new_, tag_old_;
  using Base_t::S_;
  using Base_t::vo_;

};
