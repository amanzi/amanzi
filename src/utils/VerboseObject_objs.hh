#include "VerboseObject.hh"

// The default global verbosity level.
Teuchos::EVerbosityLevel Amanzi::VerboseObject::global_default_level = Teuchos::VERB_MEDIUM;

// Show or hide line prefixes
bool Amanzi::VerboseObject::global_hide_line_prefix = false;

// Size of the left column of names.
unsigned int Amanzi::VerboseObject::global_line_prefix_size = 18;

// rank on which to write
unsigned int Amanzi::VerboseObject::global_writing_rank = 0;
