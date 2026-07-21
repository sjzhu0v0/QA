// Build this source instead of the unit-weight executable to use the
// standalone inverse-variance fitter.  The original source stays unchanged.
#include "MFit.h"

// MFit.h has already been included above, so the include guard prevents the
// following source include from renaming the original MFitterPoly definition.
#define MFitterPoly MFitterPolyInvSigma2
#include "AssoYieldFit_noScale_smeared_2template_abs_ptV2_123_fitRange.cpp"
#undef MFitterPoly
