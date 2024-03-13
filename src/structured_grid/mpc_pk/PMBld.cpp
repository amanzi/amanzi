/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <winstd.H>

#include <PorousMedia.H>

// --------------------------------------------------------------------
// -----   PMBld class instantiation
// --------------------------------------------------------------------

PMBld nsbld;

LevelBld* getLevelBld()
{
    return &nsbld;
}

// --------------------------------------------------------------------
// -----   PMBld class implementation
// --------------------------------------------------------------------

void
PMBld::variableSetUp()
{
    PorousMedia::variableSetUp();
}

void
PMBld::variableCleanUp()
{
    PorousMedia::variableCleanUp();
}

AmrLevel*
PMBld::operator()()
{
    return new PorousMedia;
}

AmrLevel*
PMBld::operator()(Amr &papa, int lev, const Geometry &level_geom,
                  const BoxArray &ba, Real time)
{
    return new PorousMedia(papa, lev, level_geom, ba, time);
}
