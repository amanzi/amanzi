
//
// $Id: PMBld.cpp,v 1.2 2011-08-10 17:27:52 lijewski Exp $
//

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
