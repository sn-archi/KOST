/*************************************************
*                                                *
* Project: SpaceFuckery                          *
*                                                *
*   Copyright 2015 - 2016 - Marc-Olivier Barre   *
*              All rights reserved               *
*                                                *
**************************************************/

#include "Elements.h"

namespace mKOST
{
  Elements::Elements()
  : a   (SIMD_INFINITY),
    Ecc (SIMD_INFINITY),
    i   (SIMD_INFINITY),
    LAN (SIMD_INFINITY),
    LoP (SIMD_INFINITY)
  {}

  Elements::~Elements()
  {
  }

  std::ostream &operator << (std::ostream &o, const Elements &v)
  {
    o     <<
    v.a   << "," <<
    v.Ecc << "," <<
    v.i   << "," <<
    v.LAN << "," <<
    v.LoP;
    return o;
  }
}