/*************************************************
*                                                *
* Project: SpaceFuckery                          *
*                                                *
*   Copyright 2015 - 2016 - Marc-Olivier Barre   *
*              All rights reserved               *
*                                                *
**************************************************/

#include "StateVectors.h"

namespace mKOST
{
  StateVectors::StateVectors()
    : pos (0.0, 0.0, 0.0),
      vel (0.0, 0.0, 0.0),
      MeL (SIMD_INFINITY)
  {}

  StateVectors::~StateVectors()
  {}

  std::ostream &operator << (std::ostream &o, const StateVectors &v)
  {
    o     <<
    v.pos.getX() << "," <<
    v.pos.getY() << "," <<
    v.pos.getZ() << "," <<
    v.vel.getX() << "," <<
    v.vel.getY() << "," <<
    v.vel.getZ() << "," <<
    v.MeL;
    return o;
  }
}