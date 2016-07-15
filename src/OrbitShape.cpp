/*************************************************
*                                                *
* Project: SpaceFuckery                          *
*                                                *
*   Copyright 2015 - 2016 - Marc-Olivier Barre   *
*              All rights reserved               *
*                                                *
**************************************************/

#include "OrbitShape.h"

namespace mKOST
{
  OrbitShape::OrbitShape()
    : pe (0.0, 0.0, 0.0),
      ap (0.0, 0.0, 0.0),
      dn (0.0, 0.0, 0.0),
      an (0.0, 0.0, 0.0),
      points (nullptr),
      numPoints (0)
  {
  }

  OrbitShape::~OrbitShape()
  {
  }

  std::ostream &operator << (std::ostream &o, const OrbitShape &v)
  {
    o     <<
    v.pe.getX() << "," <<
    v.pe.getY() << "," <<
    v.pe.getZ() << "," <<
    v.ap.getX() << "," <<
    v.ap.getY() << "," <<
    v.ap.getZ() << "," <<
    v.dn.getX() << "," <<
    v.dn.getY() << "," <<
    v.dn.getZ() << "," <<
    v.an.getX() << "," <<
    v.an.getY() << "," <<
    v.an.getZ() << "," <<
    v.points << "," <<
    v.numPoints;
    return o;
  }
}