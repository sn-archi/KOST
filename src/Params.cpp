/*************************************************
*                                                *
* Project: SpaceFuckery                          *
*                                                *
*   Copyright 2015 - 2016 - Marc-Olivier Barre   *
*              All rights reserved               *
*                                                *
**************************************************/

#include "Params.h"

namespace mKOST
{
  Params::Params()
    : SMi (SIMD_INFINITY),
      PeD (SIMD_INFINITY),
      ApD (SIMD_INFINITY),
      MnA (SIMD_INFINITY),
      TrA (SIMD_INFINITY),
      TrL (SIMD_INFINITY),
      EcA (SIMD_INFINITY),
      Lec (SIMD_INFINITY),
      T   (SIMD_INFINITY),
      PeT (SIMD_INFINITY),
      ApT (SIMD_INFINITY),
      AgP (SIMD_INFINITY)
  {}

  Params::~Params()
  {}

  std::ostream &operator << (std::ostream &o, const Params &v)
  {
    o     <<
    v.PeD << "," <<
    v.ApD << "," <<
    v.MnA << "," <<
    v.TrA << "," <<
    v.EcA << "," <<
    v.T   << "," <<
    v.AgP << "," <<
    v.TrL << "," <<
    v.Lec << "," <<
    v.PeT << "," <<
    v.ApT << "," <<
    v.SMi;
    return o;
  }
}