/***************************************************************************
 *   Copyright (C) 2008 by C J Plooy                                       *
 *   cornwarecjp@lavabit.com                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2.1 of the  *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU Lesser General Public License for more details.                   *
 *                                                                         *
 *   You should have received a copy of the GNU Lesser General Public      *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef TYPES_H
#define TYPES_H

#include "LinearMath/btScalar.h"
#include "LinearMath/btVector3.h"
#include "LinearMath/btMatrix3x3.h"

namespace mKOST
{
  ATTRIBUTE_ALIGNED16(struct) Elements
  {
    btScalar a;    /**< Semi-major axis */
    btScalar Ecc;  /**< Eccentricity */
    btScalar i;    /**< Inclination */
    btScalar LAN;  /**< Longitude of ascending node */
    btScalar LoP;  /**< Longitude of periapsis */
    btScalar L;    /**< Mean longitude at epoch */

    BT_DECLARE_ALIGNED_ALLOCATOR();

    Elements()
    : a (SIMD_INFINITY),
    Ecc (SIMD_INFINITY),
    i (SIMD_INFINITY),
    LAN (SIMD_INFINITY),
    LoP (SIMD_INFINITY),
    L (SIMD_INFINITY)
    {}
  };

  ATTRIBUTE_ALIGNED16(struct) Params
  {
    btScalar SMi;  /**< semi-minor axis */
    btScalar PeD;  /**< periapsis distance */
    btScalar ApD;  /**< apoapsis distance */
    btScalar MnA;  /**< mean anomaly */
    btScalar TrA;  /**< true anomaly */
    btScalar MnL;  /**< mean longitude */
    btScalar TrL;  /**< true longitude */
    btScalar EcA;  /**< eccentric anomaly */
    btScalar Lec;  /**< linear eccentricity */
    btScalar T;    /**< orbit period */
    btScalar PeT;  /**< time to next periapsis passage */
    btScalar ApT;  /**< time to next apoapsis passage */
    btScalar AgP;  /**< argument of periapsis */

    BT_DECLARE_ALIGNED_ALLOCATOR();

    Params()
    : SMi (SIMD_INFINITY),
    PeD (SIMD_INFINITY),
    ApD (SIMD_INFINITY),
    MnA (SIMD_INFINITY),
    TrA (SIMD_INFINITY),
    MnL (SIMD_INFINITY),
    TrL (SIMD_INFINITY),
    EcA (SIMD_INFINITY),
    Lec (SIMD_INFINITY),
    T (SIMD_INFINITY),
    PeT (SIMD_INFINITY),
    ApT (SIMD_INFINITY),
    AgP (SIMD_INFINITY)
    {}
  };

  ATTRIBUTE_ALIGNED16(struct) StateVectors
  {
    btVector3 pos;  /**< Position */
    btVector3 vel;  /**< Velocity */
    btScalar MeL;   /**< Mean longitude at epoch */

    BT_DECLARE_ALIGNED_ALLOCATOR();

    StateVectors()
    : pos (0.0, 0.0, 0.0),
    vel (0.0, 0.0, 0.0),
    MeL (SIMD_INFINITY)
    {}
  };

  ATTRIBUTE_ALIGNED16(struct) OrbitShape
  {
    btVector3 pe;            /**< Periaspsis position */
    btVector3 ap;            /**< Apoapsis position */
    btVector3 dn;            /**< Descending node position */
    btVector3 an;            /**< Ascending node position */
    btVector3* points;       /**< points */
    unsigned int numPoints;  /**< numPoints */

    BT_DECLARE_ALIGNED_ALLOCATOR();

    OrbitShape()
    : pe (0.0, 0.0, 0.0),
    ap (0.0, 0.0, 0.0),
    dn (0.0, 0.0, 0.0),
    an (0.0, 0.0, 0.0),
    points (nullptr),
    numPoints (0)
    {}
  };
}
#endif // TYPES_H