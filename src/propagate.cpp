/***************************************************************************
 *   Copyright (C) 2008 by C J Plooy / Tim Blaxland                        *
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

/*
Changelog:

CJP       = C J Plooy
TBlaxland = Tim Blaxland

Format:
-----------------
date (dd-mm-yyyy)
author
changes

-----------------
27-02-09
CJP
Merged kostGetElementsAtTime and kostGetElementsAtTime2
Renamed these files to kost_propagate.*

-----------------
03-02-09
TBlaxland
Finished kostGetElementsAtTime
Added kostGetElementsAtTime2
Modified kostElements2StateVectorAtTime to use kostGetElementsAtTime2

-----------------
01-02-09
TBlaxland
Fixed bugs in kostGetMeanAnomalyAtTime

-----------------
24-01-2009
CJP
Initial version
Moved functions from kost_elements.c to here
Renamed several functions
*/

/***************************************************************************
 * Usage notes - kostElements2StateVectorX and related functions:
 *
 * Parabolic orbits are NOT currently supported.
 *
 * Position as a function of time can be found by either:
 *
 * 1. A call to kostElements2StateVector1. Depending on settings for
 *    maxIterations and maxRelativeError, this may adversly affect frame
 *    rates in graphical applications.
 *
 * 2. To minimise impact on frame rates:
 *    2.1. Call kostGetMeanAnomaly
 *    2.2. Call kostGetEccentricAnomaly on successive time steps with a
 *         small number of iterations in each call. Each call takes the
 *         result of the previous call as its eccentricAnomalyEstimate.
 *         Repeat until kostGetEccentricAnomaly returns > 0.
 *    2.3. Call kostGetTrueAnomaly2.
 *    2.4. Call kostElements2StateVector2.
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "elements.h"

#include "propagate.h"

namespace mKOST
{
  btScalar getMeanAnomalyAtTime (
    btScalar mu,                  /* standard gravitational parameter */
    const sElements* elements, /* pointer to orbital elements at epoch */
    btScalar timeSinceEpoch)      /* time since epoch in seconds */
  {
    /* Pseudocode
     *
     * calc mean motion
     * calc change in mean anomaly
     * calc mean anomaly */

    btScalar meanMotion, deltaMeanAnomaly, meanAnomaly;

    /* calc mean motion */
    meanMotion = std::sqrt (mu / std::pow (std::fabs (elements->a), 3.0) );

    /* calc change in mean anomaly */
    deltaMeanAnomaly = timeSinceEpoch * meanMotion;

    /* calc mean anomaly */
    /* fmod takes care of overflow in the case where the mean anomaly exceeds one revolution */
    meanAnomaly = std::fmod (elements->L - elements->omegab + deltaMeanAnomaly, SIMD_2_PI);

    if (meanAnomaly < 0) meanAnomaly += SIMD_2_PI;

    return meanAnomaly;
  }

  int getTrueAnomalyAtTime (
    btScalar mu,                  /* standard gravitational parameter */
    const sElements* elements, /* pointer to orbital elements at epoch */
    btScalar* trueAnomaly,        /* location where result will be stored */
    btScalar timeSinceEpoch,      /* time since epoch in seconds */
    btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
    int maxIterations)            /* max number of iterations for calculating eccentric anomaly */
  {
    /* Returns number of iterations if successful, returns 0 otherwise. */

    /* Pseudocode
     *
     * get mean anomaly
     * get eccentric anomaly
     * calc true anomaly */

    int ret;
    btScalar meanAnomaly, eccentricAnomaly;

    /* get mean anomaly */
    meanAnomaly = getMeanAnomalyAtTime (mu, elements, timeSinceEpoch);

    /* get eccentric anomaly */
    ret = getEccentricAnomaly (elements, &eccentricAnomaly, meanAnomaly, meanAnomaly, maxRelativeError, maxIterations);

    /* calc true anomaly */
    *trueAnomaly = getTrueAnomaly2 (mu, elements, eccentricAnomaly);

    return ret;
  }

  btScalar getLANAtTime (
    btScalar mu,                  /* standard gravitational parameter */
    const sElements* elements, /* pointer to orbital elements at epoch */
    btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
    btScalar jTwoCoeff,           /* J2 coefficient of the non-spherical body being orbited */
    btScalar timeSinceEpoch)      /* time since epoch in seconds */
  {
    btScalar meanMotion;

    if (elements->e < 1.0) /* elliptical orbit */
      {
        meanMotion = std::sqrt (mu / pow (elements->a, 3.0) );
        return ( elements->theta + timeSinceEpoch * (-3.0 * meanMotion / 2.0) * pow (bodyRadius / elements->a, 2.0) * (cos (elements->i) / pow (1.0 - pow (elements->e, 2.0), 2.0) ) * jTwoCoeff );
      }
    else /* hyperbolic orbit - non spherical effect is negligible */
      return elements->theta;
  }

  btScalar getArgPeAtTime (
    btScalar mu,                  /* standard gravitational parameter */
    const sElements* elements, /* pointer to orbital elements at epoch */
    btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
    btScalar jTwoCoeff,           /* J2 coefficient of the non-spherical body being orbited */
    btScalar timeSinceEpoch)      /* time since epoch in seconds */
  {
    btScalar meanMotion;

    if (elements->e < 1.0) /* elliptical orbit */
      {
        meanMotion = std::sqrt (mu / std::pow (elements->a, 3.0) );
        return ( elements->omegab - elements->theta + timeSinceEpoch * (3.0 * meanMotion / 4.0) * std::pow (bodyRadius / elements->a, 2.0) * ( (5.0 * std::pow (std::cos (elements->i),
                 2.0) - 1.0) / std::pow (1.0 - std::pow (elements->e, 2.0), 2.0) ) * jTwoCoeff );
      }
    else /* hyperbolic orbit - non spherical effect is negligible */
      return ( elements->omegab - elements->theta );
  }

  int elements2StateVectorAtTime (
    btScalar mu,                  /* standard gravitational parameter */
    const sElements* elements, /* pointer to orbital elements at epoch */
    sStateVector* state,       /* pointer to location where state vectors at epoch+timeSinceEpoch will be stored */
    btScalar timeSinceEpoch,      /* time since epoch in seconds */
    btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
    int maxIterations,            /* max number of iterations for calculating eccentric anomaly */
    btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
    btScalar jTwoCoeff)           /* J2 coefficient of the non-spherical body being orbited */
  {
    /* Returns number of iterations if successful, returns 0 otherwise. */

    /* Pseudocode
     *
     * get true anomaly
     * get longitude of ascending node and argument of periapsis
     * calc state vectors */

    int ret;
    btScalar trueAnomaly;
    sElements updatedElements;

    /* get true anomaly */
    ret = getTrueAnomalyAtTime (mu, elements, &trueAnomaly, timeSinceEpoch, maxRelativeError, maxIterations);

    /* update elements for new epoch */
    getElementsAtTime (mu, elements, &updatedElements, timeSinceEpoch, bodyRadius, jTwoCoeff);

    /* calc state vectors */
    elements2StateVector2 (mu, &updatedElements, state, trueAnomaly);

    return ret;
  }

  void getElementsAtTime (
    btScalar mu,                     /* standard gravitational parameter */
    const sElements* elements,    /* pointer to orbital elements at epoch */
    sElements* newElements,       /* pointer to location where elements at epoch+timeSinceEpoch will be stored */
    btScalar timeSinceEpoch,         /* time since epoch in seconds */
    btScalar bodyRadius,             /* mean radius of the non-spherical body being orbited */
    btScalar jTwoCoeff)              /* J2 coefficient of the non-spherical body being orbited */
  {
    *newElements = *elements;

    /* Mean longitude: */
    newElements->L = getMeanAnomalyAtTime (mu, newElements, timeSinceEpoch) + newElements->omegab;

    if (bodyRadius > SIMD_EPSILON)
      {
        /* longitude of ascending node */
        newElements->theta =
          getLANAtTime (mu, newElements, bodyRadius, jTwoCoeff, timeSinceEpoch);

        /* argument of periapsis */
        newElements->omegab =
          newElements->theta +
          getArgPeAtTime (mu, newElements, bodyRadius, jTwoCoeff, timeSinceEpoch);
      }
  }
}