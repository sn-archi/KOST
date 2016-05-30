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
17-04-2010
CJP
Fixed the size of a small addition to the velocity
vector for radial orbits, so that extreme cases are
better supported.
Various fixes to make these functions pass the
accuracy test

-----------------
19-03-2010
CJP
Applied some of the patches as suggested by
TBlaxland to fix a bug reported by Pagnatious
and to avoid having an undefined velocity in
kostElements2StateVector2.

Made some other changes to how velocity calculation
is done in kostElements2StateVector2.

Changed kostStateVector2Elements, so that it approximates
parabolic orbits by ellipses in some cases.

Changed the output range of kostGetMeanAnomaly, the meanAnomaly
input range of kostGetEccentricAnomaly, and the output range
of kostGetEccentricAnomaly, to fix a bug for hyperbolic
orbits.

Added support for radial orbits

-----------------
01-03-2009
CJP
Replaced while loops with fmod, as suggested by TBlaxland

-----------------
01-02-2009
TBlaxland
Fixed bugs in kostGetMeanAnomaly and kostElements2StateVector2

-----------------
24-01-2009
CJP
Moved some functions to kost_time.c

-----------------
19-01-2009
TBlaxland
Fixed some minor bugs
Changed kostElements2StateVector2 from left-hand to right-hand coordinates

-----------------
17-01-2009
CJP
Changed some whitespace formatting
Made TBlaxland's functions ANSI compatible

-----------------
14-01-2009
TBlaxland
Added implementation of:
	kostGetMeanAnomaly
	kostGetEccentricAnomaly
	kostGetTrueAnomaly1
	kostGetTrueAnomaly2
	kostGetLAN
	kostGetArgPe
	kostElements2StateVector1
	kostElements2StateVector2
Fixed division by zero at e==1.0
-----------------
16-11-2008
CJP
Initial version
Added implementation of:
	kostStateVector2Elements
*/

/***************************************************************************
 * Usage notes - kostElements2StateVectorX and related functions:
 *
 * Parabolic orbits are NOT currently supported.
 *
 * Position as a function of time can be found by either:
 *
 * 1. A call to mKOST::elements2StateVector1. Depending on settings for
 *    maxIterations and maxRelativeError, this may adversly affect frame
 *    rates in graphical applications.
 *
 * 2. To minimise impact on frame rates:
 *    2.1. Call mKOST::getMeanAnomaly
 *    2.2. Call mKOST::getEccentricAnomaly on successive time steps with a
 *         small number of iterations in each call. Each call takes the
 *         result of the previous call as its eccentricAnomalyEstimate.
 *         Repeat until kostGetEccentricAnomaly returns > 0.
 *    2.3. Call mKOST::getTrueAnomaly2.
 *    2.4. Call mKOST::elements2StateVector2.
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "elements.h"

namespace mKOST
{
  btVector3 gete(btScalar mu, sStateVector state)
  {
    btVector3 h (state.pos.cross (state.vel));
    return (state.vel.cross(h) / mu - state.pos / state.pos.length());
  }

  btVector3 geth(sStateVector state)
  {
    return state.pos.cross (state.vel);
  }

  btVector3 getn(btVector3 h)
  {
    return btVector3 (0.0, 1.0, 0.0).cross (h);
  }

  btVector3 getn(btScalar LaN)
  {
    return btVector3 (std::cos (LaN), 0.0, -std::sin (LaN));
  }

  btScalar getLaN(btVector3 n)
  {
    if (n.length() < SIMD_EPSILON)
    {
      return 0.0;
    }
    else
    {
      btScalar LaN = std::acos (n.getX() / n.length());
      return (n.getZ() > 0.0)?SIMD_2_PI - LaN:LaN;
    }
  }

  btScalar getMeanAnomaly (
    btScalar mu,                   /* standard gravitational parameter */
    const sElements* elements)   /* pointer to orbital elements at epoch */
  {
    btScalar meanAnomaly;
    /* calc mean anomaly */
    meanAnomaly = elements->L - elements->LoP;

    if (elements->e < 1.0)
      {
        /* check range is in 0 to 2π */
        if (meanAnomaly < 0.0) meanAnomaly += SIMD_2_PI;
        if (meanAnomaly >  SIMD_2_PI) meanAnomaly -= SIMD_2_PI;
      }

    return meanAnomaly;
  }

  int getEccentricAnomaly (
    const sElements* elements,      /* pointer to orbital elements at epoch */
    btScalar* eccentricAnomaly,        /* location where result will be stored */
    btScalar meanAnomaly,              /* mean anomaly */
    btScalar eccentricAnomalyEstimate, /* initial estimate of eccentric anomaly, start with mean anomaly if no better estimate available */
    btScalar maxRelativeError,         /* maximum relative error in eccentric anomaly */
    int maxIterations)                 /* max number of iterations for calculating eccentric anomaly */
  {
    /* Code will terminate when either maxIterations or maxRelativeError is reached.
     * Returns number of iterations if relativeError < maxRelativeError, returns 0 otherwise. */

    /* Pseudocode
     *
     * do
     *  calculate next estimate of the root of Kepler's equation using Newton's method
     *  calculate estimate of mean anomaly from estimate of eccentric anomaly
     *  calculate relativeError
     * while ((iterations<=maxIterations)&&(relativeError<maxRelativeError))
     * if iterations<=maxIterations return iterations, else return 0 */

    int i (0);
    btScalar relativeError, meanAnomalyEstimate, e;

    if (elements->e == 0.0)   /* circular orbit */
      {
        *eccentricAnomaly = meanAnomaly;
        return 1;
      }

    if (elements->e == 1.0) /* parabolic orbit - approximate to hyperbolic */
      e = elements->e + SIMD_EPSILON;
    else
      e = elements->e;

    do
      {
        if (elements->e < 1.0)   /* elliptical orbit */
          {
            /* calculate next estimate of the root of Kepler's equation using Newton's method */
            eccentricAnomalyEstimate = eccentricAnomalyEstimate - (eccentricAnomalyEstimate - e * std::sin (eccentricAnomalyEstimate) - meanAnomaly) / (1.0 - e * std::cos (eccentricAnomalyEstimate) );
            /* calculate estimate of mean anomaly from estimate of eccentric anomaly */
            meanAnomalyEstimate = eccentricAnomalyEstimate - e * std::sin (eccentricAnomalyEstimate);
          }
        else   /* hyperbolic orbit */
          {
            /* calculate next estimate of the root of Kepler's equation using Newton's method */
            eccentricAnomalyEstimate = eccentricAnomalyEstimate - (e * std::sinh (eccentricAnomalyEstimate) - eccentricAnomalyEstimate - meanAnomaly) / (e * std::cosh (eccentricAnomalyEstimate) - 1.0);
            /* calculate estimate of mean anomaly from estimate of eccentric anomaly */
            meanAnomalyEstimate = e * std::sinh (eccentricAnomalyEstimate) - eccentricAnomalyEstimate;
          }
        /* calculate relativeError */
        relativeError = 1.0 - meanAnomalyEstimate / meanAnomaly;
        ++i;
      }
    while ( (i < maxIterations) && (fabs (relativeError) > fabs (maxRelativeError) ) );

    if (elements->e < 1.0)
      {
        /* check range is in 0 to 2π */
        eccentricAnomalyEstimate = fmod (eccentricAnomalyEstimate, SIMD_2_PI);
        if (eccentricAnomalyEstimate < 0.0) eccentricAnomalyEstimate += SIMD_2_PI;
        if (eccentricAnomalyEstimate > SIMD_2_PI) eccentricAnomalyEstimate -= SIMD_2_PI;
      }

    *eccentricAnomaly = eccentricAnomalyEstimate;

    if ( fabs (relativeError) < fabs (maxRelativeError) )
      return i;
    else
      return 0;
  }

  int getTrueAnomaly (
    btScalar mu,                  /* standard gravitational parameter */
    const sElements* elements,    /* pointer to orbital elements at epoch */
    btScalar* trueAnomaly,        /* location where result will be stored */
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
    meanAnomaly = getMeanAnomaly (mu, elements);

    /* get eccentric anomaly */
    if (elements->e < 1.0)
      {
        ret = getEccentricAnomaly (
                elements, &eccentricAnomaly,
                meanAnomaly, meanAnomaly,
                maxRelativeError, maxIterations);
      }
    else
      {
        ret = getEccentricAnomaly (
                elements, &eccentricAnomaly,
                meanAnomaly, std::log (2.0 * meanAnomaly / elements->e + 1.8),
                maxRelativeError, maxIterations);
      }

    /* calc true anomaly */
    *trueAnomaly = getTrueAnomaly2 (mu, elements, eccentricAnomaly);

    return ret;
  }

  btScalar getTrueAnomaly2 (
    btScalar mu,                  /* standard gravitational parameter */
    const sElements* elements, /* pointer to orbital elements at epoch */
    btScalar eccentricAnomaly)    /* eccentric anomaly */
  {
    btScalar ret;
    if (elements->e < 1.0)   /* elliptical orbit */
      {
        ret = 2.0 * std::atan (std::sqrt ( (1.0 + elements->e) / (1.0 - elements->e) ) * std::tan (eccentricAnomaly / 2.0) );
      }
    else   /* hyperbolic orbit */
      {
        ret = std::acos ( (std::cosh (eccentricAnomaly) - elements->e) / (1 - elements->e * std::cosh (eccentricAnomaly) ) );
        if (eccentricAnomaly < 0.0) ret = -ret; /* Always the same sign */
      }

    return ret;
  }

  btScalar getAgP (btVector3 e, btVector3 n, btVector3 h, bool isCircular, bool isEquatorial)
  {
    if (isCircular)
      {
        return 0.0;
      }
    else if (isEquatorial)
      {
        btScalar AgP (atan2 (e.getZ(), e.getX()));
        return (h.getY() < 0.0)?-AgP:AgP;
      }
    else
      {
        btScalar AgP (std::acos (n.dot (e) / (n.length() * e.length())));
        return (e.getY() < 0.0)?SIMD_2_PI - AgP:AgP;
      }
  }

  btScalar getTrAFromState(btVector3 pos, btVector3 vel, btVector3 n, btVector3 e, bool isCircular, bool isEquatorial, bool* sin_TrA_isNegative)
  {
    btScalar TrA (0.0);
    if (isCircular)
    {
      if (isEquatorial)
      {
        TrA = std::acos (pos.getX() / pos.length());
        if (vel.getX() > 0.0)
        {
          *sin_TrA_isNegative = true;
          TrA = SIMD_2_PI - TrA;
        }
      }
      else
      {
        TrA = std::acos (n.dot (pos) ) / (n.length() * pos.length());
        if (n.dot (vel) > 0.0)
        {
          *sin_TrA_isNegative = true;
          TrA = SIMD_2_PI - TrA;
        }
      }
    }
    else
    {
      btScalar tmp = e.dot (pos) / (e.length() * pos.length());

      /* Avoid acos out of range. 1 and -1 are included cause we know the result. */
      if (tmp <= -1.0)
      {
        TrA = SIMD_PI;
      }
      else if (tmp >= 1.0)
      {
        TrA = 0.0;
      }
      else
      {
        TrA = std::acos (tmp);
      }

      /* Changing sin_TrA_isNegative on a negative epsilon value doesn't seem right */
      if (pos.dot (vel) + SIMD_EPSILON < 0.0)
      {
        *sin_TrA_isNegative = true;
        TrA = SIMD_2_PI - TrA;
      }
    }
    return TrA;
  }

  void elements2StateVector2 (
    btScalar mu,                  /* standard gravitational parameter */
    const sElements* elements, /* pointer to orbital elements at epoch */
    sStateVector* state,       /* pointer to location where state vector at epoch will be stored */
    btScalar trueAnomaly)         /* true anomaly */
  {
    /* Pseudocode
     *
     * calc nodal vector, n
     * calc angular momentum vector, h
     * calc position vector
     * calc argument of position
     * calc direction of position vector
     * calc length of position vector
     * calc velocity vector
     * calculate magnitude of velocity vector
     * calc components of velocity vector perpendicular and parallel to radius vector
     * add to get velocity vector */

    /* eccentricity */
    btScalar e (elements->e);
    if (e == 1.0) /* parabolic orbit - approximate to hyperbolic orbit */
    {
      e += SIMD_EPSILON;
    }

    /* unit vector in direction of ascending node */
    btVector3 n (getn (elements->LaN));

    /* unit vector pointing ecliptic north */
    btVector3 north (btVector3 (0.0, 1.0, 0.0));

    /* calc angular momentum vector, h */
    /* projection of h in ecliptic (xz) plane */
    btVector3 h (n.cross (north));
    h *= std::sin(elements->i);
    /* elevation of h */
    h.setY (std::cos(elements->i));
    h.normalize();
    /* calc magnitude of h */

    /* calc radius and velocity at periapsis */
    btScalar rPe, vPe;
    if (e < 1.0)   /* elliptical orbit */
      {
        rPe = elements->a * (1.0 - e * e) / (1.0 + e);
        vPe = std::sqrt (mu * (2.0 / rPe - 1.0 / elements->a) );
      }
    else   /* hyperbolic orbit */
      {
        rPe = std::fabs (elements->a) * (e * e - 1.0) / (1.0 + e);
        vPe = std::sqrt (mu * (2.0 / rPe + 1.0 / std::fabs (elements->a) ) );
      }
    /* calc h */
    h *= rPe * vPe;

    /* calc position vector */
    /* argument of position, measured from the longitude of ascending node */
    btScalar argPos (elements->LoP - elements->LaN + trueAnomaly);

    /*
    calc direction of position vector:
    r/|r| = sin(ArgPos) * ((h / |h|) x n) + cos(argPos) * n
    */
    state->pos = std::sin(argPos) * h.normalized().cross(n) + std::cos(argPos) * n;
    //btVector3 tmpv (std::cos(argPos) * n);
    //state->pos = h.normalized();
    //state->pos = state->pos.cross (n);
    //state->pos *= std::sin(argPos);
    //state->pos += tmpv;

    /* calc length of position vector */
    if (e < 1.0)                  /* elliptical orbit */
      state->pos = elements->a * (1.0 - e * e) / (1.0 + e * std::cos(trueAnomaly)) * state->pos;
    else                          /* hyperbolic orbit */
      state->pos = std::fabs (elements->a) * (e * e - 1.0) / (1.0 + e * std::cos(trueAnomaly)) * state->pos;

    /* calc velocity vector */
    /* calculate squared magnitude of velocity vector */
    btScalar v2;
    if (e < 1.0)   /* elliptical orbit */
      {
        v2 = mu * (2.0 / state->pos.length() - 1.0 / elements->a);
      }
    else   /* hyperbolic orbit */
      {
        v2 = mu * (2.0 / state->pos.length() + 1.0 / std::fabs (elements->a));
      }

    /* calc components of velocity vector perpendicular and parallel to radius vector:

    perpendicular:
    vPro = (|h|/|pos|) * normal(h x pos)

    parallel:
    vO = √(v² - |vPro|²) * sign(sin(trueAnomaly)) * normal(pos)
    */
    btVector3 vPro (h.length() / state->pos.length() * h.cross (state->pos).normalized());
    //btVector3 vPro (h.cross (state->pos));
    //vPro.normalize();
    //vPro *= h.length() / state->pos.length();

    btScalar tmpr (std::sin(trueAnomaly));
    btVector3 vO;
    if (tmpr == 0.0)   /* check for apsis condition to avoid divide by zero */
      {
        vO = btVector3 (0.0, 0.0, 0.0);
      }
    else
      {
        btScalar signSinTrueAnomaly (tmpr / std::fabs (tmpr));

        btScalar v0_sq (v2 - vPro.length2());
        /* check for small negative numbers resulting from rounding */
        if (v0_sq < 0.0) v0_sq = 0.0;

        vO = state->pos.normalized();
        vO *= (std::sqrt (v0_sq) * signSinTrueAnomaly);
      }

    /* add to get velocity vector */
    state->vel = vPro + vO;
  }

  int elements2StateVector (
    btScalar mu,                  /* standard gravitational parameter */
    const sElements* elements,    /* pointer to orbital elements at epoch */
    sStateVector* state,          /* pointer to location where state vector will be stored */
    btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
    int maxIterations)            /* max number of iterations for calculating eccentric anomaly */
  {
    /*
    Code based on kostElements2StateVector1, made by TBlaxland
    */

    /* Returns number of iterations if successful, returns 0 otherwise. */

    /* Pseudocode
     *
     * get true anomaly
     * get longitude of ascending node and argument of periapsis
     * calc state vectors */

    /* get true anomaly */
    btScalar trueAnomaly;
    int ret (getTrueAnomaly (mu, elements, &trueAnomaly, maxRelativeError, maxIterations));

    /* calc state vectors */
    elements2StateVector2 (mu, elements, state, trueAnomaly);

    return ret;
  }

  void stateVector2Elements (
    btScalar mu,
    const sStateVector* state,
    sElements* elements,
    sOrbitParam* params)
  {
    if (params == NULL)
    {
      static sOrbitParam dummy;
      params = &dummy;
    }

    btVector3 vel (state->vel);

    btVector3 h (geth(*state));

    /*
    Radial orbits are not supported.
    e is not significantly different from 1.0
    if |h| < √(ε x µ x |r|)
    */
    if (h.length2() < SIMD_EPSILON * mu * state->pos.length())
    {
      /*
      We assume that the position is non-zero.
      Otherwise, we are in a singularity anyway,
      and no formula will get us out of there.
      */

      /* component of v parallel to pos */
      btVector3 v_parallel ((state->pos.dot (vel) / state->pos.length2()) * state->pos);

      /*
      Calculate how large the orthogonal component
      should be to make e significantly different
      from 1.0:

      |v_ortho| = √(ε x μ / |r|)
      */
      btScalar v_ortho_size (std::sqrt (SIMD_EPSILON * mu / state->pos.length()));

      /* New orthogonal component */
      btVector3 v_ortho = btVector3 (0.0, 1.0, 0.0);
      v_ortho = state->pos.cross (v_ortho);
      v_ortho = (v_ortho_size / v_ortho.length()) * v_ortho;

      /* replace the old orthogonal part */
      vel = v_parallel + v_ortho;

      h = state->pos.cross (vel);
    }

    btVector3 n (getn(h));

    btScalar E (vel.length2() / 2 - mu / state->pos.length());
    if (E == 0.0)
      E = SIMD_EPSILON;

    /*
    Alternative formula for e:
    e = (v x h) / μ - r / |r|
    */
    btVector3 e (gete(mu, *state));

    btScalar abse (e.length());
    /* parabolic orbit are not supported */
    if (abse > 1.0 - SIMD_EPSILON && abse < 1.0 + SIMD_EPSILON)
    {
      if (E >= 0.0)
      {
        abse = 1.0 + SIMD_EPSILON; /* Approximate with hyperbolic */
      }
      else
      {
        abse = 1.0 - SIMD_EPSILON; /* Approximate with elliptic */
      }
    }

    bool isEquatorial (n.length() < SIMD_EPSILON);
    bool isCircular (abse < SIMD_EPSILON);
    bool isHyperbola (abse >= 1.0);

    /*Ecc*/
    elements->e = abse;

    /*
      SMa
      dp = a(1 - e)
      da = a(1 + e)
    */
    elements->a = -mu / (2.0 * E);
    if (isHyperbola)
    {
      params->PeD = h.length2() / mu;
      params->ApD = INFINITY;
    }
    else
    {
      params->ApD = elements->a * (1.0 + elements->e);
      params->PeD = elements->a * (1.0 - elements->e);
    }

    /*Inc*/
    if (h.length() == 0.0)
      {
        /*
        Avoid division by zero absh
        By convention, take the smallest possible i,
        which is the angle between r and the
        equatorial plane.
        */
        elements->i = std::fmod (std::asin (state->pos.getY() / state->pos.length()), SIMD_PI);
      }
    else
      {
        elements->i = std::acos (h.getY() / h.length());
      }

    /* Longitude of ascending node */
    elements->LaN = getLaN(n);

    /* Argument of Periapsis*/
    params->AgP = getAgP (e, n, h, isCircular, isEquatorial);

    /*TrA*/
    bool sin_TrA_isNegative (false);
    getTrAFromState(state->pos, vel, n, e, isCircular, isEquatorial, &sin_TrA_isNegative);

    /*Lec*/
    params->Lec = elements->a * elements->e;

    /*
    SMi
    b² = a²(1 - e²)
    */
    if (isHyperbola)
      {
        params->SMi =
          std::sqrt (elements->a * elements->a * (elements->e * elements->e - 1.0) );
      }
    else
      {
        params->SMi =
          std::sqrt (elements->a * elements->a * (1.0 - elements->e * elements->e) );
      }

    /*LPe*/
    elements->LoP = std::fmod (elements->LaN + params->AgP, SIMD_2_PI);

    /*EcA*/
    if (isHyperbola)
      {
        btScalar tmp = (1.0 - state->pos.length() / elements->a) / elements->e;

        /*Avoid acosh out of range:*/
        if (tmp <= 1.0)
          {
            params->EcA = 0.0;
          }
        else
          {
            params->EcA = std::acosh (tmp);
          }
      }
    else if (isCircular)
      {
        params->EcA = 0.0;
      }
    else
      {
        btScalar tmp = (1.0 - state->pos.length() / elements->a) / elements->e;

        /* Avoid acos out of range. 1 and -1 are included cause we now the result. */
        if (tmp <= -1.0)
          {
            params->EcA = SIMD_PI;
          }
        else if (tmp >= 1.0)
          {
            params->EcA = 0.0;
          }
        else
          {
            params->EcA = std::acos (tmp);
          }
      }

    if (isHyperbola)
      {
        /*Copy sign from sin(TrA)*/
        if (sin_TrA_isNegative != (params->EcA < 0.0) )
          params->EcA = -params->EcA;
      }
    else
      {
        /*Same rule basically, but with EcA in 0..2π range*/
        if (sin_TrA_isNegative)
          params->EcA = SIMD_2_PI - params->EcA;
      }

    /*MnA*/
    if (isHyperbola)
      {
        params->MnA = elements->e * std::sinh (params->EcA) - params->EcA;
      }
    else
      {
        params->MnA = params->EcA - elements->e * std::sin (params->EcA);
      }

    /*MnL*/
    elements->L = params->MnA + elements->LoP;
    if (!isHyperbola)
      elements->L = std::fmod (elements->L, SIMD_2_PI);

    /*TrL*/
    params->TrL = std::fmod (elements->LoP + params->TrA, SIMD_2_PI);

    /*
    T = 2π√(a³/μ)

    fabs is for supporting hyperbola
    */
    params->T = SIMD_2_PI * std::sqrt (std::fabs(std::pow(elements->a, 3) / mu));

    /*
    Calculating PeT and ApT:
    */
    btScalar tPe (params->MnA * params->T / SIMD_2_PI); /*Time since last Pe*/

    if (isHyperbola)
      {
        params->PeT = -tPe;
      }
    else
      {
        params->PeT = params->T - tPe;
      }

    params->ApT = 0.5 * params->T - tPe;
    if (params->ApT < 0.0) params->ApT += params->T;
  }
}
