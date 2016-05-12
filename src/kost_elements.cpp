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

#include "kost_linalg.h"
#include "kost_constants.h"

#include "kost_elements.h"

namespace mKOST
{
  btScalar getMeanAnomaly (
    btScalar mu,                   /* standard gravitational parameter */
    const sElements* elements)   /* pointer to orbital elements at epoch */
  {
    btScalar meanAnomaly;
    /* calc mean anomaly */
    meanAnomaly = elements->L - elements->omegab;

    if (elements->e < 1.0)
      {
        /* check range is in 0 to 2*pi */
        if (meanAnomaly < 0.0) meanAnomaly += M_TWOPI;
        if (meanAnomaly >  M_TWOPI) meanAnomaly -= M_TWOPI;
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

    int i;
    btScalar relativeError, meanAnomalyEstimate, e;

    if (elements->e == 0.0)   /* circular orbit */
      {
        *eccentricAnomaly = meanAnomaly;
        return 1;
      }

    if (elements->e == 1.0) /* parabolic orbit - approximate to hyperbolic */
      e = elements->e + KOST_VERYSMALL;
    else
      e = elements->e;

    i = 0;
    do
      {
        if (elements->e < 1.0)   /* elliptical orbit */
          {
            /* calculate next estimate of the root of Kepler's equation using Newton's method */
            eccentricAnomalyEstimate = eccentricAnomalyEstimate - (eccentricAnomalyEstimate - e * sin (eccentricAnomalyEstimate) - meanAnomaly) / (1 - e * cos (eccentricAnomalyEstimate) );
            /* calculate estimate of mean anomaly from estimate of eccentric anomaly */
            meanAnomalyEstimate = eccentricAnomalyEstimate - e * sin (eccentricAnomalyEstimate);
          }
        else   /* hyperbolic orbit */
          {
            /* calculate next estimate of the root of Kepler's equation using Newton's method */
            eccentricAnomalyEstimate = eccentricAnomalyEstimate - (e * sinh (eccentricAnomalyEstimate) - eccentricAnomalyEstimate - meanAnomaly) / (e * cosh (eccentricAnomalyEstimate) - 1.0);
            /* calculate estimate of mean anomaly from estimate of eccentric anomaly */
            meanAnomalyEstimate = e * sinh (eccentricAnomalyEstimate) - eccentricAnomalyEstimate;
          }
        /* calculate relativeError */
        relativeError = 1.0 - meanAnomalyEstimate / meanAnomaly;
        i++;
      }
    while ( (i < maxIterations) && (fabs (relativeError) > fabs (maxRelativeError) ) );

    if (elements->e < 1.0)
      {
        /* check range is in 0 to 2*pi */
        eccentricAnomalyEstimate = fmod (eccentricAnomalyEstimate, M_TWOPI);
        if (eccentricAnomalyEstimate < 0.0) eccentricAnomalyEstimate += M_TWOPI;
        if (eccentricAnomalyEstimate > M_TWOPI) eccentricAnomalyEstimate -= M_TWOPI;
      }

    *eccentricAnomaly = eccentricAnomalyEstimate;

    if ( fabs (relativeError) < fabs (maxRelativeError) )
      return i;
    else
      return 0;
  }

  int getTrueAnomaly (
    btScalar mu,                  /* standard gravitational parameter */
    const sElements* elements, /* pointer to orbital elements at epoch */
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
                meanAnomaly, log (2.0 * meanAnomaly / elements->e + 1.8),
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
        ret = 2.0 * atan (sqrt ( (1.0 + elements->e) / (1.0 - elements->e) ) * tan (eccentricAnomaly / 2.0) );
      }
    else   /* hyperbolic orbit */
      {
        ret = acos ( (cosh (eccentricAnomaly) - elements->e) / (1 - elements->e * cosh (eccentricAnomaly) ) );
        if (eccentricAnomaly < 0.0) ret = -ret; /* Always the same sign */
      }

    return ret;
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

    btVector3 n; /* unit vector in direction of ascending node */
    btVector3 h; /* angular momentum vector */
    btVector3 north; /* unit vector pointing ecliptic north */
    btVector3 vPro, vO; /* prograde and outward components of velocity vector */
    btVector3 tmpv; /* temporary vector value */
    btScalar argPos; /* argument of position, measured from the longitude of ascending node */
    btScalar rPe, vPe; /* radius and velocity at periapsis */
    btScalar v2; /* magnitude squared of velocity vector */
    btScalar e; /* eccentricity */
    btScalar tmpr; /* temporary real value */

    e = elements->e;
    if (e == 1.0) /* parabolic orbit - approximate to hyperbolic orbit */
      e += KOST_VERYSMALL;

    /* calc nodal vector */
    n = btVector3 (std::cos (elements->theta), 0.0, std::sin (elements->theta) );

    /* equatorial north vector */
    north = btVector3 (0.0, 1.0, 0.0);

    /* calc angular momentum vector, h */
    /* projection of h in ecliptic (xz) plane */
    h = n.cross (north);
    h = std::sin (elements->i) * h;

    /* elevation of h */
    h.setY (std::cos (elements->i) );
    h.normalize() ;

    /* calc magnitude of h */
    /* calc radius and velocity at periapsis */
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
    h = (rPe * vPe) * h;

    /* calc position vector */
    /* calc argument of position */
    argPos = elements->omegab - elements->theta + trueAnomaly;

    /*
    calc direction of position vector:
    r/|r| = sin(ArgPos) * ((h / |h|) x n) + cos(argPos) * n
    */
    tmpv = std::cos (argPos) * n;
    state->pos = h.normalized();
    state->pos = state->pos.cross (n);
    state->pos = std::sin (argPos) * state->pos;
    state->pos = state->pos + tmpv;

    /* calc length of position vector */
    if (e < 1.0) /* elliptical orbit */
      state->pos = elements->a * (1.0 - e * e) / (1.0 + e * std::cos (trueAnomaly) ) * state->pos;
    else /* hyperbolic orbit */
      state->pos = std::fabs (elements->a) * (e * e - 1.0) / (1.0 + e * std::cos (trueAnomaly) ) * state->pos;

    /* calc velocity vector */
    /*  calculate magnitude of velocity vector */
    if (e < 1.0)   /* elliptical orbit */
      {
        v2 = mu * (2.0 / state->pos.length() ) - 1.0 / elements->a;
      }
    else   /* hyperbolic orbit */
      {
        v2 = mu * (2.0 / state->pos.length() ) + 1.0 / std::fabs (elements->a);
      }

    /* calc components of velocity vector perpendicular and parallel to radius vector:

    perpendicular:
    vPro = (|h|/|pos|) * normal(h x pos)

    parallel:
    vO = sqrt(v^2 - |vPro|^2) * sign(sin(trueAnomaly)) * normal(pos)
    */
    vPro = h.cross (state->pos);
    vPro.normalize();
    vPro = (h.length() / state->pos.length() ) * vPro;

    tmpr = std::sin (trueAnomaly);
    if (tmpr == 0.0)   /* check for apsis condition to avoid divide by zero */
      {
        vO = btVector3 (0.0, 0.0, 0.0);
      }
    else
      {
        btScalar signSinTrueAnomaly = tmpr / std::fabs (tmpr);

        btScalar v0_sq = v2 - vPro.length2();
        /* check for small negative numbers resulting from rounding */
        if (v0_sq < 0.0) v0_sq = 0.0;

        vO = state->pos.normalized();
        vO = (std::sqrt (v0_sq) * signSinTrueAnomaly) * vO;
      }

    /* add to get velocity vector */
    state->vel = vPro + vO;
  }

  int elements2StateVector (
    btScalar mu,                  /* standard gravitational parameter */
    const sElements* elements, /* pointer to orbital elements at epoch */
    sStateVector* state,       /* pointer to location where state vector will be stored */
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

    int ret;
    btScalar trueAnomaly;

    /* get true anomaly */
    ret = getTrueAnomaly (mu, elements, &trueAnomaly, maxRelativeError, maxIterations);

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
    /*
    See appendix C in orbiter.pdf
    */

    btVector3 vel, h, n, e;
    btScalar absh, absn, absr, abse, E, tPe;
    bool isEquatorial, isCircular, isHyperbola;
    bool sin_TrA_isNegative = false;

    if (params == NULL)
      {
        static sOrbitParam dummy;
        params = &dummy;
      }

    vel = state->vel;

    absr  = state->pos.length();

    h = state->pos.cross (vel);
    absh = h.length();

    /*
    Radial orbits are not supported.
    e is not significantly different from 1.0
    if |h| < sqrt(epsilon * mu * |r|)
    */
    if (absh * absh < KOST_VERYSMALL * mu * absr)
      {
        /*
        We assume that the position is non-zero.
        Otherwise, we are in a singularity anyway,
        and no formula will get us out of there.
        */
        btVector3 v_ortho, v_parallel;
        btScalar v_ortho_size;

        /* component of v parallel to pos */
        v_parallel = (state->pos.dot (vel) / state->pos.length2() ) * state->pos;

        /*
        Calculate how large the orthogonal component
        should be to make e significantly different
        from 1.0:

        |v_ortho| = sqrt(epsilon*mu / |r|)
        */
        v_ortho_size = sqrt (KOST_VERYSMALL * mu / absr);

        /* New orthogonal component */
        v_ortho = btVector3 (0.0, 1.0, 0.0);
        v_ortho = state->pos.cross (v_ortho);
        v_ortho = (v_ortho_size / v_ortho.length() ) * v_ortho;

        /* replace the old orthogonal part */
        vel = v_parallel + v_ortho;

        h = state->pos.cross (vel);
        absh = h.length();
      }

    n = btVector3 (-h.getZ(), h.getX(), 0.0);
    absn = n.length();

    E = 0.5 * vel.length2() - mu / absr;
    if (E == 0.0)
      E = KOST_VERYSMALL;

    /*
    Alternative formula for e:
    e = (v x h) / mu - r / |r|
    */
    e = vel.cross (h);
    e = (absr / mu) * e;
    e = e - state->pos;
    e = (1.0 / absr) * e;

    abse = e.length();

    /* parabolic orbit are not supported */
    if (abse > 1.0 - KOST_VERYSMALL && abse < 1.0 + KOST_VERYSMALL)
      {
        if (E >= 0.0)
          {
            abse = 1.0 + KOST_VERYSMALL; /* Approximate with hyperbolic */
          }
        else
          {
            abse = 1.0 - KOST_VERYSMALL; /* Approximate with elliptic */
          }
      }

    isEquatorial = absn < KOST_VERYSMALL;
    isCircular   = abse < KOST_VERYSMALL;
    isHyperbola  = abse >= 1.0;

    /*Ecc*/
    elements->e = abse;

    /*
      SMa
      dp = a * (1-e)
      da = a * (1+e)
    */
    if (isHyperbola)
      {
        std::cout << absh << std::endl;
        params->PeD = h.length2() / mu;
        elements->a = INFINITY;
        params->ApD = INFINITY;
      }
    else
      {
        elements->a = -mu / (2.0 * E);
        params->ApD = elements->a * (1.0 + elements->e);
        params->PeD = elements->a * (1.0 - elements->e);
      }

    /*Inc*/
    if (absh == 0.0)
      {
        /*
        Avoid division by zero absh
        By convention, take the smallest possible i,
        which is the angle between r and the
        equatorial plane.
        */
        elements->i = std::fmod (std::asin (state->pos.getY() / absr), M_PI);
      }
    else
      {
        elements->i = std::fmod (std::acos (h.getY() / absh), M_PI);
      }

    /*LAN*/
    if (isEquatorial)
      {
        elements->theta = 0.0;
      }
    else
      {
        elements->theta = std::acos (n.getX() / absn);
        if (n.getZ() < 0.0) elements->theta = M_TWOPI - elements->theta;
      }

    /*AgP*/
    params->AgP = 0.0;
    if (isCircular)
      {
        params->AgP = 0.0;
      }
    else if (isEquatorial)
      {
        params->AgP = atan2 (e.getZ(), e.getX() );
        if (h.getY() < 0.0) params->AgP = -params->AgP;
      }
    else
      {
        params->AgP = std::acos (n.dot (e) / (absn * abse));
        if (e.getY() < 0.0) params->AgP = M_TWOPI - params->AgP;
      }

    /*TrA*/
    sin_TrA_isNegative = false;
    if (isCircular)
      {
        if (isEquatorial)
          {
            params->TrA = std::acos (state->pos.getX() / absr);
            if (vel.getX() > 0.0)
              {
                sin_TrA_isNegative = true;
                params->TrA = M_TWOPI - params->TrA;
              }
          }
        else
          {
            params->TrA = std::acos (n.dot (state->pos) ) / (absn * absr);
            if (n.dot (vel) > 0.0)
              {
                sin_TrA_isNegative = true;
                params->TrA = M_TWOPI - params->TrA;
              }
          }
      }
    else
      {
        btScalar tmp = e.dot (state->pos) / (abse * absr);

        /*Avoid acos out of range:*/
        if (tmp <= -1.0)
          {
            params->TrA = M_PI;
          }
        else if (tmp >= 1.0)
          {
            params->TrA = 0.0;
          }
        else
          {
            params->TrA = std::acos (tmp);
          }

        if (state->pos.dot (vel) < 0.0)
          {
            sin_TrA_isNegative = true;
            params->TrA = M_TWOPI - params->TrA;
          }
      }

    /*Lec*/
    params->Lec = elements->a * elements->e;

    /*
    SMi
    b^2 = a^2 * (1 - e^2)
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
    elements->omegab = std::fmod (elements->theta + params->AgP, M_TWOPI);

    /*EcA*/
    if (isHyperbola)
      {
        btScalar tmp = (1.0 - absr / elements->a) / elements->e;

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
        btScalar tmp = (1.0 - absr / elements->a) / elements->e;

        /*Avoid acos out of range:*/
        if (tmp <= -1.0)
          {
            params->EcA = M_PI;
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
        /*Same rule basically, but with EcA in 0..2pi range*/
        if (sin_TrA_isNegative)
          params->EcA = M_TWOPI - params->EcA;
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
    elements->L = params->MnA + elements->omegab;
    if (!isHyperbola)
      elements->L = std::fmod (elements->L, M_TWOPI);

    /*TrL*/
    params->TrL = std::fmod (elements->omegab + params->TrA, M_TWOPI);

    /*
    T = 2*pi*sqrt(a^3 / mu)

    fabs is for supporting hyperbola
    */
    params->T =
      M_TWOPI * std::sqrt (fabs (elements->a * elements->a * elements->a / mu) );

    /*
    Calculating PeT and ApT:
    */
    tPe = params->MnA * params->T / M_TWOPI; /*Time since last Pe*/

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