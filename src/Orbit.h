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
mKOST is the modified Kepler Orbital Simulation Toolkit.
This header file contains the Orbit object
*/

#ifndef ELEMENTS_H
#define ELEMENTS_H

#include "types.h"

namespace mKOST
{
  class Orbit
  {
    public:
      btVector3 calc_e (btScalar mu, btVector3 pos, btVector3 vel);
      btVector3 geth (sStateVector state);
      btVector3 getn (btVector3 h);
      btVector3 getn (btScalar LaN);
      btScalar getLaN (btVector3 n);
      btScalar getAgP (btVector3 e, btVector3 n, btVector3 h, bool isCircular, bool isEquatorial);
      btScalar getTrAFromState (btVector3 pos, btVector3 vel, btVector3 n, btVector3 e, bool isCircular, bool isEquatorial, bool* sin_TrA_isNegative);

      btScalar getMeanAnomaly (
        btScalar mu,                   /* standard gravitational parameter */
        const sElements* elements); /* pointer to orbital elements at epoch */

      int getEccentricAnomaly (
        const sElements* elements,      /* pointer to orbital elements at epoch */
        btScalar* eccentricAnomaly,        /* location where result will be stored */
        btScalar meanAnomaly,              /* mean anomaly */
        btScalar eccentricAnomalyEstimate, /* initial estimate of eccentric anomaly, start with mean anomaly if no better estimate available */
        btScalar maxRelativeError,         /* maximum relative error in eccentric anomaly */
        int maxIterations);                /* max number of iterations for calculating eccentric anomaly */

      int getTrueAnomaly (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar* trueAnomaly,        /* location where result will be stored */
        btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
        int maxIterations);           /* max number of iterations for calculating eccentric anomaly */

      btScalar getTrueAnomaly2 (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar eccentricAnomaly);   /* eccentric anomaly */

      void elements2StateVector2 (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        sStateVector* state,       /* pointer to location where state vector at epoch will be stored */
        btScalar trueAnomaly);        /* true anomaly */

      int elements2StateVector (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        sStateVector* state,       /* pointer to location where state vector will be stored */
        btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
        int maxIterations);           /* max number of iterations for calculating eccentric anomaly */

      int stateVector2Elements (
        btScalar mu,                  /* standard gravitational parameter */
        const sStateVector* state, /* pointer to state vector at epoch */
        sElements* elements,       /* pointer to location where orbital elements at epoch will be stored */
        sOrbitParam* params);      /* pointer to location where extra orbital parameters will be stored */

      void elements2Shape (const sElements* elements, sOrbitShape* shape);

      btScalar getMeanAnomalyAtTime (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar timeSinceEpoch);     /* time since epoch in seconds */

      int getTrueAnomalyAtTime (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar* trueAnomaly,        /* location where result will be stored */
        btScalar timeSinceEpoch,      /* time since epoch in seconds */
        btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
        int maxIterations);           /* max number of iterations for calculating eccentric anomaly */

      btScalar getLANAtTime (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
        btScalar jTwoCoeff,           /* J2 coefficient of the non-spherical body being orbited */
        btScalar timeSinceEpoch);     /* time since epoch in seconds */

      btScalar getArgPeAtTime (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
        btScalar jTwoCoeff,           /* J2 coefficient of the non-spherical body being orbited */
        btScalar timeSinceEpoch);     /* time since epoch in seconds */

      int elements2StateVectorAtTime (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        sStateVector* state,       /* pointer to where state vectors at epoch+timeSinceEpoch will be stored */
        btScalar timeSinceEpoch,      /* time since epoch in seconds */
        btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
        int maxIterations,            /* max number of iterations for calculating eccentric anomaly */
        btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
        btScalar jTwoCoeff);          /* J2 coefficient of the non-spherical body being orbited */

      void getElementsAtTime (
        btScalar mu,                     /* standard gravitational parameter */
        const sElements* elements,    /* pointer to orbital elements at epoch */
        sElements* newElements,       /* pointer to where elements at epoch+timeSinceEpoch will be stored */
        btScalar timeSinceEpoch,         /* time since epoch in seconds */
        btScalar bodyRadius,             /* mean radius of the non-spherical body being orbited */
        btScalar jTwoCoeff);             /* J2 coefficient of the non-spherical body being orbited */
    private:
      btScalar a;      /* Semi-major axis */
      btScalar e;      /* Eccentricity */
      btScalar i;      /* Inclination */
      btScalar LaN;    /* Longitude of ascending node */
      btScalar LoP;    /* Longitude of periapsis */
      btScalar L;      /* Mean longitude at epoch */
      btScalar SMi;    /* semi-minor axis */
      btScalar PeD;    /* periapsis distance */
      btScalar ApD;    /* apoapsis distance */
      btScalar MnA;    /* mean anomaly */
      btScalar TrA;    /* true anomaly */
      btScalar MnL;    /* mean longitude */
      btScalar TrL;    /* true longitude */
      btScalar EcA;    /* eccentric anomaly */
      btScalar Lec;    /* linear eccentricity */
      btScalar T;      /* orbital period */
      btScalar PeT;    /* time to next periapsis passage */
      btScalar ApT;    /* time to next apoapsis passage */
      btScalar AgP;    /* argument of periapsis */
      btVector3 pos;   /* Position at epoch */
      btVector3 vel;   /* Velocity at epoch */
      btVector3 pe, ap, dn, an; /* node positions */
      btVector3* points;
      unsigned int numPoints;
      bool isCircular;
      bool isEquatorial;
  };
}
#endif // ELEMENTS_H
