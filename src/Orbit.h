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

#ifndef ELEMENTS_H
#define ELEMENTS_H

#include "types.h"

namespace mKOST
{
  //! This object describes an orbit and allows to compute positions and velocities at a chosen time

  /*!
    ### Usage notes - mKOST::Orbit::elements2StateVectorX and related methods:

    Parabolic orbits are NOT currently supported.

    Position as a function of time can be found by either:
    - A call to mKOST::Orbit::elements2StateVector. Depending on settings for
      maxIterations and maxRelativeError, this may adversly affect frame
      rates in graphical applications.
    - To minimise impact on frame rates:
      - Call mKOST::Orbit::getMeanAnomaly
      - Call mKOST::Orbit::getEccentricAnomaly on successive time steps with a
        small number of iterations in each call. Each call takes the
        result of the previous call as its eccentricAnomalyEstimate.
        Repeat until mKOST::Orbit::GetEccentricAnomaly returns > 0.
      - Call mKOST::Orbit::getTrueAnomaly2.
      - Call mKOST::Orbit::elements2StateVector2.
   */
  class Orbit
  {
    public:

      /** \brief Calculates the eccentricity vector
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param vel Velocity vector
       * \param pos Position vector
       * \return Eccentricity vector
       *
       */
      btVector3 calc_e (btScalar mu, btVector3 pos, btVector3 vel);

      /** \brief Calculate the angular momentum vector
       *
       * \param state State vector structure
       * \return The angular momentum as a btVector3
       *
       */
      btVector3 geth (sStateVector state);

      /** \brief Calculate the nodal vector n that points to the ascending node
       *
       * \param h Angular momentum vector
       * \return Nodal vector n
       *
       */
      btVector3 getn (btVector3 h);

      /** \brief Calculate the nodal vector n that points to the ascending node
       *
       * \param LaN Longitude of ascending node in radians
       * \return Nodal vector n
       *
       */
      btVector3 getn (btScalar LaN);

      /** \brief Calculate the longitude of ascending node
       *
       * \param n Nodal vector pointing to the ascending node
       * \return Longitude of ascending node in radians
       *
       */
      btScalar getLaN (btVector3 n);

      /** \brief Calculates the argument of periapsis
       *
       * \param e Eccentricity vector
       * \param n Nodal vector
       * \param h Angular momentum vector
       * \param isCircular
       * \param isEquatorial
       * \return Argument of periapsis in radians
       *
       */
      btScalar getAgP (btVector3 e, btVector3 n, btVector3 h, bool isCircular, bool isEquatorial);

      /** \brief Calculate true anomaly from state vectors and a shitload of other things
       *
       * \param pos Position vector
       * \param vel Velocity vector
       * \param n Nodal vector
       * \param e Eccentricity vector
       * \param isCircular
       * \param isEquatorial
       * \param sin_TrA_isNegative
       * \return True anomaly in radians at current epoch
       *
       */

      btScalar getTrAFromState (btVector3 pos, btVector3 vel, btVector3 n, btVector3 e, bool isCircular, bool isEquatorial, bool* sin_TrA_isNegative);

      /** \brief
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param elements
       * \return
       *
       */
      btScalar getMeanAnomaly (
        btScalar mu,                   /* standard gravitational parameter */
        const sElements* elements); /* pointer to orbital elements at epoch */

      /** \brief Determine eccentric anomaly from an estimate using a Newton-Raphson root finding algorythm
       *
       * \param elements
       * \param eccentricAnomaly
       * \param meanAnomaly
       * \param eccentricAnomalyEstimate
       * \param maxRelativeError
       * \param maxIterations
       * \return number of iterations or 0 if failed to find a solution
       *
       */
      int getEccentricAnomaly (
        const sElements* elements,      /* pointer to orbital elements at epoch */
        btScalar* eccentricAnomaly,        /* location where result will be stored */
        btScalar meanAnomaly,              /* mean anomaly */
        btScalar eccentricAnomalyEstimate, /* initial estimate of eccentric anomaly, start with mean anomaly if no better estimate available */
        btScalar maxRelativeError,         /* maximum relative error in eccentric anomaly */
        int maxIterations);                /* max number of iterations for calculating eccentric anomaly */

      /** \brief
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param elements
       * \param trueAnomaly
       * \param maxRelativeError
       * \param maxIterations
       * \return number of iterations or 0 if failed to find a solution
       *
       */
      int getTrueAnomaly (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar* trueAnomaly,        /* location where result will be stored */
        btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
        int maxIterations);           /* max number of iterations for calculating eccentric anomaly */

      /** \brief
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param elements
       * \param eccentricAnomaly
       * \return True anomaly in radians
       *
       */
      btScalar getTrueAnomaly2 (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar eccentricAnomaly);   /* eccentric anomaly */

      /** \brief Deduce state vectors from orbital elements. true anomaly will be used to determine position on the orbit
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param elements
       * \param state
       * \param trueAnomaly
       *
       */
      void elements2StateVector2 (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        sStateVector* state,       /* pointer to location where state vector at epoch will be stored */
        btScalar trueAnomaly);        /* true anomaly */

      /** \brief Deduce state vectors from orbital elements. L will be used to determine true anomaly
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param elements
       * \param state
       * \param maxRelativeError
       * \param maxIterations
       * \return number of iterations or 0 if failed to find a solution
       *
       */
      int elements2StateVector (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        sStateVector* state,       /* pointer to location where state vector will be stored */
        btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
        int maxIterations);           /* max number of iterations for calculating eccentric anomaly */

      /** \brief Deduce orbital elements from state vectors at epoch
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param state
       * \param elements
       * \param params
       * \return number of iterations or 0 if failed to find a solution
       *
       */
      int stateVector2Elements (
        btScalar mu,                  /* standard gravitational parameter */
        const sStateVector* state, /* pointer to state vector at epoch */
        sElements* elements,       /* pointer to location where orbital elements at epoch will be stored */
        sOrbitParam* params);      /* pointer to location where extra orbital parameters will be stored */

      /** \brief
       *
       * \param elements
       * \param shape
       *
       */
      void elements2Shape (const sElements* elements, sOrbitShape* shape);

      /** \brief
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param elements
       * \param timeSinceEpoch
       * \return
       *
       */
      btScalar getMeanAnomalyAtTime (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar timeSinceEpoch);     /* time since epoch in seconds */

      /** \brief
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param elements
       * \param trueAnomaly
       * \param timeSinceEpoch
       * \param maxRelativeError
       * \param maxIterations
       * \return number of iterations or 0 if failed to find a solution
       *
       */
      int getTrueAnomalyAtTime (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar* trueAnomaly,        /* location where result will be stored */
        btScalar timeSinceEpoch,      /* time since epoch in seconds */
        btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
        int maxIterations);           /* max number of iterations for calculating eccentric anomaly */

      /** \brief
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param elements
       * \param bodyRadius
       * \param jTwoCoeff
       * \param timeSinceEpoch
       * \return Longitude of ascending node at the specified time
       *
       */
      btScalar getLANAtTime (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
        btScalar jTwoCoeff,           /* J2 coefficient of the non-spherical body being orbited */
        btScalar timeSinceEpoch);     /* time since epoch in seconds */

      /** \brief
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param elements
       * \param bodyRadius
       * \param jTwoCoeff
       * \param timeSinceEpoch
       * \return Argument of periapsis at the specified time
       *
       */
      btScalar getArgPeAtTime (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
        btScalar jTwoCoeff,           /* J2 coefficient of the non-spherical body being orbited */
        btScalar timeSinceEpoch);     /* time since epoch in seconds */

      /** \brief
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param elements
       * \param state
       * \param timeSinceEpoch
       * \param maxRelativeError
       * \param maxIterations
       * \param bodyRadius
       * \param jTwoCoeff
       * \return number of iterations or 0 if failed to find a solution
       *
       */
      int elements2StateVectorAtTime (
        btScalar mu,                  /* standard gravitational parameter */
        const sElements* elements, /* pointer to orbital elements at epoch */
        sStateVector* state,       /* pointer to where state vectors at epoch+timeSinceEpoch will be stored */
        btScalar timeSinceEpoch,      /* time since epoch in seconds */
        btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
        int maxIterations,            /* max number of iterations for calculating eccentric anomaly */
        btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
        btScalar jTwoCoeff);          /* J2 coefficient of the non-spherical body being orbited */

      /** \brief
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param elements
       * \param newElements
       * \param timeSinceEpoch
       * \param bodyRadius
       * \param jTwoCoeff
       *
       */
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
