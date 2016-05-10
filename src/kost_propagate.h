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
KOST is the Kepler Orbital Simulation Toolkit.
This header file contains function for determining past/future states
*/

#ifndef KOST_TIME_H
#define KOST_TIME_H

#include "kost_types.h"

namespace mKOST
{
#ifdef __cplusplus
  extern "C" {
#endif

  btScalar getMeanAnomalyAtTime (
    btScalar mu,                  /* standard gravitational parameter */
    const kostElements* elements, /* pointer to orbital elements at epoch */
    btScalar timeSinceEpoch);     /* time since epoch in seconds */

  int getTrueAnomalyAtTime (
    btScalar mu,                  /* standard gravitational parameter */
    const kostElements* elements, /* pointer to orbital elements at epoch */
    btScalar* trueAnomaly,        /* location where result will be stored */
    btScalar timeSinceEpoch,      /* time since epoch in seconds */
    btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
    int maxIterations);           /* max number of iterations for calculating eccentric anomaly */

  btScalar getLANAtTime (
    btScalar mu,                  /* standard gravitational parameter */
    const kostElements* elements, /* pointer to orbital elements at epoch */
    btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
    btScalar jTwoCoeff,           /* J2 coefficient of the non-spherical body being orbited */
    btScalar timeSinceEpoch);     /* time since epoch in seconds */

  btScalar getArgPeAtTime (
    btScalar mu,                  /* standard gravitational parameter */
    const kostElements* elements, /* pointer to orbital elements at epoch */
    btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
    btScalar jTwoCoeff,           /* J2 coefficient of the non-spherical body being orbited */
    btScalar timeSinceEpoch);     /* time since epoch in seconds */

  int elements2StateVectorAtTime (
    btScalar mu,                  /* standard gravitational parameter */
    const kostElements* elements, /* pointer to orbital elements at epoch */
    kostStateVector* state,       /* pointer to where state vectors at epoch+timeSinceEpoch will be stored */
    btScalar timeSinceEpoch,      /* time since epoch in seconds */
    btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
    int maxIterations,            /* max number of iterations for calculating eccentric anomaly */
    btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
    btScalar jTwoCoeff);          /* J2 coefficient of the non-spherical body being orbited */

  void getElementsAtTime (
    btScalar mu,                     /* standard gravitational parameter */
    const kostElements* elements,    /* pointer to orbital elements at epoch */
    kostElements* newElements,       /* pointer to where elements at epoch+timeSinceEpoch will be stored */
    btScalar timeSinceEpoch,         /* time since epoch in seconds */
    btScalar bodyRadius,             /* mean radius of the non-spherical body being orbited */
    btScalar jTwoCoeff);             /* J2 coefficient of the non-spherical body being orbited */

#ifdef __cplusplus
}
#endif
}
#endif
