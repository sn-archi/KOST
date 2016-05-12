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
This header file contains orbital element tools
*/

#ifndef KOST_ELEMENTS_H
#define KOST_ELEMENTS_H

#include "kost_types.h"

namespace mKOST
{
#ifdef __cplusplus
  extern "C" {
#endif

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

  void stateVector2Elements (
    btScalar mu,                  /* standard gravitational parameter */
    const sStateVector* state, /* pointer to state vector at epoch */
    sElements* elements,       /* pointer to location where orbital elements at epoch will be stored */
    sOrbitParam* params);      /* pointer to location where extra orbital parameters will be stored */

#ifdef __cplusplus
}
#endif
}
#endif