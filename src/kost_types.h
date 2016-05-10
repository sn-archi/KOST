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

/*
KOST is the Kepler Orbit Simulation Toolkit.
This header defines the data types of KOST.
*/

#ifndef KOST_TYPES_H
#define KOST_TYPES_H

#include "kost_settings.h"
#include "LinearMath/btVector3.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef kostMatrix3
typedef struct
{

	btScalar
		m11, m12, m13,
		m21, m22, m23,
		m31, m32, m33;

} kostMatrix3;
#endif

#ifndef kostElements
typedef struct
{
	btScalar a;      /*Semi-major axis*/
	btScalar e;      /*Eccentricity*/
	btScalar i;      /*Inclination*/
	btScalar theta;  /*Longitude of ascending node*/
	btScalar omegab; /*Longitude of periapsis*/
	btScalar L;      /*Mean longitude at epoch*/
} kostElements;
#endif

#ifndef kostOrbitParam
typedef struct
{
	/*Same as ORBITPARAM*/
	btScalar SMi;  /*semi-minor axis*/
	btScalar PeD;  /*periapsis distance*/
	btScalar ApD;  /*apoapsis distance*/
	btScalar MnA;  /*mean anomaly*/
	btScalar TrA;  /*true anomaly*/
	btScalar MnL;  /*mean longitude*/
	btScalar TrL;  /*true longitude*/
	btScalar EcA;  /*eccentric anomaly*/
	btScalar Lec;  /*linear eccentricity*/
	btScalar T;    /*orbit period*/
	btScalar PeT;  /*time to next periapsis passage*/
	btScalar ApT;  /*time to next apoapsis passage*/

	/*Additional*/
	btScalar AgP;  /*argument of periapsis*/
} kostOrbitParam;
#endif

typedef struct
{
	btVector3 pos;
	btVector3 vel;
} kostStateVector;


typedef struct
{
	btVector3 pe, ap, dn, an;

	btVector3 *points;
	unsigned int numPoints;

} kostOrbitShape;

#ifdef __cplusplus
}
#endif

#endif


