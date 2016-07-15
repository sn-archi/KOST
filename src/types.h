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
#include "Elements.h"
#include "Params.h"
#include "StateVectors.h"
#include "OrbitShape.h"

/* Data about central body (earth) */
#define MU 3.986004418e14
#define R 6378100.0
#define EPSILON 12*SIMD_EPSILON

#endif // TYPES_H