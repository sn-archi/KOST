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

#include <cmath>

#include "shape.h"

namespace mKOST
{
  void elements2Shape (const sElements* elements, sOrbitShape* shape)
  {
    unsigned int i = 0;

    /*Some utility values: */
    btScalar multiplier = elements->a * (1.0 - elements->e * elements->e);
    btScalar AgP = elements->LoP - elements->LaN;

    /*
    First: Orbit in its own coordinate system:
    */

    /*Pe, Ap*/
    shape->pe = btVector3 ( elements->a * (1.0 - elements->e), 0.0, 0.0);
    shape->ap = btVector3 (-elements->a * (1.0 + elements->e), 0.0, 0.0);

    /*Points*/
    if (shape->numPoints == 1)
      {
        shape->points[0] = shape->pe;
      }
    else if (shape->numPoints > 1)
      {
        btScalar maxTrA, dTrA, TrA;

        /*Range of angles*/
        maxTrA = SIMD_PI;
        if (elements->e >= 1.0)
          {
            maxTrA = std::acos (-1.0 / elements->e);

            /*Make it a bit smaller to avoid division by zero:*/
            maxTrA *= ( ( (btScalar) shape->numPoints) / (shape->numPoints + 1) );
          }

        /*Angle change per segment*/
        dTrA = (2 * maxTrA) / (shape->numPoints - 1);

        TrA = -maxTrA;
        for (i = 0; i < shape->numPoints; i++)
          {
            btScalar absr = std::fabs (multiplier / (1.0 + elements->e * std::cos (TrA) ) );

            btVector3 direction = btVector3 (std::cos (TrA), std::sin (TrA), 0.0);
            shape->points[i] = absr * direction;

            TrA += dTrA;
          }
      }


    /*AN*/
    {
      btScalar TrA = -AgP;
      btScalar absr = multiplier / (1.0 + elements->e * std::cos (TrA) );

      if (absr <= 0.0)
        {
          shape->an = btVector3 (0.0, 0.0, 0.0);
        }
      else
        {
          btVector3 direction = btVector3 (std::cos (TrA), std::sin (TrA), 0.0);
          shape->an = absr * direction;
        }
    }

    /*DN*/
    {
      btScalar TrA = SIMD_PI - AgP;
      btScalar absr = multiplier / (1.0 + elements->e * std::cos (TrA) );

      if (absr <= 0.0)
        {
          shape->dn = btVector3 (0.0, 0.0, 0.0);
        }
      else
        {
          btVector3 direction = btVector3 (std::cos (TrA), std::sin (TrA), 0.0);
          shape->dn = absr * direction;
        }
    }



    /*
    Then: rotate the coordinates:
    */
    {
      btMatrix3x3 AgPMat, LANMat, IncMat, transform;

      AgPMat.setEulerZYX (0.0, AgP, 0.0);
      IncMat.setEulerZYX (elements->i, 0.0, 0.0);
      LANMat.setEulerZYX (0.0, elements->LaN, 0.0);

      /* Now, global = LANMat * IncMat * AgPMat * local: */
      transform = LANMat * IncMat;
      transform *= AgPMat;

      shape->pe = transform * shape->pe;
      shape->ap = transform * shape->ap;
      shape->an = transform * shape->an;
      shape->dn = transform * shape->dn;

      if (shape->numPoints != 0)
        for (i = 0; i < shape->numPoints; i++)
          shape->points[i] = transform * shape->points[i];
    }
  }
}
