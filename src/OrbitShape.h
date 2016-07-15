/******************************************
*                                         *
* Project: SpaceFuckery                   *
* Function:                               *
*                                         *
*   Copyright 2015 - Marc-Olivier Barre   *
*           All rights reserved           *
*                                         *
******************************************/

#ifndef ORBITSHAPE_H
#define ORBITSHAPE_H

#include <iostream>
#include "LinearMath/btScalar.h"
#include "LinearMath/btVector3.h"

namespace mKOST
{
  ATTRIBUTE_ALIGNED16(class) OrbitShape
  {
    public:
      BT_DECLARE_ALIGNED_ALLOCATOR();

      btVector3 pe;            /**< Periaspsis position */
      btVector3 ap;            /**< Apoapsis position */
      btVector3 dn;            /**< Descending node position */
      btVector3 an;            /**< Ascending node position */
      btVector3* points;       /**< points */
      unsigned int numPoints;  /**< numPoints */

      /** Default constructor */
      OrbitShape();

      /** Default destructor */
      virtual ~OrbitShape();
  };
}
#endif // ORBITSHAPE_H
