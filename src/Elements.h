/******************************************
*                                         *
* Project: SpaceFuckery                   *
* Function:                               *
*                                         *
*   Copyright 2015 - Marc-Olivier Barre   *
*           All rights reserved           *
*                                         *
******************************************/

#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <iostream>
#include "LinearMath/btScalar.h"
#include "LinearMath/btVector3.h"
#include "LinearMath/btMatrix3x3.h"

namespace mKOST
{
  ATTRIBUTE_ALIGNED16(class) Elements
  {
    public:
      BT_DECLARE_ALIGNED_ALLOCATOR();

      btScalar a;    /**< Semi-major axis */
      btScalar Ecc;  /**< Eccentricity */
      btScalar i;    /**< Inclination */
      btScalar LAN;  /**< Longitude of ascending node */
      btScalar LoP;  /**< Longitude of periapsis */

      /** Default constructor */
      Elements();

      /** Default destructor */
      virtual ~Elements();

      /** \brief Writes data as CSV to a stream
       *
       * \param o The stream to write to
       * \param v The object which content is to be writen
       * \return Pointer to a stream
       *
       */
      friend std::ostream &operator << (std::ostream &o, const Elements &v);
  };
}
#endif // ELEMENTS_H