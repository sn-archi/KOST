/******************************************
*                                         *
* Project: SpaceFuckery                   *
* Function:                               *
*                                         *
*   Copyright 2015 - Marc-Olivier Barre   *
*           All rights reserved           *
*                                         *
******************************************/
#ifndef PARAMS_H
#define PARAMS_H

#include <iostream>
#include "LinearMath/btScalar.h"
#include "LinearMath/btVector3.h"
#include "LinearMath/btMatrix3x3.h"

namespace mKOST
{
  ATTRIBUTE_ALIGNED16(class) Params
  {
    public:
      BT_DECLARE_ALIGNED_ALLOCATOR();

      btScalar SMi;  /**< semi-minor axis */
      btScalar PeD;  /**< periapsis distance */
      btScalar ApD;  /**< apoapsis distance */
      btScalar MnA;  /**< mean anomaly */
      btScalar TrA;  /**< true anomaly */
      btScalar TrL;  /**< true longitude */
      btScalar EcA;  /**< eccentric anomaly */
      btScalar Lec;  /**< linear eccentricity */
      btScalar T;    /**< orbit period */
      btScalar PeT;  /**< time to next periapsis passage */
      btScalar ApT;  /**< time to next apoapsis passage */
      btScalar AgP;  /**< argument of periapsis */

      /** Default constructor */
      Params();

      /** Default destructor */
      virtual ~Params();

      /** \brief Writes data as CSV to a stream
       *
       * \param o The stream to write to
       * \param v The object which content is to be writen
       * \return Pointer to a stream
       *
       */
      friend std::ostream &operator << (std::ostream &o, const Params &v);
  };
}
#endif // PARAMS_H