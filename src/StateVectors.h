/******************************************
*                                         *
* Project: SpaceFuckery                   *
* Function:                               *
*                                         *
*   Copyright 2015 - Marc-Olivier Barre   *
*           All rights reserved           *
*                                         *
******************************************/

#ifndef STATEVECTORS_H
#define STATEVECTORS_H

#include <iostream>
#include "LinearMath/btScalar.h"
#include "LinearMath/btVector3.h"

namespace mKOST
{
  //! This object describes a set of state vectors
  ATTRIBUTE_ALIGNED16(class) StateVectors
  {
    public:
      BT_DECLARE_ALIGNED_ALLOCATOR();

      btVector3 pos;  /**< Position */
      btVector3 vel;  /**< Velocity */
      btScalar MeL;   /**< Mean longitude at epoch */

      /** Default constructor */
      StateVectors();

      /** Default destructor */
      virtual ~StateVectors();

      /** \brief Writes data as CSV to a stream
       *
       * \param o The stream to write to
       * \param v The object which content is to be writen
       * \return Pointer to a stream
       *
       */
      friend std::ostream &operator << (std::ostream &o, const StateVectors &v);
  };
}

#endif // STATEVECTORS_H
