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

  /*!Parabolic orbits are NOT currently supported.

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

    The positions are all relative to the central body. If you are orbiting
    earth and your reference is the sun, you'll have to do the translation
    to get something that makes sense !
   */
  ATTRIBUTE_ALIGNED16(class) Orbit
  {
    private:
      btScalar     mu;         /**< Standard gravitationnal parameter for the central body */
      Elements     *mElements;  /**< Orbital elements */
      Params       *mParams;    /**< Suplemental orbital parameters */
      OrbitShape   *mShape;     /**< Orbital shape */
      btVector3    *n;          /**< Vector pointing to the AN */
      btVector3    *e;          /**< Eccentricity vector */
      bool         Circular;   /**< True for circular orbits */
      bool         Hyperbola;  /**< True for hyperbolic orbits */
      bool         Equatorial; /**< True for equatorial orbits */

      /** \brief Calculates the eccentricity vector
       *
       * \param mu µ (standard gravitational parameter) value for the central body
       * \param vel Velocity vector
       * \param pos Position vector
       * \return Eccentricity vector
       *
       */
      void calcE (StateVectors* state);

      /** \brief Calculate the angular momentum vector
       *
       * Formula used: e = (v x h) / μ - r / |r|
       *
       * \param state State vector structure
       * \return The angular momentum vector
       *
       */
      btVector3 calcH (StateVectors* state);

      /** \brief Calculate the nodal vector n that points to the ascending node from the angular momentum
       *
       * \param h Angular momentum vector
       * \return Nodal vector
       *
       */
      void calcN (btVector3 *h);

      /** \brief Calculate the nodal vector n that points to the ascending node from the longitude of ascending node
       *
       * \param LaN Longitude of ascending node in radians
       * \return Nodal vector
       *
       */
      void calcN (void);

      /** \brief Calculate the longitude of ascending node
       *
       * \return Longitude of ascending node in radians
       *
       */
      btScalar calcLAN (void);

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
      btScalar calcAgP (btVector3 *h);

      /** \brief Calculate true anomaly from state vectors
       *
       * \return True anomaly in radians at current epoch
       *
       */
      void calcTrA (StateVectors* state);

      /** \brief
       *
       * *Pseudocode*
       *
       * - get mean anomaly
       * - get eccentric anomaly
       * - calc true anomaly
       *
       * \param MeL Mean Longitude at epoch
       * \param maxRelativeError Your tolerance for error in the approximation. Influences the number
       *        of runs you'll have to make to get an appropriate result
       * \param maxIterations the maximum number of iterations before stoping to throw the result not found exception
       * \return True anomaly aproximation
       *
       */
      btScalar calcTrA (btScalar MeL, btScalar maxRelativeError, int maxIterations);

      /** \brief Calculate true anomaly from eccentric anomaly
       *
       * \param eccentricAnomaly
       * \return True anomaly in radians
       *
       */
      btScalar calcTrA (btScalar eccentricAnomaly);

      /** \brief Calculate eccentric anomaly
       *
       * \return Eccentric anomaly in radians at current epoch
       *
       */
      btScalar calcEcA (StateVectors* state);

      /** \brief Calculate mean anomaly
       *
       * \param MeL Mean longitude at epoch
       * \return Mean anomaly in radians at current epoch
       *
       */
      btScalar calcMnA (btScalar MeL);

      /** \brief
       *
       * \return Mean anomaly in radians at current epoch
       *
       */
      btScalar calcMnA (void);

      /** \brief Determine eccentric anomaly from an estimate using a Newton-Raphson root finding algorythm
       * Code will terminate when either maxIterations or maxRelativeError is reached.
       * Returns number of iterations if relativeError < maxRelativeError, returns 0 otherwise.
       * *Pseudocode*
       *
       * do
       *  calculate next estimate of the root of Kepler's equation using Newton's method
       *  calculate estimate of mean anomaly from estimate of eccentric anomaly
       *  calculate relativeError
       * while ((iterations<=maxIterations)&&(relativeError<maxRelativeError))
       * if iterations<=maxIterations return iterations, else return 0
       * \param meanAnomaly
       * \param eccentricAnomalyEstimate
       * \param maxRelativeError
       * \param maxIterations
       * \return Mean anomaly estimate
       *
       */
      btScalar calcEcA (btScalar meanAnomaly,
        btScalar ecaEstimate,
        btScalar maxRelativeError,
        int maxIterations);

    public:
      BT_DECLARE_ALIGNED_ALLOCATOR();

      /** \brief Orbit no initialization constructor */
      Orbit(void);

      /** \brief Orbit constructor from state vectors
       *
       * \param mu
       * \param statevectors
       *
       */
      Orbit(btScalar mu, StateVectors* state);

      /** \brief Orbit constructor from orbital elements
       *
       * \param mu
       * \param elements
       */
      Orbit(btScalar mu, Elements* elements);

      /** \brief Orbit destructor
       *
       *
       */
      ~Orbit();

      /** \brief Deduce state vectors from orbital elements. True anomaly will be used to determine position on the orbit
       *
       * Pseudocode
       * - calc nodal vector, n
       * - calc angular momentum vector, h
       * - calc position vector
       * - calc argument of position
       * - calc direction of position vector
       * - calc length of position vector
       * - calc velocity vector
       * - calculate magnitude of velocity vector
       * - calc components of velocity vector perpendicular and parallel to radius vector
       * - add to get velocity vector
       * \param trueAnomaly
       *
       */
      StateVectors elements2StateVector (btScalar trueAnomaly);

      /** \brief Deduce state vectors from orbital elements. L will be used to determine true anomaly
       *
       * *Pseudocode*
       *
       * - get true anomaly
       * - get longitude of ascending node and argument of periapsis
       * - calc state vectors
       * - get true anomaly
       *
       * \param MeL Mean longitude
       * \param maxRelativeError maximum relative error in eccentric anomaly
       * \param maxIterations max number of iterations for calculating eccentric anomaly
       * \return number of iterations or 0 if failed to find a solution
       *
       */
      StateVectors elements2StateVector (btScalar MeL, btScalar maxRelativeError, int maxIterations);

      /** \brief Refresh the object's orbital elements from state vectors
       *
       * \param state The state vectors to refresh from
       *
       */
      void refreshFromStateVectors (StateVectors* state);

      /** \brief
       *
       * \param elements
       * \param shape
       *
       */
      void elements2Shape (void);

      /** \brief
       *
       * \param mu Gravitationnal parameter
       *
       */
      void setMu(btScalar val);

      /** \brief
       *
       * \return Current elements of the Orbit object
       *
       */
      Elements getElements(void);

      /** \brief
       *
       * \return Current state vectors of the Orbit object
       *
       */
      StateVectors getStateVectors(void);

      /** \brief
       *
       * \return Current parameters of the Orbit object
       *
       */
      Params getParams(void);

      /** \brief
       *
       * \return Current shape of the Orbit object
       *
       */
      OrbitShape getShape(void);

      /** \brief
       * Pseudocode
       * - calc mean motion
       * - calc change in mean anomaly
       * - calc mean anomaly
       * \param timeSinceEpoch Time since epoch in seconds
       * \return A mean anomaly in radians
       *
       */
      btScalar getMeanAnomalyAtTime (btScalar timeSinceEpoch);     /* time since epoch in seconds */

      /** \brief
       * Pseudocode
       * - get mean anomaly
       * - get eccentric anomaly
       * - calc true anomaly
       * \param timeSinceEpoch
       * \param maxRelativeError
       * \param maxIterations
       * \return number of iterations or 0 if failed to find a solution
       *
       */
      btScalar getTrueAnomalyAtTime (
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
      StateVectors elements2StateVectorAtTime (
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
      Elements getElementsAtTime (
        btScalar timeSinceEpoch,         /* time since epoch in seconds */
        btScalar bodyRadius,             /* mean radius of the non-spherical body being orbited */
        btScalar jTwoCoeff);             /* J2 coefficient of the non-spherical body being orbited */
  };
}
#endif // ELEMENTS_H
