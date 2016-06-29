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

#include "Orbit.h"

namespace mKOST
{
  Orbit::Orbit (void)
    : mElements  {Elements()},
      mParams    {Params()},
      mShape     {OrbitShape()},
      n          {btVector3()},
      e          {btVector3()}
  {
  }

  Orbit::Orbit (btScalar mu, StateVectors* state)
    : mElements  {Elements()},
      mParams    {Params()},
      mShape     {OrbitShape()},
      n          {btVector3()},
      e          {btVector3()},
      mu         {mu}
  {
    refreshFromStateVectors(state);
  }

  Orbit::Orbit (btScalar mu, Elements elements)
    : mElements  {elements},
      mParams    {Params()},
      mShape     {OrbitShape()},
      n          {btVector3()},
      e          {btVector3()},
      mu         {mu}
  {
    /** parabolic orbit - approximate to hyperbolic orbit */
    if (mElements.Ecc == 1.0)
    {
      mElements.Ecc += SIMD_EPSILON;
    }
    calcN();
    refreshParams();
  }

  Orbit::~Orbit ()
  {
  }

  void Orbit::refreshFromStateVectors (StateVectors* state)
  {
    btVector3 h (calcH(state));

    /** If velocity and position vectors are parallel and aligned, throw an ecception */
    if (h.isZero())
      throw "Angular momentum is null !";

    calcN(&h);
    if (n.length() < SIMD_EPSILON)
    {
      n = btVector3 (1.0, 0.0, 0.0);
      Equatorial = true;
    }
    else Equatorial = false;

    btScalar E (state->vel.length2() / 2 - mu / state->pos.length());
    if (E == 0.0)
      E = SIMD_EPSILON;

    /** Set the e vector from state vectors */
    calcE(state);

    mElements.Ecc = e.length();

    /** parabolic orbit - approximate to hyperbolic orbit */
    if (mElements.Ecc == 1.0)
    {
      mElements.Ecc += SIMD_EPSILON;
    }

    Circular = mElements.Ecc < SIMD_EPSILON;
    Hyperbola = mElements.Ecc >= 1.0;

    /** parabolic orbit are not supported */
    if (mElements.Ecc > 1.0 - SIMD_EPSILON && mElements.Ecc < 1.0 + SIMD_EPSILON)
    {
      if (E >= 0.0)
      {
        mElements.Ecc = 1.0 + SIMD_EPSILON; /* Approximate with hyperbolic */
      }
      else
      {
        mElements.Ecc = 1.0 - SIMD_EPSILON; /* Approximate with elliptic */
      }
    }

    /*
      SMa
      dp = a(1 - e)
      da = a(1 + e)
    */
    mElements.a = -mu / (2.0 * E);

    mParams.PeD = calcPeD(&h);
    mParams.ApD = calcApD();

    /** Inc */
    mElements.i = std::acos (h.getY() / h.length());

    /** Longitude of ascending node */
    mElements.LAN = calcLAN();

    /** Argument of Periapsis */
    mParams.AgP = calcAgP(&h);

    /** EcA */
    mParams.EcA = calcEcA(state);

    /** MnA */
    mParams.MnA = calcMnA();

    /** TrA */
    mParams.TrA = calcTrA(mParams.EcA);

    /** Lec */
    mParams.Lec = mElements.a * mElements.Ecc;

    /** SMi
     * b² = a²(1 - e²) */
    mParams.SMi = calcSMi();

    /** LoP */
    mElements.LoP = std::fmod (mElements.LAN + mParams.AgP, SIMD_2_PI);

    /** Mean Longitude at epoch*/
    state->MeL = mParams.MnA + mElements.LoP;
    if (!Hyperbola) state->MeL = std::fmod (state->MeL, SIMD_2_PI);

    /** TrL */
    mParams.TrL = std::fmod (mElements.LoP + mParams.TrA, SIMD_2_PI);

    /** Orbital period */
    mParams.T = calcT();

    /** Calculating PeT and ApT: */
    btScalar tPe (mParams.MnA * mParams.T / SIMD_2_PI); /*Time since last Pe*/

    if (Hyperbola)
      {
        mParams.PeT = -tPe;
      }
    else
      {
        mParams.PeT = mParams.T - tPe;
      }

    mParams.ApT = 0.5 * mParams.T - tPe;
    if (mParams.ApT < 0.0) mParams.ApT += mParams.T;
  }

  StateVectors Orbit::elements2StateVector (btScalar trueAnomaly) const
  {
    btVector3 h (calcH());

    /** calc position vector */
    /** argument of position, measured from the longitude of ascending node */
    btScalar argPos (mElements.LoP - mElements.LAN + trueAnomaly);

    /** calc direction of position vector:
     *  r/|r| = sin(ArgPos) * ((h / |h|) x n) + cos(argPos) * n */
    StateVectors state;
    state.pos = std::sin(argPos) * h.normalized().cross(n) + std::cos(argPos) * n;

    /** calc length of position vector */
    if (mElements.Ecc < 1.0)                  /** elliptical orbit */
      state.pos = mElements.a * (1.0 - mElements.Ecc * mElements.Ecc) / (1.0 + mElements.Ecc * std::cos(trueAnomaly)) * state.pos;
    else                          /** hyperbolic orbit */
      state.pos = std::fabs (mElements.a) * (mElements.Ecc * mElements.Ecc - 1.0) / (1.0 + mElements.Ecc * std::cos(trueAnomaly)) * state.pos;

    /** calc velocity vector */
    /** calculate squared magnitude of velocity vector */
    btScalar v2;
    if (mElements.Ecc < 1.0)   /** elliptical orbit */
      {
        v2 = mu * (2.0 / state.pos.length() - 1.0 / mElements.a);
      }
    else   /** hyperbolic orbit */
      {
        v2 = mu * (2.0 / state.pos.length() + 1.0 / std::fabs (mElements.a));
      }

    /** calc components of velocity vector perpendicular and parallel to radius vector:
    - perpendicular:
      vPro = (|h|/|pos|) * normal(h x pos)
    - parallel:
      vO = √(v² - |vPro|²) * sign(sin(trueAnomaly)) * normal(pos) */
    btVector3 vPro (h.length() / state.pos.length() * h.cross (state.pos).normalized());

    btScalar tmpr (std::sin(trueAnomaly));
    btVector3 vO (btVector3(0.0, 0.0, 0.0));
    if (tmpr != 0.0) /** check for apsis condition to avoid divide by zero */
      {
        btScalar signSinTrueAnomaly (tmpr / std::fabs (tmpr));

        btScalar v0_sq (v2 - vPro.length2());
        /** check for small negative numbers resulting from rounding */
        if (v0_sq < 0.0) v0_sq = 0.0;

        vO = state.pos.normalized();
        vO *= (std::sqrt (v0_sq) * signSinTrueAnomaly);
      }

    /** add to get velocity vector */
    state.vel = vPro + vO;

    return state;
  }

  StateVectors Orbit::elements2StateVector (btScalar MeL, btScalar maxRelativeError, int maxIterations) const
  {
    btScalar trueAnomaly (calcTrA (MeL, maxRelativeError, maxIterations));
    return elements2StateVector (trueAnomaly);
  }

  void Orbit::elements2Shape ()
  {
    /** Some utility values: */
    btScalar multiplier (mElements.a * (1.0 - mElements.Ecc * mElements.Ecc));

      /** First: Orbit in its own coordinate system: */
    /** periapsis and apoapsis */
    mShape.pe = btVector3 ( mElements.a * (1.0 - mElements.Ecc), 0.0, 0.0);
    mShape.ap = btVector3 (-mElements.a * (1.0 + mElements.Ecc), 0.0, 0.0);

    /** Points */
    if (mShape.numPoints == 1)
    {
      mShape.points[0] = mShape.pe;
    }
    else if (mShape.numPoints > 1)
    {
      btScalar maxTrA, dTrA, currentTrA;

      /** Range of angles */
      maxTrA = SIMD_PI;
      if (mElements.Ecc >= 1.0)
        {
          maxTrA = std::acos (-1.0 / mElements.Ecc);

          /*Make it a bit smaller to avoid division by zero:*/
          maxTrA *= (((btScalar) mShape.numPoints) / (mShape.numPoints + 1));
        }

      /** Angle change per segment */
      dTrA = (2 * maxTrA) / (mShape.numPoints - 1);

      currentTrA = -maxTrA;
      for (unsigned int i (0); i < mShape.numPoints; ++i)
      {
        btScalar absr (std::fabs (multiplier / (1.0 + mElements.Ecc * std::cos (currentTrA))));
        btVector3 direction (btVector3 (std::cos (currentTrA), 0.0, -std::sin (currentTrA)));
        mShape.points[i] = absr * direction;
        currentTrA += dTrA;
      }
    }

    /** AN */
    {
      btScalar currentTrA (-mParams.AgP);
      btScalar absr (multiplier / (1.0 + mElements.Ecc * std::cos (currentTrA)));

      if (absr <= 0.0)
      {
        mShape.an = btVector3 (0.0, 0.0, 0.0);
      }
      else
      {
        btVector3 direction (btVector3 (std::cos (currentTrA), 0.0, -std::sin (currentTrA)));
        mShape.an = absr * direction;
      }
    }

    /** DN */
    {
      btScalar currentTrA (SIMD_PI - mParams.AgP);
      btScalar absr (multiplier / (1.0 + mElements.Ecc * std::cos (currentTrA)));

      if (absr <= 0.0)
      {
        mShape.dn = btVector3 (0.0, 0.0, 0.0);
      }
      else
      {
        btVector3 direction (btVector3 (std::cos (currentTrA), std::sin (currentTrA), 0.0));
        mShape.dn = absr * direction;
      }
    }

      /** Then: rotate the coordinates: */
    {
      btMatrix3x3 AgPMat, LANMat, IncMat, transform;

      AgPMat.setEulerZYX (0.0, mParams.AgP, 0.0);
      IncMat.setEulerZYX (mElements.i, 0.0, 0.0);
      LANMat.setEulerZYX (0.0, mElements.LAN, 0.0);

      /* Now, global = LANMat * IncMat * AgPMat * local: */
      transform = LANMat * IncMat;
      transform *= AgPMat;

      mShape.pe = transform * mShape.pe;
      mShape.ap = transform * mShape.ap;
      mShape.an = transform * mShape.an;
      mShape.dn = transform * mShape.dn;

      if (mShape.numPoints != 0)
        for (unsigned int i (0); i < mShape.numPoints; ++i)
          mShape.points[i] = transform * mShape.points[i];
    }
  }

  btScalar Orbit::getMeanAnomalyAtTime (btScalar MeL, btScalar timeSinceEpoch) const
  {
    /* calc mean motion */
    btScalar meanMotion (std::sqrt (mu / std::pow (std::fabs (mElements.a), 3.0)));

    /* calc change in mean anomaly */
    btScalar deltaMeanAnomaly (timeSinceEpoch * meanMotion);

    /* calc mean anomaly */
    /* fmod takes care of overflow in the case where the mean anomaly exceeds one revolution */
    btScalar meanAnomaly (std::fmod (MeL - mElements.LoP + deltaMeanAnomaly, SIMD_2_PI));

    if (meanAnomaly < 0) meanAnomaly += SIMD_2_PI;

    return meanAnomaly;
  }

  btScalar Orbit::getTrueAnomalyAtTime (
    btScalar MeL,
    btScalar timeSinceEpoch,
    btScalar maxRelativeError,
    int maxIterations)
  {
    int ret;
    btScalar meanAnomaly, eccentricAnomaly;

    /* get mean anomaly */
    meanAnomaly = getMeanAnomalyAtTime (MeL, timeSinceEpoch);

    /* get eccentric anomaly */
    eccentricAnomaly = calcEcA (meanAnomaly, meanAnomaly, maxRelativeError, maxIterations);

    /* calc true anomaly */
    return calcTrA (eccentricAnomaly);
  }

  btScalar Orbit::getLANAtTime (
    btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
    btScalar jTwoCoeff,           /* J2 coefficient of the non-spherical body being orbited */
    btScalar timeSinceEpoch)      /* time since epoch in seconds */
  {
    btScalar meanMotion, LANAtTime;

    if (mElements.Ecc < 1.0) /* elliptical orbit */
      {
        meanMotion = std::sqrt (mu / pow (mElements.a, 3.0) );
        return ( mElements.LAN + timeSinceEpoch * (-3.0 * meanMotion / 2.0) * pow (bodyRadius / mElements.a, 2.0) * (cos (mElements.i) / pow (1.0 - pow (mElements.Ecc, 2.0), 2.0) ) * jTwoCoeff );
      }
    else /* hyperbolic orbit - non spherical effect is negligible */
      return mElements.LAN;
  }

  btScalar Orbit::getArgPeAtTime (
    btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
    btScalar jTwoCoeff,           /* J2 coefficient of the non-spherical body being orbited */
    btScalar timeSinceEpoch)      /* time since epoch in seconds */
  {
    btScalar meanMotion;

    if (mElements.Ecc < 1.0) /* elliptical orbit */
      {
        meanMotion = std::sqrt (mu / std::pow (mElements.a, 3.0) );
        return ( mElements.LoP - mElements.LAN + timeSinceEpoch * (3.0 * meanMotion / 4.0) * std::pow (bodyRadius / mElements.a, 2.0) * ( (5.0 * std::pow (std::cos (mElements.i),
                 2.0) - 1.0) / std::pow (1.0 - std::pow (mElements.Ecc, 2.0), 2.0) ) * jTwoCoeff );
      }
    else /* hyperbolic orbit - non spherical effect is negligible */
      return ( mElements.LoP - mElements.LAN );
  }

  StateVectors Orbit::elements2StateVectorAtTime (
    btScalar MeL,
    btScalar timeSinceEpoch,      /* time since epoch in seconds */
    btScalar maxRelativeError,    /* maximum relative error in eccentric anomaly */
    int maxIterations,            /* max number of iterations for calculating eccentric anomaly */
    btScalar bodyRadius,          /* mean radius of the non-spherical body being orbited */
    btScalar jTwoCoeff)           /* J2 coefficient of the non-spherical body being orbited */
  {
    /* Returns number of iterations if successful, returns 0 otherwise. */

    /* Pseudocode
     *
     * get true anomaly
     * get longitude of ascending node and argument of periapsis
     * calc state vectors */

//    Elements updatedElements;

    /* get true anomaly */
    btScalar trueAnomaly (getTrueAnomalyAtTime (MeL, timeSinceEpoch, maxRelativeError, maxIterations));

    /* update elements for new epoch */
//    Elements elementsAtTime (getElementsAtTime (timeSinceEpoch, bodyRadius, jTwoCoeff));

    /* calc state vectors */
    return elements2StateVector (trueAnomaly);
  }

  void Orbit::calcE (StateVectors* state)
  {
    btVector3 h (calcH(state));
    e = state->vel.cross(h) / mu - state->pos / state->pos.length();
  }

  btVector3 Orbit::calcH (StateVectors* state) const
  {
    return (state->pos.cross (state->vel));
  }

  btVector3 Orbit::calcH() const
  {
    btVector3 h (n.cross (north) * std::sin(mElements.i));
    h.setY (std::cos(mElements.i));
    h.normalize();

    /** calc radius and velocity at periapsis */
    btScalar rPe, vPe;
    if (mElements.Ecc < 1.0)   /** elliptical orbit */
      {
        rPe = mElements.a * (1.0 - std::pow(mElements.Ecc,2)) / (1.0 + mElements.Ecc);
        vPe = std::sqrt (mu * (2.0 / rPe - 1.0 / mElements.a));
      }
    else   /** hyperbolic orbit */
      {
        rPe = std::fabs (mElements.a) * (mElements.Ecc * mElements.Ecc - 1.0) / (1.0 + mElements.Ecc);
        vPe = std::sqrt (mu * (2.0 / rPe + 1.0 / std::fabs (mElements.a)));
      }
    /** calc h */
    h *= rPe * vPe;

    return h;
  }

  void Orbit::calcN(btVector3 *h)
  {
    n = btVector3(h->getZ(), 0.0, -h->getX());
    n.normalize();
  }

  void Orbit::calcN ()
  {
    n = btVector3(std::cos (mElements.LAN), 0.0, -std::sin (mElements.LAN));
  }

  btScalar Orbit::calcT()
  {
    /** fabs is for supporting hyperbola */
    return SIMD_2_PI * std::sqrt (std::fabs(std::pow(mElements.a, 3) / mu));
  }

  btScalar Orbit::calcSMi() const
  {
    if (Hyperbola) return std::sqrt (mElements.a * mElements.a * (mElements.Ecc * mElements.Ecc - 1.0));
    else return std::sqrt (mElements.a * mElements.a * (1.0 - mElements.Ecc * mElements.Ecc));
  }

  btScalar Orbit::calcPeD(btVector3 *h) const
  {
    if (Hyperbola)
    {
      return h->length2() / mu;
    }
    else
    {
      return mElements.a * (1.0 - mElements.Ecc);
    }
  }

  btScalar Orbit::calcApD() const
  {
    if (Hyperbola)
    {
      return BT_INFINITY;
    }
    else
    {
      return mElements.a * (1.0 + mElements.Ecc);
    }
  }

  btScalar Orbit::calcLAN() const
  {
    if (n.length() < SIMD_EPSILON)
    {
      return 0.0;
    }
    else
    {
      btScalar LAN (std::acos (n.getX() / n.length()));
      return (n.getZ() > 0.0)?SIMD_2_PI - LAN:LAN;
    }
  }

  btScalar Orbit::calcAgP (btVector3 *h) const
  {
    if (Circular)
    {
      return 0.0;
    }
    else if (Equatorial)
    {
      btScalar AgP (atan2 (-e.getZ(), e.getX()));
      return (h->getY() < 0.0)?-AgP:AgP;
    }
    else
    {
      btScalar AgP (std::acos (n.dot (e) / (n.length() * mElements.Ecc)));
      return (e.getY() < 0.0)?SIMD_2_PI - AgP:AgP;
    }
  }

  btScalar Orbit::calcTrA(StateVectors* state) const
  {
    btScalar trueAnomaly;
    if (Circular)
    {
      if (Equatorial)
      {
        trueAnomaly = std::acos (state->pos.getX() / state->pos.length());
        if (state->vel.getX() > 0.0) trueAnomaly = SIMD_2_PI - trueAnomaly;
      }
      else
      {
        trueAnomaly = std::acos (n.dot (state->pos)) / (n.length() * state->pos.length());
        if (n.dot (state->vel) > 0.0) trueAnomaly = SIMD_2_PI - trueAnomaly;
      }
    }
    else
    {
      btScalar tmp = e.dot (state->pos) / (mElements.Ecc * state->pos.length());

      /* Avoid acos out of range. 1 and -1 are included cause we know the result. */
      if (tmp <= -1.0) trueAnomaly = SIMD_PI;
      else if (tmp >= 1.0) trueAnomaly = 0.0;
      else trueAnomaly = std::acos (tmp);

      /* Changing sin_TrA_isNegative on a negative epsilon value doesn't seem right */
      if (state->pos.dot (state->vel) < 0.0) trueAnomaly = SIMD_2_PI - trueAnomaly;
    }
    return trueAnomaly;
  }

  btScalar Orbit::calcTrA (btScalar MeL, btScalar maxRelativeError, int maxIterations) const
  {
    btScalar meanAnomaly, eccentricAnomaly;

    /* get mean anomaly */
    meanAnomaly = calcMnA(MeL);

    /* get eccentric anomaly */
    if (mElements.Ecc < 1.0)
    {
      eccentricAnomaly = calcEcA (meanAnomaly,
               meanAnomaly,
               maxRelativeError,
               maxIterations);
    }
    else
    {
      eccentricAnomaly = calcEcA (meanAnomaly,
               std::log (2.0 * meanAnomaly / mElements.Ecc + 1.8),
               maxRelativeError,
               maxIterations);
    }

    /* calc true anomaly */
    return calcTrA (eccentricAnomaly);
  }

  btScalar Orbit::calcTrA (btScalar eccentricAnomaly) const
  {
    btScalar ret;
    if (mElements.Ecc < 1.0)   /* elliptical orbit */
    {
      ret = 2.0 * std::atan (std::sqrt ( (1.0 + mElements.Ecc) / (1.0 - mElements.Ecc) ) * std::tan (eccentricAnomaly / 2.0) );
    }
    else   /* hyperbolic orbit */
    {
      ret = std::acos ( (std::cosh (eccentricAnomaly) - mElements.Ecc) / (1 - mElements.Ecc * std::cosh (eccentricAnomaly) ) );
      if (eccentricAnomaly < 0.0) ret = -ret; /* Always the same sign */
    }
    return ret;
  }

  btScalar Orbit::calcEcA (btScalar meanAnomaly, btScalar ecaEstimate, btScalar maxRelativeError, int maxIterations) const
  {
    if (Circular)
      {
        return meanAnomaly;
      }

    btScalar relativeError, mnaEstimate;

    if (meanAnomaly == 0.0) return 0.0;
    else
    {
      int i (0);
      do
        {
          if (!Circular && !Hyperbola)   /** elliptical orbit */
            {
              /** calculate next estimate of the root of Kepler's equation using Newton's method */
              ecaEstimate = ecaEstimate - (ecaEstimate - mElements.Ecc * std::sin (ecaEstimate) - meanAnomaly) / (1.0 - mElements.Ecc * std::cos (ecaEstimate));
              /** calculate estimate of mean anomaly from estimate of eccentric anomaly */
              mnaEstimate = ecaEstimate - mElements.Ecc * std::sin (ecaEstimate);
            }
          else   /** hyperbolic orbit */
            {
              /** calculate next estimate of the root of Kepler's equation using Newton's method */
              ecaEstimate = ecaEstimate - (mElements.Ecc * std::sinh (ecaEstimate) - ecaEstimate - meanAnomaly) / (mElements.Ecc * std::cosh (ecaEstimate) - 1.0);
              /** calculate estimate of mean anomaly from estimate of eccentric anomaly */
              mnaEstimate = mElements.Ecc * std::sinh (ecaEstimate) - ecaEstimate;
            }
          /** calculate relativeError */
          relativeError = 1.0 - mnaEstimate / meanAnomaly;
          ++i;
        }
      while ( (i < maxIterations) && (std::fabs (relativeError) > std::fabs (maxRelativeError) ) );
    }

    if (!Hyperbola)
      {
        /** check range is in 0 to 2π */
        ecaEstimate = std::fmod (ecaEstimate, SIMD_2_PI);
        if (ecaEstimate < 0.0) ecaEstimate += SIMD_2_PI;
        if (ecaEstimate >= SIMD_2_PI) ecaEstimate -= SIMD_2_PI;
      }

    if ( std::fabs (relativeError) <= std::fabs (maxRelativeError) )
      return ecaEstimate;
    else
      throw "No acceptable solution found";
  }

  btScalar Orbit::calcEcA(StateVectors* state) const
  {
    btScalar eccentricAnomaly (0.0);
    btScalar cosEccentricAnomaly ((1.0 - state->pos.length() / mElements.a) / mElements.Ecc);
    if (Hyperbola && cosEccentricAnomaly > 1.0)
    {
      eccentricAnomaly = std::acosh (cosEccentricAnomaly);
    }
    else
    {
      /* Avoid acos out of range. 1 and -1 are included cause we now the result. */
      if (cosEccentricAnomaly <= -1.0)
      {
        eccentricAnomaly = SIMD_PI;
      }
      else if (cosEccentricAnomaly >= 1.0)
      {
        eccentricAnomaly = 0.0;
      }
      else
      {
        eccentricAnomaly = std::acos (cosEccentricAnomaly);
      }
    }

    /** Copy sign from sin(TrA) */
    if (Hyperbola && ((std::sin(mParams.TrA) < 0.0) != (eccentricAnomaly < 0.0))) eccentricAnomaly = -eccentricAnomaly;
    /** Same rule basically, but with eccentricAnomaly in 0..2π range */
    else if ((std::sin(mParams.TrA) < 0.0)) eccentricAnomaly = SIMD_2_PI - eccentricAnomaly;
    return eccentricAnomaly;
  }

  btScalar Orbit::calcMnA (btScalar MeL) const
  {
    btScalar meanAnomaly (MeL - mElements.LoP);

    if (mElements.Ecc < 1.0)
    {
      /* check range is in 0 to 2π */
      if (meanAnomaly < 0.0) meanAnomaly += SIMD_2_PI;
      if (meanAnomaly >= SIMD_2_PI) meanAnomaly -= SIMD_2_PI;
    }
    return meanAnomaly;
  }

  btScalar Orbit::calcMnA () const
  {
    if (Hyperbola)
    {
      btScalar meanAnomaly (mElements.Ecc * std::sinh(mParams.EcA) - mParams.EcA);
      return meanAnomaly;
    }
    else
    {
      btScalar meanAnomaly (mParams.EcA - mElements.Ecc * std::sin(mParams.EcA));
      return meanAnomaly;
    }
  }

  void Orbit::refreshParams()
  {
    btVector3 h (calcH());
    Equatorial = (mElements.i == 0.0 || mElements.i == SIMD_PI);
    Circular = mElements.Ecc == 0.0;
    Hyperbola = mElements.Ecc >= 1.0;
    mParams.SMi = calcSMi();
    mParams.PeD = calcPeD(&h);
    mParams.ApD = calcApD();
    mParams.Lec = mElements.a * mElements.Ecc;
    mParams.T = calcT();
  }
}
