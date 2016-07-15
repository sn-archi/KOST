#include <iostream>
#include "../../src/mKOST.h"

/** TODO: Globals are evil */
mKOST::StateVectors sv_maxVerror, sv_maxRerror;
btScalar maxVerror (-1.0), maxRerror (-1.0);

void testState (mKOST::StateVectors* sv)
{
  /** Convert to orbital elements */
  mKOST::Orbit orbit (MU, sv);

  /** Convert back to state vector */
  mKOST::StateVectors out (orbit.elements2StateVector (sv->MeL, EPSILON, 1000000));

  /** Compute the difference between the two state vectors and get an error ratio based on the vector length */
  btVector3 diff (sv->pos - out.pos);
  btScalar error (diff.length() / sv->pos.length());

  if (error > maxRerror)
  {
    maxRerror = error;
    sv_maxRerror = *sv;
    std::cout << "New pos error max:\n" << *sv << "\nerror = " << maxRerror << std::endl;
  }

  diff = sv->vel - out.vel;
  error = (diff.length() / sv->vel.length());

  if (error > maxVerror)
  {
    maxVerror = error;
    sv_maxVerror = *sv;
    std::cout << "New vel error max:\n" << *sv << "\nerror = " << maxRerror << std::endl;
  }
}

int main (int argc, char* argv[])
{
  int rx, ry, rz, vx, vy, vz;
  mKOST::StateVectors sv;
  int rmin = 6, rmax = 12, vmin = 3, vmax = 7;
  int counter (0);
  int nosol (0);
  int mnaisnan (0);
  int hnull (0);
  int other (0);

  /** Arbitrary 6D positions */
  for (rx = rmin; rx <= rmax; ++rx)
    for (rz = rmin; rz <= rmax; ++rz)
      for (ry = rmin; ry <= rmax; ++ry)
        for (vx = vmin; vx <= vmax; ++vx)
          for (vz = vmin; vz <= vmax; ++vz)
            for (vy = vmin; vy <= vmax; ++vy)
            {
              double rxf = std::pow (10.0, rx);
              double ryf = std::pow (10.0, ry);
              double rzf = std::pow (10.0, rz);
              double vxf = std::pow (10.0, vx);
              double vyf = std::pow (10.0, vy);
              double vzf = std::pow (10.0, vz);

              int srx, sry, srz, svx, svy, svz;
              for (srx = -1; srx <= 1; ++srx)
                for (sry = -1; sry <= 1; ++sry)
                  for (srz = -1; srz <= 1; ++srz)
                    for (svx = -1; svx <= 1; ++svx)
                      for (svy = -1; svy <= 1; ++svy)
                        for (svz = -1; svz <= 1; ++svz)
                        {
                          if ((srx == 0) || (sry == 0) || (srz == 0))
                            continue;
                          if ((svx == 0) || (svy == 0) || (svz == 0))
                            continue;

                          sv.pos = btVector3 (srx * rxf, sry * ryf, srz * rzf);
                          sv.vel = btVector3 (svx * vxf, svy * vyf, svz * vzf);
                          try
                          {
                            testState (&sv);
                          }
                          catch (const int errMsg)
                          {
                            if (errMsg == 1)
                              ++nosol;
                            if (errMsg == 2)
                              ++mnaisnan;
                            if (errMsg == 3)
                              ++hnull;
                            else ++other;

                          }
                          ++counter;
                        }
            }
  std::cout << "maxRerror = " << maxRerror << std::endl;
  std::cout << "maxVerror = " << maxVerror << std::endl;
  std::cout << "tests done: " << counter << std::endl;
  std::cout << "No acceptable solution found: " << nosol << std::endl;
  std::cout << "Fucked up Mean anomaly estimate: " << mnaisnan << std::endl;
  std::cout << "Null angular momentum: " << hnull << std::endl;
  std::cout << "Other issues: " << other << std::endl;

  return 0;
}