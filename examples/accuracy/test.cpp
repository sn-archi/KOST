#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "../../src/mKOST.h"

mKOST::sStateVector sv_maxVerror, sv_maxRerror;
btScalar maxVerror = -1.0, maxRerror = -1.0;

int testState (const mKOST::sStateVector* sv)
{
  mKOST::sElements elements;
  mKOST::sStateVector out;
  btVector3 diff;
  btScalar error;
  int ret;
  mKOST::Orbit* mOrbit = new mKOST::Orbit;

  /*Convert to orbital elements*/
  if (mOrbit->stateVector2Elements (MU, sv, &elements, NULL))
    return 1;

  /*Convert back to state vector*/
  ret = mOrbit->elements2StateVector (MU, &elements, &out, 10 * SIMD_EPSILON, 1000000);

  diff = sv->pos - out.pos;
  error = diff.length();

  if (error > maxRerror)
    {
      maxRerror = error;
      sv_maxRerror = *sv;

      printf ("New pos error max:\n"
              "     pos  = %e, %e, %e\n"
              "     vel  = %e, %e, %e\n"
              "     error = %f\n",
              sv->pos.getX(), sv->pos.getY(), sv->pos.getZ(),
              sv->vel.getX(), sv->vel.getY(), sv->vel.getZ(),
              maxRerror);
    }

  diff = sv->vel - out.vel;
  error = diff.length();

  if (error > maxVerror)
    {
      maxVerror = error;
      sv_maxVerror = *sv;

      printf ("New vel error max:\n"
              "     pos = %e, %e, %e\n"
              "     vel = %e, %e, %e\n"
              "     error = %f\n",
              sv->pos.getX(), sv->pos.getY(), sv->pos.getZ(),
              sv->vel.getX(), sv->vel.getY(), sv->vel.getZ(),
              maxVerror
             );
    }
  return 0;
}

int main (int argc, char* argv[])
{
  int rx, ry, rz, vx, vy, vz;
  mKOST::sStateVector sv;
  int rmin = 4, rmax = 10, vmin = 0, vmax = 4;
  int counter (0);
  int skipped (0);

  /*Arbitrary 6D positions*/
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
                              if (testState (&sv))
                                ++skipped;
                              ++counter;
                            }
              }
  printf ("maxRerror = %e\n", maxRerror);
  printf ("maxVerror = %e\n", maxVerror);
  printf ("tests done: %i\n", counter);
  printf ("tests skipped: %i\n", skipped);

  return 0;
}
