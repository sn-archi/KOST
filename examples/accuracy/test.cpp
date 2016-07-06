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
  mKOST::StateVectors out (orbit.elements2StateVector (sv->MeL, SIMD_EPSILON, 1000000));

  /** Compute the difference between the two state vectors and get an error ratio based on the vector length */
  btVector3 diff (sv->pos - out.pos);
  btScalar error (diff.length() / sv->pos.length());

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
  error = (diff.length() / sv->vel.length());

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
}

int main (int argc, char* argv[])
{
  int rx, ry, rz, vx, vy, vz;
  mKOST::StateVectors sv;
  int rmin = 4, rmax = 10, vmin = 0, vmax = 4;
  int counter (0);
  int skipped (0);

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
                          catch (const char * errMsg)
                          {
                            //std::cout << errMsg << std::endl;
                            ++skipped;
                          }
                          ++counter;
                        }
            }
  printf ("maxRerror = %e\n", maxRerror);
  printf ("maxVerror = %e\n", maxVerror);
  printf ("tests done: %i\n", counter);
  printf ("tests skipped: %i\n", skipped);

  return 0;
}