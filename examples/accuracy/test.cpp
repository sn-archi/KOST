#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "../../src/kost.h"

/*Data about central body (earth)*/
#define R 6378100.0
#define mu 3.986004418e14


mKOST::sStateVector sv_maxVerror, sv_maxRerror;
btScalar maxVerror = -1.0, maxRerror = -1.0;

void testState (const mKOST::sStateVector* sv)
{
  mKOST::sElements elements;
  mKOST::sStateVector out;
  btVector3 diff;
  btScalar error;

  /*Convert to orbital elements*/
  mKOST::stateVector2Elements (mu, sv, &elements, NULL);

  /*Convert back to state vector*/
  mKOST::elements2StateVector (mu, &elements, &out, KOST_VERYSMALL, 1000);

  diff = sv->pos - out.pos;
  error = (diff.length() / sv->pos.length());
  //if (error > maxRerror && btFuzzyZero(elements.e - 1.0))
  if (error > maxRerror)
    {
      maxRerror = error;
      sv_maxRerror = *sv;

      printf ("New pos error max:\n"
              "     pos  = %e, %e, %e\n"
              "     diff = %e, %e, %e\n"
              "     vel  = %e, %e, %e\n"
              "     error = %f\n",
              sv_maxRerror.pos.getX(), sv_maxRerror.pos.getY(), sv_maxRerror.pos.getZ(),
              diff.getX(), diff.getY(), diff.getZ(),
              sv_maxRerror.vel.getX(), sv_maxRerror.vel.getY(), sv_maxRerror.vel.getZ(),
              maxRerror);
    }

  diff = sv->vel - out.vel;
  error = (diff.length() / sv->vel.length());
  //if (error > maxVerror && btFuzzyZero(elements.e - 1.0))
  if (error > maxVerror)
    {
      maxVerror = error;
      sv_maxVerror = *sv;

      printf ("New vel error max:\n"
              "     pos = %e, %e, %e\n"
              "     vel = %e, %e, %e\n"
              "     diff = %e, %e, %e\n"
              "     error = %f\n",
              sv_maxVerror.pos.getX(), sv_maxVerror.pos.getY(), sv_maxVerror.pos.getZ(),
              sv_maxVerror.vel.getX(), sv_maxVerror.vel.getY(), sv_maxVerror.vel.getZ(),
              diff.getX(), diff.getY(), diff.getZ(),
              maxVerror
             );
    }
}

int main (int argc, char* argv[])
{
  int rx, ry, rz, vx, vy, vz;
  mKOST::sStateVector sv;
  int rmin = 5, rmax = 12, vmin = 0, vmax = 4;
  int counter = 0;

  /*Arbitrary 6D positions*/
  for (rx = rmin; rx <= rmax; rx++)
    for (rz = rmin; rz <= rmax; rz++)
      for (ry = rmin; ry <= rmax; ry++)
        for (vx = vmin; vx <= vmax; vx++)
          for (vz = vmin; vz <= vmax; vz++)
            for (vy = vmin; vy <= vmax; vy++)
              {
                double rxf = std::pow (10.0, rx);
                double ryf = std::pow (10.0, ry);
                double rzf = std::pow (10.0, rz);
                double vxf = std::pow (10.0, vx);
                double vyf = std::pow (10.0, vy);
                double vzf = std::pow (10.0, vz);

                int srx, sry, srz, svx, svy, svz;
                for (srx = -1; srx <= 1; srx++)
                  for (sry = -1; sry <= 1; sry++)
                    for (srz = -1; srz <= 1; srz++)
                      for (svx = -1; svx <= 1; svx++)
                        for (svy = -1; svy <= 1; svy++)
                          for (svz = -1; svz <= 1; svz++)
                            {
                              if ((srx == 0) || (sry == 0) || (srz == 0))
                                continue;
                              if ((svx == 0) || (svy == 0) || (svz == 0))
                                continue;

                              sv.pos = btVector3 (srx * rxf, sry * ryf, srz * rzf);
                              sv.vel = btVector3 (svx * vxf, svy * vyf, svz * vzf);
                              testState (&sv);
                              ++counter;
                            }
              }
  printf ("maxRerror = %e\n", maxRerror);
  printf ("maxVerror = %e\n", maxVerror);
  printf ("tests done: %i\n", counter);

  return 0;
}
