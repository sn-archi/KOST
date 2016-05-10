#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "../../src/kost.h"

/*Data about central body (earth)*/
#define R 6378100.0
#define M 5.9736e24
#define mu (KOST_GRAVITATIONAL_CONSTANT * M)

kostStateVector sv_maxVerror, sv_maxRerror;
btScalar maxVerror = -1.0, maxRerror = -1.0;

void testState (const kostStateVector* sv)
{
  kostElements elements;
  kostStateVector out;
  btVector3 diff;
  btScalar error;

  /*Convert to orbital elements*/
  kostStateVector2Elements (mu, sv, &elements, NULL);

  /*Convert back to state vector*/
  kostElements2StateVector (mu, &elements, &out, KOST_VERYSMALL, 1000000);

  diff = sv->pos - out.pos;
  error = diff.length() / sv->pos.length();
  if (error > maxRerror)
    {
      maxRerror = error;
      sv_maxRerror = *sv;

      printf ("New pos error max:\n"
              "     pos = %e, %e, %e\n"
              "     vel = %e, %e, %e\n"
              "     error = %e\n",
              sv_maxRerror.pos.getX(), sv_maxRerror.pos.getY(), sv_maxRerror.pos.getZ(),
              sv_maxRerror.vel.getX(), sv_maxRerror.vel.getY(), sv_maxRerror.vel.getZ(),
              maxRerror
             );
    }

  diff = sv->vel - out.vel;
  error = diff.length() / sv->vel.length();
  if (error > maxVerror)
    {
      maxVerror = error;
      sv_maxVerror = *sv;

      printf ("New vel error max:\n"
              "     pos = %e, %e, %e\n"
              "     vel = %e, %e, %e\n"
              "     error = %e\n",
              sv_maxVerror.pos.getX(), sv_maxVerror.pos.getY(), sv_maxVerror.pos.getZ(),
              sv_maxVerror.vel.getX(), sv_maxVerror.vel.getY(), sv_maxVerror.vel.getZ(),
              maxVerror
             );
    }
}

int main (int argc, char* argv[])
{
  int rx, ry, rz, vx, vy, vz;
  kostStateVector sv;
  int rmin = 5, rmax = 12, vmin = 0, vmax = 4;

  /*Arbitrary 6D positions*/
  for (rx = rmin; rx <= rmax; rx++)
    for (ry = rmin; ry <= rmax; ry++)
      for (rz = rmin; rz <= rmax; rz++)
        for (vx = vmin; vx <= vmax; vx++)
          for (vy = vmin; vy <= vmax; vy++)
            for (vz = vmin; vz <= vmax; vz++)
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
                              if (srx == sry == srz == 0)
                                continue;
                              if (svx == svy == svz == 0)
                                continue;

                              sv.pos = btVector3 (srx * rxf, sry * ryf, srz * rzf);
                              sv.vel = btVector3 (svx * vxf, svy * vyf, svz * vzf);
                              testState (&sv);
                            }
              }

  printf ("maxRerror = %e\n", maxRerror);
  printf ("maxVerror = %e\n", maxVerror);

  return 0;
}
