#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "../../src/mKOST.h"

/* Number of output samples */
#define N 500

/* Central body params */
#define R 6378100.0
#define mu 3.986004418e14

int main (int argc, char* argv[])
{
  mKOST::sStateVector initial, out;
  mKOST::sElements elements, elements2;
  mKOST::sOrbitParam params;
  double maxt;

  btVector3 output[N];
  unsigned int i = 0;

  std::cout << std::endl << std::endl << std::endl;

  /*Initial state at t=0*/
  initial.pos = btVector3 (-1.0e8, -1.0e5, 1.0e5);
  initial.vel = btVector3 (1000.0, 1.0, -1.0);

  /*Convert to orbital elements*/
  mKOST::stateVector2Elements (mu, &initial, &elements, &params);
  mKOST::elements2StateVector (mu, &elements, &out, SIMD_EPSILON, 1000000);
  mKOST::stateVector2Elements (mu, &out, &elements2, &params);

  printf ("Orbital elements:\n"
          "     a = %f m\n"
          "     e = %f\n"
          "     i = %f\n"
          " theta = %f\n"
          "omegab = %f\n"
          "     L = %f\n",
          elements.a, elements.e, elements.i, elements.theta, elements.omegab, elements.L
         );

  printf ("Additional parameters:\n"
          "   PeD = %f m\n"
          "   ApD = %f m\n"
          "   MnA = %f\n"
          "   TrA = %f\n"
          "   EcA = %f\n"
          "     T = %f s\n"
          "   AgP = %f\n",
          params.PeD,
          params.ApD,
          params.MnA,
          params.TrA,
          params.EcA,
          params.T,
          params.AgP
         );
  printf ("initial:\n"
          "  position: %f, %f, %f\n"
          "  velocity: %f, %f, %f\n",
          initial.pos.getX(),
          initial.pos.getY(),
          initial.pos.getZ(),
          initial.vel.getX(),
          initial.vel.getY(),
          initial.vel.getZ()
          );
  printf ("reversed:\n"
          "  position: %f, %f, %f\n"
          "  velocity: %f, %f, %f\n",
          out.pos.getX(),
          out.pos.getY(),
          out.pos.getZ(),
          out.vel.getX(),
          out.vel.getY(),
          out.vel.getZ()
          );

  maxt = 1.0 * params.T;

  for (i = 0; i < N; i++)
    {
      double t = (i * maxt) / N;

      mKOST::sStateVector stateNow;

      mKOST::elements2StateVectorAtTime (mu, &elements, &stateNow, t, SIMD_EPSILON, 1000, 0.0, 0.0);

      output[i] = stateNow.pos;
    }

  {
    FILE* fp = fopen ("orbit.dat", "w");

    // GNUPlot uses Y fr depth and Z for height. This'll make things looks as expected
    // Just don't pay attention to the axix names.
    for (i = 0; i < N; i++)
      fprintf (fp, "%f\t%f\t%f\n", output[i].getX(), output[i].getZ(), output[i].getY() );

    fclose (fp);
  }

  return system ("./test.plot");
}
