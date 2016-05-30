#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "../../src/mKOST.h"

/* Number of output samples */
#define N 100

int main (int argc, char* argv[])
{
  mKOST::sStateVector initial, out;
  mKOST::sElements elements, elements2;
  mKOST::sOrbitParam params, params2;
  double maxt;
  int foundroot;

  btVector3 output[N];
  unsigned int i = 0;

  /*Initial state at t=0*/
  initial.pos = btVector3 (-1.0e4, 1.0e8, 1.0e4);
  initial.vel = btVector3 (-1.0, 1.0e4, 1.0);

  /*Convert to orbital elements*/
  mKOST::stateVector2Elements (MU, &initial, &elements, &params);

    printf ("initial:\n"
          "  position: %e, %e, %e\n"
          "  velocity: %e, %e, %e\n-----\n",
          initial.pos.getX(),
          initial.pos.getY(),
          initial.pos.getZ(),
          initial.vel.getX(),
          initial.vel.getY(),
          initial.vel.getZ()
          );

  printf ("Orbital elements:\n"
          "     a = %e m\n"
          "     e = %e\n"
          "     i = %e\n"
          "   LaN = %e\n"
          "   LoP = %e\n"
          "     L = %e\n",
          elements.a, elements.e, elements.i, elements.LaN, elements.LoP, elements.L
         );
  printf ("Additional parameters:\n"
          "   PeD = %e m\n"
          "   ApD = %e m\n"
          "   MnA = %e\n"
          "   TrA = %e\n"
          "   EcA = %e\n"
          "     T = %e s\n"
          "   AgP = %e\n-----\n",
          params.PeD,
          params.ApD,
          params.MnA,
          params.TrA,
          params.EcA,
          params.T,
          params.AgP
         );

  if (!mKOST::elements2StateVector (MU, &elements, &out, SIMD_EPSILON, 1000000))
    printf ("Couldn't find a root\n");
  printf ("reversed:\n"
          "  position: %e, %e, %e\n"
          "  velocity: %e, %e, %e\n-----\n",
          out.pos.getX(),
          out.pos.getY(),
          out.pos.getZ(),
          out.vel.getX(),
          out.vel.getY(),
          out.vel.getZ()
          );

  mKOST::stateVector2Elements (MU, &out, &elements2, &params2);
  printf ("Orbital elements:\n"
          "     a = %e m\n"
          "     e = %e\n"
          "     i = %e\n"
          "   LaN = %e\n"
          "   LoP = %e\n"
          "     L = %e\n",
          elements2.a,
          elements2.e,
          elements2.i,
          elements2.LaN,
          elements2.LoP,
          elements2.L
         );
  printf ("Additional parameters:\n"
          "   PeD = %e m\n"
          "   ApD = %e m\n"
          "   MnA = %e\n"
          "   TrA = %e\n"
          "   EcA = %e\n"
          "     T = %e s\n"
          "   AgP = %e\n",
          params2.PeD,
          params2.ApD,
          params2.MnA,
          params2.TrA,
          params2.EcA,
          params2.T,
          params2.AgP
         );


  maxt = 1.0 * params.T;

  for (i = 0; i < N; i++)
    {
      double t = (i * maxt) / N;

      mKOST::sStateVector stateNow;

      mKOST::elements2StateVectorAtTime (MU, &elements, &stateNow, t, SIMD_EPSILON, 1000000, 0.0, 0.0);

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
