#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "../../src/kost.h"

/*Number of output samples*/
#define N 100

int main (int argc, char* argv[])
{
  /*Data about central body (earth)*/
  double R = 6378100.0;
  double M = 5.97237e24;
  double mu = 3.986004418e14;

  mKOST::sStateVector initial, state2;
  mKOST::sElements elements;
  mKOST::sOrbitParam params;
  double maxt;

  btVector3 output[N];
  unsigned int i = 0;

  printf ("mu = %e\n", mu);

  /*Initial state at t=0*/
  initial.pos = btVector3 (R + 200000.0, 0.0, 0.0);
  initial.vel = btVector3 (0.0, 0.0, 10000.0);

  /*Convert to orbital elements*/
  mKOST::stateVector2Elements (mu, &initial, &elements, &params);

  printf ("Orbital elements:\n"
          "     a = %f m\n"
          "     e = %f\n"
          "     i = %f°\n"
          " theta = %f°\n"
          "omegab = %f°\n"
          "     L = %f°\n",
          elements.a, elements.e, elements.i*180/M_PI, elements.theta*180/M_PI, elements.omegab*180/M_PI, elements.L*180/M_PI
         );

  printf ("Additional parameters:\n"
          "   PeD = %f m\n"
          "   ApD = %f m\n"
          "   MnA = %f°\n"
          "   TrA = %f°\n"
          "   EcA = %f°\n"
          "     T = %f s\n"
          "   AgP = %f°\n",
          params.PeD,
          params.ApD,
          params.MnA*180/M_PI,
          params.TrA*180/M_PI,
          params.EcA*180/M_PI,
          params.T,
          params.AgP
         );

  maxt = 1.0 * params.T;

  for (i = 0; i < N; i++)
    {
      double t = (i * maxt) / N;

      mKOST::sStateVector stateNow;

      mKOST::elements2StateVectorAtTime (mu, &elements, &stateNow, t, KOST_VERYSMALL, 20, 0.0, 0.0);

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
