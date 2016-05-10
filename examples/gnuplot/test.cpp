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
  double M = 5.9736e24;
  double mu = KOST_GRAVITATIONAL_CONSTANT * M;

  mKOST::kostStateVector initial, state2;
  mKOST::kostElements elements;
  mKOST::kostOrbitParam params;
  double maxt;

  btVector3 output[N];
  unsigned int i = 0;

  printf ("mu = %e\n", mu);

  /*Initial state at t=0*/
  initial.pos = btVector3 (R + 200000.0, 0.0, 0.0);
  initial.vel = btVector3 (0.0, 10000.0, 0.0);

  /*Convert to orbital elements*/
  mKOST::stateVector2Elements (mu, &initial, &elements, &params);

  printf ("Orbital elements:\n"
          "     a = %e\n"
          "     e = %e\n"
          "     i = %e\n"
          " theta = %e\n"
          "omegab = %e\n"
          "     L = %e\n",
          elements.a, elements.e, elements.i, elements.theta, elements.omegab, elements.L
         );

  printf ("Additional parameters:\n"
          "   PeD = %e\n"
          "   ApD = %e\n"
          "   MnA = %e\n"
          "   TrA = %e\n"
          "   EcA = %e\n"
          "     T = %e\n"
          "   AgP = %e\n",
          params.PeD,
          params.ApD,
          params.MnA,
          params.TrA,
          params.EcA,
          params.T,
          params.AgP
         );

  /*
  kostElements2StateVector(mu, &elements, &state2, KOST_VERYSMALL, 1000000);

  printf("Returned state:\n"
  	"     pos = %e, %e, %e\n"
  	"     vel = %e, %e, %e\n",
  	state2.pos.x, state2.pos.y, state2.pos.z,
  	state2.vel.x, state2.vel.y, state2.vel.z
  );

  printf("Returned state error:\n"
  	"     pos = %e, %e, %e\n"
  	"     vel = %e, %e, %e\n",
  	state2.pos.x - initial.pos.x, state2.pos.y - initial.pos.y, state2.pos.z - initial.pos.z,
  	state2.vel.x - initial.vel.x, state2.vel.y - initial.vel.y, state2.vel.z - initial.vel.z
  );

  return 0;
  */

  maxt = 1.0 * params.T;

  for (i = 0; i < N; i++)
    {
      double t = (i * maxt) / N;

      mKOST::kostStateVector stateNow;

      mKOST::elements2StateVectorAtTime (mu, &elements, &stateNow, t, KOST_VERYSMALL, 20, 0.0, 0.0);

      output[i] = stateNow.pos;
    }

  {
    FILE* fp = fopen ("orbit.dat", "w");

    for (i = 0; i < N; i++)
      fprintf (fp, "%f\t%f\t%f\n", output[i].getX(), output[i].getY(), output[i].getZ() );

    fclose (fp);
  }

  return system ("./test.plot");
}
