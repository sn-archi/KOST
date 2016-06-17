#include <fstream>
#include <iostream>

#include "../../src/mKOST.h"

/* Number of output samples */
#define N 100

int main (int argc, char* argv[])
{
  /** Initial state at t=0 */
  mKOST::StateVectors initial;
  initial.pos = btVector3 (-1.000000e+04, -1.000000e+09, -1.000000e+04);
  initial.vel = btVector3 (-1.000000e+00, -1.000000e+00, -1.000000e+03);
  mKOST::Orbit inOrbit;

  /** Convert to orbital elements */
  try
  {
    inOrbit.setMu (MU);
    inOrbit.refreshFromStateVectors(&initial);
  }
  catch (const char*)
  {
    return 1;
  }

  printf ("initial:\n"
        "  position: %e, %e, %e\n"
        "  velocity: %e, %e, %e\n"
        "  Mean Longitude: %e\n-----\n",
        initial.pos.getX(),
        initial.pos.getY(),
        initial.pos.getZ(),
        initial.vel.getX(),
        initial.vel.getY(),
        initial.vel.getZ(),
        initial.MeL
        );

  mKOST::Elements elements (inOrbit.getElements());
  printf ("Orbital elements:\n"
          "     a = %e m\n"
          "     e = %e\n"
          "     i = %e\n"
          "   LaN = %e\n"
          "   LoP = %e\n"
          "   MeL = %e\n",
          elements.a, elements.Ecc, elements.i, elements.LAN, elements.LoP, initial.MeL
         );

  mKOST::Params params (inOrbit.getParams());
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

  mKOST::StateVectors out;
  try
  {
    out = inOrbit.elements2StateVector (initial.MeL, SIMD_EPSILON, 1000000);
  }
  catch (const char*)
  {
    printf ("Couldn't find a root\n");
  }

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

  mKOST::Orbit outOrbit (MU, &out);
  mKOST::Elements elements2 (outOrbit.getElements());
  printf ("Orbital elements:\n"
          "     a = %e m\n"
          "     e = %e\n"
          "     i = %e\n"
          "   LaN = %e\n"
          "   LoP = %e\n"
          "   MeL = %e\n",
          elements2.a,
          elements2.Ecc,
          elements2.i,
          elements2.LAN,
          elements2.LoP,
          out.MeL
          );

  mKOST::Params params2 (outOrbit.getParams());
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

  btScalar maxt (1.0 * params.T);

  btVector3 output[N];
  for (int i (0); i < N; ++i)
  {
    double t ((i * maxt) / N);
    mKOST::StateVectors stateNow;
    try
    {
      stateNow = inOrbit.elements2StateVectorAtTime (t, SIMD_EPSILON, 1000000, 0.0, 0.0);
    }
    catch (const char*)
    {
      printf ("Couldn't find a root\n");
    }

    output[i] = stateNow.pos;
  }

  {
    //FILE* fp = fopen ("orbit.dat", "w");
    std::ofstream fs("orbit.dat");
    // GNUPlot uses Y for depth and Z for height. This'll make things looks as expected
    // Just don't pay attention to the axix names.
    for (int i (0); i < N; ++i)
      //fprintf (fp, "%f\t%f\t%f\n", output[i].getX(), output[i].getZ(), -output[i].getY() );
      fs << output[i].getX() << "\t" << output[i].getZ() << "\t" << -output[i].getY() << std::endl;
    //fclose (fp);
  }

  return system ("./test.plot");
}
