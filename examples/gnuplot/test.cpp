#include <fstream>
#include <iostream>

#include "../../src/mKOST.h"

/* Number of output samples */
#define N 100

int main (int argc, char* argv[])
{
  /** Initial state at t=0 */
  mKOST::StateVectors initial;
  initial.pos = btVector3 (0.0, 0.0, -1.000000e+04);
  initial.vel = btVector3 (1.000000e+04, 1.0, 0.0);
  mKOST::Orbit orbit;

  /** Convert to orbital elements */
  try
  {
    orbit.setMu (MU);
    orbit.refreshFromStateVectors(&initial);
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

  mKOST::Elements elements (orbit.getElements());
  printf ("Orbital elements:\n"
          "     a = %e m\n"
          "     e = %e\n"
          "     i = %e\n"
          "   LaN = %e\n"
          "   LoP = %e\n"
          "   MeL = %e\n",
          elements.a, elements.Ecc, elements.i, elements.LAN, elements.LoP, initial.MeL
         );

  mKOST::Params params (orbit.getParams());
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
    out = orbit.elements2StateVector (initial.MeL, SIMD_EPSILON, 1000);
  }
  catch (const char * errMsg)
  {
    std::cout << errMsg << std::endl;
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

  orbit.refreshFromStateVectors (&out);
  elements = orbit.getElements();
  printf ("Orbital elements:\n"
          "     a = %e m\n"
          "     e = %e\n"
          "     i = %e\n"
          "   LaN = %e\n"
          "   LoP = %e\n"
          "   MeL = %e\n",
          elements.a,
          elements.Ecc,
          elements.i,
          elements.LAN,
          elements.LoP,
          out.MeL
          );

  params = orbit.getParams();
  printf ("Additional parameters:\n"
          "   PeD = %e m\n"
          "   ApD = %e m\n"
          "   MnA = %e\n"
          "   TrA = %e\n"
          "   EcA = %e\n"
          "     T = %e s\n"
          "   AgP = %e\n",
          params.PeD,
          params.ApD,
          params.MnA,
          params.TrA,
          params.EcA,
          params.T,
          params.AgP
         );

  btScalar maxt (1.0 * params.T);

  btVector3 output[N];
  for (int i (0); i < N; ++i)
  {
    double t ((i * maxt) / N);
    mKOST::StateVectors stateNow;
    try
    {
      stateNow = orbit.elements2StateVectorAtTime (initial.MeL, t, SIMD_EPSILON, 1000, 0.0, 0.0);
    }
    catch (const char*)
    {
      printf ("Couldn't find a root\n");
    }
    output[i] = stateNow.pos;
  }

  std::ofstream fs("orbit.dat");
  for (int i (0); i < N; ++i)
    fs << output[i].getX() << "\t" << output[i].getZ() << "\t" << -output[i].getY() << std::endl;

  return system ("./test.plot");
}
