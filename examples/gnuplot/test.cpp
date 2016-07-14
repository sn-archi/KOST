#include <fstream>
#include <iostream>

#include "../../src/mKOST.h"

/* Number of output samples */
#define N 100

int main (int argc, char* argv[])
{
  /** Initial state at t=0 */
  mKOST::StateVectors initial;
  initial.pos = btVector3 (-1.000000e+05, -1.000000e+04, -1.000000e+05);
  initial.vel = btVector3 (1.000000e+01, -1.000000e+00, 1.000000e+01);
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
  std::cout << elements << std::endl;

  mKOST::Params params (orbit.getParams());
  std::cout << params << std::endl;

  mKOST::StateVectors out;
  try
  {
    out = orbit.elements2StateVector (initial.MeL, EPSILON, 1000000);
  }
  catch (int errMsg)
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
  std::cout << elements << std::endl;

  params = orbit.getParams();
  std::cout << params << std::endl;

  btScalar maxt (1.0 * params.T);

  btVector3 output[N];
  for (int i (0); i < N; ++i)
  {
    double t ((i * maxt) / N);
    mKOST::StateVectors stateNow;
    try
    {
      stateNow = orbit.elements2StateVectorAtTime (initial.MeL, t, EPSILON, 1000000, 0.0, 0.0);
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
