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
  catch (int)
  {
    return 1;
  }

  std::cout << "initial:\n" << initial << std::endl;

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

  std::cout << "reversed:\n" << out << std::endl;

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
    catch (int)
    {
      std::cout << "Couldn't find a root\n" << std::endl;
    }
    output[i] = stateNow.pos;
  }

  std::ofstream fs("orbit.dat");
  for (int i (0); i < N; ++i)
    fs << output[i].getX() << "\t" << output[i].getZ() << "\t" << -output[i].getY() << std::endl;

  return system ("./test.plot");
}
