#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "../../src/mKOST.h"

int main (int argc, char* argv[])
{
  /** Alternative test, with an initial param list at t0 */
  mKOST::Elements elements;
  elements.a = -1926942043.120084;
  elements.Ecc = 1.0078467921690348;
  elements.i = 1.5716962815553079;
  elements.LAN = 4.7223886470713552;
  elements.LoP = 0.18601293306578182;

  mKOST::Orbit orbit (MU, elements);

  mKOST::StateVectors out = orbit.elements2StateVector (0.35485236693591626, 12 * SIMD_EPSILON, 1000000);

  std::cout << "initial:\n" << out << std::endl;

  std::cout << elements << std::endl;
  mKOST::Params params (orbit.getParams());
  std::cout << params << std::endl;

  orbit.refreshFromStateVectors (&out);
  elements = orbit.getElements();
  std::cout << elements << std::endl;
  params = orbit.getParams();
  std::cout << params << std::endl;

  out = orbit.elements2StateVector(0.35485236693591626, 12 * SIMD_EPSILON, 1000000);

  std::cout << "reversed:\n" << out << std::endl;

  return 0;
}
