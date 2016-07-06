#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "../../src/mKOST.h"

int main (int argc, char* argv[])
{
  /** Alternative test, with an initial param list at t0 */
  mKOST::Elements elements;
  elements.a = 8660.442;
  elements.Ecc = 0.9999710e-01;
  elements.i = 2.356194;
  elements.LoP = 2.526092;
  elements.LAN = 1.570796;

  mKOST::Orbit orbit (MU, elements);

  mKOST::StateVectors out = orbit.elements2StateVector (5.678450, SIMD_EPSILON, 1000000);
  printf ("state:\n"
          "  position: %f, %f, %f\n"
          "  velocity: %f, %f, %f\n",
          out.pos.getX(),
          out.pos.getY(),
          out.pos.getZ(),
          out.vel.getX(),
          out.vel.getY(),
          out.vel.getZ()
          );
  std::cout << elements << std::endl;
  mKOST::Params params (orbit.getParams());
  std::cout << params << std::endl;

  orbit.refreshFromStateVectors (&out);
  elements = orbit.getElements();
  std::cout << elements << std::endl;
  params = orbit.getParams();
  std::cout << params << std::endl;

  out = orbit.elements2StateVector(5.678450, SIMD_EPSILON, 1000);
  printf ("state:\n"
          "  position: %f, %f, %f\n"
          "  velocity: %f, %f, %f\n",
          out.pos.getX(),
          out.pos.getY(),
          out.pos.getZ(),
          out.vel.getX(),
          out.vel.getY(),
          out.vel.getZ()
          );
  return 0;
}
