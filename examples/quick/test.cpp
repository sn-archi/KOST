#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "../../src/mKOST.h"

int main (int argc, char* argv[])
{
  /** Alternative test, with an initial param list at t0 */
  mKOST::Elements elements;
  elements.a = 5.799119e5;
  elements.Ecc = 0.5;
  elements.i = SIMD_PI / 4;
  elements.LoP = SIMD_PI / 4;
  elements.LAN = 0.0;
  elements.L = SIMD_HALF_PI;

  mKOST::Orbit mOrbit (MU, &elements);

  /** Convert to orbital elements */
  mOrbit.elements2StateVector (SIMD_EPSILON, 1000000);

  printf ("Orbital elements:\n"
          "     a = %f m\n"
          "     e = %f\n"
          "     i = %f\n"
          " theta = %f\n"
          "omegab = %f\n"
          "     L = %f\n",
          elements.a, elements.Ecc, elements.i, elements.LAN, elements.LoP, elements.L
         );

  mKOST::Params params (mOrbit.getParams());
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

  mOrbit.stateVector2Elements ();
  mKOST::Elements elements2 (mOrbit.getElements());
  printf ("Orbital elements reversed:\n"
          "     a = %f m\n"
          "     e = %f\n"
          "     i = %f\n"
          " theta = %f\n"
          "omegab = %f\n"
          "     L = %f\n",
          elements2.a, elements2.Ecc, elements2.i, elements2.LAN, elements2.LoP, elements2.L
         );

  params = mOrbit.getParams();
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

  mKOST::StateVectors out (mOrbit.getStateVectors());
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
